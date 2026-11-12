"""Web/queue regression tests, including real worker subprocess lifecycle tests."""
import json
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import time
import unittest
import re
from html import unescape
from http.cookies import SimpleCookie
from urllib.parse import parse_qs, urlsplit
from unittest.mock import patch

from fastapi.testclient import TestClient
from main import create_app
from web.jobs import JobSpec, JobStore, Settings, command_for, prepare_assembly, source_catalog, validate_spec
from web_worker import run_job

BASE_URL = 'http://127.0.0.1:8000'


def fake_ccp4(root):
    programs=root/'bin'; programs.mkdir(parents=True)
    (programs/'ccp4.setup-sh').write_text('#!/bin/sh\n')
    aimless=programs/'aimless'; aimless.write_text('#!/bin/sh\nexit 0\n'); aimless.chmod(0o755)
    return root


def inspector(paths, settings):
    return [{'path':str(p),'name':p.name,'ready':True,'frames':180,'companions':[],'error':''} for p in paths]


def processing_statistics():
    return {'resolution':[2.5,10,2.5],'isigi':[10,20,2],'cc_half':[.99,.99,.7],'completeness':[95,99,90]}


def await_status(client, identifier, statuses, timeout=15):
    deadline = time.monotonic()+timeout
    while time.monotonic()<deadline:
        row = client.get('/api/jobs/'+identifier).json()
        if row['status'] in statuses:
            return row
        time.sleep(.1)
    raise AssertionError('Job did not reach %s: %s' % (statuses,row))


class WebTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        root = Path(self.temp.name)
        self.data = root/'data'; self.data.mkdir()
        self.master = self.data/'deep'/'sample_master.h5'
        self.master.parent.mkdir(); self.master.touch()
        self.output = root/'processed'; self.output.mkdir()
        self.settings = Settings(self.data,self.output,root/'jobs',ccp4=fake_ccp4(root/'default ccp4'),testing=True)
        self.client_context = TestClient(create_app(self.settings,start_worker=False,inspector=inspector),base_url=BASE_URL)
        self.client = self.client_context.__enter__()

    def tearDown(self):
        self.client_context.__exit__(None,None,None)
        self.temp.cleanup()

    def payload(self, **extra):
        return {'mode':'process','input_dir':str(self.data),'master_files':[str(self.master)],**extra}

    def processed_result(self, source, identifier='single', master=None, last_frame=1800,
                         frames=None, resolution=1.8, content=b'MTZ source'):
        directory=source/'crystals'/identifier
        mtz=directory/f'prefix_{last_frame:04d}'/'integrated.mtz'
        mtz.parent.mkdir(parents=True); mtz.write_bytes(content)
        prefix={'mtz':str(mtz),'last_frame':last_frame,'statistics':processing_statistics()}
        manifest=directory/'manifest.json'
        manifest.write_text(json.dumps({'config':{'master':str(master or self.master),'frames':frames,'resolution':resolution},
                                       'prefixes':[prefix],'selected':prefix}))
        return manifest, mtz

    def test_python_html_and_no_script_or_templates(self):
        response = self.client.get('/')
        self.assertEqual(response.status_code,200)
        self.assertIn('Process diffraction data',response.text)
        self.assertIn('sample_master.h5',response.text)
        self.assertIn('Assemble processed datasets',response.text)
        self.assertNotIn('<script',response.text.lower())
        self.assertEqual(self.client.get('/docs').status_code,404)

    def test_ccp4_default_is_detected_from_aimless_on_path(self):
        detected=fake_ccp4(self.output.parent/'detected ccp4')
        result=subprocess.CompletedProcess(['which','aimless'],0,str(detected/'bin/aimless')+'\n','')
        with patch.dict(os.environ,{'PYMDA_CCP4':''}),patch('autoDials.subprocess.run',return_value=result) as which:
            settings=Settings(self.data,self.output,self.settings.jobs_root,testing=True)
            self.assertEqual(settings.ccp4,detected)
            with TestClient(create_app(settings,start_worker=False,inspector=inspector),base_url=BASE_URL) as client:
                page=client.get('/')
                self.assertIn(f'id="ccp4" name="ccp4" type="text" value="{detected}"',page.text)
                self.assertIn('Choose CCP4 directory',page.text)
                self.assertIn('Detect from PATH',page.text)
                self.assertNotIn('notice error',page.text)
            which.assert_called_once_with(['which','aimless'],capture_output=True,text=True,timeout=5)

    def test_select_ccp4_directory_controls_discovery_jobs_and_remembered_location(self):
        chosen=fake_ccp4(self.output.parent/'another CCP4')
        inspected=[]
        def record(paths,settings):
            inspected.append(settings.ccp4)
            return inspector(paths,settings)
        with TestClient(create_app(self.settings,start_worker=False,inspector=record),base_url=BASE_URL) as client:
            response=client.get('/browse',params={'kind':'ccp4','path':str(chosen/'bin'),
                'input_dir':str(self.data),'source_dir':str(self.output),'output_dir':str(self.output)})
            self.assertEqual(response.status_code,200,response.text)
            selected=unescape(re.search(r'href="([^"]+)">Use this directory</a>',response.text).group(1))
            self.assertEqual(parse_qs(urlsplit(selected).query)['ccp4'],[str(chosen)])
            page=client.get(selected)
            self.assertEqual(page.status_code,200,page.text)
            self.assertEqual(inspected[-1],chosen)
            self.assertIn(f'id="ccp4" name="ccp4" type="text" value="{chosen}"',page.text)
            self.assertEqual(page.text.count(f'type="hidden" name="ccp4" value="{chosen}"'),4)
            output_link=unescape(re.search(r'href="([^"]+)">Choose output directory</a>',page.text).group(1))
            self.assertEqual(parse_qs(urlsplit(output_link).query)['ccp4'],[str(chosen)])
            response=client.get('/api/datasets',params={'input_dir':str(self.data),'ccp4':str(chosen/'bin')})
            self.assertEqual(response.status_code,200,response.text)
            self.assertEqual(inspected[-1],chosen)
            response=client.post('/jobs',data={'mode':'process','input_dir':str(self.data),
                'master_files':str(self.master),'output_dir':str(self.output),'ccp4':str(chosen/'bin')},follow_redirects=False)
            self.assertEqual(response.status_code,303,response.text)
            self.assertEqual(inspected[-1],chosen)
            job=client.get('/api/jobs/'+response.headers['location'].split('/')[-1]).json()
            self.assertEqual(job['params']['ccp4'],str(chosen))
            args=command_for(JobSpec(**job['params']),Path(job['output_dir']),self.settings)
            self.assertEqual(args[args.index('--ccp4')+1],str(chosen))
            self.assertIn(f'value="{chosen}"',client.get('/').text)
            client.get('/',params={'ccp4':str(self.settings.ccp4)})
            self.assertEqual(client.get('/api/jobs/'+job['id']).json()['params']['ccp4'],str(chosen))
            self.assertNotEqual(self.settings.ccp4,chosen)

    def test_detect_from_path_replaces_saved_ccp4_selection(self):
        detected=fake_ccp4(self.output.parent/'detected ccp4')
        self.client.get('/')
        result=subprocess.CompletedProcess(['which','aimless'],0,str(detected/'bin/aimless')+'\n','')
        with patch('autoDials.subprocess.run',return_value=result):
            page=self.client.get('/',params={'detect_ccp4':'true','ccp4':str(self.settings.ccp4)})
        self.assertNotIn('notice error',page.text)
        self.assertIn(f'id="ccp4" name="ccp4" type="text" value="{detected}"',page.text)
        self.assertIn(f'value="{detected}"',self.client.get('/').text)

    def test_invalid_ccp4_is_rejected_before_output_creation_and_saved_choice_falls_back(self):
        destination=self.output.parent/'unused output'
        response=self.client.post('/api/jobs',json=self.payload(ccp4=str(self.data),output_dir=str(destination)))
        self.assertEqual(response.status_code,422,response.text)
        self.assertIn('No CCP4 setup script',response.text)
        self.assertFalse(destination.exists())
        self.assertEqual(self.client.get('/api/datasets',params={'ccp4':str(self.data)}).status_code,422)
        page=self.client.get('/',params={'ccp4':str(self.data)})
        self.assertIn('No CCP4 setup script',page.text)
        page=self.client.get('/browse',params={'kind':'ccp4','path':str(self.data)})
        self.assertNotIn('Use this directory',page.text)
        self.assertIn('No CCP4 setup script',page.text)
        self.client.cookies.set('pymda_ccp4',str(self.data),domain='127.0.0.1',path='/')
        page=self.client.get('/')
        self.assertNotIn('notice error',page.text)
        self.assertIn(f'id="ccp4" name="ccp4" type="text" value="{self.settings.ccp4}"',page.text)

    def test_app_and_worker_support_manual_ccp4_when_aimless_is_not_on_path(self):
        pipeline=self.output.parent/'fake_pipeline.py'; pipeline.write_text(FAKE_PIPELINE)
        settings=Settings(self.data,self.output,self.output.parent/'manual-ccp4-jobs',ccp4=None,pipeline=pipeline,testing=True)
        with patch.dict(os.environ,{'PATH':'/usr/bin:/bin','PYMDA_CCP4':''}):
            with TestClient(create_app(settings,inspector=inspector),base_url=BASE_URL) as client:
                page=client.get('/')
                self.assertIn('CCP4 was not found by which aimless',page.text)
                self.assertIn('Choose CCP4 directory',page.text)
                response=client.post('/api/jobs',json=self.payload())
                self.assertEqual(response.status_code,422,response.text)
                response=client.post('/api/jobs',json=self.payload(ccp4=str(self.settings.ccp4)))
                self.assertEqual(response.status_code,201,response.text)
                job=await_status(client,response.json()['id'],{'succeeded'})
                self.assertEqual(job['params']['ccp4'],str(self.settings.ccp4))
                self.assertTrue(client.get('/health').json()['worker_running'])

    def test_validation_and_unknown_fields(self):
        for extra in [{'thread':20},{'reso':float('inf')},{'resolutions':[-2]}, {'resolutions':[2.261]}, {'unit_cell':'1,2,3,10,10,160'}, {'command':'touch injected'}]:
            # JSON cannot encode inf; send a string, parsed/validated as a float.
            if extra.get('reso') == float('inf'): extra['reso']='inf'
            self.assertEqual(self.client.post('/api/jobs',json=self.payload(**extra)).status_code,422)
        self.assertEqual(self.client.get('/api/jobs').json(),[])

    def test_reject_unselected_or_undiscovered_files(self):
        other=self.output.parent/'other input'; other.mkdir()
        for extra in [{'master_files':[]},{'input_dir':str(other),'master_files':[str(self.master)]},{'master_files':['/etc/passwd']}]:
            self.assertEqual(self.client.post('/api/jobs',json=self.payload(**extra)).status_code,422)

    def test_external_input_directory_discovery_browsing_and_submission(self):
        directory=self.output.parent/'external input'; directory.mkdir()
        master=directory/'another_master.h5'; master.touch()
        alias=self.data/'linked input'; alias.symlink_to(directory,target_is_directory=True)
        for path in [directory,alias]:
            with self.subTest(path=path):
                response=self.client.get('/api/datasets',params={'input_dir':str(path)})
                self.assertEqual(response.status_code,200,response.text)
                self.assertEqual([item['path'] for item in response.json()],[str(master)])
                page=self.client.get('/',params={'input_dir':str(path)})
                self.assertIn('another_master.h5',page.text)
                self.assertNotIn('notice error',page.text)
                response=self.client.get('/browse',params={'kind':'input_dir','path':str(path),'source_dir':str(self.output),'output_dir':str(self.output)})
                self.assertEqual(response.status_code,200,response.text)
                for label in ['Home folder','Filesystem root','Parent directory']:
                    self.assertIn(label,response.text)
                selected=unescape(re.search(r'href="([^"]+)">Use this directory</a>',response.text).group(1))
                self.assertIn('another_master.h5',self.client.get(selected).text)
                response=self.client.post('/api/jobs',json=self.payload(input_dir=str(path),master_files=[str(master)]))
                self.assertEqual(response.status_code,201,response.text)
                spec=JobSpec(**response.json()['params'])
                args=command_for(spec,Path(response.json()['output_dir']),self.settings)
                self.assertEqual(args[args.index('--hdf5_path')+1],str(directory))
                self.assertEqual(args[args.index('--master_files')+1],str(master))
        response=self.client.get('/browse',params={'kind':'input_dir','path':str(self.data)})
        self.assertIn('linked input',response.text)
        response=self.client.get('/browse',params={'kind':'input_dir','path':str(self.data.parent)})
        self.assertEqual(response.status_code,200,response.text)
        self.assertIn('external input',response.text)

    def test_readable_master_symlink_outside_default_input_root(self):
        directory=self.output.parent/'master targets'; directory.mkdir()
        master=directory/'external_master.h5'; master.touch()
        link=self.data/'external_master.h5'; link.symlink_to(master)
        datasets=self.client.get('/api/datasets',params={'input_dir':str(self.data)}).json()
        self.assertIn(str(master),[item['path'] for item in datasets])
        response=self.client.post('/api/jobs',json=self.payload(master_files=[str(link)]))
        self.assertEqual(response.status_code,201,response.text)
        self.assertEqual(response.json()['params']['master_files'],[str(master)])

    def test_missing_input_and_source_directories_still_report_errors(self):
        missing=self.output.parent/'missing input'
        for kind,api in [('input_dir','/api/datasets'),('source_dir','/api/processed-datasets')]:
            with self.subTest(kind=kind):
                response=self.client.get(api,params={kind:str(missing)})
                self.assertEqual(response.status_code,422,response.text)
                self.assertIn('Directory does not exist',response.text)
                self.assertEqual(self.client.get('/browse',params={'kind':kind,'path':str(missing)}).status_code,422)
                page=self.client.get('/',params={kind:str(missing)})
                self.assertIn('Directory does not exist',page.text)
        self.assertEqual(self.client.post('/api/jobs',json=self.payload(input_dir=str(missing))).status_code,422)
        self.assertEqual(self.client.post('/api/jobs',json={'mode':'assemble','source_dir':str(missing)}).status_code,422)

    def test_input_and_processed_directories_require_read_and_traversal_access(self):
        source=self.output/'source'; self.processed_result(source)
        real_access=os.access
        for directory,kind,api in [(self.data,'input_dir','/api/datasets'),(source,'source_dir','/api/processed-datasets')]:
            for denied in [os.R_OK,os.X_OK]:
                with self.subTest(directory=directory,denied=denied):
                    def access(path,mode):
                        return False if Path(path)==directory and mode & denied else real_access(path,mode)
                    with patch('web.jobs.os.access',side_effect=access):
                        response=self.client.get(api,params={kind:str(directory)})
                        self.assertEqual(response.status_code,422,response.text)
                        self.assertIn('not readable',response.text)
                        self.assertEqual(self.client.get('/browse',params={'kind':kind,'path':str(directory)}).status_code,422)
                        self.assertIn('not readable',self.client.get('/',params={kind:str(directory)}).text)
                        payload=self.payload() if kind=='input_dir' else {'mode':'assemble','source_dir':str(directory)}
                        self.assertEqual(self.client.post('/api/jobs',json=payload).status_code,422)

    def test_selected_master_requires_read_access(self):
        real_access=os.access
        with patch('web.jobs.os.access',side_effect=lambda path,mode: False if Path(path)==self.master and mode & os.R_OK else real_access(path,mode)):
            response=self.client.post('/api/jobs',json=self.payload())
            self.assertEqual(response.status_code,422,response.text)
            self.assertIn('Master file is not readable',response.text)

    def test_single_dataset_finalization_and_directory_buttons(self):
        result = self.client.post('/api/jobs',json=self.payload(mode='both', output_dir=str(self.output)))
        self.assertEqual(result.status_code,201)
        self.assertTrue(Path(result.json()['output_dir']).is_relative_to(self.output))
        response = self.client.get('/browse',params={'kind':'output_dir','path':str(self.output)})
        self.assertEqual(response.status_code,200)
        self.assertIn('Use this directory',response.text)
        self.assertEqual(self.client.get('/browse',params={'kind':'output_dir','path':str(self.master)}).status_code,422)
        self.assertEqual(self.client.post('/directories',data={'path':str(self.output),'name':'../escape'}).status_code,422)
        self.assertEqual(self.client.post('/directories',data={'path':str(self.output),'name':'new'},follow_redirects=False).status_code,303)
        self.assertTrue((self.output/'new').is_dir())
        self.assertIn('Choose output directory',self.client.get('/').text)

    def test_select_external_output_directory_and_remember_choice(self):
        directory=self.output.parent/'chosen results'; directory.mkdir()
        response=self.client.get('/browse',params={'kind':'output_dir','path':str(directory),'input_dir':str(self.data),'source_dir':str(self.output)})
        self.assertEqual(response.status_code,200)
        self.assertIn('Folder path',response.text)
        self.assertIn('Home folder',response.text)
        selected=unescape(re.search(r'href="([^"]+)">Use this directory</a>',response.text).group(1))
        page=self.client.get(selected)
        self.assertEqual(page.status_code,200)
        self.assertIn(f'id="process_output_dir" name="output_dir" type="text" value="{directory}"',page.text)
        self.assertIn(f'id="assembly_output_dir" name="output_dir" type="text" value="{directory}"',page.text)
        self.assertEqual(page.text.count('Choose output directory'),3)
        self.assertIn(str(directory),self.client.get('/').text)
        response=self.client.post('/jobs',data={'mode':'process','input_dir':str(self.data),'master_files':str(self.master),'output_dir':str(directory)},follow_redirects=False)
        self.assertEqual(response.status_code,303)
        job=self.client.get('/api/jobs/'+response.headers['location'].split('/')[-1]).json()
        self.assertEqual(Path(job['output_dir']).parent.parent,directory)
        self.assertEqual(job['params']['output_dir'],str(directory))

    def test_create_output_folder_outside_default_root_and_reject_invalid_locations(self):
        parent=self.output.parent/'external'; parent.mkdir()
        response=self.client.post('/directories',data={'path':str(parent),'name':'new results','input_dir':str(self.data),'source_dir':str(self.output)},follow_redirects=False)
        self.assertEqual(response.status_code,303)
        self.assertTrue((parent/'new results').is_dir())
        self.assertIn('source_dir=',response.headers['location'])
        self.assertIn('input_dir=',response.headers['location'])
        for path in [str(self.master),str(self.master/'child')]:
            self.assertEqual(self.client.post('/api/jobs',json=self.payload(output_dir=path)).status_code,422)
        real_access=os.access
        with patch('web.jobs.os.access',side_effect=lambda path,mode: False if Path(path)==parent and mode & os.W_OK else real_access(path,mode)):
            response=self.client.post('/api/jobs',json=self.payload(output_dir=str(parent)))
            self.assertEqual(response.status_code,422)
            self.assertIn('not writable',response.text)
            response=self.client.post('/api/jobs',json=self.payload(output_dir=str(parent/'missing'/'nested')))
            self.assertEqual(response.status_code,422,response.text)
            self.assertFalse((parent/'missing').exists())
            page=self.client.get('/browse',params={'kind':'output_dir','path':str(parent)})
            self.assertEqual(page.status_code,200,page.text)
            self.assertNotIn('Use this directory',page.text)
            self.assertIn('not writable',page.text)

    def test_new_output_path_is_created_only_for_a_valid_submission(self):
        directory=self.output.parent/'new output'/'nested results'
        page=self.client.get('/',params={'output_dir':str(directory)})
        self.assertEqual(page.status_code,200,page.text)
        self.assertIn(f'value="{directory}"',page.text)
        self.assertNotIn('notice error',page.text)
        self.assertFalse(directory.parent.exists())
        response=self.client.get('/browse',params={'kind':'output_dir','path':str(directory)})
        self.assertEqual(response.status_code,200,response.text)
        self.assertIn('does not exist yet',response.text)
        self.assertFalse(directory.parent.exists())
        response=self.client.post('/api/jobs',json=self.payload(output_dir=str(directory),split_wedges=3,frames=2))
        self.assertEqual(response.status_code,422,response.text)
        self.assertFalse(directory.parent.exists())
        unreadable=lambda paths,settings:[dict(d,ready=False) for d in inspector(paths,settings)]
        with TestClient(create_app(self.settings,start_worker=False,inspector=unreadable),base_url=BASE_URL) as client:
            response=client.post('/api/jobs',json=self.payload(output_dir=str(directory)))
            self.assertEqual(response.status_code,422,response.text)
            self.assertFalse(directory.parent.exists())
        response=self.client.post('/jobs',data={'mode':'process','input_dir':str(self.data),'master_files':str(self.master),'output_dir':str(directory)},follow_redirects=False)
        self.assertEqual(response.status_code,303,response.text)
        self.assertTrue(directory.is_dir())
        job=self.client.get('/api/jobs/'+response.headers['location'].split('/')[-1]).json()
        self.assertEqual(job['params']['output_dir'],str(directory))
        self.assertEqual(Path(job['output_dir']).parent.parent,directory)
        self.assertIn(f'value="{directory}"',self.client.get('/').text)

    def test_deleted_or_inaccessible_remembered_directories_fall_back(self):
        missing=self.output.parent/'test1'
        blocked=self.output.parent/'inaccessible'; blocked.mkdir()
        real_access=os.access
        for directory in [missing,blocked]:
            with self.subTest(directory=directory):
                self.client.cookies.set('pymda_output_dir',str(directory),domain='127.0.0.1',path='/')
                self.client.cookies.set('pymda_source_dir',str(directory),domain='127.0.0.1',path='/')
                with patch('web.jobs.os.access',side_effect=lambda path,mode: False if Path(path)==blocked else real_access(path,mode)):
                    page=self.client.get('/')
                self.assertEqual(page.status_code,200,page.text)
                self.assertNotIn('notice error',page.text)
                self.assertIn(f'id="source_dir" name="source_dir" value="{self.output}"',page.text)
                self.assertIn(f'id="process_output_dir" name="output_dir" type="text" value="{self.output}"',page.text)
                cookies=SimpleCookie()
                for header in page.headers.get_list('set-cookie'):
                    cookies.load(header)
                self.assertEqual(cookies['pymda_output_dir'].value,str(self.output))
                self.assertEqual(cookies['pymda_source_dir'].value,str(self.output))

    def test_missing_startup_locations_fall_back_for_workspace_browser_api_and_jobs(self):
        missing=self.output.parent/'removed defaults'
        settings=Settings(self.data,missing,self.settings.jobs_root,ccp4=self.settings.ccp4,testing=True)
        self.processed_result(settings.jobs_root/'previous job'/'output')
        with patch.dict(os.environ,{'PYMDA_START_INPUT_DIR':str(missing),'PYMDA_START_SOURCE_DIR':str(missing)}):
            with TestClient(create_app(settings,start_worker=False,inspector=inspector),base_url=BASE_URL) as client:
                page=client.get('/')
                self.assertEqual(page.status_code,200,page.text)
                self.assertNotIn('notice error',page.text)
                self.assertIn('sample_master.h5',page.text)
                self.assertIn('Open this result directory',page.text)
                self.assertIn(f'id="process_output_dir" name="output_dir" type="text" value="{settings.jobs_root}"',page.text)
                for kind in ['input_dir','source_dir','output_dir']:
                    self.assertEqual(client.get('/browse',params={'kind':kind}).status_code,200)
                self.assertEqual(client.get('/api/datasets').json()[0]['path'],str(self.master))
                self.assertTrue(client.get('/api/processed-datasets').json()[0]['ready'])
                response=client.post('/api/jobs',json={'mode':'process','master_files':[str(self.master)]})
                self.assertEqual(response.status_code,201,response.text)
                self.assertEqual(response.json()['params']['input_dir'],str(self.data))
                response=client.post('/api/jobs',json={'mode':'assemble'})
                self.assertEqual(response.status_code,201,response.text)
                self.assertEqual(response.json()['params']['source_dir'],str(settings.jobs_root))
                self.assertFalse(missing.exists())

    def test_worker_reports_output_creation_failure_in_job_state(self):
        job=self.client.post('/api/jobs',json=self.payload()).json()
        store=self.client.app.state.store
        job=store.claim()
        with patch('web_worker.Path.mkdir',side_effect=PermissionError('Cannot create job output')):
            run_job(store,job)
        result=store.get(job['id'])
        self.assertEqual(result['status'],'failed')
        self.assertIn('Cannot create job output',result['message'])

    def test_typed_output_directory_is_remembered_after_submission(self):
        directory=self.output.parent/'typed output'; directory.mkdir()
        self.client.get('/')
        response=self.client.post('/jobs',data={'mode':'process','input_dir':str(self.data),'master_files':str(self.master),'output_dir':str(directory)},follow_redirects=False)
        self.assertEqual(response.status_code,303)
        self.assertIn(f'value="{directory}"',self.client.get('/').text)

    def test_processed_data_can_be_reused_from_selected_output_location(self):
        source=self.output.parent/'custom outputs'/'job'/'output'
        directory=source/'crystals'/'single'; directory.mkdir(parents=True)
        mtz=directory/'integrated.mtz'; mtz.write_bytes(b'MTZ source')
        prefix={'mtz':str(mtz),'last_frame':180,'statistics':processing_statistics()}
        (directory/'manifest.json').write_text(json.dumps({'config':{'master':str(self.master),'resolution':1.8},'prefixes':[prefix],'selected':prefix}))
        page=self.client.get('/',params={'source_dir':str(source)})
        self.assertEqual(page.status_code,200)
        self.assertIn('value="single"',page.text)
        response=self.client.post('/api/jobs',json={'mode':'assemble','source_dir':str(source),'dataset_ids':['single']})
        self.assertEqual(response.status_code,201,response.text)
        job=response.json(); prepare_assembly(job,JobSpec(**job['params']))
        copied=Path(job['output_dir'])/'crystals'/'single'/'integrated.mtz'
        self.assertEqual(copied.read_bytes(),b'MTZ source')

    def test_parent_output_and_job_folders_discover_processed_results(self):
        source=self.output/'chosen results'
        result=source/'job-one'/'output'
        self.processed_result(result)
        for folder in [source,source/'job-one',result]:
            with self.subTest(folder=folder):
                response=self.client.get('/api/processed-datasets',params={'source_dir':str(folder)})
                self.assertEqual(response.status_code,200,response.text)
                entries=response.json(); self.assertEqual(len(entries),1)
                self.assertTrue(entries[0]['ready'])
                self.assertEqual(entries[0]['source_dir'],str(result))
                page=self.client.get('/',params={'source_dir':str(folder),'output_dir':str(source)})
                self.assertIn('sample_master.h5',page.text)
                self.assertIn(str(result),page.text)
                self.assertIn('Open this result directory',page.text)
        self.assertEqual(entries[0]['id'],'single')
        parent_page=self.client.get('/',params={'source_dir':str(source),'output_dir':str(source)})
        link=unescape(re.search(r'href="([^"]+)">Find processed datasets in this output folder</a>',parent_page.text).group(1))
        self.assertIn('Open this result directory',self.client.get(link).text)
        self.assertIn(f'id="source_dir" name="source_dir" value="{source}"',self.client.get('/').text)

    def test_repeated_runs_have_stable_ids_and_one_default_selection_per_scan(self):
        source=self.output/'runs'
        first,_=self.processed_result(source/'job-one'/'output',content=b'MTZ first')
        second,_=self.processed_result(source/'job-two'/'output',content=b'MTZ second')
        os.utime(first,ns=(1_000_000_000,1_000_000_000))
        os.utime(second,ns=(2_000_000_000,2_000_000_000))
        before={entry['source_dir']:entry['id'] for entry in source_catalog(source)}
        self.assertEqual(len(set(before.values())),2)
        self.processed_result(source/'pilot'/'output',frames=180,last_frame=180)
        other=self.data/'other_master.h5'
        self.processed_result(source/'job-three'/'output',identifier='other',master=other)
        response=self.client.get('/api/processed-datasets',params={'source_dir':str(source)})
        entries=response.json()
        after={entry['source_dir']:entry['id'] for entry in entries}
        self.assertTrue(all(after[key]==value for key,value in before.items()))
        preferred=[entry for entry in entries if entry['default_selected']]
        self.assertEqual(len(preferred),2)
        self.assertIn(before[str(second.parent.parent.parent)],[entry['id'] for entry in preferred])
        page=self.client.get('/',params={'source_dir':str(source)})
        checkboxes=re.findall(r'<input[^>]*name="dataset_ids"[^>]*>',page.text)
        self.assertEqual(len(checkboxes),4)
        self.assertEqual(sum('checked' in tag for tag in checkboxes),2)
        for mode in ['assemble','compare']:
            response=self.client.post('/api/jobs',json={'mode':mode,'source_dir':str(source)})
            self.assertEqual(response.status_code,201,response.text)
            self.assertEqual(set(response.json()['params']['dataset_ids']),{entry['id'] for entry in preferred})

    def test_crystal_and_export_folders_discover_the_containing_dataset(self):
        result=self.output/'job'/'output'
        manifest,mtz=self.processed_result(result)
        self.processed_result(result,identifier='other',master=self.data/'other_master.h5')
        for folder,expected in [(result/'crystals',{'single','other'}),(manifest.parent,{'single'}),(mtz.parent,{'single'})]:
            with self.subTest(folder=folder):
                response=self.client.get('/api/processed-datasets',params={'source_dir':str(folder)})
                self.assertEqual(response.status_code,200,response.text)
                self.assertEqual({entry['id'] for entry in response.json()},expected)
                response=self.client.post('/api/jobs',json={'mode':'assemble','source_dir':str(folder),'dataset_ids':['single']})
                self.assertEqual(response.status_code,201,response.text)
                job=response.json(); prepare_assembly(job,JobSpec(**job['params']))
                copied=json.loads((Path(job['output_dir'])/'crystals/single/manifest.json').read_text())
                self.assertEqual(Path(copied['selected']['mtz']).read_bytes(),b'MTZ source')

    def test_assembly_form_defaults_respect_the_selected_integration_limit(self):
        for limit,cutoff in [(3.0,'3.0'),(3.1416,'3.15')]:
            with self.subTest(limit=limit):
                source=self.output/f'limit_{limit}'
                self.processed_result(source,resolution=limit)
                page=self.client.get('/',params={'source_dir':str(source)})
                section=page.text.split('id="processed-datasets"',1)[1]
                merging=re.search(r'name="resolutions"[^>]*value="([^"]+)"',section).group(1)
                comparison=re.search(r'name="reso_cchalf"[^>]*value="([^"]+)"',section).group(1)
                self.assertEqual(merging,cutoff); self.assertEqual(comparison,cutoff)
                response=self.client.post('/jobs',data={'mode':'assemble','source_dir':str(source),'dataset_ids':['single'],
                    'resolutions':merging,'reso_cchalf':comparison,'output_dir':str(self.output)},follow_redirects=False)
                self.assertEqual(response.status_code,303,response.text)

    def test_selecting_nested_runs_copies_exact_exports_and_preserves_sources(self):
        source=self.output/'runs'
        first,_=self.processed_result(source/'job-one'/'output',content=b'MTZ first')
        second,_=self.processed_result(source/'job-two'/'output',content=b'MTZ second')
        originals={path:path.read_bytes() for path in [first,second]}
        entries=source_catalog(source)
        ids=[entry['id'] for entry in entries]
        response=self.client.post('/api/jobs',json={'mode':'assemble','source_dir':str(source),'dataset_ids':ids})
        self.assertEqual(response.status_code,201,response.text)
        job=response.json(); prepare_assembly(job,JobSpec(**job['params']))
        output=Path(job['output_dir'])
        self.assertEqual(len(list(output.glob('crystals/*/manifest.json'))),2)
        for entry,content in zip(entries,[b'MTZ first',b'MTZ second']):
            manifest=output/'crystals'/entry['id']/'manifest.json'
            copied=json.loads(manifest.read_text())
            self.assertEqual(copied['dataset_id'],entry['id'])
            self.assertEqual(copied['source_manifest'],str(entry['manifest']))
            self.assertEqual(Path(copied['selected']['mtz']).read_bytes(),content)
            self.assertEqual(copied['selected']['mtz'],copied['prefixes'][0]['mtz'])
            self.assertTrue(Path(copied['selected']['mtz']).is_relative_to(manifest.parent))
        selection=json.loads((output/'selection.json').read_text())
        self.assertEqual(selection['source_manifests'],{entry['id']:str(entry['manifest']) for entry in entries})
        self.assertTrue(all(path.read_bytes()==content for path,content in originals.items()))

    def test_incomplete_results_do_not_hide_or_select_valid_results(self):
        source=self.output/'runs'
        self.processed_result(source/'good'/'output',identifier='good')
        missing,mtz=self.processed_result(source/'missing'/'output',identifier='missing'); mtz.unlink()
        invalid,_=self.processed_result(source/'invalid'/'output',identifier='invalid'); invalid.write_text('{')
        malformed,_=self.processed_result(source/'malformed'/'output',identifier='malformed')
        record=json.loads(malformed.read_text()); record['prefixes'][0]['last_frame']='unfinished'
        malformed.write_text(json.dumps(record))
        entries=self.client.get('/api/processed-datasets',params={'source_dir':str(source)}).json()
        self.assertEqual(len(entries),4)
        ready=[entry for entry in entries if entry['ready']]
        self.assertEqual(len(ready),1); self.assertTrue(ready[0]['default_selected'])
        page=self.client.get('/',params={'source_dir':str(source)})
        checkboxes=re.findall(r'<input[^>]*name="dataset_ids"[^>]*>',page.text)
        self.assertEqual(sum('disabled' in tag for tag in checkboxes),3)
        self.assertIn('Incomplete or unavailable processing result',page.text)
        response=self.client.post('/api/jobs',json={'mode':'assemble','source_dir':str(source)})
        self.assertEqual(response.status_code,201,response.text)
        self.assertEqual(response.json()['params']['dataset_ids'],[ready[0]['id']])
        unavailable=next(entry for entry in entries if entry['dataset_id']=='missing')
        self.assertEqual(self.client.post('/api/jobs',json={'mode':'assemble','source_dir':str(source),'dataset_ids':[unavailable['id']]}).status_code,422)

    def test_moved_and_relative_export_paths_can_be_reused(self):
        source=self.output/'moved results'
        absolute,_=self.processed_result(source/'absolute'/'output',identifier='absolute',content=b'MTZ moved')
        relative,_=self.processed_result(source/'relative'/'output',identifier='relative',content=b'MTZ relative')
        for manifest in [absolute,relative]:
            record=json.loads(manifest.read_text())
            path=Path(record['selected']['mtz']).relative_to(manifest.parent)
            reference=str(Path('/unavailable/old/output/crystals')/manifest.parent.name/path) if manifest==absolute else str(path)
            record['selected']['mtz']=reference; record['prefixes'][0]['mtz']=reference
            manifest.write_text(json.dumps(record))
        entries=source_catalog(source)
        self.assertTrue(all(entry['ready'] for entry in entries))
        response=self.client.post('/api/jobs',json={'mode':'assemble','source_dir':str(source),'dataset_ids':[entry['id'] for entry in entries]})
        self.assertEqual(response.status_code,201,response.text)
        job=response.json(); prepare_assembly(job,JobSpec(**job['params']))
        contents=[]
        for manifest in Path(job['output_dir']).glob('crystals/*/manifest.json'):
            record=json.loads(manifest.read_text())
            contents.append(Path(record['selected']['mtz']).read_bytes())
        self.assertCountEqual(contents,[b'MTZ moved',b'MTZ relative'])

    def test_symlinked_crystal_directory_can_be_copied(self):
        original,_=self.processed_result(self.output/'original')
        source=self.output/'linked result'; (source/'crystals').mkdir(parents=True)
        (source/'crystals'/'single').symlink_to(original.parent,target_is_directory=True)
        response=self.client.post('/api/jobs',json={'mode':'assemble','source_dir':str(source),'dataset_ids':['single']})
        self.assertEqual(response.status_code,201,response.text)
        job=response.json(); prepare_assembly(job,JobSpec(**job['params']))
        record=json.loads((Path(job['output_dir'])/'crystals/single/manifest.json').read_text())
        self.assertEqual(Path(record['selected']['mtz']).read_bytes(),b'MTZ source')

    def test_consecutive_wedge_field_is_validated_and_passed_to_pipeline(self):
        self.assertIn('Consecutive wedges per scan',self.client.get('/').text)
        response=self.client.post('/jobs',data={'mode':'process','input_dir':str(self.data),'master_files':str(self.master),'split_wedges':'3'},follow_redirects=False)
        self.assertEqual(response.status_code,303)
        job=self.client.get('/api/jobs/'+response.headers['location'].split('/')[-1]).json()
        spec=JobSpec(**job['params'])
        args=command_for(spec,Path(job['output_dir']),self.settings)
        self.assertEqual(args[args.index('--split_wedges')+1],'3')
        for extra in [{'split_wedges':0},{'split_wedges':101},{'split_wedges':2.5},{'split_wedges':3,'frames':2}]:
            self.assertEqual(self.client.post('/api/jobs',json=self.payload(**extra)).status_code,422)
        small=lambda paths,settings:[dict(d,frames=2) for d in inspector(paths,settings)]
        with TestClient(create_app(self.settings,start_worker=False,inspector=small),base_url=BASE_URL) as client:
            self.assertEqual(client.post('/api/jobs',json=self.payload(split_wedges=3)).status_code,422)

    def test_wedge_exports_are_copied_displayed_and_downloadable(self):
        source=self.output/'source'; directory=source/'crystals'/'single'; directory.mkdir(parents=True)
        mtz=directory/'full.mtz'; mtz.write_bytes(b'MTZ full scan')
        stats={'resolution':[2.5,10,2.5],'isigi':[10,20,2],'cc_half':[.99,.99,.7],'completeness':[95,99,90]}
        prefix={'mtz':str(mtz),'last_frame':1800,'statistics':stats}
        wedges=[]
        for number,(first,last) in enumerate([(1,600),(601,1200),(1201,1800)],1):
            target=directory/f'wedge_{number}.mtz'; target.write_bytes(b'MTZ wedge')
            wedges.append({'number':number,'first_frame':first,'last_frame':last,'used_first_frame':max(2,first),'used_last_frame':last,'mtz':str(target),'statistics':stats,'error':''})
        original={'config':{'master':str(self.master),'resolution':1.8},'prefixes':[prefix],'selected':prefix,'split_wedges':3,'scan_wedges':wedges}
        manifest=directory/'manifest.json'; manifest.write_text(json.dumps(original))
        response=self.client.post('/api/jobs',json={'mode':'assemble','source_dir':str(source),'dataset_ids':['single']})
        self.assertEqual(response.status_code,201,response.text)
        job=response.json(); prepare_assembly(job,JobSpec(**job['params']))
        output=Path(job['output_dir'])
        record=json.loads((output/'crystals/single/manifest.json').read_text())
        self.assertEqual(json.loads(manifest.read_text()),original)
        self.assertTrue(all(Path(w['mtz']).is_relative_to(output) and Path(w['mtz']).is_file() for w in record['scan_wedges']))
        page=self.client.get('/jobs/'+job['id']).text
        self.assertIn('Consecutive scan wedges',page)
        for frames in ['1–600','601–1200','1201–1800']:
            self.assertIn(frames,page)
        self.assertIn('Integrated frames 2–600',page)
        relative=Path(record['scan_wedges'][0]['mtz']).relative_to(output)
        download=self.client.get('/jobs/'+job['id']+'/file',params={'path':str(relative)})
        self.assertEqual(download.content,b'MTZ wedge')

    def test_more_than_ten_datasets_can_process_compare_and_assemble(self):
        masters=[]
        for index in range(12):
            path=self.data/f'other_{index}_master.h5'; path.touch(); masters.append(str(path))
        response=self.client.post('/api/jobs',json=self.payload(mode='both',master_files=masters,max_subsets=8))
        self.assertEqual(response.status_code,201,response.text)
        spec=JobSpec(**response.json()['params'])
        args=command_for(spec,Path(response.json()['output_dir']),self.settings)
        self.assertEqual(args[args.index('--subset_search')+1],'auto')
        self.assertEqual(args[args.index('--max_subsets')+1],'8')
        source=self.output/'other_collection'
        for index in range(12):
            directory=source/'crystals'/f'other_{index}'; directory.mkdir(parents=True)
            mtz=directory/'integrated.mtz'; mtz.write_bytes(b'MTZ source')
            prefix={'mtz':str(mtz),'last_frame':180,'statistics':processing_statistics()}
            (directory/'manifest.json').write_text(json.dumps({'config':{'master':masters[index],'resolution':1.8},'prefixes':[prefix],'selected':prefix}))
        for mode in ['compare','assemble']:
            response=self.client.post('/api/jobs',json={'mode':mode,'source_dir':str(source)})
            self.assertEqual(response.status_code,201,response.text)
            self.assertEqual(len(response.json()['params']['dataset_ids']),12)
        self.assertEqual(self.client.post('/api/jobs',json={'mode':'assemble','source_dir':str(source),'subset_search':'exhaustive'}).status_code,422)

    def test_failed_quality_screen_exposes_scaling_statistics(self):
        job=self.client.post('/api/jobs',json=self.payload()).json()
        output=Path(job['output_dir']); output.mkdir(parents=True)
        (output/'trials.json').write_text(json.dumps([{'prefix_policy':'full','eligible':False,'statistics':{'resolution':[2.5,10,2.5],'isigi':[3,6,2],'cc_half':[.95,.95,.9],'completeness':[61.5,70,63]}}]))
        self.client.app.state.store.update(job['id'],status='failed')
        page=self.client.get('/jobs/'+job['id']).text
        self.assertIn('No trial met the configured quality thresholds',page)
        self.assertIn('61.5%',page)
        self.assertIn('Trial statistics JSON',page)

    def test_local_host_and_origin(self):
        self.assertEqual(self.client.get('/',headers={'host':'evil.example'}).status_code,403)
        self.assertEqual(self.client.post('/api/jobs',json=self.payload(),headers={'origin':'http://evil.example'}).status_code,403)
        self.assertEqual(self.client.post('/api/jobs',json=self.payload(),headers={'origin':BASE_URL}).status_code,201)

    def test_queued_cancel_resume_and_isolated_paths(self):
        first = self.client.post('/api/jobs',json=self.payload()).json()
        second = self.client.post('/api/jobs',json=self.payload()).json()
        self.assertNotEqual(first['output_dir'],second['output_dir'])
        self.assertEqual(self.client.post('/api/jobs/'+first['id']+'/cancel').json()['status'],'cancelled')
        self.assertEqual(self.client.post('/api/jobs/'+first['id']+'/resume').json()['status'],'queued')
        store = self.client.app.state.store
        claimed = store.claim()
        self.assertEqual(claimed['id'],first['id'])
        self.assertIsNone(store.claim())

    def test_forms_redirect_without_javascript(self):
        result = self.client.post('/jobs',data={'mode':'process','input_dir':str(self.data),'master_files':str(self.master),'resolutions':'2.5, 3.0','frames':'','nres':''},follow_redirects=False)
        self.assertEqual(result.status_code,303)
        detail = self.client.get(result.headers['location'])
        self.assertIn('http-equiv="refresh"',detail.text)
        self.assertNotIn('<script',detail.text.lower())
        self.assertEqual(self.client.post('/jobs',data={'mode':'assemble','dataset_ids':[]}).status_code,422)
        self.assertEqual(self.client.post('/jobs',data={'mode':'compare','dataset_ids':[]}).status_code,422)

    def test_log_escaping_and_artifact_confinement(self):
        job = self.client.post('/api/jobs',json=self.payload()).json()
        output = Path(job['output_dir']);output.mkdir(parents=True)
        (output/'job.stdout.log').write_text('<script>alert(1)</script>')
        (output/'escaped.log').symlink_to('/etc/passwd')
        detail = self.client.get('/jobs/'+job['id']).text
        self.assertIn('&lt;script&gt;',detail)
        self.assertNotIn('<script>',detail)
        base = '/jobs/'+job['id']+'/file'
        self.assertEqual(self.client.get(base,params={'path':'../../../../etc/passwd'}).status_code,404)
        self.assertEqual(self.client.get(base,params={'path':'escaped.log'}).status_code,404)
        self.assertEqual(self.client.get(base,params={'path':'job.stdout.log'}).text,'<script>alert(1)</script>')
        self.assertEqual(self.client.get('/jobs/'+job['id'],params={'log_file':'/etc/passwd'}).status_code,404)

    def test_command_arguments_are_not_shell_text(self):
        spec = validate_spec(JobSpec(**self.payload(spg='P 1',unit_cell='35,48,77,85,89,90')),self.settings)
        args = command_for(spec,self.settings.jobs_root/'job',self.settings)
        self.assertIn('P 1',args)
        self.assertIn('--master_files',args)
        self.assertNotIn('shell=True',' '.join(args))


FAKE_PIPELINE = '''import argparse,json,os,subprocess,sys,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--dir_output');p.add_argument('--frames',type=int,default=2);a,_=p.parse_known_args();o=Path(a.dir_output);o.mkdir(parents=True,exist_ok=True)
print('<script>tool output</script>',flush=True)
if a.frames==1:sys.exit(7)
if a.frames==99:
 child=subprocess.Popen([sys.executable,'-c','import signal,time;signal.signal(signal.SIGTERM,signal.SIG_IGN);time.sleep(120)'])
 (o/'child.pid').write_text(str(child.pid));time.sleep(120)
else:
 time.sleep(.3);(o/'crystals'/'fake').mkdir(parents=True,exist_ok=True);(o/'crystals'/'fake'/'manifest.json').write_text('{}');sys.exit(0)
'''


class WorkerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory();root=Path(self.temp.name)
        self.data=root/'data';self.data.mkdir();self.master=self.data/'sample_master.h5';self.master.touch()
        pipeline=root/'fake_pipeline.py';pipeline.write_text(FAKE_PIPELINE)
        self.settings=Settings(self.data,root/'processed',root/'jobs',ccp4=fake_ccp4(root/'default ccp4'),pipeline=pipeline,testing=True)
        self.context=TestClient(create_app(self.settings,inspector=inspector),base_url=BASE_URL)
        self.client=self.context.__enter__()

    def tearDown(self):
        self.context.__exit__(None,None,None);self.temp.cleanup()

    def submit(self, frames=2):
        response=self.client.post('/api/jobs',json={'mode':'process','input_dir':str(self.data),'master_files':[str(self.master)],'frames':frames})
        self.assertEqual(response.status_code,201,response.text)
        return response.json()

    def test_worker_success_and_failure(self):
        job=self.submit();self.assertEqual(await_status(self.client,job['id'],{'succeeded'})['status'],'succeeded')
        failed=self.submit(1);row=await_status(self.client,failed['id'],{'failed'})
        self.assertIn('7',row['message'])

    def test_worker_processes_external_input_into_new_output_directory(self):
        directory=self.data.parent/'external data'; directory.mkdir()
        master=directory/'external_master.h5'; master.touch()
        output=self.data.parent/'external results'/'new collection'
        chosen=fake_ccp4(self.data.parent/'chosen ccp4')
        response=self.client.post('/api/jobs',json={'mode':'process','input_dir':str(directory),'master_files':[str(master)],'output_dir':str(output),'ccp4':str(chosen/'bin')})
        self.assertEqual(response.status_code,201,response.text)
        job=await_status(self.client,response.json()['id'],{'succeeded'})
        self.assertTrue(Path(job['output_dir']).is_relative_to(output))
        self.assertTrue((Path(job['output_dir'])/'crystals/fake/manifest.json').is_file())
        args=json.loads((Path(job['output_dir'])/'web_command.json').read_text())
        self.assertEqual(args[args.index('--hdf5_path')+1],str(directory))
        self.assertEqual(args[args.index('--master_files')+1],str(master))
        self.assertEqual(args[args.index('--ccp4')+1],str(chosen))
        self.assertEqual(job['params']['ccp4'],str(chosen))

    def child_for(self, job):
        path=Path(job['output_dir'])/'child.pid';deadline=time.monotonic()+15
        while not path.is_file() and time.monotonic()<deadline:time.sleep(.1)
        self.assertTrue(path.is_file())
        return int(path.read_text())

    def assert_not_running(self,pid):
        deadline=time.monotonic()+5
        while time.monotonic()<deadline:
            path=Path('/proc')/str(pid)/'stat'
            if not path.exists() or path.read_text().rsplit(')',1)[1].split()[0]=='Z':return
            time.sleep(.1)
        self.fail('Descendant process is still running')

    def test_cancellation_kills_descendants_and_preserves_queue(self):
        first=self.submit(99);child=self.child_for(first);second=self.submit()
        self.assertEqual(self.client.get('/api/jobs/'+second['id']).json()['status'],'queued')
        self.client.post('/api/jobs/'+first['id']+'/cancel')
        await_status(self.client,first['id'],{'cancelled'});self.assert_not_running(child)
        await_status(self.client,second['id'],{'succeeded'})

    def test_worker_crash_recovery_and_resume(self):
        job=self.submit(99);child=self.child_for(job)
        worker=self.client.app.state.worker
        worker.kill();worker.wait()
        recovered=await_status(self.client,job['id'],{'interrupted'},timeout=20)
        self.assertIn('restarted',recovered['message'])
        self.assert_not_running(child)
        self.assertEqual(self.client.post('/api/jobs/'+job['id']+'/resume').json()['status'],'queued')
        await_status(self.client,job['id'],{'running'})
        self.client.post('/api/jobs/'+job['id']+'/cancel');await_status(self.client,job['id'],{'cancelled'})


if __name__=='__main__':
    unittest.main()
