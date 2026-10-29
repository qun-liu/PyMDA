"""Web/queue regression tests, including real worker subprocess lifecycle tests."""
import json
import os
from pathlib import Path
import signal
import tempfile
import time
import unittest

from fastapi.testclient import TestClient
from main import create_app
from web.jobs import JobSpec, JobStore, Settings, command_for, validate_spec

BASE_URL = 'http://127.0.0.1:8000'


def inspector(paths, settings):
    return [{'path':str(p),'name':p.name,'ready':True,'frames':180,'companions':[],'error':''} for p in paths]


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
        self.settings = Settings(self.data,self.output,root/'jobs',testing=True)
        self.client_context = TestClient(create_app(self.settings,start_worker=False,inspector=inspector),base_url=BASE_URL)
        self.client = self.client_context.__enter__()

    def tearDown(self):
        self.client_context.__exit__(None,None,None)
        self.temp.cleanup()

    def payload(self, **extra):
        return {'mode':'process','input_dir':str(self.data),'master_files':[str(self.master)],**extra}

    def test_python_html_and_no_script_or_templates(self):
        response = self.client.get('/')
        self.assertEqual(response.status_code,200)
        self.assertIn('Process diffraction data',response.text)
        self.assertIn('sample_master.h5',response.text)
        self.assertIn('Assemble processed datasets',response.text)
        self.assertNotIn('<script',response.text.lower())
        self.assertEqual(self.client.get('/docs').status_code,404)

    def test_validation_and_unknown_fields(self):
        for extra in [{'thread':20},{'reso':float('inf')},{'resolutions':[-2]}, {'resolutions':[2.261]}, {'unit_cell':'1,2,3,10,10,160'}, {'command':'touch injected'}]:
            # JSON cannot encode inf; send a string, parsed/validated as a float.
            if extra.get('reso') == float('inf'): extra['reso']='inf'
            self.assertEqual(self.client.post('/api/jobs',json=self.payload(**extra)).status_code,422)
        self.assertEqual(self.client.get('/api/jobs').json(),[])

    def test_reject_unselected_or_external_files(self):
        for extra in [{'master_files':[]},{'input_dir':'/tmp','master_files':[str(self.master)]},{'master_files':['/etc/passwd']}]:
            self.assertEqual(self.client.post('/api/jobs',json=self.payload(**extra)).status_code,422)

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
        self.settings=Settings(self.data,root/'processed',root/'jobs',pipeline=pipeline,testing=True)
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
