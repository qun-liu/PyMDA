"""Regression checks for real failure modes in diffraction discovery/statistics."""
import tempfile
import unittest
import os
import json
from types import SimpleNamespace
from unittest.mock import patch
from pathlib import Path
from helpers import discover_masters, dataset_id, frame_wedges, aimless_stats, admissible
from autoDials import Toolchain


class PipelineTests(unittest.TestCase):
    def test_repeated_master_versions_keep_distinct_selection_and_search_paths(self):
        from pymda import merge_subsets, record_id, subset_label
        records=[{'dataset_id':f'run_{i}', 'config':{'master':'repeated_master.h5' if i<2 else f'crystal_{i}_master.h5'},
                  'selected':{'statistics':{'cc_half':[.99,.99,.8],'isigi':[100-i,10,2],'completeness':[95,99,90]}}} for i in range(12)]
        args=SimpleNamespace(subset_search='auto',max_subsets=64)
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)
            (out/'comparison.json').write_text(json.dumps({'groups':[['run_1','run_5','run_7']]}))
            subsets,search=merge_subsets(records,out,args)
            self.assertIn((1,5,7),subsets)
            self.assertEqual(search['dataset_subsets'][0],[f'run_{i}' for i in range(12)])
            self.assertEqual(record_id(records[0]),'run_0')
            self.assertNotEqual(subset_label(records,(0,5)),subset_label(records,(1,5)))
            self.assertEqual(subset_label(records,(0,5)),subset_label(list(reversed(records)),(6,11)))
        legacy={'config':{'master':'legacy_master.h5'}}
        self.assertEqual(record_id(legacy),dataset_id(Path('legacy_master.h5')))

    def test_single_small_and_large_collections_choose_finite_subsets(self):
        from pymda import merge_subsets, subset_label
        records=[{'config':{'master':f'crystal_{i:03d}_master.h5'},'selected':{'statistics':{'cc_half':[.99,.99,.8],'isigi':[100-i,10,2],'completeness':[95,99,90]}}} for i in range(80)]
        args=SimpleNamespace(subset_search='auto',max_subsets=64)
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)
            single,search=merge_subsets(records[:1],out,args)
            self.assertEqual(single,[(0,)])
            self.assertEqual(search['mode'],'single')
            small,search=merge_subsets(records[:3],out,args)
            self.assertEqual(small,[(0,1),(0,2),(1,2),(0,1,2)])
            self.assertEqual(search['mode'],'exhaustive')
            group=[dataset_id(Path(r['config']['master'])) for r in records[20:25]]
            (out/'comparison.json').write_text(json.dumps({'groups':[group]}))
            subsets,search=merge_subsets(records,out,args)
            self.assertEqual(search['mode'],'bounded')
            self.assertLessEqual(len(subsets),64)
            self.assertEqual(len(set(subsets)),len(subsets))
            self.assertEqual(subsets[0],tuple(range(80)))
            self.assertIn(tuple(range(20,25)),subsets)
            self.assertIn(tuple(range(79)),subsets)
            self.assertIn((0,1),subsets)
            self.assertEqual(search,json.loads((out/'search.json').read_text()))
            self.assertLess(len(subset_label(records,range(80))),80)
            self.assertEqual(subset_label(records,range(80)),subset_label(list(reversed(records)),range(80)))
            args.max_subsets=1
            self.assertEqual(merge_subsets(records,out,args)[0],[tuple(range(80))])
            args.subset_search='exhaustive'
            with self.assertRaisesRegex(ValueError,'choose auto or bounded'):
                merge_subsets(records,out,args)

    def test_completed_pilot_can_resume_exports_without_reading_images(self):
        from pymda import process
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp); master=out/'pilot_master.h5'
            directory=out/'crystals'/dataset_id(master); directory.mkdir(parents=True)
            (directory/'imported.expt').write_text(json.dumps({'scan':[{'image_range':[1,1800]}]}))
            (directory/'integrated.expt').write_text(json.dumps({'scan':[{'image_range':[2,180]}]}))
            (directory/'integrated.refl').write_bytes(b'completed pilot')
            (directory/'dials.import.command.json').write_text(json.dumps({'argv':[str(master)]}))
            (directory/'dials.integrate.command.json').write_text(json.dumps({'argv':['prediction.d_min=1.8','scan_range=1,180']}))
            (directory/'dials.integrate.console.log').write_text('Saving the experiments to integrated.expt')
            config={'master':str(master),'frames':180,'wedges':2,'resolution':1.8,'space_group':None,'unit_cell':None,'fix_detector':False}
            (directory/'processing_config.json').write_text(json.dumps(config))
            class Tools:
                def run(self,name,args,directory,stdin=None,expected=()):
                    if name not in {'dials.slice_sequence','dials.export'}:
                        raise AssertionError('Resume reread images: '+name)
                    directory.mkdir(parents=True,exist_ok=True)
                    for filename in expected:
                        (directory/filename).write_bytes(b'export')
            stats={'resolution':[2.5,10,2.5],'isigi':[10,20,2],'cc_half':[.99,.99,.7],'completeness':[95,99,90]}
            args=SimpleNamespace(frames=180,wedges=2,split_wedges=3,reso=1.8,reso_cchalf=2.5,spg=None,unit_cell=None,fix_detector=False,exports_only=True,thread=1)
            with patch('pymda.scale',return_value=stats):
                record=process(Tools(),master,out,args)
            self.assertEqual([(w['first_frame'],w['last_frame']) for w in record['scan_wedges']],[(1,60),(61,120),(121,180)])
            self.assertEqual(record['config']['frames'],180)
            args.frames=90
            with self.assertRaisesRegex(ValueError,'Processing settings changed'):
                process(Tools(),master,out,args)

    def test_large_collection_can_finalize_with_bounded_search(self):
        from pymda import assemble
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)
            stats={'resolution':[2.5,10,2.5],'isigi':[10,20,2],'cc_half':[.99,.99,.7],'completeness':[95,99,90]}
            records=[]
            for index in range(12):
                mtz=out/f'input_{index}.mtz'; mtz.write_bytes(b'intensities')
                prefix={'mtz':str(mtz),'last_frame':1800,'statistics':stats}
                records.append({'config':{'master':f'crystal_{index}_master.h5'},'selected':prefix,'prefixes':[prefix]})
            calls=[]
            def scaling(tools,paths,directory,resolution):
                calls.append(len(paths)); directory.mkdir(parents=True,exist_ok=True)
                for name in ('scaled.mtz','sorted.mtz','aimless.console.log','pointless.console.log','statistics.json','inputs.json'):
                    (directory/name).write_text('data')
                return stats
            class Tools:
                def run(self,name,args,directory,stdin=None,expected=()):
                    for filename in expected:
                        (directory/filename).write_text('generated')
            args=SimpleNamespace(subset_search='auto',max_subsets=1,prefix_policy='all',resolutions=[2.5],reso=1.8,thread=1,merge_workers=1,min_completeness=90,min_outer_completeness=80,min_cc_half=.3,min_isigi=1,nres=None,crystal_identity='unknown')
            with patch('pymda.scale',side_effect=scaling),patch('builtins.print'):
                assemble(Tools(),records,out,args)
            self.assertEqual(calls,[12])
            self.assertTrue((out/'best/final.mtz').is_file())
            best=json.loads((out/'best.json').read_text())
            self.assertEqual(len(best['crystals']),12)
            self.assertEqual(best['subset_search']['mode'],'bounded')
            self.assertIn('bounded',(out/'processing_report.md').read_text())

    def test_consecutive_wedges_cover_each_frame_once(self):
        self.assertEqual(frame_wedges(1,1800,3), [(1,600),(601,1200),(1201,1800)])
        self.assertEqual(frame_wedges(10,20,3), [(10,13),(14,17),(18,20)])
        self.assertEqual(frame_wedges(1,1800,1), [(1,1800)])
        self.assertEqual(frame_wedges(5,7,3), [(5,5),(6,6),(7,7)])
        for first,last,count in [(1,2,3),(1,1800,0),(3,1,1)]:
            with self.assertRaises(ValueError):
                frame_wedges(first,last,count)

    def test_cached_integration_exports_disjoint_wedges_with_original_bounds(self):
        from pymda import process
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp); master=out/'single_master.h5'
            directory=out/'crystals'/dataset_id(master); directory.mkdir(parents=True)
            (directory/'imported.expt').write_text(json.dumps({'scan':[{'image_range':[1,1800]}]}))
            (directory/'integrated.expt').write_text(json.dumps({'scan':[{'image_range':[2,1799]}]}))
            (directory/'integrated.refl').write_bytes(b'reflections')
            mtz=directory/'full.mtz'; mtz.write_bytes(b'full scan')
            stats={'resolution':[2.5,10,2.5],'isigi':[10,20,2],'cc_half':[.99,.99,.7],'completeness':[95,99,90]}
            prefix={'mtz':str(mtz),'last_frame':1799,'statistics':stats}
            config={'master':str(master),'frames':None,'wedges':4,'resolution':1.8,'space_group':None,'unit_cell':None,'fix_detector':False}
            (directory/'manifest.json').write_text(json.dumps({'config':config,'prefixes':[prefix],'selected':prefix}))
            calls=[]
            class Tools:
                def run(self,name,args,directory,stdin=None,expected=()):
                    calls.append((name,args))
                    directory.mkdir(parents=True,exist_ok=True)
                    for filename in expected:
                        (directory/filename).write_bytes(b'wedge output')
            args=SimpleNamespace(frames=None,wedges=4,split_wedges=3,reso=1.8,reso_cchalf=2.5,spg=None,unit_cell=None,fix_detector=False)
            with patch('pymda.scale',return_value=stats):
                record=process(Tools(),master,out,args)
            self.assertEqual([(w['first_frame'],w['last_frame']) for w in record['scan_wedges']],[(1,600),(601,1200),(1201,1800)])
            self.assertEqual([next(a for a in args if str(a).startswith('image_range=')) for name,args in calls if name=='dials.slice_sequence'],['image_range=2,600','image_range=601,1200','image_range=1201,1799'])
            self.assertEqual(record['split_wedges'],3)
            self.assertTrue(all(w['statistics'] and Path(w['mtz']).is_file() for w in record['scan_wedges']))
            self.assertEqual(record['prefixes'],[prefix])
            self.assertEqual({name for name,_ in calls},{'dials.slice_sequence','dials.export'})

    def test_failed_wedge_assessment_is_not_reported_as_success(self):
        from pymda import record_scan_wedges
        with tempfile.TemporaryDirectory() as temp:
            directory=Path(temp)
            for name in ('imported.expt','integrated.expt'):
                (directory/name).write_text(json.dumps({'scan':[{'image_range':[1,6]}]}))
            class Tools:
                def run(self,name,args,directory,stdin=None,expected=()):
                    directory.mkdir(parents=True,exist_ok=True)
                    for filename in expected:
                        (directory/filename).write_bytes(b'export')
            args=SimpleNamespace(frames=None,split_wedges=3,reso_cchalf=2.5)
            with patch('pymda.scale',side_effect=ValueError('No scaling summary')):
                with self.assertRaisesRegex(RuntimeError,'No consecutive wedge processed'):
                    record_scan_wedges(Tools(),directory,{},args)
            record=json.loads((directory/'manifest.json').read_text())
            self.assertTrue(all(w['mtz'] is None and w['statistics'] is None and w['error'] for w in record['scan_wedges']))

    def test_discovery_handles_nested_data_and_duplicate_broken_links(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            folder = root/'deep'/'rotation'
            folder.mkdir(parents=True)
            master = folder/'crystal_master.h5'
            master.touch()
            (root/'duplicate_master.h5').symlink_to(master)
            (root/'missing_master.h5').symlink_to(root/'absent')
            (root/'Raster_master.h5').touch()
            self.assertEqual(discover_masters(root), [master])

    def test_summary_uses_merged_isigi_and_exact_completeness_label(self):
        with tempfile.TemporaryDirectory() as temp:
            log = Path(temp)/'aimless.log'
            log.write_text('High resolution limit 2.10 8.90 2.10\n'
                           'Mean((I)/sd(I)) 15.2 25.4 1.2\n'
                           'Mn(I) half-set correlation CC(1/2) 0.99 0.99 0.45\n'
                           'Completeness 98.5 99.1 94.2\n'
                           'Anomalous completeness 85.0 89.0 70.0\n')
            stats=aimless_stats(log)
            self.assertEqual(stats['completeness'], [98.5,99.1,94.2])
            self.assertTrue(admissible(stats))
            stats['cc_half'][2]=0.2
            self.assertFalse(admissible(stats))

    def test_one_crystal_can_produce_final_amplitudes(self):
        from pymda import assemble
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            mtz = out/'integrated.mtz'; mtz.write_bytes(b'intensities')
            stats = {'resolution':[2.5,10,2.5], 'isigi':[10,20,2], 'cc_half':[.99,.99,.7], 'completeness':[95,99,90]}
            prefix = {'mtz':str(mtz), 'last_frame':180,'statistics':stats}
            record = {'config':{'master':'one_master.h5'},'selected':prefix,'prefixes':[prefix]}
            def scaling(tools, paths, directory, resolution):
                self.assertEqual(len(paths),1)
                directory.mkdir(parents=True,exist_ok=True)
                for name in ('scaled.mtz','sorted.mtz','aimless.console.log','pointless.console.log','statistics.json','inputs.json'):
                    (directory/name).write_text('data')
                return stats
            class Tools:
                def run(self, name, args, directory, stdin=None, expected=()):
                    for target in expected:
                        (directory/target).write_text('generated')
            args = SimpleNamespace(prefix_policy='all',resolutions=[2.5],reso=1.8,thread=1,merge_workers=1,min_completeness=90,min_outer_completeness=80,min_cc_half=.3,min_isigi=1,nres=None,crystal_identity='unknown')
            with patch('pymda.scale',side_effect=scaling),patch('builtins.print'):
                assemble(Tools(),[record],out,args)
            self.assertTrue((out/'best/final.mtz').is_file())
            self.assertEqual(len(json.loads((out/'best.json').read_text())['crystals']),1)
            self.assertIn('one rotation dataset',(out/'processing_report.md').read_text())

    def test_incomplete_summary_is_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            log=Path(temp)/'aimless.log'
            log.write_text('AIMLESS failed before summary\n')
            with self.assertRaises(ValueError):
                aimless_stats(log)

    def test_failed_tool_cannot_leave_an_old_mtz_as_success(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'bin').mkdir()
            executable = root/'bin'/'failing_tool'
            executable.write_text('#!/bin/sh\necho integration_failed >&2\nexit 1\n')
            executable.chmod(0o755)
            work = root/'work'
            work.mkdir()
            (work/'integrated.mtz').write_text('stale reflections')
            tools = Toolchain.__new__(Toolchain)
            tools.root = root
            tools.env = dict(os.environ)
            with self.assertRaises(RuntimeError):
                tools.run('failing_tool', [], work, expected=['integrated.mtz'])
            self.assertFalse((work/'integrated.mtz').exists())
            self.assertIn('integration_failed', (work/'failing_tool.console.log').read_text())


if __name__ == '__main__':
    unittest.main()
