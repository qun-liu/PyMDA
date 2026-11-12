"""CCP4-runtime checks. Set PYMDA_TEST_MTZ to a measured exported intensity MTZ."""
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

from compare_crystals import compare, intensity


@unittest.skipUnless(os.getenv('PYMDA_TEST_MTZ'), 'Set PYMDA_TEST_MTZ to an exported MTZ')
class ComparisonTests(unittest.TestCase):
    def test_common_intensities_and_complete_link_groups(self):
        original = Path(os.environ['PYMDA_TEST_MTZ'])
        x = intensity(original, 2.5)
        self.assertAlmostEqual(x.correlation(x).coefficient(),1)
        with tempfile.TemporaryDirectory(prefix='pymda comparison with spaces ') as temp:
            root = Path(temp)
            for name in ('A','B','C'):
                directory = root/'crystals'/name
                directory.mkdir(parents=True)
                target = directory/'integrated.mtz'
                shutil.copy2(original,target)
                (directory/'manifest.json').write_text(json.dumps({'config':{'master':name+'_master.h5'},'selected':{'mtz':str(target)}}))
            r = compare(root,os.getenv('PYMDA_CCP4','/home/aa1/ccp47'),2.5,2,2,.99,100)
            self.assertEqual(r['groups'],[['A','B','C']])
            self.assertTrue(all(p['correlation']>.99 and p['common_reflections']>=100 for p in r['pairs']))
            r = compare(root,os.getenv('PYMDA_CCP4','/home/aa1/ccp47'),2.5,2,2,.99,10000000)
            self.assertEqual(r['groups'],[['A'],['B'],['C']])
            self.assertTrue(all(p['correlation'] is None and 'Insufficient' in p['error'] for p in r['pairs']))


if __name__ == '__main__':
    unittest.main()
