"""Regression checks for real failure modes in diffraction discovery/statistics."""
import tempfile
import unittest
import os
from pathlib import Path
from helpers import discover_masters, aimless_stats, admissible
from autoDials import Toolchain


class PipelineTests(unittest.TestCase):
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
