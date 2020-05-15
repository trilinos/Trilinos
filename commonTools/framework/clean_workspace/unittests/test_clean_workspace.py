#!/usr/bin/env python3
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""Implements tests for the clean_sentinel script."""
from __future__ import print_function
import sys
sys.dont_write_bytecode = True

import os
sys.path.insert(1, os.path.dirname(os.path.dirname(__file__)))

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

import unittest
try:
    import mock
except ImportError:
    import unittest.mock as mock

from argparse import Namespace
from datetime import datetime

# pylint: disable=import-error
from clean_workspace import Cleaner

class TestRun(unittest.TestCase):

    def test_no_directory_raises(self):
        """If there is no dir attribute in the args raise SystemExit"""
        test_args = Namespace()
        setattr(test_args, 'dir', None)
        setattr(test_args, 'force_clean', False)
        with mock.patch.object(Cleaner, 'parse_args', return_value=test_args):
            cleanerInst = Cleaner()
            if sys.version_info.major is not 3:
                with self.assertRaisesRegexp(SystemExit, "No directory passed - exiting!"):
                    cleanerInst.run()
            else:
                with self.assertRaisesRegex(SystemExit, "No directory passed - exiting!"):
                    cleanerInst.run()


    def test_force_calls_clean(self):
        """If force is passed go straight to cleanup"""
        test_args = Namespace()
        setattr(test_args, 'dir', os.path.join(os.path.sep, 'dev', 'null', 'force_cleaned'))
        setattr(test_args, 'force_clean', True)
        with mock.patch.object(Cleaner, 'parse_args', return_value=test_args):
            cleanerInst = Cleaner()
            with mock.patch.object(Cleaner, 'force_clean_space') as force_clean, \
                 mock.patch('clean_workspace.print') as m_print:
                cleanerInst.run()
            force_clean.assert_called_once_with()
            m_print.assert_called_once_with("Cleaning directory /dev/null/force_cleaned "
                                            "due to command line option")

    def test_dir_calls_clean_by_date(self):
        """If force is passed go straight to cleanup"""
        test_args = Namespace()
        setattr(test_args, 'dir', os.path.join(os.path.sep, 'dev', 'null'))
        setattr(test_args, 'force_clean', False)
        with mock.patch.object(Cleaner, 'parse_args', return_value=test_args):
            cleanerInst = Cleaner()
            with mock.patch.object(Cleaner, 'clean_space_by_date') as force_clean:
                cleanerInst.run()
            force_clean.assert_called_once_with()


class TestParseArgs(unittest.TestCase):

    def test_no_args_gives_help_and_exits(self):
        """Test that the function does the right thing when given no arguments"""
        if sys.version_info.major is not 3:
            usage_message = ('usage: programName [-h] [--force-clean] dir\n'
                             'programName: error: too few arguments\n')
        else:
            usage_message = ('usage: programName [-h] [--force-clean] dir\n'
                             'programName: error: the following arguments are required: dir\n')
        with self.assertRaises(SystemExit), \
             mock.patch.object(sys, 'argv', ['programName']), \
             mock.patch('sys.stderr', new_callable=StringIO) as cleanOut:
            cleanerInst = Cleaner()
            cleanerInst.parse_args()
        self.assertEqual(usage_message, cleanOut.getvalue())

    def test_dir_only(self):
        """Passing a directory will eliminate the SystemExit and set args.dir"""
        with mock.patch.object(sys, 'argv', ['programName',
                                             os.path.join(os.path.sep, 'scratch',
                                                          'trilinos', 'workspace')]):
            cleanerInst = Cleaner()
            args = cleanerInst.parse_args()
        self.assertEqual(os.path.join(os.path.sep, 'scratch',
                                      'trilinos', 'workspace'),
                         args.dir)
        self.assertEqual(False, args.force_clean)

    def test_force(self):
        """Adding the --force option must set the action to True"""
        with mock.patch.object(sys, 'argv', ['programName',
                                             os.path.join(os.path.sep, 'scratch',
                                                          'trilinos', 'workspace'),
                                             '--force']):
            cleanerInst = Cleaner()
            args = cleanerInst.parse_args()
        self.assertEqual(os.path.join(os.path.sep, 'scratch',
                                      'trilinos', 'workspace'),
                         args.dir)
        self.assertEqual(True, args.force_clean)


class TestForceCleanSpace(unittest.TestCase):

    def test_calls_with_dir(self):
        """This function does the final cleanup  so its just os.unlink"""
        test_args = Namespace()
        setattr(test_args, 'dir', os.path.join(os.path.sep, 'dev', 'null'))
        setattr(test_args, 'force_clean', True)
        cleanerInst = Cleaner()
        cleanerInst.args = test_args
        with mock.patch('shutil.rmtree') as m_unlink, \
             mock.patch('os.path.isdir', return_value=True):
            cleanerInst.force_clean_space()
        m_unlink.assert_called_once_with(test_args.dir)


    def test_no_call_without_dir(self):
        """This function does the final cleanup  so its just os.unlink"""
        test_args = Namespace()
        setattr(test_args, 'dir', os.path.join(os.path.sep, 'dev', 'null'))
        setattr(test_args, 'force_clean', True)
        cleanerInst = Cleaner()
        cleanerInst.args = test_args
        with mock.patch('shutil.rmtree') as m_unlink, \
             mock.patch('os.path.isdir', return_value=False):
            cleanerInst.force_clean_space()
        m_unlink.assert_not_called()


class TestCleanSpaceByDate(unittest.TestCase):

    def test_defaults_do_not_clean(self):
        """The default dates should result in no clean action"""
        cleanerInst = Cleaner()
        with mock.patch.dict('os.environ',
                             {'WORKSPACE': '/scratch/Trilinos/foo/bar'}):
            with mock.patch('clean_workspace.Cleaner.force_clean_space') as force_clean, \
                 mock.patch('clean_workspace.update_last_clean_date') as update, \
                 mock.patch("clean_sentinel.open", side_effect=IOError):
                cleanerInst.clean_space_by_date()
            force_clean.assert_not_called()
            update.assert_not_called()


    def test_newer_reference_does_clean(self):
        """The default dates should result in no clean action"""
        testDate = datetime(2019, 2, 4, hour=10, minute=48)
        test_args = Namespace()
        setattr(test_args, "dir", os.path.join(os.path.sep, "dev", "null", "fake_directory"))
        with mock.patch('clean_workspace.clean_reference_date', return_value=testDate):
            cleanerInst = Cleaner()
            cleanerInst.args = test_args
            with mock.patch.dict('os.environ',
                                 {'WORKSPACE': '/scratch/Trilinos/foo/bar'}):
                with mock.patch('clean_workspace.Cleaner.force_clean_space') as force_clean, \
                     mock.patch('clean_workspace.update_last_clean_date') as update, \
                     mock.patch('clean_workspace.print') as m_print:
                    cleanerInst.clean_space_by_date()
                force_clean.assert_called_once_with()
                update.assert_called_once_with()
                m_print.assert_called_once_with("Cleaning directory /dev/null/fake_directory "
                                                "due to newer reference date")


    def test_older_date_does_clean(self):
        """The default dates should result in no clean action"""
        testDate = datetime(2019, 2, 4, hour=10, minute=0)
        test_args = Namespace()
        setattr(test_args, "dir", os.path.join(os.path.sep, "dev", "null", "will_clean"))
        with mock.patch('clean_workspace.last_clean_date', return_value=testDate):
            cleanerInst = Cleaner()
            cleanerInst.args = test_args
            with mock.patch('clean_workspace.Cleaner.force_clean_space') as force_clean, \
                 mock.patch('clean_workspace.update_last_clean_date') as update, \
                 mock.patch('clean_workspace.print') as m_print, \
                 mock.patch("clean_sentinel.open", side_effect=IOError):
                cleanerInst.clean_space_by_date()
            force_clean.assert_called_once_with()
            update.assert_called_once_with()
            m_print.assert_called_once_with("Cleaning directory /dev/null/will_clean "
                                            "due to newer reference date")


if __name__ == '__main__':
    unittest.main()
