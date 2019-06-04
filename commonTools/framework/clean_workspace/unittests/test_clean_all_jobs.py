#!/usr/bin/env python
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

import pickle

# pylint: disable=import-error
from clean_all_jobs import CleanReference

class TestRun(unittest.TestCase):

    def test_clean_run(self):
        """A clean run should call set_reference_date"""
        cleanRefInst = CleanReference()
        with mock.patch('clean_all_jobs.set_reference_date') as refDate:
            cleanRefInst.run()
        refDate.assert_called_once_with()

    def test_pickling_exception(self):
        """A Pickling exception should give an error and exit"""
        cleanRefInst = CleanReference()
        with mock.patch('clean_all_jobs.set_reference_date',
                        side_effect=pickle.PicklingError()), \
             mock.patch('sys.stderr', new_callable=StringIO):
            with self.assertRaises(SystemExit):
                cleanRefInst.run()

    def test_io_exception(self):
        """A Pickling exception should give an error and exit"""
        cleanRefInst = CleanReference()
        with mock.patch('clean_all_jobs.set_reference_date',
                        side_effect=IOError()), \
             mock.patch('sys.stderr', new_callable=StringIO):
            with self.assertRaises(SystemExit):
                cleanRefInst.run()

if __name__ == '__main__':
    unittest.main()
