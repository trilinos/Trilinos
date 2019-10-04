#!/usr/bin/env python
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""Implements tests for the clean_sentinel script."""
# pylint: disable=wrong-import-position
from __future__ import print_function
import sys
sys.dont_write_bytecode = True

import os
sys.path.insert(1, os.path.dirname(os.path.dirname(__file__)))


import unittest
try:
    import mock
except ImportError:
    import unittest.mock as mock

from datetime import datetime
import pickle

# pylint: disable=import-error
from clean_sentinel import clean_reference_date
# pylint: disable=import-error
from clean_sentinel import set_reference_date
# pylint: disable=import-error
from clean_sentinel import last_clean_date
# pylint: disable=import-error
from clean_sentinel import update_last_clean_date


class TestCleanReferenceDate(unittest.TestCase):
    """Function to return the date of the last requested clean"""

    def test_default_date(self):
        """Test that the function returns the expected date"""
        testDate = datetime(2019, 2, 4, hour=10, minute=45,
                            second=0, microsecond=0, tzinfo=None)
        with mock.patch("clean_sentinel.open", side_effect=IOError):
            self.assertEqual(testDate, clean_reference_date())

    def test_stored_date(self):
        """If a different date has been stored make sure it gets retrieved"""
        testDate = datetime(2021, 3, 31, hour=3, minute=4,
                            second=32)
        with mock.patch('clean_sentinel.pickle.load',
                        return_value=testDate), \
             mock.patch('clean_sentinel.open'):
                self.assertEqual(testDate, clean_reference_date())

    def test_missing_reference_file(self):
        """No file present should result in the default result"""
        testDate = datetime(2021, 3, 31, hour=3, minute=4,
                            second=32)
        defaultDate = datetime(2019, 2, 4, hour=10, minute=45,
                            second=0, microsecond=0, tzinfo=None)

        with mock.patch('clean_sentinel.pickle.load',
                        return_value=testDate), \
             mock.patch('clean_sentinel.open',
                        side_effect=IOError()):
            self.assertEqual(defaultDate, clean_reference_date())


class TestSetReferenceDate(unittest.TestCase):
    """Function to set the date of the last requested clean"""

    def test_store_date(self):
        """If store the current date"""
        testDate = datetime.now()
        with mock.patch('clean_sentinel.pickle.dump') as pDump, \
             mock.patch('clean_sentinel.open'):
            set_reference_date()
        self.assertTrue(testDate < pDump.call_args[0][0])


class TestLastCleanDate(unittest.TestCase):
    """Function to return the date of the last completed clean"""

    def test_default_date(self):
        """Test that the default previous clean date is built
           - currently set a year back"""
        testDate = datetime(2019, 2, 4, hour=10, minute=46,
                            second=0, microsecond=0, tzinfo=None)
        with mock.patch.dict('os.environ',
                             {'WORKSPACE': '/scratch/Trilinos/foo/bar'}):
            self.assertEqual(testDate, last_clean_date())

    def test_stored_date(self):
        """Test that a stored date will be properly retrived"""
        testDate = datetime(2016, 12, 4, hour=1, minute=5,
                            second=0, microsecond=0)

        with mock.patch('clean_sentinel.pickle.load',
                        return_value=testDate), \
             mock.patch('clean_sentinel.open'), \
             mock.patch.dict('os.environ',
                             {'WORKSPACE': '/scratch/trilinos/foo/bar'}):
            self.assertEqual(testDate, last_clean_date())

    def test_missing_reference_file(self):
        """No file present should result in the default result"""
        testDate = datetime(2021, 3, 29, hour=3, minute=4,
                            second=32)
        defaultDate = datetime(2019, 2, 4, hour=10, minute=45,
                               second=0, microsecond=0, tzinfo=None)

        with mock.patch('clean_sentinel.pickle.load',
                        return_value=testDate), \
             mock.patch('clean_sentinel.open',
                        side_effect=IOError()):
            self.assertEqual(defaultDate, clean_reference_date())


class TestUpdateLastCleanDate(unittest.TestCase):
    """Function to update the date of the last completed clean"""

    def test_update_last_clean_date_writes(self):
        """Test that the current date is stored"""
        with mock.patch('clean_sentinel.pickle.dump') as p_save, \
             mock.patch('clean_sentinel.datetime') as dt, \
             mock.patch('clean_sentinel.open') as pFile, \
             mock.patch.dict('os.environ',
                             {'WORKSPACE': '/scratch/trilinos/foo/bar'}):
            self.assertEqual(None, update_last_clean_date())
        p_save.assert_called_once_with(dt.now(),
                                       pFile.return_value.__enter__())
        pFile.assert_called_once_with('/scratch/trilinos/foo/bar/lastCleanDate', 'w')

    def test_update_last_clean_date_raises_with_no_file(self):
        """Test that the current date is stored"""
        with mock.patch.dict('os.environ',
                             {'WORKSPACE': '/scratch/trilinos/foo/bar'}):
            with self.assertRaises(IOError):
                update_last_clean_date()


if __name__ == '__main__':
    unittest.main()
