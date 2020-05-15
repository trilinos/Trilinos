#!/usr/bin/env python3
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
# Change above to '/usr/bin/python -3' for python 3.x porting warnings
"""
This contains abstractions to handle the reference date for the last required
full clean.
"""
import sys
sys.dont_write_bytecode = True

import os

from datetime import datetime
import pickle


def clean_reference_date():
    """Returns a datetime object of the last requested clean state"""
    defaultReturnDate = datetime(2019, 2, 4, hour=10, minute=45,
                                 second=0, microsecond=0, tzinfo=None)
    date_file_name = os.path.join(os.path.sep, 'ascldap', 'users',
                                  'trilinos', '.cleanAllWorkspacesDatefile')
    return read_date(date_file_name, defaultReturnDate)


def set_reference_date():
    """Sets the requested clean state to the current time"""
    newReferenceDate = datetime.now()
    date_file_name = os.path.join(os.path.sep, 'ascldap', 'users',
                                  'trilinos', '.cleanAllWorkspacesDatefile')
    with open(date_file_name, "w") as pickleFile:
        pickle.dump(newReferenceDate, pickleFile)


def last_clean_date():
    """Returns a datetime object of the last requested clean state"""
    defaultReturnDate = datetime(2019, 2, 4, hour=10, minute=46,
                                 second=0, microsecond=0, tzinfo=None)
    date_file_name = os.path.join(os.environ['WORKSPACE'], 'lastCleanDate')
    return read_date(date_file_name, defaultReturnDate)


def update_last_clean_date():
    """stores a datetime object of the current date"""
    currentDate = datetime.now()
    date_file_name = os.path.join(os.environ['WORKSPACE'], 'lastCleanDate')
    with open(date_file_name, 'w') as sFile:
        pickle.dump(currentDate, sFile)
    return


def read_date(date_file, defaultReturnDate):
    returnValue = defaultReturnDate
    try:
        with open(date_file) as pickleFile:
            try:
                returnValue = pickle.load(pickleFile)
            except pickle.UnpicklingError:
                pass
    except IOError:
        pass
    return returnValue
