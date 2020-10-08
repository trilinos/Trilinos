#!/usr/bin/env python
# -*- coding: utf-8; mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""
"""
from __future__ import print_function
import sys
sys.dont_write_bytecode = True

import os
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
from unittest import TestCase
try:                                                      # pragma: no cover
    import unittest.mock as mock                          # pragma: no cover
except:                                                   # pragma: no cover
    import mock                                           # pragma: no cover
from mock import Mock
from mock import MagicMock
from mock import patch

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from trilinosprhelpers.setenvironment import ModuleHelper



class mock_popen(object):
    """
    Abstract base class for popen mock
    """
    #__metaclass__ = abc.ABCMeta
    def __init__(self, cmd, stdout=None, stderr=None):
        print("mock_popen> {}".format(cmd))
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = None

    #@abc.abstractmethod
    def communicate(self):
        print("mock_popen> communicate()")
        stdout = b"stdout=1"
        stderr = b"stderr=2"
        self.returncode = 0
        return (stdout,stderr)



class mock_popen_status_ok(mock_popen):
    """
    Specialization of popen mock that will return with success.
    """
    def __init__(self, cmd, stdout=None, stderr=None):
        super(mock_popen_status_ok, self).__init__(cmd,stdout,stderr)



class mock_popen_status_error(mock_popen):
    """
    Specialization of popen mock that will return with error.
    """
    def __init__(self, cmd, stdout=None, stderr=None):
        super(mock_popen_status_error, self).__init__(cmd,stdout,stderr)

    def communicate(self):
        print("mock_popen> communicate()")
        stdout = b"stdout=1"
        stderr = b"ERROR: Something wrong happened."
        self.returncode = 1
        return (stdout,stderr)


class mock_popen_status_error_mlstatus(mock_popen):
    """
    Specialization of popen mock that will return with error (status==1)
    and stderr containing '_mlstatus = False' which happens on some systems.
    """
    def __init__(self, cmd, stdout=None, stderr=None):
        super(mock_popen_status_error_mlstatus, self).__init__(cmd,stdout,stderr)

    def communicate(self):
        print("mock_popen> communicate()")
        stdout = b"stdout=1"
        stderr = b"_mlstatus = False"
        self.returncode = 1
        return (stdout,stderr)



class moduleHelperTest(TestCase):
    """
    Main test driver for the module() function provided by the
    ModuleHelper.py file
    """
    def setUp(self):
        pass


    @patch('subprocess.Popen', side_effect=mock_popen_status_ok)
    def test_module_load_status_ok(self, arg_popen):
        r = ModuleHelper.module("load", "sems-gcc/4.8.4")
        print("result = {}".format(r))
        self.assertEqual(0, r)


    @patch('subprocess.Popen', side_effect=mock_popen_status_error)
    def test_module_load_status_error(self, arg_popen):
        r = ModuleHelper.module("load", "sems-gcc/4.8.4")
        print("result = {}".format(r))
        self.assertEqual(1, r)


    @patch('subprocess.Popen', side_effect=mock_popen_status_ok)
    def test_module_swap_status_ok(self, arg_popen):
        r = ModuleHelper.module("swap", "sems-gcc/4.8.4", "sems-gcc/7.3.0")
        print("result = {}".format(r))
        self.assertEqual(0, r)


    #@patch('subprocess.Popen', side_effect=mock_popen_status_ok)
    def test_module_unload_status_ok(self):
        with patch('subprocess.Popen', side_effect=mock_popen_status_ok):
            r = ModuleHelper.module("unload", "sems-gcc/4.8.4")
        print("result = {}".format(r))
        self.assertEqual(0, r)


    def test_module_load_args_as_list(self):
        """
        The `module()` function can take arguments in as a list.
        This tests that module works when the parameter is a list of arguments.
        """
        with patch('subprocess.Popen', side_effect=mock_popen_status_ok):
            r = ModuleHelper.module( [ "unload", "sems-gcc/4.8.4" ] )
        print("result = {}".format(r))
        self.assertEqual(0, r)


    def test_module_load_error_by_mlstatus(self):
        with patch('subprocess.Popen', side_effect=mock_popen_status_error_mlstatus):
            r = ModuleHelper.module("load", "sems-gcc/4.8.4")
        print("result = {}".format(r))
        self.assertEqual(1, r)


    def test_module_load_error_no_modulecmd(self):
        with patch('distutils.spawn.find_executable', side_effect=Exception("mock side-effect error")):
            with patch('subprocess.Popen', side_effect=mock_popen_status_error_mlstatus):
                r = ModuleHelper.module("load", "sems-gcc/4.8.4")
        print("result = {}".format(r))
        self.assertEqual(1, r)

