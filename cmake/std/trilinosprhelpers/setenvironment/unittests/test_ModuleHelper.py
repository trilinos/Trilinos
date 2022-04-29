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
    from unittest.mock import Mock
    from unittest.mock import MagicMock
    from unittest.mock import patch
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
    def __init__(self, cmd, stdout=None, stderr=None):
        print("mock_popen> {}".format(cmd))
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = None

    def communicate(self):
        print("mock_popen> communicate()")
        stdout = b"os.environ['__foobar__'] ='baz'\ndel os.environ['__foobar__']"
        stderr = b"stderr=1"
        self.returncode = 0
        return (stdout,stderr)



class mock_popen_status_ok(mock_popen):
    """
    Specialization of popen mock that will return with success.
    """
    def __init__(self, cmd, stdout=None, stderr=None):
        super(mock_popen_status_ok, self).__init__(cmd,stdout,stderr)



class mock_popen_status_error_return_nonetype(mock_popen):
    """
    Specialization of popen mock that will return with error.
    """
    def __init__(self, cmd, stdout=None, stderr=None):
        super(mock_popen_status_error_return_nonetype, self).__init__(cmd,stdout,stderr)

    def communicate(self):
        print("mock_popen> communicate()")
        stdout = b""
        stderr = b""
        self.returncode = None
        return (stdout,stderr)



class mock_popen_status_error_rc1(mock_popen):
    """
    Specialization of popen mock that will return with error.

    Test the condition where modulecmd returned a status of 1 and
    has `ERROR:` in its stderr field.
    """
    def __init__(self, cmd, stdout=None, stderr=None):
        super(mock_popen_status_error_rc1, self).__init__(cmd,stdout,stderr)

    def communicate(self):
        print("mock_popen> communicate()")
        stdout = b"stdout=1"
        stderr = b"ERROR: Something wrong happened."
        self.returncode = 1
        return (stdout,stderr)



class mock_popen_status_error_rc0(mock_popen):
    """
    Specialization of popen mock.

    Simulates the results from a modulecmd operation that had
    an error loading a module (maybe not found). Modulecmd will tend
    to have a message like "ERROR: could not load module" in its stderr
    field but it will generally return an exit status of 0.
    """
    def __init__(self, cmd, stdout=None, stderr=None):
        super(mock_popen_status_error_rc0, self).__init__(cmd,stdout,stderr)

    def communicate(self):
        print("mock_popen> communicate()")
        stdout = b"os.environ['__foobar__'] ='baz'\ndel os.environ['__foobar__']"
        stderr = b"ERROR: Something wrong happened."
        self.returncode = 0
        return (stdout,stderr)



class mock_popen_status_mlstatus_success(mock_popen):
    """
    Specialization of popen mock that will return with error (status==1)
    and stderr containing '_mlstatus = True' which happens on some systems.
    """
    def __init__(self, cmd, stdout=None, stderr=None):
        super(mock_popen_status_mlstatus_success, self).__init__(cmd,stdout,stderr)

    def communicate(self):
        print("mock_popen> communicate()")
        stdout = b"_mlstatus = True"
        stderr = b""
        self.returncode = 0
        return (stdout,stderr)



class mock_popen_status_mlstatus_error(mock_popen):
    """
    Specialization of popen mock that will return with error (status==1)
    and stderr containing '_mlstatus = False' which happens on some systems.
    """
    def __init__(self, cmd, stdout=None, stderr=None):
        super(mock_popen_status_mlstatus_error, self).__init__(cmd,stdout,stderr)

    def communicate(self):
        print("mock_popen> communicate()")
        stdout = b"_mlstatus = False"
        stderr = b""
        self.returncode = 0
        return (stdout,stderr)



class mock_popen_stdout_throws_on_exec(mock_popen):
    """
    This is a special mock function to replace the Popen instruction
    in ModuleHelper with one that returns a python command in its stdout
    which will raise an exception when it is run through `exec()`
    """
    def __init__(self, cmd, stdout=None, stderr=None):
        super(mock_popen_stdout_throws_on_exec, self).__init__(cmd,stdout,stderr)

    def communicate(self):
        print("mock_popen> communicate()")
        stdout = b"raise Exception('Exception generated by mock popen')"
        stderr = b""
        self.returncode = 0
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
        r = ModuleHelper.module("load", "sems-archive-gcc/4.8.4")
        print("result = {}".format(r))
        self.assertEqual(0, r)


    @patch('subprocess.Popen', side_effect=mock_popen_status_error_rc1)
    def test_module_load_status_error(self, arg_popen):
        r = ModuleHelper.module("load", "sems-archive-gcc/4.8.4")
        print("result = {}".format(r))
        self.assertEqual(1, r)


    @patch('subprocess.Popen', side_effect=mock_popen_status_error_rc0)
    def test_module_load_status_error(self, arg_popen):
        r = ModuleHelper.module("load", "sems-archive-gcc/4.8.4")
        print("result = {}".format(r))
        self.assertEqual(1, r)


    @patch('subprocess.Popen', side_effect=mock_popen_status_ok)
    def test_module_swap_status_ok(self, arg_popen):
        r = ModuleHelper.module("swap", "sems-archive-gcc/4.8.4", "sems-archive-gcc/7.3.0")
        print("result = {}".format(r))
        self.assertEqual(0, r)


    #@patch('subprocess.Popen', side_effect=mock_popen_status_ok)
    def test_module_unload_status_ok(self):
        with patch('subprocess.Popen', side_effect=mock_popen_status_ok):
            r = ModuleHelper.module("unload", "sems-archive-gcc/4.8.4")
        print("result = {}".format(r))
        self.assertEqual(0, r)


    def test_module_load_args_as_list(self):
        """
        The `module()` function can take arguments in as a list.
        This tests that module works when the parameter is a list of arguments.
        """
        with patch('subprocess.Popen', side_effect=mock_popen_status_ok):
            r = ModuleHelper.module( [ "unload", "sems-archive-gcc/4.8.4" ] )
        print("result = {}".format(r))
        self.assertEqual(0, r)


    def test_module_load_success_by_mlstatus(self):
        with patch('subprocess.Popen', side_effect=mock_popen_status_mlstatus_success):
            r = ModuleHelper.module("load", "sems-archive-gcc/4.8.4")
        print("result = {}".format(r))
        self.assertEqual(0, r)


    def test_module_load_error_by_mlstatus(self):
        with patch('subprocess.Popen', side_effect=mock_popen_status_mlstatus_error):
            r = ModuleHelper.module("load", "sems-archive-gcc/4.8.4")
        print("result = {}".format(r))
        self.assertEqual(1, r)


    def test_module_load_error_no_modulecmd(self):
        with patch('distutils.spawn.find_executable', side_effect=Exception("mock side-effect error")):
            with patch('subprocess.Popen', side_effect=mock_popen_status_mlstatus_error):
                r = ModuleHelper.module("load", "sems-archive-gcc/4.8.4")
        print("result = {}".format(r))
        self.assertEqual(1, r)


    def test_module_load_error_exception(self):
        """
        This tests the execution path in ModuleHelper in which the exec() command would
        throw an exception.
        In the Module call, we execute `modulehelper` so that it generates a list of
        Python commands. If one of them throws an exception we catch it, print a log message,
        and re-raise a BaseException.
        """
        with patch('subprocess.Popen', side_effect=mock_popen_stdout_throws_on_exec):
            with self.assertRaises(BaseException):
                ModuleHelper.module("load", "gcc/4.8.4")


    def test_module_load_error_module_returns_nonetype(self):
        """
        This tests a failure when the `module()` function returns a NoneType
        object (i.e., like what happens with LMOD)
        """
        with patch('subprocess.Popen', side_effect=mock_popen_status_error_return_nonetype):
            with self.assertRaises(TypeError):
                ModuleHelper.module("load", "gcc/4.8.4")
