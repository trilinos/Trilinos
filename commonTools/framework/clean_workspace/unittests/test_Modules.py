#!/usr/bin/env python
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""Implements tests for Modules.py."""
# pylint: disable=wrong-import-position
# pylint: disable=invalid-name
from __future__ import print_function
import sys
sys.dont_write_bytecode = True

import os
import unittest
import mock
from cStringIO import StringIO

sys.path.insert(1, os.path.join(os.path.dirname(__file__), os.pardir))
import Modules


class TestModules(unittest.TestCase):
    """Test Modules script"""

    def setUp(self):
        modules_home = {'MODULESHOME': '/dummy/path/1'}
        which_side_effects = ['/path/to/modulecmd', None, None]
        find_side_effects = [None, '/fake/path/modules/init/python.py']
        with mock.patch.dict(os.environ, modules_home), \
                mock.patch('Modules.which',
                           side_effect=which_side_effects), \
                mock.patch('Modules.find_first_binary',
                           return_value='/fake/path/modulecmd'), \
                mock.patch('Modules.find_file_in_list',
                           side_effect=find_side_effects), \
                mock.patch('__builtin__.execfile'):
            self.module_obj = Modules.Module()

    def test_module_setup(self):
        """Test ability to instantiate the class"""
        modules_home = {'MODULESHOME': '/dummy/path/1'}
        which_side_effects = ['/path/to/modulecmd', None, None]
        find_side_effects = [None, '/fake/path/modules/init/python.py']
        with mock.patch.dict(os.environ, modules_home), \
                mock.patch('Modules.which',
                           side_effect=which_side_effects), \
                mock.patch('Modules.find_first_binary',
                           return_value='/fake/path/modulecmd'), \
                mock.patch('Modules.find_file_in_list',
                           side_effect=find_side_effects), \
                mock.patch('__builtin__.execfile') as mock_exec:
            result = Modules.Module()
        mock_exec.assert_called_once_with(find_side_effects[-1])
        self.assertEqual('/fake/path/modulecmd', result.command)
        self.assertEqual('modulecmd', result.command_name)
        self.assertEqual('/fake/path/modules/init/python.py', result.init_file)

    def test_module_setup_lmod(self):
        """Test ability to instantiate the class if on a system using lmod"""
        modules_env = {'MODULESHOME': '/dummy/path/1', 'LMOD_CMD': 'lmodcmd'}
        which_side_effects = ['/path/to/lmodcmd', None, None]
        find_side_effects = [None, '/fake/path/modules/init/python.py']
        with mock.patch.dict(os.environ, modules_env), \
                mock.patch('Modules.which',
                           side_effect=which_side_effects), \
                mock.patch('Modules.find_first_binary',
                           return_value='/fake/path/lmodcmd'), \
                mock.patch('Modules.find_file_in_list',
                           side_effect=find_side_effects), \
                mock.patch('__builtin__.execfile') as mock_exec:
            result = Modules.Module()
        mock_exec.assert_not_called()
        self.assertEqual('lmodcmd', result.command)
        self.assertEqual('lmodcmd', result.command_name)
        self.assertEqual(None, result.init_file)

    def test_module_list(self):
        """Test the module_list method"""
        loaded_modules = '''Currently Loaded Modulefiles:
sierra-python/2.7
sierra-git/2.6.1
'''
        expected = [loaded_modules.splitlines()[0],
                    loaded_modules.splitlines()[2],
                    loaded_modules.splitlines()[1]]

        with mock.patch('sys.stderr', new_callable=StringIO) as err_output, \
             mock.patch('sys.stdout', new_callable=StringIO) as std_output:
            with mock.patch('Modules.Module.module',
                            return_value=loaded_modules):
                self.module_obj.module_list()
            stderr = err_output.getvalue()
            stdout = std_output.getvalue()
        self.assertEqual(expected, stderr.splitlines())
        self.assertEqual('', stdout)

    def test_module(self):
        """Test module method"""
        with mock.patch('sys.stderr', new_callable=StringIO) as err_output, \
             mock.patch('sys.stdout', new_callable=StringIO) as std_output:
            with mock.patch.dict('os.environ', {}), \
                    mock.patch('os.path.exists', return_value=True), \
                    mock.patch('subprocess.Popen') as mock_subproc_popen:
                process_mock = mock.Mock()
                attrs = {'communicate.return_value': ('print("exec_cmd")',
                                                      'error output'),
                         'wait.return_value': 0}
                process_mock.configure_mock(**attrs)
                mock_subproc_popen.return_value = process_mock
                result = self.module_obj.module('load', 'sierra')
            stdout = std_output.getvalue()
            stderr = err_output.getvalue()
        self.assertEqual('error output', result)
        self.assertEqual('exec_cmd\n', stdout)
        self.assertEqual('error output\n', stderr)

    def test_module_no_command(self):
        """Test module method when unable to find the module command in the
        environment"""
        with mock.patch('sys.stderr', new_callable=StringIO) as err_output, \
             mock.patch('sys.stdout', new_callable=StringIO) as std_output:
            with mock.patch.dict('os.environ', {}), \
                    mock.patch('os.path.exists', return_value=False), \
                    mock.patch('subprocess.Popen') as mock_subproc_popen:
                process_mock = mock.Mock()
                attrs = {'communicate.return_value': ('print("exec_cmd")',
                                                      'error output'),
                         'wait.return_value': 0}
                process_mock.configure_mock(**attrs)
                mock_subproc_popen.return_value = process_mock
                result = self.module_obj.module('load', 'sierra')
            stdout = std_output.getvalue()
            stderr = err_output.getvalue()
        self.assertEqual('', result)
        self.assertEqual('', stdout)
        self.assertEqual('Unable to load modules; no modulecmd found.\n',
                         stderr)

    def test_module_error(self):
        """Test module command when there is an error loading the module"""
        with mock.patch('sys.stdout', new_callable=StringIO) as std_output:
            with mock.patch.dict('os.environ', {}), \
                    mock.patch('os.path.exists', return_value=True), \
                    mock.patch('subprocess.Popen') as mock_subproc_popen:
                process_mock = mock.Mock()
                attrs = {'communicate.return_value': ('print("exec_cmd")',
                                                      'module:ERROR: error '
                                                      'encountered'),
                         'wait.return_value': 0}
                process_mock.configure_mock(**attrs)
                mock_subproc_popen.return_value = process_mock
                error_msg = 'module:ERROR: error encountered'
                with self.assertRaisesRegexp(RuntimeError, error_msg):
                    self.module_obj.module('load', 'sierra')
            stdout = std_output.getvalue()
        self.assertEqual('exec_cmd\n', stdout)

    def test_module_no_stderr_in_args(self):
        """Test the module method when passed no-stderr"""
        with mock.patch('sys.stderr', new_callable=StringIO) as err_output, \
             mock.patch('sys.stdout', new_callable=StringIO) as std_output:
            with mock.patch.dict('os.environ', {}), \
                    mock.patch('os.path.exists', return_value=True), \
                    mock.patch('subprocess.Popen') as mock_subproc_popen:
                process_mock = mock.Mock()
                attrs = {'communicate.return_value': ('print("exec_cmd")',
                                                      'error_output'),
                         'wait.return_value': 0}
                process_mock.configure_mock(**attrs)
                mock_subproc_popen.return_value = process_mock
                result = self.module_obj.module('load', 'sierra', 'no-stderr')
            stdout = std_output.getvalue()
            stderr = err_output.getvalue()
        self.assertEqual('error_output', result)
        self.assertEqual('exec_cmd\n', stdout)
        self.assertEqual('', stderr)


class TestModuleCommands(unittest.TestCase):
    """Test the module commands present for legacy support"""
    @staticmethod
    def test_module():
        """Test the module function"""
        with mock.patch('Modules.Module.module') as mock_module:
            Modules.module('load', 'sierra/version')
        mock_module.assert_called_once_with('load', 'sierra/version')

    @staticmethod
    def test_module_list():
        """Test the module_list function"""
        with mock.patch('Modules.Module.module_list') as mock_module:
            Modules.module_list()
        mock_module.assert_called_once_with()


class TestGetInitFile(unittest.TestCase):
    """Test the function that gets called when script is run directly"""
    def test_get_init_file(self):
        """Test the method that gets called when script is run directly"""

        modules_home = {'MODULESHOME': '/dummy/path/1'}
        which_side_effects = ['/path/to/modulecmd', None, None]
        find_side_effects = [None, '/fake/path/modules/init/python.py']
        with mock.patch('sys.stdout', new_callable=StringIO) as output:
            with mock.patch.object(sys, 'argv', ['Modules', 'python']):
                with mock.patch.dict(os.environ, modules_home), \
                        mock.patch('Modules.which',
                                   side_effect=which_side_effects), \
                        mock.patch('Modules.find_first_binary',
                                   return_value='/fake/path/modulecmd'), \
                        mock.patch('Modules.find_file_in_list',
                                   side_effect=find_side_effects), \
                        mock.patch('__builtin__.execfile') as exec_file:
                    Modules.get_init_file()
            stdout = output.getvalue()
        self.assertEqual('/fake/path/modules/init/python.py\n', stdout)
        exec_file.assert_not_called()

    def test_get_init_file_no_file_found(self):
        """Test the function that gets called when script is run directly
        when no init file is returned"""
        modules_home = {'MODULESHOME': '/dummy/path/1'}
        which_side_effects = ['/path/to/modulecmd', None, None]
        find_side_effects = [None, None]
        with mock.patch.object(sys, 'argv', ['Modules', 'python']):
            with mock.patch.dict(os.environ, modules_home), \
                    mock.patch('Modules.which',
                               side_effect=which_side_effects), \
                    mock.patch('Modules.find_first_binary',
                               return_value=None), \
                    mock.patch('Modules.find_file_in_list',
                               side_effect=find_side_effects), \
                    mock.patch('__builtin__.execfile') as exec_file:
                error = 'Unable to determine init file for python in ' \
                        'Modules.py'
                with self.assertRaisesRegexp(RuntimeError, error):
                    Modules.get_init_file()
        exec_file.assert_not_called()

    def test_get_init_file_no_param_passed(self):
        """Test the method that gets called when script is run directly when
        no params are passed"""
        error = 'A single init file name parameter is required by ' \
                'Modules.py when invoked directly'
        with mock.patch.object(sys, 'argv', ['Modules']):
            with self.assertRaisesRegexp(RuntimeError, error):
                Modules.get_init_file()


if __name__ == '__main__':
    unittest.main()
