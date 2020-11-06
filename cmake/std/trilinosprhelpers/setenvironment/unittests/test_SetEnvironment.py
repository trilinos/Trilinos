#!/usr/bin/env python
# -*- coding: utf-8; mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""
"""
from __future__ import print_function
import sys
sys.dont_write_bytecode = True

import os
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pprint import pprint

import unittest
from unittest import TestCase

# Coverage will always miss one of these depending on the system
# and what is available.
try:                                               # pragma: no cover
    import unittest.mock as mock                   # pragma: no cover
except:                                            # pragma: no cover
    import mock                                    # pragma: no cover

from mock import Mock
from mock import MagicMock
from mock import patch

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from trilinosprhelpers.setenvironment import SetEnvironment


#===============================================================================
#
# Mock Helpers
#
#===============================================================================


def mock_module_pass(*args):
    """
    Mock the module() command that 'passes', returning a 0.
    """
    print("\nmock> module({}) ==> 0".format(args))
    return 0


def mock_module_fail(*args):
    """
    Mock the module() command that 'fails', returning a 1.
    """
    print("\nmock> module({}) ==> 1".format(args))
    return 1


#===============================================================================
#
# Tests
#
#===============================================================================

class SetEnvironmentTest(TestCase):
    """
    Main test driver for the SetEnvironment class
    """
    def setUp(self):
        print("")
        self.maxDiff = None

        self._filename = self.find_config_ini(filename="setenvironment_config.ini")

        self._test_profile_001_truth = {
            'module-list': {
                'sems-gcc': True,
                'sems-boost': True,
                'sems-cmake': True,
                'sems-python': True
            },
            'module-op': [
                ['purge', ''],
                ['use', '/projects/sems/modulefiles/projects'],
                ['load', 'sems-env'],
                ['load', 'sems-gcc/7.3.0'],
                ['load', 'sems-boost/1.63.0/base'],
                ['load', 'sems-cmake/3.10.3'],
                ['load', 'sems-python/3.5.2']
            ],
            'setenv': [   {'key': 'OMP_NUM_THREADS', 'value': '2'},
                          {'key': 'CC', 'value': 'gcc'},
                          {'key': 'CXX', 'value': 'g++'},
                          {'key': 'FC', 'value': 'gfortran'},
                          {'key': 'TEST_ENVVAR_001', 'value': 'TEST_ENVVAR_001_VALUE'}
            ],
            'unsetenv': []
        }

        self._test_profile_002_truth = {
            'module-list': {
                'sems-boost': True,
                'sems-cmake': True,
                'sems-gcc': True,
                'sems-python': True
            },
            'module-op': [
                ['purge', ''],
                ['use', '/projects/sems/modulefiles/projects'],
                ['load', 'sems-env'],
                ['load', 'sems-gcc/7.3.0'],
                ['load', 'sems-boost/1.63.0/base'],
                ['load', 'sems-cmake/3.10.3'],
                ['load', 'sems-python/3.5.2']
            ],
            'setenv': [
                {'key': 'OMP_NUM_THREADS', 'value': '2'},
                {'key': 'CC', 'value': 'gcc'},
                {'key': 'CXX', 'value': 'g++'},
                {'key': 'FC', 'value': 'gfortran'},
                {'key': 'TEST_ENVVAR_001', 'value': 'TEST_ENVVAR_001_VALUE'},
                {'key': 'CC', 'value': 'g++'}
            ],
            'unsetenv': ['CXX']
        }

        self._test_profile_003_truth = {
            'module-list': {
                'sems-boost': True,
                'sems-cmake': True,
                'sems-gcc': True
            },
            'module-op': [
                ['purge', ''],
                ['use', '/projects/sems/modulefiles/projects'],
                ['load', 'sems-env'],
                ['load', 'sems-gcc/7.3.0'],
                ['load', 'sems-boost/1.63.0/base'],
                ['load', 'sems-cmake/3.10.3']
            ],
            'setenv': [
                {'key': 'OMP_NUM_THREADS', 'value': '2'},
                {'key': 'CC', 'value': 'gcc'},
                {'key': 'CXX', 'value': 'g++'},
                {'key': 'FC', 'value': 'gfortran'},
                {'key': 'TEST_ENVVAR_001', 'value': 'TEST_ENVVAR_001_VALUE'}
            ],
            'unsetenv': []
        }

        self._test_profile_004_truth = {
            'module-list': {
                'sems-boost': True,
                'sems-cmake': True,
                'sems-gcc': True
            },
            'module-op': [
                ['purge', ''],
                ['use', '/projects/sems/modulefiles/projects'],
                ['load', 'sems-env'],
                ['load', 'sems-gcc/7.3.0'],
                ['load', 'sems-boost/1.63.0/base'],
                ['load', 'sems-cmake/3.10.3'],
                ['load', 'sems-python/3.5.2'],
                ['unload', 'sems-python']
            ],
            'setenv': [
                {'key': 'OMP_NUM_THREADS', 'value': '2'},
                {'key': 'CC', 'value': 'gcc'},
                {'key': 'CXX', 'value': 'g++'},
                {'key': 'FC', 'value': 'gfortran'},
                {'key': 'TEST_ENVVAR_001', 'value': 'TEST_ENVVAR_001_VALUE'}
            ],
            'unsetenv': []
        }

        self._test_profile_005_truth = {
            'module-list': {
                'sems-boost': True,
                'sems-cmake': True,
                'sems-gcc': True,
                'sems-python': True
            },
            'module-op': [
                ['purge', ''],
                ['use', '/projects/sems/modulefiles/projects'],
                ['load', 'sems-env'],
                ['load', 'sems-gcc/7.3.0'],
                ['load', 'sems-boost/1.63.0/base'],
                ['load', 'sems-cmake/3.10.3'],
                ['load', 'sems-python/3.5.2']
            ],
            'setenv': [
                {'key': 'OMP_NUM_THREADS', 'value': '2'},
                {'key': 'CC', 'value': 'gcc'},
                {'key': 'CXX', 'value': 'g++'},
                {'key': 'FC', 'value': 'gfortran'},
                {'key': 'TEST_ENVVAR_001', 'value': 'TEST_ENVVAR_001_VALUE'},
                {'key': 'TEST_ENVVAR_002', 'value': 'TEST_ENVVAR_002_VALUE/${TEST_ENVVAR_001}'}
            ],
            'unsetenv': []
        }

        self._test_profile_006_truth = {
            'module-list': {
                'sems-boost': True,
                'sems-cmake': True,
                'sems-gcc': True,
                'sems-python': True
            },
            'module-op': [
                ['purge', ''],
                ['use', '/projects/sems/modulefiles/projects'],
                ['load', 'sems-env'],
                ['load', 'sems-gcc/7.3.0'],
                ['load', 'sems-boost/1.63.0/base'],
                ['load', 'sems-cmake/3.10.3'],
                ['load', 'sems-python/3.5.2'],
                ['swap', 'sems-python/3.5.2', 'sems-python/3.8.0']
            ],
            'setenv': [
                {'key': 'OMP_NUM_THREADS', 'value': '2'},
                {'key': 'CC', 'value': 'gcc'},
                {'key': 'CXX', 'value': 'g++'},
                {'key': 'FC', 'value': 'gfortran'},
                {'key': 'TEST_ENVVAR_001', 'value': 'TEST_ENVVAR_001_VALUE'}
            ],
            'unsetenv': []
        }

        self._test_profile_007_truth = {}


    def find_config_ini(self, filename="config.ini"):
        rootpath = "."
        output = None
        for dirpath,dirnames,filename_list in os.walk(rootpath):
            if filename in filename_list:
                output = os.path.join(dirpath, filename)
                break
        if output is None:
            raise FileNotFoundError("Unable to find {} in {}".format(filename, os.getcwd()))
        return output


    def test_SetEnvironment_Profile001_P(self):
        """
        Test loading Profile 001
        - the module() command(s) will be mocked and set to pass (return 0) during the apply() step.
        """
        filename    = self._filename
        profile     = "TEST_PROFILE_001"
        truth       = self._test_profile_001_truth
        module_fail = False
        self._setEnv_test(filename, profile, truth=truth, module_fail=module_fail)


    def test_SetEnvironment_Profile001_F(self):
        """
        Test loading Profile 001
        - the module() command(s) will be mocked and set to fail (return 1) during the apply() step.
        """
        filename    = self._filename
        profile     = "TEST_PROFILE_001"
        truth       = self._test_profile_001_truth
        module_fail = True
        self._setEnv_test(filename, profile, truth=truth, module_fail=module_fail)


    def test_SetEnvironment_Profile002_P(self):
        """
        Test loading Profile 002
        - the module() command(s) will be mocked and set to pass (return 0) during the apply() step.
        - Section name is CaMel CaSe
        - imports TEST_PROFILE_001
        - unsetenv CXX
        - overwrite CC to be g++
        """
        filename    = self._filename
        profile     = "Test_Profile_002"
        truth       = self._test_profile_002_truth
        module_fail = False
        self._setEnv_test(filename, profile, truth=truth, module_fail=module_fail)


    def test_SetEnvironment_Profile003_P(self):
        """
        Test loading Profile 003
        - the module() command(s) will be mocked and set to pass (return 0) during the apply() step.
        - Loads Profile_001
        - module-remove sems-python
        """
        filename    = self._filename
        profile     = "TEST_PROFILE_003"
        truth       = self._test_profile_003_truth
        module_fail = False
        self._setEnv_test(filename, profile, truth=truth, module_fail=module_fail)


    def test_SetEnvironment_Profile004_P(self):
        """
        Test loading Profile 004
        - the module() command(s) will be mocked and set to pass (return 0) during the apply() step.
        - Loads Profile_001
        - add a module-unload of sems-python
        """
        filename    = self._filename
        profile     = "TEST_PROFILE_004"
        truth       = self._test_profile_004_truth
        module_fail = False
        self._setEnv_test(filename, profile, truth=truth, module_fail=module_fail)


    def test_SetEnvironment_Profile005_P(self):
        """
        Test loading Profile 005
        - the module() command(s) will be mocked and set to pass (return 0) during the apply() step.
        - Loads Profile_001
        - add a setenv with expansion of a `${envvar}`
        """
        filename    = self._filename
        profile     = "TEST_PROFILE_005"
        truth       = self._test_profile_005_truth
        module_fail = False
        self._setEnv_test(filename, profile, truth=truth, module_fail=module_fail)

        # Verify that TEST_ENVVAR_002 was expanded properly
        expected_test_envvar_002 = "TEST_ENVVAR_002_VALUE/TEST_ENVVAR_001_VALUE"
        actual_test_envvar_002   = os.environ['TEST_ENVVAR_002']
        print("")
        print("--- {}".format(actual_test_envvar_002))
        print("")
        self.assertEqual(expected_test_envvar_002, actual_test_envvar_002)


    def test_SetEnvironment_Profile006_P(self):
        """
        Test loading Profile 006
        - the module() command(s) will be mocked and set to pass (return 0) during the apply() step.
        - Loads Profile_001
        - adds a module swap command
        """
        filename    = self._filename
        profile     = "TEST_PROFILE_006"
        truth       = self._test_profile_006_truth
        module_fail = False
        self._setEnv_test(filename, profile, truth=truth, module_fail=module_fail)


    def test_SetEnvironment_ModuleRemoveMulti(self):
        """
        Test Module Remove.

        This test loads a sequence of sections that adds a module-load gcc/1.0
        and a module-load gcc/2.0 so that there are two 'gcc' modules being loaded.
        The Module_Remove section then issues a module-remove gcc.

        The expected result of this is that we have an empty set of actions
        since module-remove removes ALL occurrences of the module load operation.
        """
        filename    = self._filename
        profile     = "Module_Remove"
        module_fail = False

        truth = {
            'module-list': {},
            'module-op': [],
            'setenv': [],
            'unsetenv': []
        }

        self._setEnv_test(filename, profile, truth=truth, module_fail=module_fail)


    def test_SetEnvironment_ModuleUnloadMulti(self):
        """
        Test Module Unload.

        This test loads a sequence of sections that adds a module-load gcc/1.0
        and a module-load gcc/2.0 so that there are two 'gcc' modules being loaded.
        The Module_Unload section then issues a `module-unload gcc`.

        The expected result of this is that we will have the following items in the
        module operation list:
        - ['load', 'gcc/1.0']
        - ['load', 'gcc/2.0']
        - ['unload', 'gcc']
        In practice, this wouldn't be likely to happen but one might have a project
        which needs two versions of Python so they might have a module-load python/2.7.3
        but that might be something that environment modules don't allow (I'd need to
        check) since two modules of the same kind could cause problems w/rt to envvars
        such as PYTHON_HOME in the python exampe.

        One thing that the module-unload does do though is that it will remove the
        'gcc' entry from the modle-list dictionary, which is just used for bookkeeping.
        """
        filename    = self._filename
        profile     = "Module_Unload"
        module_fail = False

        truth = {
            'module-list': {},
            'module-op': [['load','gcc/1.0'],['load','gcc/2.0'],['unload','gcc']],
            'setenv': [],
            'unsetenv': []
        }

        self._setEnv_test(filename, profile, truth=truth, module_fail=module_fail)


    def test_SetEnvironment_ModulePurge(self):
        """
        Test the `module-purge` operation.
        """
        filename    = self._filename
        profile     = "Module_Purge"
        module_fail = False

        truth = {
            'module-list': {},
            'module-op': [['purge','']],
            'setenv': [],
            'unsetenv': []
        }

        self._setEnv_test(filename, profile, truth=truth, module_fail=module_fail)


    def test_SetEnvironment_actions_A(self):
        setEnv = SetEnvironment(self._filename, "TEST_PROFILE_001")
        #with self.assertRaisesRegexp(TypeError, "Invalid type provided"):
        with self.assertRaises(TypeError):
            setEnv.actions = "This should raise a TypeError"


    def test_SetEnvironment_actions_B(self):
        setEnv = SetEnvironment(self._filename, "TEST_PROFILE_001")
        #with self.assertRaisesRegexp(KeyError, "setenv"):
        with self.assertRaises(KeyError):
            setEnv.actions = {}


    def test_SetEnvironment_actions_C(self):
        setEnv = SetEnvironment(self._filename, "TEST_PROFILE_001")
        #with self.assertRaisesRegexp(KeyError, "unsetenv"):
        with self.assertRaises(KeyError):
            setEnv.actions = {'setenv': None}


    def test_SetEnvironment_actions_D(self):
        setEnv = SetEnvironment(self._filename, "TEST_PROFILE_001")
        #with self.assertRaisesRegexp(KeyError, "module-op"):
        with self.assertRaises(KeyError):
            setEnv.actions = {'setenv': None, 'unsetenv': None}


    def test_SetEnvironment_actions_E(self):
        setEnv = SetEnvironment(self._filename, "TEST_PROFILE_001")
        #with self.assertRaisesRegexp(KeyError, "module-list"):
        with self.assertRaises(KeyError):
            setEnv.actions = {'setenv': None, 'unsetenv': None, 'module-op': None}


    def test_SetEnvironment_missing_section_name(self):
        setEnv = SetEnvironment(self._filename, "MISSING SECTION NAME!")
        with self.assertRaises(KeyError):
            setEnv.pretty_print()       # this should throw a KeyError


    def test_SetEnvironment_module_too_many_params(self):
        setEnv = SetEnvironment(self._filename, "TEST_PROFILE_001")
        setEnv.actions = {'setenv': None,
                          'unsetenv': None,
                          'module-op': [ ['load','a','b','c'] ],
                          'module-list': None
                          }
        #with self.assertRaisesRegexp(IndexError, 'Invalid number of parameters'):
        with self.assertRaises(IndexError):
            setEnv.pretty_print()
        #with self.assertRaisesRegexp(IndexError, 'Invalid number of parameters'):
        with self.assertRaises(IndexError):
            setEnv.apply()


    def test_SetEnvironment_expand_envvars_error(self):
        setEnv = SetEnvironment(self._filename, "TEST_PROFILE_001")
        #with self.assertRaisesRegexp(KeyError, 'Required environment variable .+ does not exist'):
        with self.assertRaises(KeyError):
            setEnv._expand_envvars_in_string("envvar ${ZzZZLk23j45hnApDf} should not be found.")


    def test_SetEnvironment_missing_file(self):
        #with self.assertRaisesRegexp(FileNotFoundError, 'No such file or directory:'):
        setEnv = SetEnvironment("no_file", "TEST_PROFILE_001")
        with patch("sys.exit", return_value=1) as m_exit:
            setEnv.config
            m_exit.assert_called_once()


    def test_SetEnvironment_config(self):
        print("")
        setEnv = SetEnvironment(self._filename, "TEST_PROFILE_001")

        self.assertEqual(setEnv.config.has_section("SAMPLE_MAP"), True)
        self.assertEqual(setEnv.config.has_option("SAMPLE_MAP", "key1"), True)
        self.assertEqual(setEnv.config.get("SAMPLE_MAP", "key1"), "value1")

        options = setEnv.config.options("SAMPLE_MAP")
        for opt in options:
            print("---- {}".format(opt))
        for k,v in setEnv.config.items("SAMPLE_MAP"):
            print("---- {} : {}".format(k,v))


    def test_SetEnvironment_pretty_print_envvars(self):
        """

        """
        print("")
        envvars = [
            "SETENVIRONMENT_TEST_"
            ]

        os.environ["SETENVIRONMENT_TEST_ENVVAR_001"]   = "foobar"
        os.environ["SETENVIRONMENT_HIDDEN_ENVVAR_001"] = "baz"

        setEnv = SetEnvironment(self._filename, "TEST_PROFILE_001")
        with patch('sys.stdout', new=StringIO()) as m_stdout:
            setEnv.pretty_print_envvars(envvar_filter=envvars)
            self.assertIn("-- SETENVIRONMENT_TEST_ENVVAR_001 = foobar", m_stdout.getvalue())
            self.assertIn("-- SETENVIRONMENT_HIDDEN_ENVVAR_001", m_stdout.getvalue())
            self.assertNotIn("-- SETENVIRONMENT_HIDDEN_ENVVAR_001 = baz", m_stdout.getvalue())

        # Filtered + keys_only should print out only the one key
        with patch('sys.stdout', new=StringIO()) as m_stdout:
            setEnv.pretty_print_envvars(envvar_filter=envvars, filtered_keys_only=True)
            self.assertIn("-- SETENVIRONMENT_TEST_ENVVAR_001 = foobar", m_stdout.getvalue())
            self.assertNotIn("-- SETENVIRONMENT_HIDDEN", m_stdout.getvalue())

        # No options should print all envvars + values
        with patch('sys.stdout', new=StringIO()) as m_stdout:
            setEnv.pretty_print_envvars()
            self.assertIn("-- SETENVIRONMENT_TEST", m_stdout.getvalue())
            self.assertIn("-- SETENVIRONMENT_HIDDEN", m_stdout.getvalue())

        # No filter but we say show filtered keys only should result in
        # print all keys + values.
        with patch('sys.stdout', new=StringIO()) as m_stdout:
            setEnv.pretty_print_envvars(filtered_keys_only=True)
            self.assertIn("-- SETENVIRONMENT_TEST", m_stdout.getvalue())
            self.assertIn("-- SETENVIRONMENT_HIDDEN", m_stdout.getvalue())

        # cleanup
        del os.environ["SETENVIRONMENT_TEST_ENVVAR_001"]
        del os.environ["SETENVIRONMENT_HIDDEN_ENVVAR_001"]


    def _setEnv_test(self, filename, profile, truth=None, module_fail=False):
        """
        Test the instantiation of a SetEnvironment class. Loads the
        ini file and applyies its settings.
        The 'module' command is mocked out during the apply() step.

        Args:
            filename (string) : The filename of the .ini file to load.
            profile (string)  : The profile is the <section> name in the .ini file to load.
            truth (dict)      : A dict contianing the ground-truth of the `actions` entry from
                                SetEnvironment.
            module_fail (bool): If true, the module() command(s) executed by SetEnvironment.apply()
                                will fail (returning a 1). Default: False.
        """
        print("\nmodule_fail = {}".format(module_fail))
        setEnv = SetEnvironment(filename, profile)

        print("\nsetEnv.pretty_print():")
        setEnv.pretty_print()

        print("\nsetEnv.actions():")
        pprint(setEnv.actions, indent=4, width=90)

        print("\nExpected Actions:")
        pprint(truth, indent=4, width=90)
        print("")

        # Validation Checks
        self.assertIsInstance( setEnv.actions, dict )
        self.assertDictEqual( truth, setEnv.actions )

        if False == module_fail:
            # patch the 'module()' command to 'pass'
            with patch('trilinosprhelpers.setenvironment.ModuleHelper.module', side_effect=mock_module_pass):
                s = setEnv.apply()
            print("STATUS: {}".format(s))
            self.assertEqual(0, s)
        else:
            with patch('trilinosprhelpers.setenvironment.ModuleHelper.module', side_effect=mock_module_fail):
                s = setEnv.apply()
            print("STATUS: {}".format(s))
            self.assertEqual(1, s)



if __name__ == "__main__":
    unittest.main()  # pragma nocover


