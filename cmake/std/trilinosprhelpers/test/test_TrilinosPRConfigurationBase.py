#!/usr/bin/env python
# -*- coding: utf-8; mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""
"""
from __future__ import print_function

try:                                        # pragma: no cover
    import builtins                         # pragma: no cover
except ImportError:                         # pragma: no cover
    import __builtin__ as builtins          # pragma: no cover

import copy
from glob import glob
import os
import sys

if sys.version_info >= (3,0):               # pragma: no cover
    from io import StringIO                 # pragma: no cover
else:                                       # pragma: no cover
    from io import BytesIO as StringIO      # pragma: no cover

#import unittest
from unittest import TestCase

try:                                        # pragma: no cover
    import unittest.mock as mock            # pragma: no cover
except:                                     # pragma: no cover
    import mock                             # pragma: no cover

from mock import call
from mock import Mock
from mock import mock_open
from mock import MagicMock
from mock import patch

import argparse
import multiprocessing
import subprocess

import trilinosprhelpers
from trilinosprhelpers import setenvironment
from trilinosprhelpers import sysinfo




#==============================================================================
#
#                         M O C K   H E L P E R S
#
#==============================================================================
def mock_str_builtins_open():
    """
    This is a helper for the builtins issue where python3 uses 'builtins' and
    python2 uses __builtin__ :/
    """
    output = "builtins.open"
    if sys.version_info < (3,0):                                              # pragma: no cover
        output = "__builtin__.open"                                           # pragma: no cover
    print("MOCK> open()")                                                     # pragma: no cover
    return output


def mock_subprocess_check_output(*args, **kwargs):
    """
    Mock out a subprocess.check_output()
    """
    params = copy.deepcopy(args[0])
    for k,v in kwargs.items():
        params.append("{}={}".format(k,v))
    output = "--- subprocess.check_output({})".format(", ".join(params))

    print("MOCK> mock_packageEnables_check_output()")
    for k in args[0]:                                                         # pragma: no cover
        print("    - '{}'".format(k))                                         # pragma: no cover
    for k,v in kwargs.items():                                                # pragma: no cover
        print("    - {}={}".format(k,v))                                      # pragma: no cover
    print("")
    return str.encode(output)


def mock_subprocess_check_call(*args, **kwargs):
    """
    Mock out a subprocess.check_call()
    """
    params = copy.deepcopy(args[0])                                            # pragma: no cover
    for k,v in kwargs.items():                                                 # pragma: no cover
        params.append("{}={}".format(k,v))                                     # pragma: no cover
    output = "--- subprocess.check_call({})".format(", ".join(params))         # pragma: no cover

    print("MOCK> mock_packageEnables_check_call()")                            # pragma: no cover
    for k in args[0]:                                                          # pragma: no cover
        print("    - '{}'".format(k))                                          # pragma: no cover
    for k,v in kwargs.items():                                                 # pragma: no cover
        print("    - {}={}".format(k,v))                                       # pragma: no cover
    print("")                                                                  # pragma: no cover
    return str.encode(output)                                                  # pragma: no cover


def mock_subprocess_raise_CalledProcessError(*args, **kwargs):
    cmd = " ".join(args[0])                                                    # pragma: no cover
    print("MOCK: Raising subprocess.CalledProcessError")                       # pragma: no cover
    print("    : Command is `{}`".format(cmd))                                 # pragma: no cover
    raise subprocess.CalledProcessError(1, cmd)                                # pragma: no cover


def mock_modulehelper_module_ok(*args, **kwargs):
    """
    Mock replacement for a ModuleHelper.module() call that succeeds
    """
    cmd = ", ".join(["'{}'".format(x) for x in args])
    print("MOCK> module({})    QAPLA'!".format(cmd))                           # pragma: no cover
    return 0                                                                   # pragma: no cover


def mock_modulehelper_module_fail(*args, **kwargs):
    """
    Mock replacement for a ModuleHelper.module() call that fails
    """
    cmd = ", ".join(["'{}'".format(x) for x in args])
    print("MOCK> module({})  # FAILS".format(cmd))                             # pragma: no cover
    return 1                                                                   # pragma: no cover


def mock_early_return(*args, **kwargs):
    """
    Mock helper that prints out the args and quits with nonzero return value.
    """
    print("MOCK>", args[0])                                                    # pragma: no cover
    return 1                                                                   # pragma: no cover


def mock_se_apply_pass(*args, **kwargs):
    """
    Mock for SetEnvironment.apply() call that would pass (i.e., return 0)

    Returns:
        int 0 for PASS
    """
    params = ", ".join(["{}={}".format(k,v) for k,v in kwargs.items()])         # pragma: no cover
    print("MOCK> SetEnvironment.apply({})".format(params))                      # pragma: no cover
    return 0                                                                    # pragma: no cover


def mock_se_apply_fail(*args, **kwargs):
    """
    Mock for SetEnvironment.apply() call that would fail (i.e., return 1)

    Returns:
        int 1 for FAIL
    """
    params = ", ".join(["{}={}".format(k,v) for k,v in kwargs.items()])         # pragma: no cover
    print("MOCK> SetEnvironment.apply({})".format(params))                      # pragma: no cover
    return 1                                                                    # pragma: no cover


#==============================================================================
#
#                                T E S T S
#
#==============================================================================
class TrilinosPRConfigurationTest(TestCase):
    """
    Test TrilinsoPRConfigurationBase class
    """
    def setUp(self):
        #os.environ["TRILINOS_SOURCE_BRANCH"]  = "trilinos_source_branch_value"
        #os.environ["TRILINOS_SOURCE_REPO"]    = "trilinos_source_repo_value"
        #os.environ["TRILINOS_SOURCE_SHA"]     = "trilinos_source_sha_value"
        #os.environ["TRILINOS_TARGET_BRANCH"]  = "trilinos_target_branch_value"
        #os.environ["TRILINOS_TARGET_REPO"]    = "trilinos_target_repo_value"
        #os.environ["TRILINOS_TARGET_SHA"]     = "trilinos_target_sha_value"
        #os.environ["PULLREQUESTNUM"]          = "0000"
        #os.environ["PULLREQUEST_CDASH_TRACK"] = "Pull Request"

        # Find the config file
        #config_file = 'test_config.ini'
        config_file = 'trilinos_pr_test.ini'
        self._config_file = self.find_config_ini(config_file)

        # Set up dummy command line arguments
        self._args = self. dummy_args_python2()

        # Create SetEnvironment object for tests
        self._env  = setenvironment.SetEnvironment(filename=self._config_file, profile="Trilinos_pullrequest_gcc_4.8.4")

        # Set up some global mock patches
        self.patch_cpu_count = patch('multiprocessing.cpu_count', return_value=64)
        self.mock_cpu_count  = self.patch_cpu_count.start()



    def tearDown(self):
        #del os.environ["PULLREQUESTNUM"]
        #del os.environ["PULLREQUEST_CDASH_TRACK"]

        self.patch_cpu_count.stop()



    def dummy_args(self):
        """
        Generate dummy command line arguments
        """
        output = argparse.Namespace(
            source_repo_url="https://github.com/trilinos/Trilinos",
            source_branch_name="source_branch_name",
            target_repo_url="https://github.com/trilinos/Trilinos",
            target_branch_name="develop",
            pullrequest_build_name="Trilinos_pullrequest_gcc_7.2.0",
            jenkins_job_number=99,
            pullrequest_number='0000',
            pullrequest_cdash_track="Pull Request",
            pullrequest_config_file=self._config_file,
            workspace_dir=".",
            filename_packageenables="../packageEnables.cmake",
            filename_subprojects="../package_subproject_list.cmake",
            mode="standard",
            req_mem_per_core=3.0,
            max_cores_allowed=12,
            num_concurrent_tests=-1,
            dry_run=False
        )
        return output


    def dummy_args_num_concurrent_tests(self, value=4):
        """
        Testing the case where args.num_concurrent_tests
        is changed.
        """
        args = copy.deepcopy(self.dummy_args())
        args.num_concurrent_tests = value
        return args


    def dummy_args_python2(self):
        """
        Generate dummy command line arguments
        """
        args = copy.deepcopy(self.dummy_args())
        args.pullrequest_build_name = "Trilinos_pullrequest_python_2"
        return args


    def dummy_args_gcc_720(self):
        args = copy.deepcopy(self.dummy_args())
        args.pullrequest_build_name = "Trilinos_pullrequest_gcc_7.2.0"
        return args


    def dummy_args_master_pass(self):
        """
        Modify arguments to test a develop->master with a valid
        incoming branch name (master_merge_YYYYMMDD_HHMMSS)
        """
        args = copy.deepcopy(self.dummy_args())
        args.target_branch_name = "master"
        args.source_branch_name = "master_merge_20200101_000000"
        return args


    def dummy_args_master_fail(self):
        """
        Modify arguments to test a develop->master with an invalid
        incoming branch name (master_merge_YYYYMMDD_HHMMSS)
        """
        args = copy.deepcopy(self.dummy_args())
        args.target_branch_name = "master"
        args.source_branch_name = "invalid_source_branch_name"
        return args


    def find_config_ini(self, filename="config.ini"):
        rootpath = "."
        output = None
        for dirpath,dirnames,filename_list in os.walk(rootpath):
            if filename in filename_list:
                output = os.path.join(dirpath, filename)
                break
        return output


    def test_TrilinosPRConfigurationBase(self):
        """
        Tests if we can instantiate a TrilinosPRConfiguration object.
        """
        args = self.dummy_args_python2()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        #with patch('sysinfo.SysInfo.compute_num_usable_cores', return_value=6):
        #    print("pr_config.concurrency_build    = {}".format(pr_config.concurrency_build))
        print("pr_config.max_cores_allowed    = {}".format(pr_config.max_cores_allowed))
        print("pr_config.arg_req_mem_per_core = {}".format(pr_config.arg_req_mem_per_core))
        print("pr_config.max_test_parallelism = {}".format(pr_config.max_test_parallelism))
        print("pr_config.concurrency_build    = {}".format(pr_config.concurrency_build))
        print("pr_config.concurrency_test     = {}".format(pr_config.concurrency_test))

        self.assertEqual(pr_config.max_cores_allowed,          12)
        self.assertAlmostEqual(pr_config.arg_req_mem_per_core, 3.0)
        self.assertEqual(pr_config.max_test_parallelism,       4)

        # This will vary b/c it's proven tricky to mock out the system memory for this
        # test so far.
        self.assertGreaterEqual(pr_config.concurrency_build, 1)

        self.assertEqual(pr_config.concurrency_test, 3)


    def test_TrilinosPRConfigurationValidateBranchNameDevelop(self):
        args = self.dummy_args()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)
        pr_config.validate_branch_constraints()


    def test_TrilinosPRConfigurationValidateBranchNameMasterPASS(self):
        args = self.dummy_args_master_pass()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)
        pr_config.validate_branch_constraints()


    def test_TrilinosPRConfigurationValidateBranchNameMasterFAIL(self):
        args = self.dummy_args_master_fail()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)
        with patch('sys.stdout', new=StringIO()) as fake_out:
            with patch('sys.exit', side_effect=mock_early_return) as m:
                pr_config.validate_branch_constraints()
                m.assert_called_once()
                self.assertTrue( "ERROR:" in fake_out.getvalue())


    def test_TrilinosPRConfigurationCDashTrack(self):
        args = args = self.dummy_args_python2()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)
        cdash_track = pr_config.arg_pullrequest_cdash_track
        print("--- cdash_track = {}".format(cdash_track))
        self.assertEqual(cdash_track, "Pull Request")


    def test_TrilinosPRConfigurationBuildNamePython2(self):
        args = self.dummy_args_python2()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)
        build_name = pr_config.pullrequest_build_name
        print("--- build_name = {}".format(build_name))
        expected_build_name = "PR-{}-test-{}-{}".format(args.pullrequest_number, args.pullrequest_build_name, args.jenkins_job_number)
        self.assertEqual(build_name, expected_build_name)


    def test_TrilinosPRConfigurationBaseBuildNameGCC720(self):
        args = self.dummy_args_gcc_720()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)
        build_name = pr_config. pullrequest_build_name
        print("--- build_name = {}".format(build_name))
        expected_build_name = "PR-{}-test-{}-{}".format(args.pullrequest_number, args.pullrequest_build_name, args.jenkins_job_number)
        self.assertEqual(build_name, expected_build_name)


    def test_TrilinosPRConfigurationBasePackageEnablesPython2(self):
        print("")
        args = self.dummy_args_python2()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        # pre-load the property
        pr_config.get_property_from_config("ENABLE_MAP", args.pullrequest_build_name)

        with patch(mock_str_builtins_open(), new_callable=mock_open()) as m_open:
            with patch('subprocess.check_output', side_effect=mock_subprocess_check_output) as m_call:
                pr_config.create_package_enables_file(dryrun=False)
                calls = [ call('packageEnables.cmake','w'),
                          call('package_subproject_list.cmake','w')
                        ]
                m_open.assert_has_calls(calls, any_order=True)


    def test_TrilinosPRConfigurationBasePackageEnablesPython2_dryrun(self):
        """
        Test the PackageEnables generator in DryRun mode
        """
        print("")
        args = self.dummy_args_python2()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        # pre-load the property
        pr_config.get_property_from_config("ENABLE_MAP", args.pullrequest_build_name)

        with patch(mock_str_builtins_open(), new_callable=mock_open()) as m_open:
            pr_config.create_package_enables_file(dryrun=True)
            calls = [ call('packageEnables.cmake','w'),
                      call('package_subproject_list.cmake','w')
                    ]
            m_open.assert_has_calls(calls, any_order=True)
            print(".")


    def test_TrilinosPRConfigurationBasePackageEnablesGCC720(self):
        print("")
        args = self.dummy_args_gcc_720()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        # pre-load the property
        pr_config.get_property_from_config("ENABLE_MAP", args.pullrequest_build_name)

        with patch('subprocess.check_call', side_effect=mock_subprocess_check_call) as m_call:
            with patch('subprocess.check_output', side_effect=mock_subprocess_check_output) as m_output:
                pr_config.create_package_enables_file(dryrun=False)
                m_call.assert_called_once()
                m_output.assert_called_once()


    def test_TrilinosPRConfigurationBasePackageEnablesGCC720_dryrun(self):
        """
        Test the PackageEnables generator in DryRun mode
        We can remove this when we take down the SCAFFOLDING from development.
        """
        print("")
        args = self.dummy_args_gcc_720()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        # pre-load the property
        pr_config.get_property_from_config("ENABLE_MAP", args.pullrequest_build_name)
        pr_config.create_package_enables_file(dryrun=True)


    def test_TrilinosPRConfigurationBasePackageEnablesGCC720_ERROR(self):
        """
        Test the call to create_package_enables_file() where there should
        be an error.
        """
        print("")
        args = self.dummy_args_gcc_720()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        # pre-load the property
        pr_config.get_property_from_config("ENABLE_MAP", args.pullrequest_build_name)

        with patch('sys.stdout', new=StringIO()) as fake_out:
            with self.assertRaises( subprocess.CalledProcessError ) as m_err:
                with patch('subprocess.check_call', side_effect=mock_subprocess_check_call) as m_call:
                    with patch('subprocess.check_output', side_effect=mock_subprocess_raise_CalledProcessError) as m_output:
                        pr_config.create_package_enables_file(dryrun=False)

        self.assertTrue( "There was an issue generating packageEnables.cmake" in fake_out.getvalue())


    def test_TrilinosPRConfigurationBaseProperty_concurrency_test_err(self):
        """
        Test the condition where an incorrect type is passed into the config_data
        setter.

        Mock multiprocessing.cpu_count to return 64 and we use 4 for the
        num_concurrent_tests parameter. concurrency_test should return 16.
        """
        args = self.dummy_args_num_concurrent_tests(value=4)
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        concurrency_test = pr_config.concurrency_test
        print("concurrency_test = {}".format(concurrency_test))
        self.assertEqual(concurrency_test, 4)


    def test_TrilinosPRConfigurationBaseProperty_config_script(self):
        """
        Validate that the property config_script loads properly.
        Since dummy args is loading the configuration for "Trilinos_pullrequest_gcc_7.2.0"
        the mapped configuration script should be loading "PullRequestLinuxGCC7.2.0TestingSettings.cmake"
        """
        args = self.dummy_args()

        # Test the gcc 7.2.0 mapping
        args.pullrequest_build_name = "Trilinos_pullrequest_gcc_7.2.0"
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)
        self.assertEqual(pr_config.config_script, "PullRequestLinuxGCC7.2.0TestingSettings.cmake")


    def test_TrilinosPRConfigurationBaseProperty_get_property_from_config_PASS(self):
        print("")
        args = self.dummy_args()
        args.pullrequest_build_name = "Trilinos_pullrequest_python_3"
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        # Load the ENABLE_MAP job information to test `get_property_from_config`.
        package_enables = pr_config.get_property_from_config("ENABLE_MAP",
                                                             args.pullrequest_build_name)

        print("--- package_enables = {}".format(package_enables))
        print("--- expected_value  = TrilinosFrameworkTests")
        self.assertEqual(package_enables, "TrilinosFrameworkTests")


    def test_TrilinosPRConfigurationBaseProperty_get_property_from_config_FAIL(self):
        """
        Test the condition that the SECTION is missing from the .ini file.
        This should return the default value and print out a warning.
        """
        print("")
        args = self.dummy_args()
        args.pullrequest_build_name = "Trilinos_pullrequest_python_3"
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        # Try to load a nonexistent section.
        print("-----[ TEST 1 ]-----------------------")
        with patch('sys.stdout', new=StringIO()) as fake_out:
            package_enables = pr_config.get_property_from_config("NONEXISTENTSECTION", "N/A")
        print(fake_out.getvalue())
        print("--- package_enables = {}".format(package_enables))
        print("--- expected_value  = {}".format(None))
        self.assertEqual(package_enables, None)
        self.assertIn("WARNING", fake_out.getvalue())

        # Try to load nonexistent section but with a different default value
        print("-----[ TEST 2 ]-----------------------")
        default_value="OOPS!"
        with patch('sys.stdout', new=StringIO()) as fake_out:
            package_enables = pr_config.get_property_from_config("NONEXISTENTSECTION", "N/A", default=default_value)
        print(fake_out.getvalue())
        print("--- package_enables = {}".format(package_enables))
        print("--- expected_value  = {}".format(default_value))
        self.assertEqual(package_enables, default_value)
        self.assertIn("WARNING", fake_out.getvalue())


    def test_TrilinosPRConfigurationBaseProperty_get_multi_property_from_config_PASS(self):
        print("")
        args = self.dummy_args()
        args.pullrequest_build_name = "Trilinos_pullrequest_gcc_8.3.0_installation_testing"
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        # Load the ENABLE_MAP job information to test `get_multi_property_from_config`.
        print("-----[ TEST 1 ]-----------------------")
        package_enables = pr_config.get_multi_property_from_config("ENABLE_MAP",
                                                                   args.pullrequest_build_name)
        print("--- package_enables = {}".format(package_enables))
        print("--- expected_value  = Teuchos,Tpetra")
        self.assertEqual(package_enables, "Teuchos,Tpetra")

        # Change the delimiter
        print("-----[ TEST 2 ]-----------------------")
        package_enables = pr_config.get_multi_property_from_config("ENABLE_MAP",
                                                                   args.pullrequest_build_name,
                                                                   delimeter=" ")
        print("--- package_enables = {}".format(package_enables))
        print("--- expected_value  = Teuchos Tpetra")
        self.assertEqual(package_enables, "Teuchos Tpetra")


    def test_TrilinosPRConfigurationBaseProperty_get_multi_property_from_config_FAIL(self):
        """
        Test a failure condition where the section does not exist.
        """
        print("")
        args = self.dummy_args()
        args.pullrequest_build_name = "Trilinos_pullrequest_gcc_8.3.0_installation_testing"
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        # Load the ENABLE_MAP job information to test `get_multi_property_from_config`.
        print("-----[ TEST 1 ]-----------------------")
        package_enables = pr_config.get_multi_property_from_config("NONEXISTENTSECTION",
                                                                   args.pullrequest_build_name)
        print("--- package_enables = {}".format(package_enables))
        self.assertEqual(package_enables, None)

        # Change the default
        print("-----[ TEST 2 ]-----------------------")
        package_enables = pr_config.get_multi_property_from_config("NONEXISTENTSECTION",
                                                                   args.pullrequest_build_name,
                                                                   default="OOPS!")
        print("--- package_enables = {}".format(package_enables))
        self.assertEqual(package_enables, "OOPS!")

        # Test SECTION exists, OPTION does not exist
        print("-----[ TEST 3 ]-----------------------")
        package_enables = pr_config.get_multi_property_from_config("ENABLE_MAP",
                                                                   "NONEXISTENT_OPTION")
        print("--- package_enables = {}".format(package_enables))
        self.assertEqual(package_enables, None)


    def test_TrilinosPRConfigurationBaseProperty_subprojects_file(self):
        """
        Check property: filename_subprojects
        """
        print("")
        args = self.dummy_args()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        filename_subprojects = pr_config.arg_filename_subprojects
        print("--- filename_subprojects = {}".format(filename_subprojects))

        self.assertEqual(filename_subprojects, "../package_subproject_list.cmake")


    def test_TrilinosPRConfigurationBaseProperty_package_enables_file(self):
        """
        Check property: package_enables_file
        """
        print("")
        args = self.dummy_args()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        package_enables_file = pr_config.arg_filename_packageenables
        print("--- package_enables_file = {}".format(package_enables_file))
        self.assertEqual(package_enables_file, "../packageEnables.cmake")


    def test_TrilinosPRConfigurationBaseProperty_working_directory_ctest(self):
        """
        Check property: working_directory_ctest
        """
        print("")
        args = self.dummy_args()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        working_directory_ctest = pr_config.working_directory_ctest
        print("--- working_directory:       {}".format(pr_config.arg_workspace_dir))
        print("--- working_directory_ctest: {}".format(working_directory_ctest))
        self.assertIn("TFW_testing_single_configure_prototype", working_directory_ctest)


    def test_TrilinosPRConfigurationBase_prepare_test(self):
        """
        Test the prepare_test method
        """
        print("")
        args = self.dummy_args()
        args.max_cores_allowed=-1

        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        # Let's check the expected values for concurrency since this
        # test sets max_cores_allowed to the max detected.
        self.assertGreater(pr_config.concurrency_build, 1)
        self.assertEqual(pr_config.concurrency_test,   16)
        self.assertEqual(pr_config.max_cores_allowed,  64)

        with patch('trilinosprhelpers.setenvironment.ModuleHelper.module',
                   side_effect=mock_modulehelper_module_ok):
            with patch('subprocess.check_call',
                       side_effect=mock_subprocess_check_call) as m_call:
                with patch('subprocess.check_output',
                           side_effect=mock_subprocess_check_output) as m_output:
                    ret = pr_config.prepare_test()
                    self.assertEqual(ret, 0)


    def test_TrilinosPRConfigurationBase_prepare_test_FAIL(self):
        """
        Test the prepare_test method where it would fail due to
        mock_modulehelper_module_failodule loading problems in
        SetEnvironment()
        """
        print("")
        args = self.dummy_args()
        args.max_cores_allowed=-1

        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        # Let's check the expected values for concurrency since this
        # test sets max_cores_allowed to the max detected.
        self.assertGreater(pr_config.concurrency_build, 1)
        self.assertEqual(pr_config.concurrency_test,   16)
        self.assertEqual(pr_config.max_cores_allowed,  64)

        # Test the case that
        with patch('trilinosprhelpers.setenvironment.ModuleHelper.module',
                   side_effect=mock_modulehelper_module_ok):
            with patch('trilinosprhelpers.setenvironment.SetEnvironment.apply',
                       side_effect=mock_se_apply_fail):
                with self.assertRaises(Exception):
                    pr_config.prepare_test()

        # Test the case that the module() command fails down in SetEnvironment.apply()
        # when we throw an exception on failure.
        with patch('trilinosprhelpers.setenvironment.ModuleHelper.module',
                   side_effect=mock_modulehelper_module_fail):
            with self.assertRaises(Exception):
                pr_config.prepare_test()


    def test_TrilinosPRConfigurationBase_execute_test(self):
        """
        Executes the test. This method should be considered 'virtual' and should
        be overridden.
        """
        print("")
        args = self.dummy_args()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationBase(args)

        with self.assertRaises(NotImplementedError):
            pr_config.execute_test()

mock_modulehelper_module_fail
