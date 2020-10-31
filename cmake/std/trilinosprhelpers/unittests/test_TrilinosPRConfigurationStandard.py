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
#import subprocess

#from trilinosprhelpers import setenvironment
#from trilinosprhelpers import sysinfo


import trilinosprhelpers


#==============================================================================
#
#                         M O C K   H E L P E R S
#
#==============================================================================
def mock_chdir(*args, **kwargs):
    print("MOCK> os.chdir('{}')".format(args[0]))
    return 0


def mock_subprocess_check_call(*args, **kwargs):

    cmd = " ".join(args[0])
    print("MOCK> subprocess.check_call({})".format(cmd))
    return 0


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


def mock_module_apply(*args, **kwargs):
    """
    Mock handler for ModuleHelper.module() calls.
    """
    cmd = ", ".join(["'{}'".format(x) for x in args])                          # pragma: no cover
    print("MOCK> module({})".format(cmd))                                      # pragma: no cover
    return 0



#==============================================================================
#
#                                T E S T S
#
#==============================================================================
class TrilinosPRConfigurationStandardTest(TestCase):
    """
    Test TrilinsoPRConfigurationStandard class
    """
    def setUp(self):
        os.environ["PULLREQUEST_CDASH_TRACK"] = "Pull Request"

        # Find the config file
        config_file = 'trilinos_pr_test.ini'
        self._config_file = self.find_config_ini(config_file)

        # Set up dummy command line arguments
        self._args = self. dummy_args()

        # Set up some global mock patches
        self.patch_cpu_count = patch('multiprocessing.cpu_count', return_value=64)
        self.mock_cpu_count  = self.patch_cpu_count.start()

        self.patch_os_chdir = patch('os.chdir', side_effect=mock_chdir)
        self.mock_chdir = self.patch_os_chdir.start()

        self.patch_subprocess_check_call = patch('subprocess.check_call', side_effect=mock_subprocess_check_call)
        self.mock_subprocess_check_call = self.patch_subprocess_check_call.start()

        self.patch_subprocess_check_output = patch('subprocess.check_output',
                                                   side_effect=mock_subprocess_check_output)
        self.mock_subprocess_check_output = self.patch_subprocess_check_output.start()

        self.patch_modulehelper_module = patch('trilinosprhelpers.setenvironment.ModuleHelper.module',
                                               side_effect=mock_module_apply)
        self.mock_modulehelper_module  = self.patch_modulehelper_module.start()


    def tearDown(self):
        del os.environ["PULLREQUEST_CDASH_TRACK"]

        # Shut down global mock patches
        self.patch_cpu_count.stop()
        self.patch_os_chdir.stop()
        self.patch_subprocess_check_call.stop()
        self.patch_subprocess_check_output.stop()
        self.patch_modulehelper_module.stop()


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
            pullrequest_cdash_track="Pull Request",
            jenkins_job_number=99,
            pullrequest_number='0000',
            pullrequest_config_file=self._config_file,
            workspace_dir=".",
            filename_packageenables="../packageEnables.cmake",
            filename_subprojects="../package_subproject_list.cmake",
            mode="standard",
            req_mem_per_core=3.0,
            max_cores_allowed=12,
            num_concurrent_tests=-1,
            dry_run = False
        )
        return output


    def find_config_ini(self, filename="config.ini"):
        rootpath = "."
        output = None
        for dirpath,dirnames,filename_list in os.walk(rootpath):
            if filename in filename_list:
                output = os.path.join(dirpath, filename)
                break
        return output


    def test_TrilinosPRConfigurationStandardExec(self):
        """
        Test the Standard Configuration
        """
        print("")
        args = self.dummy_args()
        pr_config = trilinosprhelpers.TrilinosPRConfigurationStandard(args)

        # prepare step
        ret = pr_config.prepare_test()
        self.assertEqual(ret, 0)
        self.mock_cpu_count.assert_called()

        # execute step
        ret = pr_config.execute_test()
        self.mock_chdir.assert_called_once()
        self.mock_subprocess_check_call.assert_called()
        self.assertEqual(ret, 0)


    def test_TrilinosPRConfigurationStandardDryRun(self):
        """
        Test the Standard Configuration
        - Change args to enable dry_run mode.
        """
        args = self.dummy_args()
        args.dry_run = True
        pr_config = trilinosprhelpers.TrilinosPRConfigurationStandard(args)

        # prepare step
        ret = pr_config.prepare_test()
        self.assertEqual(ret, 0)
        self.mock_cpu_count.assert_called()

        # execute step
        ret = pr_config.execute_test()
        self.assertEqual(ret, 0)


    def test_TrilinosPRConfigurationStandardPython3(self):
        """
        Test the Standard Configuration
        - Change args to enable:
            - pullrequest_build_name = "Trilinos_pullrequest_python_3"
            - dry_run = True
        - Change args to enable dry_run mode.
        """
        args = self.dummy_args()
        args.pullrequest_build_name = "Trilinos_pullrequest_python_3"
        pr_config = trilinosprhelpers.TrilinosPRConfigurationStandard(args)

        # prepare step
        ret = pr_config.prepare_test()
        self.assertEqual(ret, 0)
        self.mock_cpu_count.assert_called()

        # execute step
        #ret = pr_config.execute_test()
        #self.assertEqual(ret, 0)
