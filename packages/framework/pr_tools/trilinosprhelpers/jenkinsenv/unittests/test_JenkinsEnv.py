#!/usr/bin/env python3
# -*- coding: utf-8; mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""
"""
from __future__ import print_function

try:                                        # pragma: no cover
    import builtins                         # pragma: no cover
except ImportError:                         # pragma: no cover
    import __builtin__ as builtins          # pragma: no cover

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

#from mock import Mock
#from mock import mock_open
#from mock import MagicMock
#from mock import patch


import trilinosprhelpers.jenkinsenv as jenkinsenv



#==============================================================================
#
#                         M O C K   H E L P E R S
#
#==============================================================================




#==============================================================================
#
#                                T E S T S
#
#==============================================================================
class JenkinsEnvTest(TestCase):
    """
    Test JenkinsEnv
    """
    def setUp(self):
        os.environ["BUILD_NUMBER"]    = "build_number_value"
        os.environ["BUILD_ID"]        = "build_id_value"
        os.environ["BUILD_URL"]       = "build_url_value"
        os.environ["NODE_NAME"]       = "node_name_value"
        os.environ["JOB_NAME"]        = "job_name_value"
        os.environ["JOB_BASE_NAME"]   = "job_base_name_value"
        os.environ["BUILD_TAG"]       = "build_tag_value"
        os.environ["JENKINS_HOME"]    = "jenkins_home_value"
        os.environ["JENKINS_URL"]     = "jenkins_url_value"
        os.environ["EXECUTOR_NUMBER"] = "executor_number_value"
        os.environ["NODE_LABELS"]     = "node_labels_value"
        os.environ["JAVA_HOME"]       = "java_home_value"
        os.environ["WORKSPACE"]       = "workspace_value"
        os.environ["SVN_REVISION"]    = "svn_revision_value"
        os.environ["CVS_BRANCH"]      = "cvs_branch_value"
        os.environ["GIT_COMMIT"]      = "git_commit_value"
        os.environ["GIT_URL"]         = "git_url_value"
        os.environ["GIT_BRANCH"]      = "git_branch_value"
        os.environ["JENKINS_JOB_WEIGHT"]  = "jenkins_job_weight_value"
        os.environ["JENKINS_TEST_WEIGHT"] = "jenkins_test_weight_value"

        self._envvar_helper = jenkinsenv.JenkinsEnv()


    def tearDown(self):
        del os.environ["BUILD_NUMBER"]
        del os.environ["BUILD_ID"]
        del os.environ["BUILD_URL"]
        del os.environ["NODE_NAME"]
        del os.environ["JOB_NAME"]
        del os.environ["JOB_BASE_NAME"]
        del os.environ["BUILD_TAG"]
        del os.environ["JENKINS_HOME"]
        del os.environ["JENKINS_URL"]
        del os.environ["EXECUTOR_NUMBER"]
        del os.environ["NODE_LABELS"]
        del os.environ["JAVA_HOME"]
        del os.environ["WORKSPACE"]
        del os.environ["SVN_REVISION"]
        del os.environ["CVS_BRANCH"]
        del os.environ["GIT_COMMIT"]
        del os.environ["GIT_URL"]
        del os.environ["GIT_BRANCH"]
        del os.environ["JENKINS_JOB_WEIGHT"]
        del os.environ["JENKINS_TEST_WEIGHT"]


    def test_build_number(self):
        value = self._envvar_helper.build_number
        self.assertEqual(value, "build_number_value")


    def test_build_id(self):
        value = self._envvar_helper.build_id
        self.assertEqual(value, "build_id_value")


    def test_build_url(self):
        value = self._envvar_helper.build_url
        self.assertEqual(value, "build_url_value")


    def test_node_name(self):
        value = self._envvar_helper.node_name
        self.assertEqual(value, "node_name_value")


    def test_job_name(self):
        value = self._envvar_helper.job_name
        self.assertEqual(value, "job_name_value")


    def test_job_base_name(self):
        value = self._envvar_helper.job_base_name
        self.assertEqual(value, "job_base_name_value")


    def test_build_tag(self):
        value = self._envvar_helper.build_tag
        self.assertEqual(value, "build_tag_value")


    def test_jenkins_home(self):
        value = self._envvar_helper.jenkins_home
        self.assertEqual(value, "jenkins_home_value")


    def test_jenkins_url(self):
        value = self._envvar_helper.jenkins_url
        self.assertEqual(value, "jenkins_url_value")


    def test_executor_number(self):
        value = self._envvar_helper.executor_number
        self.assertEqual(value, "executor_number_value")


    def test_node_labels(self):
        value = self._envvar_helper.node_labels
        self.assertEqual(value, "node_labels_value")


    def test_java_home(self):
        value = self._envvar_helper.java_home
        self.assertEqual(value, "java_home_value")


    def test_workspace(self):
        value = self._envvar_helper.workspace
        self.assertEqual(value, "workspace_value")


    def test_svn_revision(self):
        value = self._envvar_helper.svn_revision
        self.assertEqual(value, "svn_revision_value")


    def test_cvs_branch(self):
        value = self._envvar_helper.cvs_branch
        self.assertEqual(value, "cvs_branch_value")


    def test_git_commit(self):
        value = self._envvar_helper.git_commit
        self.assertEqual(value, "git_commit_value")


    def test_git_url(self):
        value = self._envvar_helper.git_url
        self.assertEqual(value, "git_url_value")


    def test_git_branch(self):
        value = self._envvar_helper.git_branch
        self.assertEqual(value, "git_branch_value")


    def test_jenkins_job_weight(self):
        value = self._envvar_helper.jenkins_job_weight
        self.assertEqual(value, "jenkins_job_weight_value")


    def test_jenkins_test_weight(self):
        value = self._envvar_helper.jenkins_test_weight
        self.assertEqual(value, "jenkins_test_weight_value")

    #def test_(self):
        #value = self._envvar_helper.
        #self.assertEqual(value, "_value")
