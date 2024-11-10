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
class TrilinosJenkinsEnvTest(TestCase):
    """
    Test TrilinosJenkinsEnv
    """
    def setUp(self):

        os.environ["TRILINOS_SOURCE_REPO"]    = "trilinos_source_repo_value"
        os.environ["TRILINOS_SOURCE_SHA"]     = "trilinos_source_sha_value"
        os.environ["TRILINOS_TARGET_BRANCH"]  = "trilinos_target_branch_value"
        os.environ["TRILINOS_TARGET_REPO"]    = "trilinos_target_repo_value"
        os.environ["TRILINOS_TARGET_SHA"]     = "trilinos_target_sha_value"
        os.environ["PULLREQUESTNUM"]          = "pullrequestnum_value"
        os.environ["PULLREQUEST_CDASH_TRACK"] = "pullrequest_cdash_track_value"
        os.environ["FORCE_CLEAN"]             = "force_clean_value"

        self._envvar_helper = jenkinsenv.TrilinosJenkinsEnv()


    def tearDown(self):
        del os.environ["TRILINOS_SOURCE_REPO"]
        del os.environ["TRILINOS_SOURCE_SHA"]
        del os.environ["TRILINOS_TARGET_BRANCH"]
        del os.environ["TRILINOS_TARGET_REPO"]
        del os.environ["TRILINOS_TARGET_SHA"]
        del os.environ["PULLREQUESTNUM"]
        del os.environ["PULLREQUEST_CDASH_TRACK"]
        del os.environ["FORCE_CLEAN"]



    def test_trilinos_source_repo(self):
        value = self._envvar_helper.trilinos_source_repo
        self.assertEqual(value, "trilinos_source_repo_value")


    def test_trilinos_source_sha(self):
        value = self._envvar_helper.trilinos_source_sha
        self.assertEqual(value, "trilinos_source_sha_value")


    def test_trilinos_target_branch(self):
        value = self._envvar_helper.trilinos_target_branch
        self.assertEqual(value, "trilinos_target_branch_value")


    def test_trilinos_target_repo(self):
        value = self._envvar_helper.trilinos_target_repo
        self.assertEqual(value, "trilinos_target_repo_value")


    def test_trilinos_target_sha(self):
        value = self._envvar_helper.trilinos_target_sha
        self.assertEqual(value, "trilinos_target_sha_value")


    def test_pullrequestnum(self):
        value = self._envvar_helper.pullrequestnum
        self.assertEqual(value, "pullrequestnum_value")


    def test_pullrequest_cdash_track(self):
        value = self._envvar_helper.pullrequest_cdash_track
        self.assertEqual(value, "pullrequest_cdash_track_value")


    def test_force_clean(self):
        value = self._envvar_helper.force_clean
        self.assertEqual(value, "force_clean_value")
