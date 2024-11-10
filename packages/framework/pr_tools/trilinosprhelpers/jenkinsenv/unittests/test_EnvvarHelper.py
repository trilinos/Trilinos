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
class EnvvarHelperTest(TestCase):
    """
    Test EnvvarHelper
    """
    def setUp(self):
        os.environ["__TEST_ENVVAR__"] = "__TEST_VALUE__"

        self._bogus_envvar = "__ASD75FO3QE134WNR2BJI12POSMBNEIDLA__"
        self._envvar_helper = jenkinsenv.EnvvarHelper()


    def tearDown(self):
        del os.environ["__TEST_ENVVAR__"]


    def test_envvar_found(self):
        #self._envvar_helper.get_envvar_str(self, envvar_name, error_if_missing=False)
        value = self._envvar_helper.get_envvar_str("__TEST_ENVVAR__")
        self.assertEqual(value, "__TEST_VALUE__")


    def test_envvar_notfound_without_error(self):
        value = self._envvar_helper.get_envvar_str(self._bogus_envvar)
        self.assertEqual(value, None)


    def test_envvar_notfound_with_error(self):
        with self.assertRaises(KeyError):
            self._envvar_helper.get_envvar_str(self._bogus_envvar, error_if_missing=True)


    def test_get_or_create_if_missing(self):
        value = self._envvar_helper.get_or_create_if_missing("__TEST_ENVVAR__", default_value="__NEW_VALUE__")
        self.assertEqual(value, "__TEST_VALUE__")

        value = self._envvar_helper.get_or_create_if_missing("__TEST_MISSING_ENVVAR__", default_value="__NEW_VALUE__")
        self.assertEqual(value, "__NEW_VALUE__")

        self.assertIn("__TEST_MISSING_ENVVAR__", os.environ)
        self.assertEqual(os.environ["__TEST_MISSING_ENVVAR__"], "__NEW_VALUE__")


