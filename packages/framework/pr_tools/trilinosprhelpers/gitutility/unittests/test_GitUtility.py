#!/usr/bin/env python3
# -*- coding: utf-8; mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""
"""
from __future__ import print_function

try:                                        # pragma: no cover
    import builtins                         # pragma: no cover
except ImportError:                         # pragma: no cover
    import __builtin__ as builtins          # pragma: no cover

import sys

if sys.version_info >= (3,0):               # pragma: no cover
    from io import StringIO                 # pragma: no cover
else:                                       # pragma: no cover
    from io import BytesIO as StringIO      # pragma: no cover

#import unittest
from unittest import TestCase

try:                                        # pragma: no cover
    import unittest.mock as mock            # pragma: no cover
    from unittest.mock import patch
except:                                     # pragma: no cover
    import mock                             # pragma: no cover
    from mock import patch

#from mock import Mock
#from mock import mock_open
#from mock import MagicMock


import trilinosprhelpers.gitutility as gitutility



#==============================================================================
#
#                         M O C K   H E L P E R S
#
#==============================================================================

def mock_subprocess_get_git_version_str(*args):
    output = None
    vstr   = "2.27.0"
    try:
        output = bytes(vstr, "utf-8")
    except:                                 # pragma: no cover
        output = bytes(vstr)                # pragma: no cover
    return output



#==============================================================================
#
#                                T E S T S
#
#==============================================================================
class GitUtilityTest(TestCase):
    """
    """
    def setUp(self):
        pass


    def test_version_str(self):
        print("")
        gitutil = gitutility.GitUtility()
        with patch('subprocess.check_output', side_effect=mock_subprocess_get_git_version_str) as m:
            version_str = gitutil.version_str
        m.assert_called_once_with(["git", "--version"])
        print("version_str = {}".format(version_str))


    def test_version(self):
        print("")
        gitutil = gitutility.GitUtility()
        with patch('subprocess.check_output', side_effect=mock_subprocess_get_git_version_str) as m:
            version = gitutil.version
        m.assert_called_once_with(["git", "--version"])

        self.assertIsInstance(version, dict)
        self.assertEqual(2, version['major'])
        self.assertEqual(27, version['minor'])
        self.assertEqual(0, version['patch'])

        print("version = {}".format(version))


    def test_check_minimum_version(self):
        print("")
        gitutil = gitutility.GitUtility()
        with patch('subprocess.check_output', side_effect=mock_subprocess_get_git_version_str) as m:
            version = gitutil.version
        m.assert_called_once_with(["git", "--version"])

        ret = gitutil.check_minimum_version(1)
        self.assertEqual(0, ret)

        ret = gitutil.check_minimum_version(1, 20)
        self.assertEqual(0, ret)

        ret = gitutil.check_minimum_version(1, 30)
        self.assertEqual(0, ret)

        ret = gitutil.check_minimum_version(2)
        self.assertEqual(0, ret)

        ret = gitutil.check_minimum_version(2, 20)
        self.assertEqual(0, ret)

        with self.assertRaises( SystemExit ) as m:
            gitutil.check_minimum_version(2, 30)

        with self.assertRaises( SystemExit ) as m:
            gitutil.check_minimum_version(3)

        with self.assertRaises( SystemExit ) as m:
            gitutil.check_minimum_version(3, 20)

        with self.assertRaises( SystemExit ) as m:
            gitutil.check_minimum_version(3, 30)

        with self.assertRaises( TypeError ) as m:
            gitutil.check_minimum_version("3")

        with self.assertRaises( TypeError ) as m:
            gitutil.check_minimum_version(2, "27")


    def test_pretty_print(self):
        print("")
        gitutil = gitutility.GitUtility()
        with patch('subprocess.check_output', side_effect=mock_subprocess_get_git_version_str) as m:
            version = gitutil.version
        m.assert_called_once_with(["git", "--version"])

        with patch('sys.stdout', new = StringIO()) as m_out:
            gitutil.pretty_print()

        self.assertEqual(m_out.getvalue().strip(), "Git Version Detected: 2.27.0")





