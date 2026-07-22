#!/usr/bin/env python3
# -*- coding: utf-8; mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""
"""
from __future__ import print_function

try:                                 # pragma: no cover
    import builtins                  # pragma: no cover
except ImportError:                  # pragma: no cover
    import __builtin__ as builtins   # pragma: no cover

import sys

#import unittest
from unittest import TestCase

try:                                    # pragma: no cover
    import unittest.mock as mock        # pragma: no cover
except:                                 # pragma: no cover
    import mock                         # pragma: no cover

#from mock import Mock
#from mock import mock_open
#from mock import MagicMock

from textwrap import dedent

import trilinosprhelpers.sysinfo as sysinfo


# TODO: It would be useful to learn how to run three scenerios:
#       1) psutil exists
#       2) psutil does not exist but /proc/meminfo does.
#       3) psutil does not exist and /proc/meminfo also does not exist.
#   MagicMock should be able to do something with this but I'm not quite
#   sure how to 'intercept' the `import psutil` command in the SysInfo.py file.

#==============================================================================
#
#                         M O C K   H E L P E R S
#
#==============================================================================

def mock_psutil_vm_total():
    return 67108864 * 1024


def mock_multiprocessing_cpu_count():
    return 20


#==============================================================================
#
#                                T E S T S
#
#==============================================================================

class SysInfoTest(TestCase):
    """
    """
    def setUp(self):
        self.maxDiff = None

        # for now we just need MemTotal. Let's set it to a 64 GB system w/ 56 GB free
        self._meminfo = dedent("""\
                MemTotal:       67108864 kB
                MemFree:        58720256 kB
                """)


    def _builtins(self):
        output = "builtins"             # pragma: no cover
        if sys.version_info < (3,0):    # pragma: no cover
            output = "__builtin__"      # pragma: no cover
        return output

    def test_SysInfo_have_psutil_error(self):
        """
        Test that a TypeError is raised if we try to assign something other than a
        bool to SysInfo.have_psutil
        """
        print("")
        si = sysinfo.SysInfo()
        with self.assertRaises(TypeError):
            si.have_psutil = "XXX"


    def test_GetMem_meminfo(self):
        """
        """
        print("")
        si = sysinfo.SysInfo()

        # Force SysInfo to use the meminfo method.
        si.have_psutil = False

        if not si.have_psutil:
            mo = mock.mock_open(read_data=self._meminfo)

            builtins_open = "{}.open".format(self._builtins())

            meminfo = None
            with mock.patch(builtins_open, mo) as m:
                meminfo = si.meminfo

            print("SysInfo.meminfo = {}".format(meminfo))
            print("- mem_kb = {}".format(meminfo['mem_kb']))
            print("- mem_gb = {}".format(meminfo['mem_gb']))

            m.assert_called_once()
            self.assertEqual( 67108864, int(meminfo['mem_kb']) )
            self.assertEqual(       64, int(meminfo['mem_gb']) )


    def test_GetMem_psutil(self):
        """
        """
        print("")
        si = sysinfo.SysInfo()

        # Force SysInfo to use the psutil method (we mock this so even if psutil isn't
        # available this should work)
        si.have_psutil = True

        if si.have_psutil:
            meminfo = None
            with mock.patch.object(si, "_get_psutil_vm_total", side_effect=mock_psutil_vm_total) as m:
                meminfo = si.meminfo

            print("SysInfo.meminfo = {}".format(meminfo))
            print("- mem_kb = {}".format(meminfo['mem_kb']))
            print("- mem_gb = {}".format(meminfo['mem_gb']))

            m.assert_called_once()
            self.assertEqual( 67108864, int(meminfo['mem_kb']) )
            self.assertEqual(       64, int(meminfo['mem_gb']) )


    def test_compute_num_usable_cores(self):
        # compute_num_usable_cores(self, req_mem_gb_per_core=3.0, max_cores_allowed=32):
        si = sysinfo.SysInfo()
        si.have_psutil = True

        if si.have_psutil:
            meminfo = None
            with mock.patch.object(si, "_get_psutil_vm_total", side_effect=mock_psutil_vm_total) as m:
                meminfo = si.meminfo

                print("SysInfo.meminfo = {}".format(meminfo))
                print("- mem_kb = {}".format(meminfo['mem_kb']))
                print("- mem_gb = {}".format(meminfo['mem_gb']))

                # For all:
                # sys memory: 64 GB  (mocked)
                # sys cpus  : 20     (mocked)

                # Test Case:
                # req memory: 6 GB / Core
                # max cores : 3
                # expected  : 3
                with mock.patch("multiprocessing.cpu_count", side_effect=mock_multiprocessing_cpu_count) as m:
                    n = si.compute_num_usable_cores(6, 3)
                    print("- n = {}".format(n))
                    self.assertEqual(3, n)
                    m.assert_called_once()

                # Test Case:
                # req memory: 1 GB / Core
                # max cores : 32
                # expected  : 20
                with mock.patch("multiprocessing.cpu_count", side_effect=mock_multiprocessing_cpu_count) as m:
                    n = si.compute_num_usable_cores(1, 32)
                    print("- n = {}".format(n))
                    self.assertEqual(20, n)
                    m.assert_called_once()

                # Test Case:
                # req memory: 5 GB / Core
                # max cores : 32
                # expected  : 12
                with mock.patch("multiprocessing.cpu_count", side_effect=mock_multiprocessing_cpu_count) as m:
                    n = si.compute_num_usable_cores(5, 32)
                    print("- n = {}".format(n))
                    self.assertEqual(12, n)
                    m.assert_called_once()

                # Test Case:
                # req memory: 5 GB / Core
                # max cores : 0
                # expected  : 1
                with mock.patch("multiprocessing.cpu_count", side_effect=mock_multiprocessing_cpu_count) as m:
                    n = si.compute_num_usable_cores(5, 0)
                    print("- n = {}".format(n))
                    self.assertEqual(1, n)
                    m.assert_called_once()



