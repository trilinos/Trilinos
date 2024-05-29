#!/usr/bin/env python3
# -*- coding: utf-8; mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""
"""
from __future__ import print_function

import os

from unittest import TestCase

try:                                    # pragma: no cover
    import unittest.mock as mock        # pragma: no cover
    from unittest.mock import patch
except:                                 # pragma: no cover
    import mock                         # pragma: no cover
    from mock import patch

import trilinosprhelpers.sysinfo as sysinfo


#==============================================================================
#
#                         M O C K   H E L P E R S
#
#==============================================================================

def mock_nvidia_smi():
    return ["GPU 0: Tesla V100S-PCIE-32GB (UUID: GPU-somehash1)",
            "GPU 1: Tesla V100S-PCIE-32GB (UUID: GPU-somehash2)",
            "GPU 2: Tesla V100S-PCIE-32GB (UUID: GPU-somehash3)",
            "GPU 3: Tesla V100S-PCIE-32GB (UUID: GPU-somehash4)"]

def mock_which(thing_to_find):
    return os.path.join(os.getcwd(), thing_to_find)


#==============================================================================
#
#                                T E S T S
#
#==============================================================================

class GpuUtilsTest(TestCase):
    """
    Tests for gpu_utils.
    """
    def setUp(self):
        self.maxDiff = None

    def test_list_nvidia_gpus(self):
        """
        Test that sane output from nvidia-smi yields a sane list of gpu indices.
        """
        print("")
        with patch("trilinosprhelpers.sysinfo.gpu_utils._nvidia_smi", side_effect=mock_nvidia_smi):
            ret = sysinfo.gpu_utils.list_nvidia_gpus()
        self.assertEqual(["0", "1", "2", "3"], ret)

    def test_has_nvidia_gpus(self):
        """
        Test that sane output from nvidia-smi yields positive for system possessing NVidia GPUs.
        """
        print("")
        with patch("trilinosprhelpers.sysinfo.gpu_utils._nvidia_smi", side_effect=mock_nvidia_smi):
            ret = sysinfo.gpu_utils.has_nvidia_gpus()
        self.assertTrue(ret)

    def test_nvidia_smi_output_without_smi(self):
        """
        Test that without nvidia-smi available the smi interface returns an empty list of output.
        """
        print("")
        with patch("trilinosprhelpers.sysinfo.gpu_utils.which", return_value=None):
            ret = sysinfo.gpu_utils._nvidia_smi()
        self.assertEqual([], ret)
