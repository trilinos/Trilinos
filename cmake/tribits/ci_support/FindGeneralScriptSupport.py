# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

import os
import sys


ciSupportDir = os.path.dirname(os.path.abspath(__file__))
tribitsDir = os.path.abspath(os.path.join(ciSupportDir, ".."))
pythonUtilsDir = os.path.join(tribitsDir, "python_utils")
defaultProjectDir = os.path.abspath(os.path.join(tribitsDir, "../.."))

sys.path = [pythonUtilsDir] + sys.path

from GeneralScriptSupport import *
