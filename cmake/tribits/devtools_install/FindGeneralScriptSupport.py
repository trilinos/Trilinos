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

pythonUtilsDir = os.path.abspath(
  os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "python_utils" )
  )
sys.path = [pythonUtilsDir] + sys.path

from GeneralScriptSupport import *
