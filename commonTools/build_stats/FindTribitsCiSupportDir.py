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

thisScriptsDir = os.path.dirname(os.path.abspath(__file__))
trilinosTriBITSDirEnv = os.environ.get('Trilinos_TRIBITS_DIR', None)
if trilinosTriBITSDirEnv:
  tribitsDir = trilinosTriBITSDirEnv
else:
  tribitsDir = os.path.abspath(
    os.path.join(thisScriptsDir, "../../cmake/tribits") )
ciSupportDir = os.path.join(tribitsDir, "ci_support")
pythonUtilsDir = os.path.join(tribitsDir, "python_utils")
#print "ciSupportDir =", ciSupportDir

sys.path = [ciSupportDir, pythonUtilsDir] + sys.path
