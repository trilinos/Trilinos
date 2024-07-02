#!/usr/bin/env python
# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


import sys
import os


#
# Get the location of the scripts directory whether from a sym link or the
# actual
#

checkinTestFilePath = os.path.abspath(sys.argv[0])
checkinTestFileRealPath = os.path.realpath(checkinTestFilePath)
scriptsDir = os.path.dirname(checkinTestFileRealPath)+"/python"
#print("scriptsDir='" + scriptsDir + "'")
sys.path.insert(0, scriptsDir)



#
# Import and run
#

import SnapshotDir

snapshotDirDummyDefaults = os.environ.get("SNAPSHOT_DIR_DUMMY_DEFAULTS","")

defaultOptions = None
if snapshotDirDummyDefaults:
  defaultOptions = SnapshotDir.DefaultOptions()
  defaultOptions.origDir = "<orig-dir>"
  defaultOptions.destDir = "<dest-dir>"

success = SnapshotDir.snapshotDirMainDriver(sys.argv[1:], defaultOptions)

if success:
  rtnCode = 0
else:
  rtnCode = 1

sys.exit(rtnCode)
