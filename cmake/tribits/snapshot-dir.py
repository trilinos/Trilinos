#!/usr/bin/env python

import sys
import os


#
# Get the location of the scripts directory whether from a sym link or the
# actual
#

checkinTestFilePath = os.path.abspath(sys.argv[0])
checkinTestFileRealPath = os.path.realpath(checkinTestFilePath)
scriptsDir = os.path.dirname(checkinTestFileRealPath)+"/python"
#print "scriptsDir='"+scriptsDir+"'"
sys.path.insert(0, scriptsDir)


#
# Import and run
#

import SnapshotDir
success = SnapshotDir.snapshotDirMainDriver(sys.argv[1:])
if success: rtnCode = 0
else: rtnCode = 1
sys.exit(rtnCode)
