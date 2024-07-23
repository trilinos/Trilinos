#!/usr/bin/env python

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# Usage: do-something ; mailmsg.py "something has happened"

# Send email to yourself about something that has happened.  This is used as a
# reminder to check on something so that you can do something else until it is
# finished.  It will send to your git-configured email by default and to
# $USER@sandia.gov as a backup.

import sys
import os

scriptsDir = os.path.abspath(os.path.dirname(sys.argv[0]))+"/cmake/python"
sys.path.insert(0, scriptsDir)

from GeneralScriptSupport import *

emailAddress = getCmndOutput("git config --get user.email", True, False)
if not emailAddress:
  emailAddress = os.environ['USER']+"@sandia.gov"

msg = sys.argv[1]
sysname = os.uname()[1]
pwd = os.getcwd()

emailBody = ""\
  +"Message:  "+msg+"\n" \
  +"Machine:  "+sysname+"\n" \
  +"PWD:      "+pwd

cmnd = "echo \""+emailBody+"\" | mailx -s \""+msg+"\" "+emailAddress
print(cmnd)
os.system(cmnd)
