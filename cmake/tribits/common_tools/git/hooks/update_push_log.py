#!/usr/bin/env python
#
# Simple python script that updates a push log file

import commands
import os
import re
import sys

if len(sys.argv) != 5:
  raise SystemExit("Syntax:\n  %s oldrev newrev refname pushdate " % sys.argv[0])

oldrev=sys.argv[1]
newrev=sys.argv[2]
refname=sys.argv[3]
pushdate=sys.argv[4]

#print "pushdate = '"+pushdate+"'"

commitsStr = commands.getoutput(
  "git log --pretty=format:'    %h \"%s\" <%ae> [%ad] (%ar)' --topo-order "+newrev+" ^"+oldrev)
#print "commitsStr:\n"+commitsStr

pushLogEntry = \
  "\nPush to "+refname+" on "+pushdate+"\n\n"+ \
  commitsStr + "\n"

#print "pushLogEntry:\n", pushLogEntry
               
logFileName = "push.log"

existingFileStr = ""
if os.path.exists(logFileName):
  existingFileStr = open(logFileName, 'r').read()

open(logFileName, 'w').write(pushLogEntry+"\n"+existingFileStr)
