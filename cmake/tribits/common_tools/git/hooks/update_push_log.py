#!/usr/bin/env python

import commands
import os
import re
import sys

if len(sys.argv) != 3:
  raise SystemExit("Syntax:\n  %s oldrev newrev" % sys.argv[0])

oldrev=sys.argv[1]
newrev=sys.argv[2]

dateStr = commands.getoutput("date")
#print "dateStr = '"+dateStr+"'"

commitsStr = commands.getoutput(
  "git log --pretty=format:'    %h \"%s\" <%ae> [%ad] (%ar)' --topo-order "+newrev+" ^"+oldrev)
#print "commitsStr:\n"+commitsStr

pushLogEntry = \
  "\nPush on "+dateStr+"\n\n"+ \
  commitsStr + "\n"

#print "pushLogEntry:\n", pushLogEntry
               
logFileName = "push.log"

existingFileStr = ""
if os.path.exists(logFileName):
  existingFileStr = open(logFileName, 'r').read()

open(logFileName, 'w').write(pushLogEntry+"\n"+existingFileStr)
