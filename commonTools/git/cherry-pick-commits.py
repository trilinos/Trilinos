#!/usr/bin/env python

import sys
import os
import traceback
import re

scriptsDir = os.path.abspath(os.path.dirname(sys.argv[0]))+"/../../cmake/python"
sys.path.insert(0, scriptsDir)

from GeneralScriptSupport import *


#
# Read in the command-line arguments
#

usageHelp = r"""cherry-pick-commits.py [OPTIONS]

Example:

eg log --oneline --topo-order --reverse \
  trilinos-release-10-0-branch ^trilinos-release-10-0-branch-init \
  -- packages/zoltan \
  | ./commonTools/git/cherry-pick-commits.py 2>&1 \
  | tee cherry-pick-commits.log

ToDo: Finish documentation!

"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--back-out-bad-merges", dest="backOutBadMerges", action="store_true",
  help="Back out cherry-picks with bad merges [default]" )
clp.add_option(
  "--die-on-bad-merge", dest="backOutBadMerges", action="store_false",
  help="Stop and die on bad merges",
  default=True )

clp.add_option(
  "--no-op", dest="noOp", action="store_true", default=False,
  help="Enable automatic check for the right versions of eg and git. [default]" )

# NOTE: Above, in the pairs of boolean options, the *last* add_option(...) 
# takes effect!  That is why the commands are ordered the way they are!

(options, args) = clp.parse_args()


#
# Echo the command-line
#

print ""
print "**************************************************************************"
print "Script: cherry-pick-commits.py \\"
if options.backOutBadMerges:
  print "  --back-out-bad-merges \\"
else:
  print "  --die-on-bad-merge \\"
if options.noOp:
  print "  --no-op \\"
print ""

#
# Helper functions
#

def readyToCherryPack():
  egStatusOutput = getCmndOutput("git status | grep 'nothing to commit'", throwOnError=False)
  #print "\n"+egStatusOutput
  if reCleanWorkingDir.match(egStatusOutput):
    return True
  return False
    

#
# Execute
#

reCleanWorkingDir = re.compile("""nothing to commit""")

oneLineOutputList = sys.stdin.readlines()

numCherryPicksToAttempt = len(oneLineOutputList)

print "Attempting "+str(numCherryPicksToAttempt)+" cherry picks:"
print "--------------------------------------------------------"
print "".join(oneLineOutputList)


for i in range(len(oneLineOutputList)):

  line = oneLineOutputList[i].strip()

  print "\n--------------------------------------------------------"
  print "\nCherry-pick #"+str(i+1)+"/"+str(numCherryPicksToAttempt)+": "+line

  lineList = line.split(' ')
  sha1 = lineList[0]
  msg = (' '.join(lineList[1:])).strip()
  #print "sha1 = "+sha1
  #print "msg = '"+msg+"'"

  #
  print "\nA) Check the status of the repository:\n"
  #

  if readyToCherryPack():
    print "\nWorking directory is clean and ready to apply cherry-pick!"
  else:
    print "\nThe working directory is *not* clean, aborting cherry-picks!"
    sys.exit(1)

  #
  print "\nB) Do the cherry-pick:\n"
  #

  cmnd = "git cherry-pick -x "+sha1

  if options.noOp:
    print "Would perform: "+cmnd
  else:
    rtn = echoRunSysCmnd(cmnd, throwExcept=False)
    if rtn == 0:
      None # The cherry pick was successful!
    elif rtn == 1:
      None # The cherry pick was not performed because the directory is clean!
    else:
      raise Exception("Error, '"+cmnd+"' return "+str(rtn)+"!" \
        "  Aborting remaining cherry picks!")

  #
  print "\nD) Back out the cherry-pick if it failed:\n"
  #
  
  if not readyToCherryPack():
    if options.backOutBadMerges:
      print "\nBacking out last bad cherry-pick on request ...\n"
      echoRunSysCmnd("git reset --hard HEAD")
    else:
      print "\nThe last cherry-pick failed, aborting cherry-picks!"
      sys.exit(2)
