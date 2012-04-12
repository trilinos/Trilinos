#!/usr/bin/env python

import sys
import commands

#
# Commandline options
#

pstreeCmnd = "pstree -plAca"

usageHelp = r"""kill-pstree.py PID

Kill a whole tree of process starting with the given root process PID.  This script
uses:

    """ \
+pstreeCmnd+ \
"""

to get the commands and their PIDs to kill the tree.
"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--show-kills", dest="showKills", action="store_true",
  help="Show the processes that are being killed in the order they are killed.",
  default=False )

clp.add_option(
  "--no-op", dest="noOp", action="store_true",
  help="Only show what would be killed, don't actually kill any processes.  Implies --show-kills",
  default=False )

(options, args) = clp.parse_args()

if options.noOp:
  options.showKills = True

if len(args) != 1:
  print "Error, must pass in PID as single argument!"
  sys.exit(1)

pid = args[0]



#
# Helper functions
#

class CmndAndPID:
  def __init__(self, shortCmnd, PID, cmndArgs):
    self.shortCmnd = shortCmnd
    self.PID = PID
    self.cmndArgs = cmndArgs


def getCommandsAndProcIDs(pstreeOut):

  rawSplitCmndAndProcIDs = pstreeOut.split("\n")
  #print "rawSplitCmndAndProcIDs =", rawSplitCmndAndProcIDs

  cmndsAndProcIDs = []

  for rawSplitItem in rawSplitCmndAndProcIDs:

    #print "rawSplitItem =", rawSplitItem

    # A) Strip off the leading formatting characters:   |   `-
    cmndAndProcIDStr = rawSplitItem
    for i in range(0,4):
      cmndAndProcIDStr = cmndAndProcIDStr.lstrip(" ") # Strip off leading spaces
      cmndAndProcIDStr = cmndAndProcIDStr.lstrip("|") # Strip off vertical connector
      cmndAndProcIDStr = cmndAndProcIDStr.lstrip(" ") # Strip off remaining leading spaces
      cmndAndProcIDStr = cmndAndProcIDStr.lstrip("`") # Strip off `
      cmndAndProcIDStr = cmndAndProcIDStr.lstrip("-") # Strip off -
    #print "cmndAndProcIDStr =", cmndAndProcIDStr
    if cmndAndProcIDStr == "":
      continue

    # B) Find the parts of 'cmnd,PID arguments'
    firstComaIdx = cmndAndProcIDStr.find(",")
    firstSpaceIDx = cmndAndProcIDStr.find(" ")
    shortCmnd = cmndAndProcIDStr[0:firstComaIdx]
    #print "shortCmnd =", shortCmnd
    PID = cmndAndProcIDStr[firstComaIdx+1:firstSpaceIDx]
    #print "PID =", PID
    cmndArgs = cmndAndProcIDStr[firstSpaceIDx+1:-1]
    #print "cmndArgs =", cmndArgs
    cmndAndPIDAndArgs = CmndAndPID(shortCmnd, PID, cmndArgs)
    
    # C) Append the item
    cmndsAndProcIDs.append(cmndAndPIDAndArgs)

  return cmndsAndProcIDs
  

#
# Executable code
#

pstreeOut = commands.getoutput("pstree -plAca "+pid)

#print "pstreeOut:\n", pstreeOut

cmndsAndProcIDs = getCommandsAndProcIDs(pstreeOut)

for cmndAndProcID in cmndsAndProcIDs:
  if options.showKills:
     print "killing PID = "+cmndAndProcID.PID+", short command = "+cmndAndProcID.shortCmnd \
       +", args: "+cmndAndProcID.cmndArgs
  if not options.noOp:
     commands.getstatusoutput("kill -s 9 "+cmndAndProcID.PID)
