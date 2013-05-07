#!/usr/bin/env python
# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER


import sys
import os
import traceback
import re

scriptsDir = os.path.abspath(os.path.dirname(sys.argv[0]))+"/../../cmake/python"
sys.path.insert(0, scriptsDir)

from GeneralScriptSupport import *

#
# Commandline options
#

usageHelp = r"""commit-summary-over-periods.py [OPTION]

This script collects statistics for the number of commit by the given authors
over a number of time periods (default months) and produces a summary table
sutable for import and plotting.
"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--period-unit", dest="periodUnit", type="string",
  default="month",
  help="The unit of time period used in --before and --after log arguments"
  );

clp.add_option(
  "--start-ago", dest="startAgo", type="int",
  default=0
  );

clp.add_option(
  "--end-ago", dest="endAgo", type="int",
  default=12
  );

clp.add_option(
  "--authors", dest="authors", type="string",
  default="",
  help="List (comma separated) of authors to get summaries for"
  );

clp.add_option(
  "--files", dest="files", type="string",
  default="",
  help="List (comma sepraged) of files and directires to include in the search"
  );

clp.add_option(
  "--exclude-files", dest="excludeFiles", type="string",
  default="",
  help="List (comma sepraged) of files and directires to exclude in the statistics"
  );

clp.add_option(
  "--show-log-cmnd", dest="showLogCmnd", action="store_true",
  help="Show the log command for each invocation.",
  default=False )

clp.add_option(
  "--show-commits", dest="showCommits", action="store_true",
  help="Show the actual commits in the output (implies --show-log-cmnd).",
  default=False )

clp.add_option(
  "--debug", dest="debug", action="store_true",
  help="Show debug output.",
  default=False )


(options, args) = clp.parse_args()


if options.showCommits:
  options.showLogCmnd = True


#
# Executable code
#

authorsList = options.authors.split(",")
#print "authorsList =", authorsList

excludeFilesList = options.excludeFiles.split(",")
#print "excludeFilesList =", excludeFilesList

# Gather the number of commits for each author over the given range of periods

authorsNumCommits = []
authorsNumLinesAdded = []

for period_i in range(options.startAgo, options.endAgo):

  #print "period_i =", period_i

  authorsNumCommitsThisPeriod = []
  authorsNumLinesAddedThisPeriod = []

  for author in authorsList:

    #print "author =", author

    cmnd = "eg log --pretty=format:'%n%h \"%s\" <%ae> [%ad] (%ar)' -w -C" + \
      " --before=\""+str(period_i)+" "+options.periodUnit+" ago\"" + \
      " --after=\""+str(period_i+1)+" "+options.periodUnit+" ago\"" + \
      " --numstat" + \
      " --author="+author
    if options.files:
      cmnd += " -- " + " ".join(options.files.split(","))

    if options.showLogCmnd:
      print "\n"+cmnd

    authorCommitsThisPeriod = getCmndOutput(cmnd, True)

    if options.showCommits:
      print authorCommitsThisPeriod

    authorNumCommitsThisPeriod = 0
    authorNumLinesAddedThisPeriod = 0;

    if authorCommitsThisPeriod != "":

      authorNumCommitsThisPeriod = 0
      authorNumLinesAddedThisPeriod = 0

      authorCommitsThisPeriodLines = authorCommitsThisPeriod.split("\n")

      # Loop through the commits

      line_i = 0

      while line_i < len(authorCommitsThisPeriodLines):

        if options.debug:
          print "\nTOP: authorCommitsThisPeriodLines[line_i] = '"+authorCommitsThisPeriodLines[line_i]+"'"

        # Skip empty lines between commits
        if authorCommitsThisPeriodLines[line_i].strip() == "":
          line_i += 1
          continue

        # Process current commit (i.e. ThisCommit)

        authorNumLinesAddedThisCommit = 0

        commitMsg = authorCommitsThisPeriodLines[line_i]
        if options.debug:
          print "\ncommitMsg =", commitMsg
        line_i += 1

        commitSHA1 = commitMsg.split(" ")[0].strip()
        if options.debug:
          print "commitSHA1 =", commitSHA1

        # Process modified files (until a newline is found)

        if options.debug:
          print "\nLoop through changed files:\n"

        while line_i < len(authorCommitsThisPeriodLines) \
          and authorCommitsThisPeriodLines[line_i].strip() != "" \
          :

          if options.debug: print ""

          if options.debug:
            print "\nauthorCommitsThisPeriodLines[line_i] = '"+authorCommitsThisPeriodLines[line_i]+"'"

          fileStatArray = authorCommitsThisPeriodLines[line_i].split("\t")
          if options.debug: print "fileStatArray =", fileStatArray

          fileName = fileStatArray[2].strip()
          if options.debug: print "fileName =", fileName

          if fileStatArray[0] == "-":
            numLinesAdded = 0
          else:
            numLinesAdded = int(fileStatArray[0])
          if options.debug: print "numLinesAdded =", numLinesAdded

          excludeThisFile = False

          if len(excludeFilesList):
            for excludeFile in excludeFilesList:
              if excludeFile and re.match(excludeFile, fileName):
                if options.showLogCmnd:
                  print "NOTE: Excluding the file "+fileName+" in commit "+commitSHA1+" because it matches "+excludeFile+"!"
                excludeThisFile = True
                break

          if not excludeThisFile:
            authorNumLinesAddedThisCommit += numLinesAdded
            
          line_i += 1

        if authorNumLinesAddedThisCommit > 0:
          authorNumCommitsThisPeriod += 1
          authorNumLinesAddedThisPeriod += authorNumLinesAddedThisCommit
        else:
          if options.showLogCmnd:
            print "NOTE: Excluding commit "+commitSHA1+" because all files matched exclude pattern(s) or is merge commit!"

    if options.showLogCmnd:
      print "authorNumCommitsThisPeriod =", authorNumCommitsThisPeriod
      print "authorNumLinesAddedThisPeriod =", authorNumLinesAddedThisPeriod

    authorsNumCommitsThisPeriod.append(authorNumCommitsThisPeriod)
    authorsNumLinesAddedThisPeriod.append(authorNumLinesAddedThisPeriod)

  authorsNumCommits.append(authorsNumCommitsThisPeriod)
  authorsNumLinesAdded.append(authorsNumLinesAddedThisPeriod)


#
# Output the statistics
#

maxAuthorName = 7
for author in authorsList:
  maxAuthorName = max(maxAuthorName, len(author))

paddingWidth = 2
periodWidth = 7
authorCountWidth = maxAuthorName

def printTableSeprators():
  periodSepStr = "-"*periodWidth
  authorSepStr = "-"*authorCountWidth
  print periodSepStr.ljust(periodWidth),
  for author_i in range(0,len(authorsList)):
    print "".rjust(paddingWidth) + authorSepStr.rjust(authorCountWidth),
  print ""
  

def printTable(tableData):
  # Top line
  print options.periodUnit.rjust(periodWidth),
  for author in authorsList:
    print "".rjust(paddingWidth) + author.ljust(authorCountWidth),
  print ""
  printTableSeprators()
  # Period lines
  totals = []
  for author_i in range(0,len(authorsList)):
    totals.append(0)
  for period_i in range(options.startAgo, options.endAgo):
    print repr(period_i).rjust(periodWidth),
    period_idx = period_i - options.startAgo
    for author_i in range(0,len(authorsList)):
      #print "period_i="+str(period_i)+", period_idx="+str(period_idx)+", author_i="+str(author_i)
      amount = tableData[period_idx][author_i]
      print "".rjust(paddingWidth) + repr(amount).rjust(authorCountWidth),
      totals[author_i] += amount
    print ""
  # Totals
  printTableSeprators()
  print "totals".rjust(periodWidth),
  for author_i in range(0,len(authorsList)):
    print "".rjust(paddingWidth) + repr(totals[author_i]).rjust(authorCountWidth),
  print ""



#
# Print a tables of commits
#

print "\n\n"+getCmndOutput("date").strip()

print "\n\nNumber of commits by each author over several "+options.periodUnit+"s ago:\n"
printTable(authorsNumCommits)

print "\n\nNumber of lines added by each author over several "+options.periodUnit+"s ago:\n"
printTable(authorsNumLinesAdded)
