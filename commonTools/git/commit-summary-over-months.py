#!/usr/bin/env python

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

usageHelp = r"""commit-summary-over-months.py [OPTION]

This script collects statistics for the number of commit by the given authors
over a number of months and produces a summary table sutable for import and
ploting.
"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--start-month-ago", dest="startMonthAgo", type="int",
  default=0
  );

clp.add_option(
  "--end-month-ago", dest="endMonthAgo", type="int",
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
  help="Show the actual commits in the output.",
  default=False )

clp.add_option(
  "--debug", dest="debug", action="store_true",
  help="Show debug output.",
  default=False )


(options, args) = clp.parse_args()


#
# Executable code
#

authorsList = options.authors.split(",")
#print "authorsList =", authorsList

excludeFilesList = options.excludeFiles.split(",")
#print "excludeFilesList =", excludeFilesList

# Gather the number of commits for each author over the given range of months

authorsNumCommits = []
authorsNumLinesAdded = []

for month_i in range(options.startMonthAgo, options.endMonthAgo):

  #print "month_i =", month_i

  authorsNumCommitsThisMonth = []
  authorsNumLinesAddedThisMonth = []

  for author in authorsList:

    #print "author =", author

    cmnd = "eg log --pretty=format:'%n%h \"%s\" <%ae> [%ad] (%ar)' -w -C" + \
      " --after=\""+str(month_i+1)+" months ago\"" + \
      " --numstat" + \
      " --before=\""+str(month_i)+" months ago\"" + \
      " --author="+author
    if options.files:
      cmnd += " -- " + " ".join(options.files.split(","))

    if options.showLogCmnd:
      print "\n"+cmnd

    authorCommitsThisMonth = getCmndOutput(cmnd, True)

    if options.showCommits:
      print authorCommitsThisMonth

    authorNumCommitsThisMonth = 0
    authorNumLinesAddedThisMonth = 0;

    if authorCommitsThisMonth != "":

      authorNumCommitsThisMonth = 0
      authorNumLinesAddedThisMonth = 0

      authorCommitsThisMonthLines = authorCommitsThisMonth.split("\n")

      # Loop through the commits

      line_i = 0

      while line_i < len(authorCommitsThisMonthLines):

        if options.debug:
          print "\nTOP: authorCommitsThisMonthLines[line_i] = '"+authorCommitsThisMonthLines[line_i]+"'"

        # Skip empty lines between commits
        if authorCommitsThisMonthLines[line_i].strip() == "":
          line_i += 1
          continue

        # Process current commit (i.e. ThisCommit)

        authorNumLinesAddedThisCommit = 0

        commitMsg = authorCommitsThisMonthLines[line_i]
        if options.debug:
          print "\ncommitMsg =", commitMsg
        line_i += 1

        commitSHA1 = commitMsg.split(" ")[0].strip()
        if options.debug:
          print "commitSHA1 =", commitSHA1

        # Process modified files (until a newline is found)

        if options.debug:
          print "\nLoop through changed files:\n"

        while line_i < len(authorCommitsThisMonthLines) \
          and authorCommitsThisMonthLines[line_i].strip() != "" \
          :

          if options.debug: print ""

          if options.debug:
            print "\nauthorCommitsThisMonthLines[line_i] = '"+authorCommitsThisMonthLines[line_i]+"'"

          fileStatArray = authorCommitsThisMonthLines[line_i].split("\t")
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
                if options.showCommits or options.showLogCmnd:
                  print "NOTE: Excluding the file "+fileName+" in commit "+commitSHA1+" because it matches "+excludeFile+"!"
                excludeThisFile = True
                break

          if not excludeThisFile:
            authorNumLinesAddedThisCommit += numLinesAdded
            
          line_i += 1

        if authorNumLinesAddedThisCommit > 0:
          authorNumCommitsThisMonth += 1
          authorNumLinesAddedThisMonth += authorNumLinesAddedThisCommit
        else:
          if options.showLogCmnd or options.showCommits:
            print "NOTE: Excluding commit "+commitSHA1+" because all files matched exclude pattern(s) or is merge commit!"

    if options.showLogCmnd:
      print "authorNumCommitsThisMonth =", authorNumCommitsThisMonth
      print "authorNumLinesAddedThisMonth =", authorNumLinesAddedThisMonth

    authorsNumCommitsThisMonth.append(authorNumCommitsThisMonth)
    authorsNumLinesAddedThisMonth.append(authorNumLinesAddedThisMonth)

  authorsNumCommits.append(authorsNumCommitsThisMonth)
  authorsNumLinesAdded.append(authorsNumLinesAddedThisMonth)


#
# Output the statistics
#

maxAuthorName = 7
for author in authorsList:
  maxAuthorName = max(maxAuthorName, len(author))

paddingWidth = 2
monthWidth = 7
authorCountWidth = maxAuthorName

def printTableSeprators():
  monthSepStr = "-"*monthWidth
  authorSepStr = "-"*authorCountWidth
  print monthSepStr.ljust(monthWidth),
  for author_i in range(0,len(authorsList)):
    print "".rjust(paddingWidth) + authorSepStr.rjust(authorCountWidth),
  print ""
  

def printTable(tableData):
  # Top line
  print "month".rjust(monthWidth),
  for author in authorsList:
    print "".rjust(paddingWidth) + author.ljust(authorCountWidth),
  print ""
  printTableSeprators()
  # Month lines
  totals = []
  for author_i in range(0,len(authorsList)):
    totals.append(0)
  for month_i in range(options.startMonthAgo, options.endMonthAgo):
    print repr(month_i).rjust(monthWidth),
    month_idx = month_i - options.startMonthAgo
    for author_i in range(0,len(authorsList)):
      #print "month_i="+str(month_i)+", month_idx="+str(month_idx)+", author_i="+str(author_i)
      amount = tableData[month_idx][author_i]
      print "".rjust(paddingWidth) + repr(amount).rjust(authorCountWidth),
      totals[author_i] += amount
    print ""
  # Totals
  printTableSeprators()
  print "totals".rjust(monthWidth),
  for author_i in range(0,len(authorsList)):
    print "".rjust(paddingWidth) + repr(totals[author_i]).rjust(authorCountWidth),
  print ""



#
# Print a tables of commits
#

print "\n\n"+getCmndOutput("date").strip()

print "\n\nNumber of commits by each author over several months ago:\n"
printTable(authorsNumCommits)

print "\n\nNumber of lines added by each author over several months ago\n"
printTable(authorsNumLinesAdded)
