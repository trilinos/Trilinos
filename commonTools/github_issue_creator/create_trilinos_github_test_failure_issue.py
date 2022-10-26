#!/usr/bin/env python


#
# Implementation code
#


from FindTribitsModuleDirs import *

import CreateIssueTrackerFromCDashQuery as CITFCQ


usageHelp = \
r"""Create text for the body for a new Trilinos GitHub Markdown-formatted Issue
for a set of nonpassing PR builds tests specified in a cdash/queryTests.php query.

Typical usage is:

  create_trilinos_github_test_failure_issue.py \
    -u "<failing-tests-cdash-url>" \
    -s "<summary-line>" \
    -i newGitHubIssueBody.md

The contents of the generated file 'newGitHubIssueBody.md' can then be used to
manually create a new Trilinos GitHub issue where the rest of the details can
be manually filled in.

This tool fills in some details that can be automatically determined just by
analyzing the downloaded list of nonpassing tests.  The rest of the
information must be provided once the text is copied and pasted into a new
GitHub issue.
"""


# The main function
def main():
  issueTrackerCreator = \
    CITFCQ.CreateIssueTrackerFromCDashQueryDriver(
      TrilinosIssueTrackerFormatter(),
      cdashProjectStartTimeUtc="4:00",
      usageHelp=usageHelp,
      issueTrackerUrlTemplate="https://github.com/trilinos/Trilinos/issues/<newissueid>",
      issueTrackerTemplate="#<newissueid>" )
  issueTrackerCreator.runDriver()


# Nonmember function to actually create the body of the new GitHub marddown text
#
# NOTE: This was made a nonmember function to put this to the bottom and not
# obscure the implementing class 'TrilinosIssueTrackerFormatter'
#
def getTrilinosGithubIssueBodyMarkdown(
    itd,  # issueTrackerData (type IssueTrackerData)
    uniqueGenConfigBuildNamesList,
  ):
  issueTrackerText = \
r"""
SUMMARY: """+itd.summaryLine+" "+itd.testingDayStartNonpassingDate+r"""

CC: @trilinos/<package-name\>, @<triage-contact\> (Trilinos <product-area-name\> Triage Contact)

<Checklist>
<???: Add label "type: bug"?>
<???: Add label "impacting: [tests, build, performance]"?>
<???: Add label "pkg: [from <TestName>]>
<???: Add label "PA: [see https://sems-atlassian-son.sandia.gov/confluence/display/TRIL/Triaging+contacts]?>


## Next Action Status

<status-and-or-first-action>


## Description

As shown in [this query]("""+itd.nonpassingTestsUrl+r""") (click "Shown Matching Output" in upper right) the tests:

""" + CITFCQ.getMarkdownListStr(itd.testnameList, '`') + \
r"""
in the unique GenConfig builds:

""" + CITFCQ.getMarkdownListStr(uniqueGenConfigBuildNamesList, '`') + \
r"""
started failing on testing day """+itd.testingDayStartNonpassingDate+r""".

The specific set of CDash builds impacted where:

""" + CITFCQ.getMarkdownListStr(itd.buildnameList, '`') + \
r"""

<Add details about what is failing and what the failures look like.  Make sure to include strings that are easy to match with GitHub Issue searches.\>


## Current Status on CDash

Run the [above query]("""+itd.nonpassingTestsUrl+r""") adjusting the "Begin" and "End" dates to match today any other date range or just click "CURRENT" in the top bar to see results for the current testing day.


## Steps to Reproduce

See:

* https://github.com/trilinos/Trilinos/wiki/Reproducing-PR-Testing-Errors

If you can't figure out what commands to run to reproduce the problem given this documentation, then please post a comment here and we will give you the exact minimal commands.
"""

# NOTE: ABOVE: It is important to keep entire paragraphs on one line.
# Otherwise, GitHub will show the line-breaks and it looks terrible.

  return issueTrackerText

# END FUNCTION: getTrilinosGithubIssueBodyMarkdown()


################################################################################
#
# EVERYTHING BELOW HERE SHOULD NEED TO BE MODIFIED.  IT IS JUST BOILERPLATE
# CODE
#
################################################################################


# Class implementation of callback to fill in the new ATDM Trilinos GitHub
# issue.
#
class TrilinosIssueTrackerFormatter:


  def createFormattedIssueTracker(self, issueTrackerData):
    uniqueGenConfigBuildNamesList = getUniqueGenConfigBuildNamesList(
      issueTrackerData.buildnameList)
    return getTrilinosGithubIssueBodyMarkdown(issueTrackerData,
      uniqueGenConfigBuildNamesList )


def getUniqueGenConfigBuildNamesList(buildnameList):
  # Get sorted list of GenConfig build names stripping of prefix and suffix
  allGenConfigBuildNamesList = []
  for buildname in buildnameList:
    allGenConfigBuildNamesList.append(stripGentConfigBuildName(buildname))
  allGenConfigBuildNamesList.sort()
  # Get unique list of GenConfig build names
  uniqueGenConfigBuildNamesList = []
  lastGenconfigBuild = None
  for genconfigbuild in allGenConfigBuildNamesList:
    if genconfigbuild != lastGenconfigBuild:
      uniqueGenConfigBuildNamesList.append(genconfigbuild)
      lastGenconfigBuild = genconfigbuild
  #
  return uniqueGenConfigBuildNamesList


def stripGentConfigBuildName(buildname):
  #print("buildname = '"+buildname+"'")
  #print("len(buildname) = '"+str(len(buildname))+"'")
  beginIdx = buildname.find("-test-", 0)+6
  #print("beginIdx = '"+str(beginIdx)+"'")
  endIdx = buildname.rfind("-", 0)
  #print("endIdx = '"+str(endIdx)+"'")
  return buildname[beginIdx:endIdx]


#
# Execute main if this is being run as a script
#

if __name__ == '__main__':
  sys.exit(main())
