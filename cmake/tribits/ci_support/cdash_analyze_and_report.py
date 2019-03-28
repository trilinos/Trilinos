#!/usr/bin/env python

# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
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
import pprint
import datetime

from FindGeneralScriptSupport import *
from GeneralScriptSupport import *
import CDashQueryAnalyzeReport as CDQAR
from gitdist import addOptionParserChoiceOption

#
# Help message
#


usageHelp = r"""cdash_analyze_and_report.py [options]

This script takes in CDash URL information and other data as command-line
arguments and then analyzes it to look for missing expected builds, failed
tests,q and various types of failures and then reports the findings as an HTML
file written to disk and/or as HTML-formatted emails sent to one or more email
addresses.

If all of the expected builds are found (and all of them have test results)
and there are no other failures found, then the script returns 0.  Otherwise
the script returns non-zero.  Therefore, this script can be used to drive
automated workflows by examining data on CDash.

ToDo: Finish documentation!
"""


#
# Helper functions
#


def injectCmndLineOptionsInParser(clp, gitoliteRootDefault=""):

  yesterday = (datetime.date.today()+datetime.timedelta(days=-1)).isoformat()

  clp.add_option(
    "--date", dest="date", type="string", default=yesterday,
    help="Date for the testing day <YYYY-MM-DD>."+\
      " [default yesterday '"+yesterday+"']" )

  clp.add_option(
    "--cdash-project-name", dest="cdashProjectName", type="string", default="",
    help="CDash project name (e.g. 'Trilinos'). [REQUIRED] [default = '']" )

  clp.add_option(
    "--build-set-name", dest="buildSetName", type="string", default="",
    help="Name for the set of builds, (e.g. 'Trilinos Nightly Builds)."+\
      "  This used in the email summary line and in the HTML file body"+\
      " to identify the set of builds and tests being examined."+\
      " [REQUIRED] [default = '']" )

  clp.add_option(
    "--cdash-site-url", dest="cdashSiteUrl", type="string", default="",
    help="Base CDash site (e.g. 'https://testing.sandia.gov/cdash')."+\
      " [REQUIRED] [default = '']" )

  clp.add_option(
    "--cdash-builds-filters", dest="cdashBuildsFilters", type="string",
    default="",
    help="Partial URL fragment for index.php making of the filters for"+\
      " the set of builds (e.g. 'filtercount=1&showfilters=1&field1=groupname&compare1=61&value1=ATDM')."+\
      " [REQUIRED] [default = '']" )

  clp.add_option(
    "--cdash-nonpassed-tests-filters", dest="cdashNonpassedTestsFilters", type="string",
    default="",
    help="Partial URL fragment for queryTests.php making of the filters for"+\
      " the set of non-passing tests matching this set of builds (e.g."+\
      " 'filtercombine=and&filtercount=1&showfilters=1&filtercombine=and&field1=groupname&compare1=61&value1=ATDM'"+\
      " [REQUIRED] [default = '']" )

  clp.add_option(
    "--expected-builds-file", dest="expectedBuildsFile", type="string",
    default="",
    help="Path to CSV file that lists the expected builds.  Each of these builds"+\
      " must have unique 'site' and 'buildname' field pairs or an error will be"+\
      " raised and the tool will abort.  [default = '']" )

  clp.add_option(
    "--tests-with-issue-trackers-file", dest="testsWithIssueTrackersFile",
    type="string",  default="",
    help="Path to CSV file that lists tests with issue trackers (and other data)."+\
    "  Each of these tests must have a unique 'site', 'buildName', and 'testname'"+\
    " sets or an error will be raised and the tool will abort.  [default = '']" )

  cdashQueriesCacheDir_default=os.getcwd()

  clp.add_option(
    "--cdash-queries-cache-dir", dest="cdashQueriesCacheDir", type="string",
    default=cdashQueriesCacheDir_default,
    help="Cache CDash query data this directory." \
      +" [Default = '"+cdashQueriesCacheDir_default+"']" )

  clp.add_option(
    "--cdash-base-cache-files-prefix", dest="cdashBaseCacheFilesPrefix", type="string",
    default="",
    help="Prefix given to the base-level cache files outside of the test_history/"+\
      " directory.   This is to allow multiple invocations of this script to share"+\
      " the same base cache directory and share the test_history/ in case there are"+\
      " overrlapping sets of tests where the test history cache could be reused."+\
      " [default is derived from the --build-set-name=<build_set_name> argument where"+\
      " spaces and punctuation in <build_set_name> are replaced with '_']" )

  addOptionParserChoiceOption(
    "--use-cached-cdash-data", "useCachedCDashDataStr",
    ("on", "off"), 1,
    "Use data downloaded from CDash already cached.  Note that this only"+\
    " impacts the reuse of base-level cache files and does not impact the usage"+\
    " of test history cache in the <cacheDir>/test_history/ directory."+\
    "  If a test history file for a given testing day exists under the test_history/"+\
    " directory it is used unconditionally.",
    clp )

  testHistoryDaysDefault= 30

  clp.add_option(
    "--limit-test-history-days", dest="testHistoryDays",
    default=testHistoryDaysDefault, type="int",
    help="Number of days to go back in history for each test."+\
      "  [default = '"+str(testHistoryDaysDefault)+"']" )

  limitTableRows = 10

  clp.add_option(
    "--limit-table-rows", dest="limitTableRows", type="int",
    default=limitTableRows,
    help="Limit to the number of rows displayed in many of"+\
      " the tables.  This impacts tables like 'twoif' and 'twoinr'"+\
      " that could have thousands of entries for some projects."+\
      "   This limits the number of tests for which detailed test history"+\
      " is downloaded from CDash and is therefore important to ensure the"+\
      " tool does not take too long to execute.  However, this does NOT"+\
      " limit the number of rows in many other tables that should be bounded like"+\
      " any of the tables related to the builds or the list of tests with"+\
      " issue trackers.  (The number of those should never be extremely high.)"+\
       "  [default '"+str(limitTableRows)+"']" )

  addOptionParserChoiceOption(
    "--print-details", "printDetailsStr",
    ("on", "off"), 1,
    "Print more info like the CDash URLs for downloaded data and the cache"+\
      " file names.",
    clp )

  clp.add_option(
    "--write-failing-tests-without-issue-trackers-to-file",
    dest="writeFailingTestsWithoutIssueTrackersToFile", type="string", default="",
    help="Write CSV file with a list of tets with issue trackers failed 'twif'." \
    +"  This is to make it easy to add new entires to the file read by" \
    +" the option --tests-with-issue-trackers-file=<file>. [default = '']" )

  clp.add_option(
    "--write-email-to-file", dest="writeEmailToFile", type="string", default="",
    help="Write the body of the HTML email to this file. [default = '']" )

  clp.add_option(
    "--email-from-address=", dest="emailFromAddress", type="string", default="",
    help="Address reported in the sent email. [default '']" )

  clp.add_option(
    "--send-email-to=", dest="sendEmailTo", type="string", default="",
    help="Send email to 'address1, address2, ...'.  [default '']" )


def validateCmndLineOptions(inOptions):
  
  if inOptions.date == "":
    print "Error, can't have empty --date, must pass in --date=YYYY-MM-DD!"
    sys.exit(1)
  else:
    CDQAR.validateAndConvertYYYYMMDD(inOptions.date)

  # ToDo: Assert more of the options to make sure they are correct!


def setExtraCmndLineOptionsAfterParse(inOptions_inout):

  if inOptions_inout.useCachedCDashDataStr == "on":
    setattr(inOptions_inout, 'useCachedCDashData', True)
  else:
    setattr(inOptions_inout, 'useCachedCDashData', False)

  if inOptions_inout.printDetailsStr == "on":
    setattr(inOptions_inout, 'printDetails', True)
  else:
    setattr(inOptions_inout, 'printDetails', False)

  if inOptions_inout.cdashBaseCacheFilesPrefix == "":
    inOptions_inout.cdashBaseCacheFilesPrefix = \
     CDQAR.getFileNameStrFromText(inOptions_inout.buildSetName)+"_"


def getCmndLineOptions():
  from optparse import OptionParser
  clp = OptionParser(usage=usageHelp)
  injectCmndLineOptionsInParser(clp)
  (options, args) = clp.parse_args()
  validateCmndLineOptions(options)
  setExtraCmndLineOptionsAfterParse(options)
  return options


def fwdCmndLineOptions(inOptions, lt=""):
  cmndLineOpts = \
    "  --date='"+inOptions.date+"'"+lt+\
    "  --cdash-project-name='"+inOptions.cdashProjectName+"'"+lt+\
    "  --build-set-name='"+inOptions.buildSetName+"'"+lt+\
    "  --cdash-site-url='"+inOptions.cdashSiteUrl+"'"+lt+\
    "  --cdash-builds-filters='"+inOptions.cdashBuildsFilters+"'"+lt+\
    "  --cdash-nonpassed-tests-filters='"+inOptions.cdashNonpassedTestsFilters+"'"+lt+\
    "  --expected-builds-file='"+inOptions.expectedBuildsFile+"'"+lt+\
    "  --tests-with-issue-trackers-file='"+inOptions.testsWithIssueTrackersFile+"'"+lt+\
    "  --cdash-queries-cache-dir='"+inOptions.cdashQueriesCacheDir+"'"+lt+\
    "  --cdash-base-cache-files-prefix='"+inOptions.cdashBaseCacheFilesPrefix+"'"+lt+\
    "  --use-cached-cdash-data='"+inOptions.useCachedCDashDataStr+"'"+lt+\
    "  --limit-test-history-days='"+str(inOptions.expectedBuildsFile)+"'"+lt+\
    "  --limit-table-rows='"+str(inOptions.limitTableRows)+"'"+lt+\
    "  --print-details='"+inOptions.printDetailsStr+"'"+lt+\
    "  --write-failing-tests-without-issue-trackers-to-file='"+inOptions.writeFailingTestsWithoutIssueTrackersToFile+"'"+lt+\
    "  --write-email-to-file='"+inOptions.writeEmailToFile+"'"+lt+\
    "  --email-from-address='"+inOptions.emailFromAddress+"'"+lt+\
    "  --send-email-to='"+inOptions.sendEmailTo+"'"+lt
  return cmndLineOpts 


def echoCmndLineOptions(inOptions):
  print(fwdCmndLineOptions(inOptions, " \\\n"))


def echoCmndLine(inOptions):

  print("")
  print("**************************************************************************")
  print("cdash_analyze_and_report.py \\")

  echoCmndLineOptions(inOptions)


# Class object to store and manipulate vars the top-level main() vars that are
# operated on by various functions.
#
# NOTE: This is put into a class object so that these vars can be updated in
# place when passed to a function.
#
class OverallVars(object):
  def __init__(self):
    # Gives the final result (assume passing by defualt)
    self.globalPass = True
    # This is the top of the body
    self.htmlEmailBodyTop = ""
    # This is the bottom of the email body
    self.htmlEmailBodyBottom = ""
    # This var will store the list of data numbers for the summary line
    self.summaryLineDataNumbersList = []


# Class to help get test history and then analyze and report for each test
# set.
#
# NOTE: The reason this is a class is that the data inOptions and overallVars
# never changes once this object is constructed in main().  This avoids having
# to pass these options in every function call for each test set.
#
class TestSetGetDataAnayzeReporter(object):


  def __init__(self, inOptions, testsSortOrder, testHistoryCacheDir, overallVars):
    self.inOptions = inOptions
    self.testsSortOrder = testsSortOrder
    self.testHistoryCacheDir = testHistoryCacheDir
    self.overallVars = overallVars


  def testSetGetDataAnalyzeReport( self,
      testSetType,
      testSetDescr, testSetAcro, testSetTotalSize, testSetLOD,
      testSetNonzeroSizeTriggerGlobalFail=True,
      colorTestSet=None,     # Change to one of the supported colors
      sortTests=True,
      limitTableRows=None,   # Change to 'int' > 0 to limit to this this
      getTestHistory=False,
    ):
  
    print("")
  
    testSetSummaryStr =  CDQAR.getCDashDataSummaryHtmlTableTitleStr(testSetDescr,
      testSetAcro, testSetTotalSize)
  
    print(testSetSummaryStr)
  
    if testSetTotalSize > 0:
  
      self.overallVars.globalPass = False
  
      self.overallVars.summaryLineDataNumbersList.append(
        testSetAcro+"="+str(testSetTotalSize))
  
      self.overallVars.htmlEmailBodyTop += \
        CDQAR.colorHtmlText(testSetSummaryStr, colorTestSet)+"<br>\n"
  
      if sortTests or limitTableRows:
        testSetSortedLimitedLOD = CDQAR.sortAndLimitListOfDicts(
          testSetLOD, self.testsSortOrder,
          limitTableRows )
      else:
        testSetSortedLimitedLOD = testSetLOD
  
      if getTestHistory:

        CDQAR.foreachTransform(
          testSetSortedLimitedLOD,
          CDQAR.AddTestHistoryToTestDictFunctor(
            self.inOptions.cdashSiteUrl,
            self.inOptions.cdashProjectName,
            self.inOptions.date,
            self.inOptions.testHistoryDays,
            self.testHistoryCacheDir,
            useCachedCDashData=self.inOptions.useCachedCDashData,
            alwaysUseCacheFileIfExists=True,
            verbose=True,
            printDetails=self.inOptions.printDetails,
            )
          )
  
      self.overallVars.htmlEmailBodyBottom += CDQAR.createCDashTestHtmlTableStr(
        testSetType,
        testSetDescr, testSetAcro, testSetTotalSize, testSetSortedLimitedLOD,
        self.inOptions.testHistoryDays, limitRowsToDisplay=limitTableRows,
        testSetColor=colorTestSet )


#
# Run the script
#

if __name__ == '__main__':

  #
  # Get commandline options
  #

  inOptions = getCmndLineOptions()
  echoCmndLine(inOptions)

  cacheDirAndBaseFilePrefix = \
    inOptions.cdashQueriesCacheDir+"/"+inOptions.cdashBaseCacheFilesPrefix

  #
  # A) Define common data, etc
  #

  tcd = CDQAR.TableColumnData
  pp = pprint.PrettyPrinter(indent=2)

  groupSiteBuildNameSortOrder = ['group', 'site', 'buildname']

  #
  # B) Sound off
  #

  print "***"
  print "*** Query and analyze CDash results for "+inOptions.buildSetName+\
        " for testing day "+inOptions.date
  print "***"

  #
  # C) Create beginning of email body (that does not require getting any data off CDash)
  #

  # Aggregation of vars that get updated in this main() body and by functions
  # called.
  overallVars = OverallVars()

  overallVars.htmlEmailBodyTop += \
   "<h2>Build and Test results for "+inOptions.buildSetName \
      +" on "+inOptions.date+"</h2>\n\n"

  #
  # D) Read data files, get data off of CDash, do analysis, and construct HTML
  # body parts
  #

  try:

    # Beginning of top full bulid and tests CDash links paragraph 
    overallVars.htmlEmailBodyTop += "<p>\n"

    #
    # D.1) Read data from input files, set up cache directories
    #
    # Assert this data is correct and abort if there is an error before we run
    # expensive CDash queries!
    #

    # Get list of expected builds from input CSV file
    expectedBuildsLOD = []
    if inOptions.expectedBuildsFile:
      expectedBuildsLOD = \
        CDQAR.getExpectedBuildsListfromCsvFile(inOptions.expectedBuildsFile)
    print("\nNum expected builds = "+str(len(expectedBuildsLOD)))

    # Create a SearchableListOfDicts that will look up an expected build given
    # just a test dict fields ['site', 'buildName']. (The list of tests with
    # issue trackers does not have 'group' since cdash/queryTests.php does not
    # give the 'group' associated with each test.  Also, note that we need
    # this special SearchableListOfDicts since the Build Name key name
    # different for a cdash/queryTests.php test dict 'buildName' and a
    # cdash/index.php build dict 'buildname'.)
    testsToExpectedBuildsSLOD = \
      CDQAR.createTestToBuildSearchableListOfDicts(expectedBuildsLOD)
    # ToDo: Put in try/except to print about error in duplicate rows in the
    # list of expected builds.

    # Get list of tests with issue trackers from the input CSV file
    testsWithIssueTrackersLOD = []
    if inOptions.testsWithIssueTrackersFile:
      testsWithIssueTrackersLOD = CDQAR.getTestsWtihIssueTrackersListFromCsvFile(
        inOptions.testsWithIssueTrackersFile)
    print("\nNum tests with issue trackers = "+str(len(testsWithIssueTrackersLOD)))

    # Get a SearchableListOfDicts for the tests with issue trackers to allow
    # them to be looked up based on matching ['site', 'buildName', 'testname']
    # key/value pairs.
    testsWithIssueTrackersSLOD = \
      CDQAR.createSearchableListOfTests(testsWithIssueTrackersLOD)
    # ToDo: Put in try/except to print about error in duplicate rows in the
    # list of tests with issue trackers.

    # Get a functor that will return True if a passed-in test dict matches a
    # test with an issue tracker for the test key/value pairs ['site',
    # 'buildName', and 'testname'].
    testsWithIssueTrackerMatchFunctor = \
      CDQAR.MatchDictKeysValuesFunctor(testsWithIssueTrackersSLOD)

    # Assert they the list of tests with issue trackers matches the list of
    # expected builds
    (allTestsMatch, errMsg) = CDQAR.testsWithIssueTrackersMatchExpectedBuilds(
      testsWithIssueTrackersLOD, testsToExpectedBuildsSLOD)
    if not allTestsMatch:
      raise Exception(errMsg)

    # Test history cache dir
    testHistoryCacheDir = inOptions.cdashQueriesCacheDir+"/test_history"
    if not os.path.exists(testHistoryCacheDir):
      print("\nCreating new test cache directory '"+testHistoryCacheDir+"'") 
      os.mkdir(testHistoryCacheDir)

    #
    # D.2) Get top-level lists of build and nonpassing tests off CDash
    #

    #
    # D.2.a) Get list of dicts of builds off cdash/index.phpp
    #

    cdashIndexBuildsBrowserUrl = CDQAR.getCDashIndexBrowserUrl(
      inOptions.cdashSiteUrl, inOptions.cdashProjectName, inOptions.date,
      inOptions.cdashBuildsFilters)

    print("\nCDash builds browser URL:\n\n  "+cdashIndexBuildsBrowserUrl+"\n")
   
    cdashIndexBuildsQueryUrl = CDQAR.getCDashIndexQueryUrl(
      inOptions.cdashSiteUrl,
      inOptions.cdashProjectName,
      inOptions.date,
      inOptions.cdashBuildsFilters )

    fullCDashIndexBuildsJsonCacheFile = \
      cacheDirAndBaseFilePrefix+"fullCDashIndexBuilds.json"

    buildsLOD = CDQAR.downloadBuildsOffCDashAndFlatten(
      cdashIndexBuildsQueryUrl,
      fullCDashIndexBuildsJsonCacheFile,
      inOptions.useCachedCDashData )
    print("\nNum builds = "+str(len(buildsLOD)))
  
    # HTML line "Builds on CDash" 
    overallVars.htmlEmailBodyTop += \
     "<a href=\""+cdashIndexBuildsBrowserUrl+"\">"+\
     "Builds on CDash</a> (num/expected="+\
     str(len(buildsLOD))+"/"+str(len(expectedBuildsLOD))+")<br>\n"

    # Create a SearchableListOfDict object to help look up builds given a
    # build dict by key/value pairs 'group', 'site', and 'buildname' (requires
    # unique builds with these key/value pairs)
    buildsSLOD = CDQAR.createSearchableListOfBuilds(buildsLOD)
    # ToDo: Add try/except to report duplicate builds in case this raises an
    # exception.

    #
    # D.2.b) Get list of dicts of all nonpassing tests off
    # cdash/queryTests.php
    #

    cdashNonpassingTestsBrowserUrl = CDQAR.getCDashQueryTestsBrowserUrl(
      inOptions.cdashSiteUrl, inOptions.cdashProjectName, inOptions.date,
      inOptions.cdashNonpassedTestsFilters)

    print("\nGetting list of nonpassing tests from CDash ...\n")

    print("\nCDash nonpassing tests browser URL:\n\n"+\
      "  "+cdashNonpassingTestsBrowserUrl+"\n")

    cdashNonpassingTestsQueryUrl = CDQAR.getCDashQueryTestsQueryUrl(
      inOptions.cdashSiteUrl, inOptions.cdashProjectName, inOptions.date,
      inOptions.cdashNonpassedTestsFilters)

    cdashNonpassingTestsQueryJsonCacheFile = \
      cacheDirAndBaseFilePrefix+"fullCDashNonpassingTests.json"

    nonpassingTestsLOD = CDQAR.downloadTestsOffCDashQueryTestsAndFlatten(
      cdashNonpassingTestsQueryUrl, cdashNonpassingTestsQueryJsonCacheFile,
      inOptions.useCachedCDashData )
    print("\nNum nonpassing tests direct from CDash query = "+\
      str(len(nonpassingTestsLOD)))
  
    # HTML line "Nonpassing Tests on CDash"
    overallVars.htmlEmailBodyTop += \
     "<a href=\""+cdashNonpassingTestsBrowserUrl+"\">"+\
     "Non-passing Tests on CDash</a> (num="+str(len(nonpassingTestsLOD))+")<br>\n"
  
    # End of full build and test link paragraph and start the next paragraph
    # for the summary of failures and other tables
    overallVars.htmlEmailBodyTop += \
      "</p>\n\n"+\
      "<p>\n"

    # Create a SearchableListOfDicts object for looking up a nonpassing test
    # given the test dict fields 'site', 'buildName', and 'testname'.
    nonpassingTestsSLOD = CDQAR.createSearchableListOfTests(
      nonpassingTestsLOD, removeExactDuplicateElements=True,
      checkDictsAreSame_in=CDQAR.checkCDashTestDictsAreSame )
    # NOTE: Above we add the option to remove exact duplicate tests since
    # cdash/queryTests.php can return duplicate tests (i.e. all test dict
    # fields are the same except and has the same buildid but could have
    # different testids!)
    # ToDo: Add try/except for above code in order to report duplicate tests
    # where the buildid (and other fields) not match.

    print("Num nonpassing tests after removing duplicate tests = "+\
      str(len(nonpassingTestsLOD)))

    # Create a functor to to see if a test dict matches one of the nonpassing
    # tests downloaded from cdash/queryTests.php.
    nonpassingTestsMatchFunctor = \
      CDQAR.MatchDictKeysValuesFunctor(nonpassingTestsSLOD)

    #
    # D.3) Partition the varous list of tests into different sets that will
    # be displayed in different tables.
    #

    # Add issue tracker info for all nonpassing tests (including adding empty
    # issue tracker fields for tests that don't have issue trackers)
    CDQAR.foreachTransform( nonpassingTestsLOD,
      CDQAR.AddIssueTrackerInfoToTestDictFunctor(testsWithIssueTrackersSLOD))

    # Split the list of nonpassing tests into those with and without issue
    # trackers
    (nonpassingTestsWithIssueTrackersLOD,nonpassingTestsWithoutIssueTrackersLOD)=\
      CDQAR.splitListOnMatch(nonpassingTestsLOD, testsWithIssueTrackerMatchFunctor)
    print("Num nonpassing tests without issue trackers = "+\
      str(len(nonpassingTestsWithoutIssueTrackersLOD)))
    print("Num nonpassing tests with issue trackers = "+\
      str(len(nonpassingTestsWithIssueTrackersLOD)))

    # Split the list nonpassing tests without issue trackers into 'twoif' and
    # 'twoinp'
    (twoifLOD, twoinrLOD) = CDQAR.splitListOnMatch(
      nonpassingTestsWithoutIssueTrackersLOD, CDQAR.isTestFailed)
    print("Num nonpassing tests without issue trackers Failed = "+str(len(twoifLOD)))
    print("Num nonpassing tests without issue trackers Not Run = "+str(len(twoinrLOD)))

    # Split the list nonpassing tests with issue trackers into 'twif' and
    # 'twinp'
    (twifLOD, twinrLOD) = CDQAR.splitListOnMatch(
      nonpassingTestsWithIssueTrackersLOD, CDQAR.isTestFailed)
    print("Num nonpassing tests with issue trackers Failed = "+str(len(twifLOD)))
    print("Num nonpassing tests with issue trackers Not Run = "+str(len(twinrLOD)))

    # Get list of tests with issue trackers that are not in the list of
    # nonpassing tests (and therefore these are passing or missing)
    testsWithIssueTrackersGrossPassingOrMissingLOD = CDQAR.getFilteredList(
      testsWithIssueTrackersSLOD,
      CDQAR.NotMatchFunctor(nonpassingTestsMatchFunctor) )
    print("Num tests with issue trackers gross passing or missing = "+\
      str(len(testsWithIssueTrackersGrossPassingOrMissingLOD)))

    #
    # D.4) Process and tabulate lists of builds
    #

    #
    # 'bm'
    #

    print("\nSearch for any missing expected builds ...\n")

    missingExpectedBuildsLOD = CDQAR.getMissingExpectedBuildsList(
      buildsSLOD, expectedBuildsLOD)
    #print("\nmissingExpectedBuildsLOD:")
    #pp.pprint(missingExpectedBuildsLOD)

    bmDescr = "Builds Missing"
    bmAcro = "bm"
    bmNum = len(missingExpectedBuildsLOD)

    bmSummaryStr = \
      CDQAR.getCDashDataSummaryHtmlTableTitleStr(bmDescr,  bmAcro, bmNum)

    print(bmSummaryStr)

    if bmNum > 0:

      overallVars.globalPass = False

      overallVars.summaryLineDataNumbersList.append(bmAcro+"="+str(bmNum))

      overallVars.htmlEmailBodyTop += CDQAR.colorHtmlText(bmSummaryStr,CDQAR.cdashColorFailed())+"<br>\n"

      bmColDataList = [
        tcd("Group", 'group'),
        tcd("Site", 'site'),
        tcd("Build Name", 'buildname'),
        tcd("Missing Status", 'status'),
        ]

      overallVars.htmlEmailBodyBottom += CDQAR.createCDashDataSummaryHtmlTableStr(
         bmDescr,  bmAcro, bmColDataList, missingExpectedBuildsLOD,
        groupSiteBuildNameSortOrder, None )
      # NOTE: Above we don't want to limit any missing builds in this table
      # because that data is not shown on CDash and that list will never be
      # super big.

    #
    # 'cf'
    #

    print("\nSearch for any builds with configure failures ...\n")

    buildsWithConfigureFailuresLOD = \
      CDQAR.getFilteredList(buildsSLOD, CDQAR.buildHasConfigureFailures)

    cDescr = "Builds with Configure Failures"
    cAcro = "cf"
    cNum = len(buildsWithConfigureFailuresLOD)

    cSummaryStr = \
      CDQAR.getCDashDataSummaryHtmlTableTitleStr(cDescr,  cAcro, cNum)

    print(cSummaryStr)

    if cNum > 0:

      overallVars.globalPass = False

      overallVars.summaryLineDataNumbersList.append(cAcro+"="+str(cNum))

      overallVars.htmlEmailBodyTop += CDQAR.colorHtmlText(cSummaryStr,CDQAR.cdashColorFailed())+"<br>\n"

      cColDataList = [
        tcd("Group", 'group'),
        tcd("Site", 'site'),
        tcd("Build Name", 'buildname'),
        ]

      overallVars.htmlEmailBodyBottom += CDQAR.createCDashDataSummaryHtmlTableStr(
        cDescr,  cAcro, cColDataList, buildsWithConfigureFailuresLOD,
        groupSiteBuildNameSortOrder, inOptions.limitTableRows )

      # ToDo: Update to show number of configure failures and the history info
      # for that build with hyperlinks and don't limit the number of builds
      # shown.

    #
    # 'bf'
    #

    print("\nSearch for any builds with compilation (build) failures ...\n")

    buildsWithBuildFailuresLOD = \
      CDQAR.getFilteredList(buildsSLOD, CDQAR.buildHasBuildFailures)

    bDescr = "Builds with Build Failures"
    bAcro = "bf"
    bNum = len(buildsWithBuildFailuresLOD)

    bSummaryStr = \
      CDQAR.getCDashDataSummaryHtmlTableTitleStr(bDescr,  bAcro, bNum)

    print(bSummaryStr)

    if bNum > 0:

      overallVars.globalPass = False

      overallVars.summaryLineDataNumbersList.append(bAcro+"="+str(bNum))

      overallVars.htmlEmailBodyTop += CDQAR.colorHtmlText(bSummaryStr,CDQAR.cdashColorFailed())+"<br>\n"

      cColDataList = [
        tcd("Group", 'group'),
        tcd("Site", 'site'),
        tcd("Build Name", 'buildname'),
        ]

      overallVars.htmlEmailBodyBottom += CDQAR.createCDashDataSummaryHtmlTableStr(
        bDescr,  bAcro, cColDataList, buildsWithBuildFailuresLOD,
        groupSiteBuildNameSortOrder, inOptions.limitTableRows )

      # ToDo: Update to show number of builds failures and the history info
      # for that build with hyperlinks and don't limit the number of builds
      # shown.

    #
    # D.5) Analyaize and report the different sets of tests
    #

    #
    # D.5.a) Final processing of lists of tests and splitting into the
    # different tests sets to report
    #

    # Sort order for tests to display in tables
    testsSortOrder = ['testname', 'buildName', 'site']

    # Object to make it easy to process the different test sets
    testSetGetDataAnayzeReporter = TestSetGetDataAnayzeReporter(inOptions,
      testsSortOrder, testHistoryCacheDir, overallVars)

    # Special functor to look up missing expected build given a test dict
    testsToMissingExpectedBuildsSLOD = \
      CDQAR.createTestToBuildSearchableListOfDicts(missingExpectedBuildsLOD)

    # Functor for matching an missing expected build given a test dict
    testMatchesMissingExpectedBuildsFunctor = CDQAR.MatchDictKeysValuesFunctor(
      testsToMissingExpectedBuildsSLOD)

    # Get list of tests with issue trackers that are not in the list of
    # nonpassing tests and don't match expected builds (and therefore will be
    # included in the sets 'twip' and 'twim').
    ( testsWithIssueTrackersMatchingMissingExpectedBuildsLOD,
      testsWithIssueTrackersPassingOrMissingLOD ) \
      = \
      CDQAR.splitListOnMatch( testsWithIssueTrackersGrossPassingOrMissingLOD,
        testMatchesMissingExpectedBuildsFunctor )
    print("\nNum tests with issue trackers passing or missing matching"+\
      " posted builds = "+str(len(testsWithIssueTrackersPassingOrMissingLOD)))
    print("\nTests with issue trackers missing that match"+\
      " missing expected builds: num="+\
      str(len(testsWithIssueTrackersMatchingMissingExpectedBuildsLOD)))
    if len(testsWithIssueTrackersMatchingMissingExpectedBuildsLOD) > 0:
      for testDict in testsWithIssueTrackersMatchingMissingExpectedBuildsLOD:
        print("  "+sorted_dict_str(testDict))
      print("\nNOTE: The above tests will NOT be listed in the set 'twim'!")

    # Get test history for all of the tests with issue trackers that are not
    # passing or missing.  These will either be tests that are passing today
    # (and therefore have history) or they will be tests that are missing.
    # (But don't get test history or list out tests with issue trackers that
    # match missing expected builds that did not submit any test data today.)

    twipLOD = []
    twimLOD = []

    if testsWithIssueTrackersPassingOrMissingLOD:

      print("\nGetting test history for tests with issue trackers"+\
        " passing or missing: num="+str(len(testsWithIssueTrackersPassingOrMissingLOD)))

      CDQAR.foreachTransform(
        testsWithIssueTrackersPassingOrMissingLOD,
        CDQAR.AddTestHistoryToTestDictFunctor(
          inOptions.cdashSiteUrl,
          inOptions.cdashProjectName,
          inOptions.date,
          inOptions.testHistoryDays,
          testHistoryCacheDir,
          useCachedCDashData=inOptions.useCachedCDashData,
          alwaysUseCacheFileIfExists=True,
          verbose=True,
          printDetails=inOptions.printDetails,
          )
        )

      # Split into lists for 'twip' and 'twim'
      (twipLOD, twimLOD) = CDQAR.splitListOnMatch(
        testsWithIssueTrackersPassingOrMissingLOD, CDQAR.isTestPassed )

    print("\nNum tests with issue trackers Passed = "+str(len(twipLOD)))
    print("Num tests with issue trackers Missing = "+str(len(twimLOD)))

    #
    # D.5.b) Report the different sets of tests
    #
    # NOTE: The order of these is chosen so those that require action of the
    # person doing the triaging are sorted to the top.
    #

    # twoif
    testSetGetDataAnayzeReporter.testSetGetDataAnalyzeReport( 'nopass',
      "Tests without issue trackers Failed",
      "twoif",
      len(twoifLOD),
      twoifLOD,
      colorTestSet=CDQAR.cdashColorFailed(),
      limitTableRows=inOptions.limitTableRows,
      getTestHistory=True,
      )

    # twoinr
    testSetGetDataAnayzeReporter.testSetGetDataAnalyzeReport( 'nopass',
      "Tests without issue trackers Not Run",
      "twoinr",
      len(twoinrLOD),
      twoinrLOD,
      colorTestSet=CDQAR.cdashColorNotRun(),
      limitTableRows=inOptions.limitTableRows,
      getTestHistory=True,
      )

    # twip
    testSetGetDataAnayzeReporter.testSetGetDataAnalyzeReport( 'pass',
      "Tests with issue trackers Passed",
      "twip",
      len(twipLOD),
      twipLOD,
      colorTestSet=CDQAR.cdashColorPassed(),
      limitTableRows=None,
      getTestHistory=False,  # Already got it above!
      )

    # twim
    testSetGetDataAnayzeReporter.testSetGetDataAnalyzeReport( 'missing',
      "Tests with issue trackers Missing",
      "twim",
      len(twimLOD),
      twimLOD,
      colorTestSet=None,
      limitTableRows=None,
      getTestHistory=False,  # Already got it above!
      )

    # twif
    testSetGetDataAnayzeReporter.testSetGetDataAnalyzeReport( 'nopass',
      "Tests with issue trackers Failed",
      "twif",
      len(twifLOD),
      twifLOD,
      colorTestSet=None,
      limitTableRows=None,
      getTestHistory=True,
      )

    # twinr
    testSetGetDataAnayzeReporter.testSetGetDataAnalyzeReport( 'nopass',
      "Tests with issue trackers Not Run",
      "twinr",
      len(twinrLOD),
      twinrLOD,
      colorTestSet=None,
      limitTableRows=None,
      getTestHistory=True,
      )

    #
    # D.6) Write out list twiof to CSV file
    #

    if inOptions.writeFailingTestsWithoutIssueTrackersToFile:

      twoifCsvFileName = inOptions.writeFailingTestsWithoutIssueTrackersToFile

      print("\nWriting list of 'twiof' to file "+twoifCsvFileName+" ...")

      CDQAR.writeTestsLODToCsvFile(twoifLOD, twoifCsvFileName)

  except Exception:
    # Traceback!
    print("")
    sys.stdout.flush()
    traceback.print_exc()
    # Report the error
    overallVars.htmlEmailBodyBottom += "\n<pre><code>\n"+\
      traceback.format_exc()+"\n</code></pre>\n"
    print("\nError, could not compute the analysis due to"+\
      " above error so return failed!")
    overallVars.globalPass = False
    overallVars.summaryLineDataNumbersList.append("SCRIPT CRASHED")

  #
  # E) Put together final email summary line
  #

  if overallVars.globalPass:
    summaryLine = "PASSED"
  else:
    summaryLine = "FAILED"

  if overallVars.summaryLineDataNumbersList:
    summaryLine += " (" + ", ".join(overallVars.summaryLineDataNumbersList) + ")"

  summaryLine += ": "+inOptions.buildSetName+" on "+inOptions.date

  #
  # F) Finish off HTML body guts and define overall HTML body style
  #

  # Finish off the top paragraph of the summary lines
  overallVars.htmlEmailBodyTop += \
    "</p>"
    
  # Construct HTML body guts without header or begin/end body.
  htmlEmaiBodyGuts = \
    overallVars.htmlEmailBodyTop+\
    "\n\n"+\
    overallVars.htmlEmailBodyBottom

  htmlHeaderAndBeginBody = \
    "<html>\n"+\
    "<head>\n"+\
    "<style>\n"+\
    "h1 {\n"+\
    "  font-size: 40px;\n"+\
    "}\n"+\
    "h2 {\n"+\
    "  font-size: 30px;\n"+\
    "}\n"+\
    "h3 {\n"+\
    "  font-size: 24px;\n"+\
    "}\n"+\
    "p {\n"+\
    "  font-size: 18px;\n"+\
    "}\n"+\
    "</style>\n"+\
    "</head>\n"+\
    "\n"+\
    "<body>\n"+\
    "\n"

  htmlEndBody = \
    "</body>\n"+\
    "</html>\n"

  #
  # G) Write HTML body file and/or send HTML email(s)
  #

  if inOptions.writeEmailToFile:
    print("\nWriting HTML file '"+inOptions.writeEmailToFile+"' ...")
    htmlEmaiBodyFileStr = \
      htmlHeaderAndBeginBody+\
      "<h2>"+summaryLine+"</h2>\n\n"+\
      htmlEmaiBodyGuts+"\n"+\
      htmlEndBody
    with open(inOptions.writeEmailToFile, 'w') as outFile:
      outFile.write(htmlEmaiBodyFileStr)

  if inOptions.sendEmailTo:
    htmlEmaiBody = \
      htmlHeaderAndBeginBody+\
      htmlEmaiBodyGuts+"\n"+\
      htmlEndBody
    for emailAddress in inOptions.sendEmailTo.split(','):
      emailAddress = emailAddress.strip()
      print("\nSending email to '"+emailAddress+"' ...")
      msg=CDQAR.createHtmlMimeEmail(
        inOptions.emailFromAddress, emailAddress, summaryLine, "", htmlEmaiBody)
      CDQAR.sendMineEmail(msg)

  #
  # H) Return final global pass/fail
  #

  print("\n"+summaryLine+"\n")

  if overallVars.globalPass:
    sys.exit(0)
  else:
    sys.exit(1)
