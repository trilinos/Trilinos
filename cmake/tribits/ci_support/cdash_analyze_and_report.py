#!/usr/bin/env python

# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

import sys
import pprint
import datetime

from FindGeneralScriptSupport import *
from GeneralScriptSupport import *
import CDashQueryAnalyzeReport as CDQAR
import cdash_build_testing_date as CBTD
from gitdist import addOptionParserChoiceOption

#
# Help message
#


usageHelp = r"""cdash_analyze_and_report.py [options]

This script takes in CDash URL information and other data as command-line
arguments and then analyzes it to look for missing expected builds, failed
tests, and various types of other failures and then reports the findings as an
HTML file written to disk and/or an HTML-formatted email sent to one or more
email addresses.  (Other types of output can be produced as well in different
files.)

If all of the expected builds are found (and all of them have test results)
and there are no other failures found, then the script returns 0.  Otherwise
the script returns non-zero.  Therefore, this script can be used to drive
automated workflows by examining data on CDash.
"""


#
# Helper functions
#


def injectCmndLineOptionsInParser(clp, gitoliteRootDefault=""):

  clp.add_option(
    "--date", dest="date", type="string", default='yesterday',
    help="Date for the testing day <YYYY-MM-DD> or special values 'today'"+\
      " or 'yesterday'. [default 'yesterday']" )

  clp.add_option(
    "--cdash-project-testing-day-start-time", dest="cdashProjectTestingDayStartTime",
    type="string", default="00:00",
    help="The CDash project testing day build star time in UTC in format '<hh>:<mm>'."+\
      " [default = '00:00']" )

  clp.add_option(
    "--cdash-project-name", dest="cdashProjectName", type="string", default="",
    help="CDash project name (e.g. 'Trilinos'). [REQUIRED]" )

  clp.add_option(
    "--build-set-name", dest="buildSetName", type="string", default="",
    help="Name for the set of builds, (e.g. 'Trilinos Nightly Builds)."+\
      "  This used in the email summary line and in the HTML file body"+\
      " to identify the set of builds and tests being examined."+\
      " [REQUIRED]" )

  clp.add_option(
    "--cdash-site-url", dest="cdashSiteUrl", type="string", default="",
    help="Base CDash site (e.g. 'https://testing.sandia.gov/cdash')."+\
      " [REQUIRED]" )

  clp.add_option(
    "--cdash-builds-filters", dest="cdashBuildsFilters", type="string",
    default="",
    help="Partial URL fragment for index.php making of the filters for"+\
      " the set of builds (e.g. 'filtercount=1&showfilters=1&field1=groupname&compare1=61&value1=ATDM')."+\
      " [REQUIRED]" )

  clp.add_option(
    "--cdash-nonpassed-tests-filters", dest="cdashNonpassedTestsFilters", type="string",
    default="",
    help="Partial URL fragment for queryTests.php making of the filters for"+\
      " the set of non-passing tests matching this set of builds (e.g."+\
      " 'filtercombine=and&filtercount=1&showfilters=1&filtercombine=and&field1=groupname&compare1=61&value1=ATDM')."+\
      "  This set of filter fields may also filter out extra nonpassing tests"+\
      " such for known random system failures to avoid flooding the output.  In this"+\
      " case, one should also set --require-test-history-match-nonpassing-tests=off."+\
      " [REQUIRED]" )

  clp.add_option(
    "--expected-builds-file", dest="expectedBuildsFile", type="string",
    default="",
    help="Path to a CSV file that lists the expected builds.  Each of these builds"+\
      " must have unique 'site' and 'buildname' field pairs or an error will be"+\
      " raised and the tool will abort.  A list of files is also allowed that are"+\
      " separated with ',' as <file1>,<file2>,... [default = '']" )

  clp.add_option(
    "--tests-with-issue-trackers-file", dest="testsWithIssueTrackersFile",
    type="string",  default="",
    help="Path to CSV file that lists tests with issue trackers (and other data)."+\
    "  Each of these tests must have a unique 'site', 'buildName', and 'testname'"+\
    " sets or an error will be raised and the tool will abort.  [default = '']" )

  addOptionParserChoiceOption(
    "--filter-out-builds-and-tests-not-matching-expected-builds",
    "filterOutBuildsAndTestsNotMatchingExpectedBuildsStr",
    ("on", "off"), 1,
    "Filter out build and test data not matching input list of expected builds."+\
    "  If set to 'on', this will filter out build and test data downloaded"+\
    " from CDash that does not match the list of expected builds provided in"+\
    " --expected-builds-file=<csv-file>.  This will also filter out any tests"+\
    " with issue trackers listed in"+\
    " --tests-with-issue-trackers-file=<csv-file>.",
    clp )

  cdashQueriesCacheDir_default=os.getcwd()

  clp.add_option(
    "--cdash-queries-cache-dir", dest="cdashQueriesCacheDir", type="string",
    default=cdashQueriesCacheDir_default,
    help="Cache CDash query data this directory." \
      +" [default = '"+cdashQueriesCacheDir_default+"']" )

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
    "--require-test-history-match-nonpassing-tests",
    "requireTestHistoryMatchNonpassingTestsStr",
    ("on", "off"), 0,
    "Require that the status for each tracked test listed in the tests with issue"\
    +" trackers CSV file match the status of that test returned from the test history"\
    +" returned from CDash.  In general, these should match up but these may not if extra"\
    +" filter criteria has been added to the list on nonpassing tests in the"\
    +" --cdash-nonpassed-tests-filters=<filters> set of filters (such as to filter out"\
    +" a large number of random system failures).  In this case, an error will be"\
    +" returned by default and the script will crash.  But this can be relaxed by"\
    +" setting this to 'off' which will result in these tracked tests being listed in"\
    +" the 'twim' table but typically with status 'Failed'.",
    clp )

  addOptionParserChoiceOption(
    "--print-details", "printDetailsStr",
    ("on", "off"), 1,
    "Print more info like the CDash URLs for downloaded data and the cache"+\
      " file names.",
    clp )

  addOptionParserChoiceOption(
    "--list-unexpected-builds",
    "listUnexpectedBuildsStr",
    ("on", "off"), 1,
    "List unexpected builds downloaded from CDash (i.e. not matching an expected build)'.",
    clp )

  clp.add_option(
    "--write-unexpected-builds-to-file",
    dest="writeUnexpectedBuildsToFile", type="string", default="",
    help="Write a CSV file with a list of unexpected builds 'bu'." \
    +"  This is to make it easy to add new entries to the file read by" \
    +" the option --expected-builds-file=<csv-file>. [default = '']" )

  clp.add_option(
    "--write-failing-tests-without-issue-trackers-to-file",
    dest="writeFailingTestsWithoutIssueTrackersToFile", type="string", default="",
    help="Write a CSV file with a list of tests with issue trackers failed 'twif'." \
    +"  This is to make it easy to add new entries to the file read by" \
    +" the option --tests-with-issue-trackers-file=<csv-file>. [default = '']" )

  clp.add_option(
    "--write-test-data-to-file",
    dest="writeTestDataToFile", type="string", default="",
    help="Write pretty-printed Python list of dictionaries for tests" \
    +" with issue trackers.  This includes the history of the tests for" \
    +" --limit-test-history-days=<days> of history.  This contains all of the" \
    +" information that appears in the generated summary tables for tests with" \
    +" associated issue trackers.  [default = '']" )

  clp.add_option(
    "--write-email-to-file", dest="writeEmailToFile", type="string", default="",
    help="Write the body of the HTML email to this file. [default = '']" )

  clp.add_option(
    "--email-from-address=", dest="emailFromAddress", type="string", default="",
    help="Address reported in the sent email. [default '']" )

  clp.add_option(
    "--send-email-to=", dest="sendEmailTo", type="string", default="",
    help="Send email to 'address1, address2, ...'.  [default '']" )

  addOptionParserChoiceOption(
    "--email-without-soft-hyphens",
    "emailWithoutSoftHyphensStr",
    ("on", "off"), 1,
    "Remove soft hyphens from emails.",
    clp )


def validateAndConvertCmndLineOptions(inOptions):

  if inOptions.date == "":
    print("Error, can't have empty --date, must pass in --date=YYYY-MM-DD"+\
      " or special values --date=today or --date=yesterday!")
    sys.exit(1)
  else:
    dateTimeObj = CDQAR.convertInputDateArgToYYYYMMDD(
      inOptions.cdashProjectTestingDayStartTime,
      inOptions.date)
    inOptions.date = CBTD.getDateStrFromDateTime(dateTimeObj)

  # ToDo: Assert more of the options to make sure they are correct!


def setExtraCmndLineOptionsAfterParse(inOptions_inout):

  setattr(inOptions_inout, 'filterOutBuildsAndTestsNotMatchingExpectedBuilds',
    inOptions_inout.filterOutBuildsAndTestsNotMatchingExpectedBuildsStr == "on")

  setattr(inOptions_inout, 'useCachedCDashData',
    inOptions_inout.useCachedCDashDataStr == "on")

  setattr(inOptions_inout, 'requireTestHistoryMatchNonpassingTests',
    inOptions_inout.requireTestHistoryMatchNonpassingTestsStr == "on")

  setattr(inOptions_inout, 'printDetails',
    inOptions_inout.printDetailsStr == "on")

  setattr(inOptions_inout, 'listUnexpectedBuilds',
    inOptions_inout.listUnexpectedBuildsStr == "on")

  if inOptions_inout.cdashBaseCacheFilesPrefix == "":
    inOptions_inout.cdashBaseCacheFilesPrefix = \
     CDQAR.getFileNameStrFromText(inOptions_inout.buildSetName)+"_"

  setattr(inOptions_inout, 'emailWithoutSoftHyphens',
    inOptions_inout.emailWithoutSoftHyphensStr == "on")


def getCmndLineOptions():
  from optparse import OptionParser
  clp = OptionParser(usage=usageHelp)
  injectCmndLineOptionsInParser(clp)
  (options, args) = clp.parse_args()
  validateAndConvertCmndLineOptions(options)
  setExtraCmndLineOptionsAfterParse(options)
  return options


def fwdCmndLineOptions(inOptions, lt=""):
  io = inOptions
  cmndLineOpts = \
    "  --date='"+io.date+"'"+lt+\
    "  --cdash-project-testing-day-start-time='"+io.cdashProjectTestingDayStartTime+"'"+lt+\
    "  --cdash-project-name='"+io.cdashProjectName+"'"+lt+\
    "  --build-set-name='"+io.buildSetName+"'"+lt+\
    "  --cdash-site-url='"+io.cdashSiteUrl+"'"+lt+\
    "  --cdash-builds-filters='"+io.cdashBuildsFilters+"'"+lt+\
    "  --cdash-nonpassed-tests-filters='"+io.cdashNonpassedTestsFilters+"'"+lt+\
    "  --expected-builds-file='"+io.expectedBuildsFile+"'"+lt+\
    "  --tests-with-issue-trackers-file='"+io.testsWithIssueTrackersFile+"'"+lt+\
    "  --filter-out-builds-and-tests-not-matching-expected-builds='"+\
      io.filterOutBuildsAndTestsNotMatchingExpectedBuildsStr+"'"+lt+\
    "  --cdash-queries-cache-dir='"+io.cdashQueriesCacheDir+"'"+lt+\
    "  --cdash-base-cache-files-prefix='"+io.cdashBaseCacheFilesPrefix+"'"+lt+\
    "  --use-cached-cdash-data='"+io.useCachedCDashDataStr+"'"+lt+\
    "  --limit-test-history-days='"+str(io.testHistoryDays)+"'"+lt+\
    "  --limit-table-rows='"+str(io.limitTableRows)+"'"+lt+\
    "  --require-test-history-match-nonpassing-tests='"+io.requireTestHistoryMatchNonpassingTestsStr+"'"+lt+\
    "  --print-details='"+io.printDetailsStr+"'"+lt+\
    "  --list-unexpected-builds='"+io.listUnexpectedBuildsStr+"'"+lt+\
    "  --write-unexpected-builds-to-fileo='"+io.writeUnexpectedBuildsToFile+"'"+lt+\
    "  --write-failing-tests-without-issue-trackers-to-file='"+io.writeFailingTestsWithoutIssueTrackersToFile+"'"+lt+\
    "  --write-test-data-to-file='"+io.writeTestDataToFile+"'"+lt+\
    "  --write-email-to-file='"+io.writeEmailToFile+"'"+lt+\
    "  --email-from-address='"+io.emailFromAddress+"'"+lt+\
    "  --send-email-to='"+io.sendEmailTo+"'"+lt+\
    "  --email-without-soft-hyphens='"+io.emailWithoutSoftHyphensStr+"'"+lt
  return cmndLineOpts


def echoCmndLineOptions(inOptions):
  print(fwdCmndLineOptions(inOptions, " \\\n"))


def echoCmndLine(inOptions):

  print("")
  print("**************************************************************************")
  print("cdash_analyze_and_report.py \\")

  echoCmndLineOptions(inOptions)


# Strategy class that can get test history for a list of tests and set them in
# the test dicts taking input from the cdash_analyze_and_report.py commandline
# arguments.
#
class AddTestHistoryStrategy(object):


  def __init__(self, inOptions, testHistoryCacheDir):
    self.inOptions = inOptions
    self.testHistoryCacheDir = testHistoryCacheDir


  def getTestHistory(self, testLOD):

    sio = self.inOptions

    CDQAR.foreachTransform(
      testLOD,
      CDQAR.AddTestHistoryToTestDictFunctor(
        cdashUrl=sio.cdashSiteUrl,
        projectName=sio.cdashProjectName,
        date=sio.date,
        testingDayStartTimeUtc=sio.cdashProjectTestingDayStartTime,
        daysOfHistory=sio.testHistoryDays,
        testCacheDir=self.testHistoryCacheDir,
        useCachedCDashData=sio.useCachedCDashData,
        alwaysUseCacheFileIfExists=True,
        verbose=True,
        printDetails=sio.printDetails,
        requireMatchTestTopTestHistory=sio.requireTestHistoryMatchNonpassingTests,
        )
      )


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

  print("***")
  print("*** Query and analyze CDash results for "+inOptions.buildSetName+\
        " for testing day "+inOptions.date)
  print("***")

  #
  # C) Create beginning of email body (that does not require getting any data off CDash)
  #

  # Aggregation of vars that get updated in this main() body and by functions
  # called.
  cdashReportData = CDQAR.CDashReportData()

  cdashReportData.htmlEmailBodyTop += \
   "<h2>Build and Test results for "+inOptions.buildSetName \
      +" on "+inOptions.date+"</h2>\n\n"

  #
  # D) Read data files, get data off of CDash, do analysis, and construct HTML
  # body parts
  #

  try:

    # Beginning of top full build and tests CDash links paragraph
    cdashReportData.htmlEmailBodyTop += "<p>\n"

    #
    # D.1) Read data from input files, set up cache directories
    #
    # Assert this data is correct and abort if there is an error before we run
    # expensive CDash queries!
    #

    # Get list of expected builds from input CSV file
    expectedBuildsLOD = CDQAR.getExpectedBuildsListOfDictsFromCsvFileArg(
      inOptions.expectedBuildsFile)
    print("\nNum expected builds = "+str(len(expectedBuildsLOD)))

    # Create a SearchableListOfDict object to help look up expected builds
    # given a build dict by key/value pairs 'group', 'site', and 'buildname'
    # (requires unique builds with these key/value pairs)
    expectedBuildsSLOD = CDQAR.createSearchableListOfBuilds(expectedBuildsLOD)

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
    if inOptions.testsWithIssueTrackersFile:
      fullTestsWithIssueTrackersLOD = CDQAR.getTestsWtihIssueTrackersListFromCsvFile(
        inOptions.testsWithIssueTrackersFile)
    else:
      fullTestsWithIssueTrackersLOD = []
    print("\nNum tests with issue trackers read from CSV file = "+\
      str(len(fullTestsWithIssueTrackersLOD)))

    if inOptions.filterOutBuildsAndTestsNotMatchingExpectedBuilds:
      (testsWithIssueTrackersLOD, testsWithIssueTrackersNotExpectedLOD) = \
        CDQAR.splitTestsOnMatchExpectedBuilds(fullTestsWithIssueTrackersLOD,
          testsToExpectedBuildsSLOD)
      print("Num tests with issue trackers matching expected builds = "+\
        str(len(testsWithIssueTrackersLOD)))
    else:
      testsWithIssueTrackersLOD = fullTestsWithIssueTrackersLOD
    print("Num tests with issue trackers = "+\
      str(len(testsWithIssueTrackersLOD)))

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

    # Assert that the list of tests with issue trackers matches the list of
    # expected builds
    (allTestsMatch, errMsg) = CDQAR.doTestsWithIssueTrackersMatchExpectedBuilds(
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

    # @arghdos: note, we do not have to normalize the URLs from the input
    # options because they are currently taken from the cdash site already
    # (i.e., they are already in normalized form).

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

    fullBuildsLOD = CDQAR.downloadBuildsOffCDashAndFlatten(
      cdashIndexBuildsQueryUrl,
      fullCDashIndexBuildsJsonCacheFile,
      inOptions.useCachedCDashData )

    print("\nNum builds downloaded from CDash = "+str(len(fullBuildsLOD)))

    (buildsExpectedLOD, buildsUnexpectedLOD) = \
      CDQAR.splitTestsOnMatchExpectedBuilds(fullBuildsLOD, expectedBuildsSLOD)

    if inOptions.filterOutBuildsAndTestsNotMatchingExpectedBuilds:
      print("Num builds matching expected builds = "+str(len(buildsExpectedLOD)))
      buildsLOD = buildsExpectedLOD
    else:
      buildsLOD = fullBuildsLOD

    if inOptions.listUnexpectedBuilds:
      print("Num builds unexpected = "+str(len(buildsUnexpectedLOD)))

    print("Num builds = "+str(len(buildsLOD)))

    # HTML line "Builds on CDash"
    cdashReportData.htmlEmailBodyTop += \
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

    fullNonpassingTestsLOD = CDQAR.downloadTestsOffCDashQueryTestsAndFlatten(
      cdashNonpassingTestsQueryUrl, cdashNonpassingTestsQueryJsonCacheFile,
      inOptions.useCachedCDashData )

    print("\nNum nonpassing tests direct from CDash query = "+\
      str(len(fullNonpassingTestsLOD)))

    if inOptions.filterOutBuildsAndTestsNotMatchingExpectedBuilds:
      (nonpassingTestsLOD, nonpassingTestsNotExpectedLOD) = \
        CDQAR.splitTestsOnMatchExpectedBuilds(fullNonpassingTestsLOD,
          testsToExpectedBuildsSLOD)
      print("Num nonpassing tests matching expected builds = "+\
       str(len(nonpassingTestsLOD)))
    else:
      nonpassingTestsLOD = fullNonpassingTestsLOD

    print("Num nonpassing tests = "+\
      str(len(nonpassingTestsLOD)))

    # HTML line "Nonpassing Tests on CDash"
    cdashReportData.htmlEmailBodyTop += \
     "<a href=\""+cdashNonpassingTestsBrowserUrl+"\">"+\
     "Non-passing Tests on CDash</a> (num="+str(len(nonpassingTestsLOD))+")<br>\n"

    # End of full build and test link paragraph and start the next paragraph
    # for the summary of failures and other tables
    cdashReportData.htmlEmailBodyTop += \
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
    # D.3) Partition the various list of tests into different sets that will
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

    buildsetReporter = CDQAR.SingleBuildsetReporter(cdashReportData)

    #
    # 'bm'
    #

    print("\nSearch for any missing expected builds ...\n")

    missingExpectedBuildsLOD = CDQAR.getMissingExpectedBuildsList(
      buildsSLOD, expectedBuildsLOD)

    buildsetReporter.reportSingleBuildset("Builds Missing", "bm",
      missingExpectedBuildsLOD,
      buildsetGlobalPass=False,
      buildsetColor=CDQAR.cdashColorFailed(),
      buildsetColDataList=[
        tcd("Group", 'group'),
        tcd("Site", 'site'),
        tcd("Build Name", 'buildname'),
        tcd("Missing Status", 'status'),
        ],
      )

    #
    # 'cf'
    #

    print("\nSearch for any builds with configure failures ...\n")

    buildsWithConfigureFailuresLOD = \
      CDQAR.getFilteredList(buildsSLOD, CDQAR.buildHasConfigureFailures)

    buildsetReporter.reportSingleBuildset("Builds with Configure Failures", "cf",
      buildsWithConfigureFailuresLOD,
      buildsetGlobalPass=False,
      buildsetColor=CDQAR.cdashColorFailed(),
      )

    #
    # 'bf'
    #

    print("\nSearch for any builds with compilation (build) failures ...\n")

    buildsWithBuildFailuresLOD = \
      CDQAR.getFilteredList(buildsSLOD, CDQAR.buildHasBuildFailures)

    buildsetReporter.reportSingleBuildset("Builds with Build Failures", "bf",
      buildsWithBuildFailuresLOD,
      buildsetGlobalPass=False,
      buildsetColor=CDQAR.cdashColorFailed(),
      )

    #
    # 'bu'
    #

    if inOptions.listUnexpectedBuilds:
      buildsetReporter.reportSingleBuildset("Builds Unexpected", "bu",
        buildsUnexpectedLOD,
        buildsetGlobalPass=True,
        buildsetColor=None,
        )

    #
    # D.5) Analyaize and report the different sets of tests
    #

    #
    # D.5.a) Final processing of lists of tests and splitting into the
    # different tests sets to report
    #

    # Object to make it easy to process the different test sets
    addTestHistoryStrategy = AddTestHistoryStrategy(inOptions, testHistoryCacheDir)
    testsetReporter = CDQAR.SingleTestsetReporter(cdashReportData,
      addTestHistoryStrategy=addTestHistoryStrategy)

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
          inOptions.cdashProjectTestingDayStartTime,
          inOptions.testHistoryDays,
          testHistoryCacheDir,
          useCachedCDashData=inOptions.useCachedCDashData,
          alwaysUseCacheFileIfExists=True,
          verbose=True,
          printDetails=inOptions.printDetails,
          requireMatchTestTopTestHistory=inOptions.requireTestHistoryMatchNonpassingTests,

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

    testsetReporter.reportSingleTestset(
      CDQAR.getStandardTestsetTypeInfo('twoif'),
      len(twoifLOD), twoifLOD,
      limitTableRows=inOptions.limitTableRows,
      getTestHistory=True,
      )

    testsetReporter.reportSingleTestset(
      CDQAR.getStandardTestsetTypeInfo('twoinr'),
      len(twoinrLOD), twoinrLOD,
      limitTableRows=inOptions.limitTableRows,
      getTestHistory=True,
      )

    testsetReporter.reportSingleTestset(
      CDQAR.getStandardTestsetTypeInfo('twip'),
      len(twipLOD), twipLOD,
      limitTableRows=None,
      getTestHistory=False,  # Already got it above!
      )

    testsetReporter.reportSingleTestset(
      CDQAR.getStandardTestsetTypeInfo('twim', ""),
      len(twimLOD), twimLOD,
      limitTableRows=None,
      getTestHistory=False,  # Already got it above!
      )

    testsetReporter.reportSingleTestset(
      CDQAR.getStandardTestsetTypeInfo('twif', ""),
      len(twifLOD), twifLOD,
      limitTableRows=None,
      getTestHistory=True,
      )

    testsetReporter.reportSingleTestset(
      CDQAR.getStandardTestsetTypeInfo('twinr', ""),
      len(twinrLOD), twinrLOD,
      limitTableRows=None,
      getTestHistory=True,
      )

    #
    # D.6) Write out list of unexpected builds to CSV file
    #

    if inOptions.writeUnexpectedBuildsToFile:
      unexpectedBuildsCsvFileName = inOptions.writeUnexpectedBuildsToFile
      print("\nWriting list of unexpected builds to file "\
        +unexpectedBuildsCsvFileName+" ...")
      CDQAR.writeExpectedBuildsListOfDictsToCsvFile(buildsUnexpectedLOD,
        unexpectedBuildsCsvFileName)

    #
    # D.7) Write out list twiof to CSV file
    #

    if inOptions.writeFailingTestsWithoutIssueTrackersToFile:
      twoifCsvFileName = inOptions.writeFailingTestsWithoutIssueTrackersToFile
      print("\nWriting list of 'twiof' to file "+twoifCsvFileName+" ...")
      CDQAR.writeTestsListOfDictsToCsvFile(twoifLOD, twoifCsvFileName)

    #
    # D.8) Write out test data to CSV file
    #

    if inOptions.writeTestDataToFile:
      testDataFileName = inOptions.writeTestDataToFile
      print("\nWriting out gathered test data to file "+testDataFileName+" ...")
      testDataLOD = twipLOD + twimLOD + twifLOD + twinrLOD
      CDQAR.foreachTransform(testDataLOD,
        CDQAR.AddCDashTestingDayFunctor(inOptions.date))
      # ToDo: Add the first inOptions.limitTableRows elements of twiofLOD and twoinrLOD?
      CDQAR.pprintPythonDataToFile(testDataLOD, testDataFileName)

  except Exception:
    # Traceback!
    print("")
    sys.stdout.flush()
    traceback.print_exc()
    # Report the error
    cdashReportData.htmlEmailBodyBottom += "\n<pre><code>\n"+\
      traceback.format_exc()+"\n</code></pre>\n"
    print("\nError, could not compute the analysis due to"+\
      " above error so return failed!")
    cdashReportData.globalPass = False
    cdashReportData.summaryLineDataNumbersList.append("SCRIPT CRASHED")

  #
  # E) Put together final email summary line
  #

  summaryLine = CDQAR.getOverallCDashReportSummaryLine(cdashReportData,
    inOptions.buildSetName, inOptions.date)

  #
  # F) Finish off HTML body guts and define overall HTML body style
  #

  # Finish off the top paragraph of the summary lines
  cdashReportData.htmlEmailBodyTop += \
    "</p>\n"

  #
  # G) Write HTML body file and/or send HTML email(s)
  #

  defaultPageStyle = CDQAR.getDefaultHtmlPageStyleStr()

  if inOptions.writeEmailToFile:
    print("\nWriting HTML file '"+inOptions.writeEmailToFile+"' ...")
    fullCDashHtmlReportPageStr = CDQAR.getFullCDashHtmlReportPageStr(cdashReportData,
      pageTitle=summaryLine, pageStyle=defaultPageStyle)
    with open(inOptions.writeEmailToFile, 'w') as outFile:
      outFile.write(fullCDashHtmlReportPageStr)

  if inOptions.sendEmailTo:
    htmlEmailBodyStr = CDQAR.getFullCDashHtmlReportPageStr(cdashReportData,
      pageStyle=defaultPageStyle)
    for emailAddress in inOptions.sendEmailTo.split(','):
      emailAddress = emailAddress.strip()
      print("\nSending email to '"+emailAddress+"' ...")
      msg=CDQAR.createHtmlMimeEmail(
        inOptions.emailFromAddress, emailAddress, summaryLine, "",
        htmlEmailBodyStr, inOptions.emailWithoutSoftHyphens)
      CDQAR.sendMineEmail(msg)

  #
  # H) Return final global pass/fail
  #

  print("\n"+summaryLine+"\n")

  if cdashReportData.globalPass:
    sys.exit(0)
  else:
    sys.exit(1)
