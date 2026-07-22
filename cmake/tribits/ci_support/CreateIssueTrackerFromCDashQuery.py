import os
import sys
import pprint

import CDashQueryAnalyzeReport as CDQAR
import cdash_build_testing_date as CDBTD
from Python2and3 import dictItems

g_pp = pprint.PrettyPrinter(indent=2)


# Class that is used to create issue trackers for nonpassing tests downloaded
# from a CDash queryTests.php query.
#
# To make this general, the user has to create a class object that is passed
# in as the argument 'issueTrackerFormatter' in this class that supports
# the following functions:
#
# issueTrackerFormatter.createFormattedIssueTracker((issueTrackerData):
#
#   This function is expected to return the created text for the issue tracker
#   given the data provided by the argument 'issueTrackerData' of type
#   'IssueTrackerData'.  The format of this text can be anything such as
#   Markdown, RestructredText, JIRA Issue Markup, etc.  The only assumption
#   about formatting is in the HTML table passed in the argument
#   'issueTrackerData.testHistoryHtmlTableText'.
#
# The object 'issueTrackerFormatter' can format the issue tracker text anyway
# it wants.  The object 'issueTrackerFormatter' is stateless w.r.t. this class
# CreateIssueTrackerFromCDashQueryDriver.
#
class CreateIssueTrackerFromCDashQueryDriver:


  def __init__(self, issueTrackerFormatter, cdashProjectStartTimeUtc=None,
      usageHelp="", issueTrackerUrlTemplate="", issueTrackerTemplate="",
    ):
    self.issueTrackerFormatter = issueTrackerFormatter
    self.cdashProjectStartTimeUtc = cdashProjectStartTimeUtc
    self.usageHelp = usageHelp
    self.issueTrackerUrlTemplate = issueTrackerUrlTemplate
    self.issueTrackerTemplate = issueTrackerTemplate
    self.options = None


  def runDriver(self):

    print("\n***")
    print("*** Getting data to create a new issue tracker")
    print("***\n")

    self.getCmndLineOptions()

    print("Downloading full list of nonpassing tests from CDash URL:\n")
    print("   "+self.options.nonpassingTestsUrl+"\n")

    nonpassingTestsLOD = self.downloadNonpassingTestsData()
    print("\nTotal number of nonpassing tests over all days = "\
      +str(len(nonpassingTestsLOD)))

    uniqNonpassingTestsLOD = getUniqueTestsListOfDicts(nonpassingTestsLOD)
    print("\nTotal number of unique nonpassing test/build pairs over all days = "\
          +str(len(uniqNonpassingTestsLOD)))

    (testnameList, buildnameList) = getTestnameAndBuildnameLists(uniqNonpassingTestsLOD)
    print("\nNumber of test names = "+str(len(testnameList)))
    print("\nNumber of build names = "+str(len(buildnameList)))

    testingDayStartNonpassingDate = getTestingDayStartNonpassingDate(
      nonpassingTestsLOD, self.options.cdashProjectStartTimeUtc)
    #print("testingDayStartNonpassingDate = "+testingDayStartNonpassingDate)

    testHistoryHtmlTableText = ""  # ToDo: Implement!

    issueTrackerText = self.issueTrackerFormatter.createFormattedIssueTracker(
      IssueTrackerData(
        summaryLine=self.options.summaryLine,
        testingDayStartNonpassingDate=testingDayStartNonpassingDate,
        nonpassingTestsUrl=self.options.nonpassingTestsUrl,
        uniqNonpassingTestsLOD=uniqNonpassingTestsLOD,
        buildnameList=buildnameList,
        testnameList=testnameList,
        testHistoryHtmlTableText=testHistoryHtmlTableText,
        )
      )

    self.writeNewIssueTrackerFile(issueTrackerText)

    self.writeNewTestsWithIssueTrackersFile(uniqNonpassingTestsLOD)


  def injectExtraCmndLineOptions(self, clp):
    clp.add_argument(
      "--nonpassing-tests-url", "-u", dest="nonpassingTestsUrl", default="",
      help="Full CDash queryTest.php URL for the list of nonpassing tests over a"\
        +" time range given as 'begin' and 'end' fields [Required]." )
    clp.add_argument(
      "--cdash-project-start-time", dest="cdashProjectStartTimeUtc", default="",
      help="Starting time for the CDash testing day in 'hh:mm' in UTC."\
        + " Check the CDash project settings for the testing day start time." )
    clp.add_argument(
      "--summary-line", "-s", dest="summaryLine", default="<summary-line>",
      help="The summary line text for the issue tracker [Default '<summary-line>']." )
    clp.add_argument(
      "--new-issue-tracker-file", "-i", dest="newIssueTrackerFile", default="",
      help="File created with the new issue tracker text if"\
        +" specified [Default empty '']." )
    clp.add_argument(
      "--new-tests-with-issue-trackers-file", "-t", dest="newTestsWithIssueTrackersFile",
      default="",
      help="CSV file created with entries for the list of nonpassing tests "\
        +" if specified [Default empty '']." )


  def getCmndLineOptions(self):
    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    clp = ArgumentParser(
      description=self.usageHelp,
      formatter_class=RawDescriptionHelpFormatter)
    self.injectExtraCmndLineOptions(clp)
    self.options = clp.parse_args(sys.argv[1:])
    self.postReadFixupCommandlineOptions()
    self.validateCommandLineOptions()


  def postReadFixupCommandlineOptions(self):
    if (self.options.cdashProjectStartTimeUtc == "") \
        and self.cdashProjectStartTimeUtc \
      :
      self.options.cdashProjectStartTimeUtc = self.cdashProjectStartTimeUtc


  def validateCommandLineOptions(self):
    if self.options.nonpassingTestsUrl == "":
      print("Error, the argument --nopassing-tests-url is required (see --help)")
      sys.exit(1)
    if self.options.cdashProjectStartTimeUtc == "":
      print("Error, the argument --cdash-project-start-time is required (see --help)")
      sys.exit(1)


  def downloadNonpassingTestsData(self):
    apiQueryTestsUrl = self.options.nonpassingTestsUrl.replace(
      "/queryTests.php", "/api/v1/queryTests.php")
    cdashDownloadFileForTesting = \
      os.environ.get('CREATE_ISSUE_TRACKER_FROM_CDASH_QUERY_FILE_FOR_UNIT_TESTING',
      '' )
    nonpassingTestsLOD = \
       CDQAR.downloadTestsOffCDashQueryTestsAndFlatten(apiQueryTestsUrl,
         fullCDashQueryTestsJsonCacheFile=cdashDownloadFileForTesting,
         alwaysUseCacheFileIfExists=True )
    return nonpassingTestsLOD


  def writeNewIssueTrackerFile(self, issueTrackerText):
    if self.options.newIssueTrackerFile:
      print("\nWriting out new issue tracker text to '"\
        +self.options.newIssueTrackerFile+"'")
      with open(self.options.newIssueTrackerFile, 'w') as fileHandle:
        fileHandle.write(issueTrackerText)


  def writeNewTestsWithIssueTrackersFile(self, uniqNonpassingTestsLOD):
    if self.options.newTestsWithIssueTrackersFile:
      print("\nWriting out list of test/build pairs for CSV file '"\
        +self.options.newTestsWithIssueTrackersFile+"'")
      csvFileStruct = CDQAR.writeTestsListOfDictsToCsvFileStructure(uniqNonpassingTestsLOD,
        self.issueTrackerUrlTemplate, self.issueTrackerTemplate )
      with open(self.options.newTestsWithIssueTrackersFile, 'w') as csvFile:
        csvFile.write(CDQAR.writeCsvFileStructureToStr(csvFileStruct))


# Class object for issue tracker data
#
class IssueTrackerData:

  def __init__(self,
      summaryLine=None,
      testingDayStartNonpassingDate=None,
      nonpassingTestsUrl=None,
      uniqNonpassingTestsLOD=None,
      buildnameList=None,
      testnameList=None,
      testHistoryHtmlTableText=None,
    ):
    self.summaryLine = summaryLine
    self.testingDayStartNonpassingDate = testingDayStartNonpassingDate
    self.nonpassingTestsUrl = nonpassingTestsUrl
    self.uniqNonpassingTestsLOD = uniqNonpassingTestsLOD
    self.buildnameList = buildnameList
    self.testnameList = testnameList
    self.testHistoryHtmlTableText = testHistoryHtmlTableText


#
# Nonmember functions
#


def getUniqueTestsListOfDicts(cdashTestsLOD): #LOD == "List of Dicts"
  # Get a dict that has just unique keys 'site', 'buildName' and 'testname' values
  uniqTestDict = {}
  for testDict in cdashTestsLOD:
    site = testDict.get('site')
    buildName = testDict.get('buildName')
    testname = testDict.get('testname')
    uniqTestDict.update( {
      site+"_"+buildName+"_"+testname :
        {'site':site, 'buildName':buildName, 'testname':testname}
      } )
  # Get a flat list of test dicts
  uniqNonpassingTestsLOD= []
  for key, value in dictItems(uniqTestDict):
    uniqNonpassingTestsLOD.append(value)
  uniqNonpassingTestsLOD = CDQAR.sortAndLimitListOfDicts(uniqNonpassingTestsLOD,
    CDQAR.getDefaultTestDictsSortKeyList())
  return uniqNonpassingTestsLOD


def getTestnameAndBuildnameLists(uniqNonpassingTestsLOD):
  testnameSet = set()
  buildnameSet = set()
  for testDict in uniqNonpassingTestsLOD:
    testnameSet.add(testDict.get('testname'))
    buildnameSet.add(testDict.get('buildName'))
  testnameList = extractListAndSort(testnameSet)
  buildnameList = extractListAndSort(buildnameSet)
  return (testnameList, buildnameList) 


def extractListAndSort(setIn):
  listOut = []
  for item in setIn: listOut.append(item)
  listOut.sort()
  return listOut


def getTestingDayStartNonpassingDate(nonpassingTestsLOD, cdashProjectStartTimeUtcStr):
  cdashProjectStartTimeUtcStrTD = CDBTD.getProjectTestingDayStartTimeDeltaFromStr(
    cdashProjectStartTimeUtcStr)
  testingDayStartNonpassingDate = None
  for testDict in nonpassingTestsLOD:
    #print("")
    #g_pp.pprint(testDict)
    buildStartTime = testDict.get('buildstarttime')
    #print("buildStartTime = "+str(buildStartTime))
    testingDate = CDBTD.getTestingDayDateFromBuildStartTimeStr(buildStartTime,
      cdashProjectStartTimeUtcStrTD)
    if (testingDayStartNonpassingDate==None) \
        or (testingDate  < testingDayStartNonpassingDate) \
      :
      testingDayStartNonpassingDate = testingDate
  #print("testingDayStartNonpassingDate = "+str(testingDayStartNonpassingDate))
  return testingDayStartNonpassingDate


#
# Helper functions to be used in concrete drivers
#


def getMarkdownListStr(inputList, formatChars):
  markdownText = ""
  for item in inputList:
    markdownText += "* "+formatChars+item+formatChars+"\n"
  return markdownText
