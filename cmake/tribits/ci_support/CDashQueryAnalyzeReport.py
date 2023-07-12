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

try:
  # Python 2
  from urllib2 import urlopen
  from urllib2 import quote as urlquote
except ImportError:
  # Python 3
  from urllib.request import urlopen
  from urllib.parse import quote as urlquote

import sys
import hashlib
import json
import datetime
import copy
import pprint
import csv

from FindGeneralScriptSupport import *
from GeneralScriptSupport import *
from Python2and3 import u, csvReaderNext

import cdash_build_testing_date as CBTD


# Accept the --date input option with values 'today', 'yesterday', or some
# 'YYYY-MM-DD' value.
#
def convertInputDateArgToYYYYMMDD(cdashProjectTestingDayStartTime, dateText,
    currentDateTimeStr=None,  # Used for unit testing only
  ):
  if dateText == "yesterday" or dateText == "today":
    if dateText == "yesterday": dayIncr = -1
    else: dayIncr = 0
    dateTime = CBTD.getRelativeCDashBuildStartTimeFromCmndLineArgs(
      currentDateTimeStr, cdashProjectTestingDayStartTime, dayIncr)
    rtnDate = CBTD.getDateOnlyFromDateTime(dateTime)
  else:
    rtnDate = validateAndConvertYYYYMMDD(dateText)
  return rtnDate


# Validate a date YYYY-MM-DD string and return a date object for the
# 'datetime' module.
#
def validateAndConvertYYYYMMDD(dateText):
  try:
    return datetime.datetime.strptime(dateText, '%Y-%m-%d')
  except ValueError:
    raise ValueError("Incorrect data format for '"+dateText+"', should be YYYY-MM-DD")


# Get a file name string from a general text string.
#
# This replaces non-alphanumeric chars with '_'.
#
def getFileNameStrFromText(inputStr):
  fileNameStr = ""
  for char in inputStr:
    if char.isalnum():
      fileNameStr += char
    else:
       fileNameStr += "_"
  return fileNameStr


# Check if the key/value pairs for two dicts are the same and if so, return an
# error message explaining how they are different.
#
# Returns tuple (hasSameKeyValuePairs, errMsg).  If
# hasSameKeyValuePairs==True, then errMsg==None.  Otherwise, if
# hasSameKeyValuePairs==False, then errMsg gives a string that explains how
# they are different.
#
# This improves on a simple check dict_1 == dict_2 in that shows exactly why
# the dicts are different for a single key/value pair.
#
def checkDictsAreSame(dict_1, dict_1_name, dict_2, dict_2_name):
  # Assume all passing unless we find a difference
  hasSameKeyValuePairs = True
  errMsg = None
  # Start with the fast internal Python check
  if dict_1 == dict_2:
    return (True, None)
  # Check if they have the same number of keys
  if hasSameKeyValuePairs and (len(dict_1.keys()) != len(dict_2.keys())):
    hasSameKeyValuePairs = False
    errMsg = "len("+dict_1_name+".keys())="+str(len(dict_1.keys()))+\
      " != len("+dict_2_name+".keys())="+str(len(dict_2.keys()))
  # Check that they have the same key/value pairs
  if hasSameKeyValuePairs:
    for key_1 in dict_1.keys():
      if not key_1 in dict_2.keys():
        hasSameKeyValuePairs = False
        errMsg = dict_1_name+"['"+key_1+"'] does not exist in "+dict_2_name
        break
      keyVal_1 = dict_1[key_1]
      keyVal_2 = dict_2[key_1]
      if keyVal_1 != keyVal_2:
        hasSameKeyValuePairs = False
        errMsg = dict_1_name+"['"+key_1+"'] = '"+str(keyVal_1)+"' != "+\
          dict_2_name+"['"+key_1+"'] = '"+str(keyVal_2)+"'"
        break
    #end for
  #end if
  # Return the final result
  return (hasSameKeyValuePairs, errMsg)


# Compress a long file name to avoid open() error
#
# If the full file name must be shorted and if prefix!="", then it is added to
# the beginning of the shortened filename.  Also, if ext!="", then "."+ext is
# added to the end of the shortened filename.  Otherwise, if inputFileName is
# not too long, then it is returned without modification (i.e. 'prefix' and
# 'ext' are ignored).  NOTE: If 'prefix' and 'ext' are too long, then the
# returned shortened filename may also be too long.
#
# This function should return a shorter unique file name that is platform
# independent.
#
def getCompressedFileNameIfTooLong(inputFileName, prefix="", ext=""):
  maxFileNameLength = 255 # ToDo: Figure out for this system?
  if len(inputFileName) > maxFileNameLength:
    hashObject = hashlib.sha1(str(inputFileName).encode('utf-8'))
    hashStr = hashObject.hexdigest()
    newFileName = prefix+hashObject.hexdigest()
    if ext: newFileName += "." + ext
    return newFileName
  return inputFileName


# Filter and input list and return a list with elements where
# matchFunctor(inputList[i])==True.
#
def getFilteredList(inputList, matchFunctor):
  filteredList = []
  for ele in inputList:
    if matchFunctor(ele): filteredList.append(ele)
  return filteredList


# Filter an input list returning a two lists (matchList, nomatchList) where
# the first list has elements where matchFunctor(inputList[i])==True and the
# second list has elements where matchFunctor(inputList[i])==False.
#
def splitListOnMatch(inputList, matchFunctor):
  #print("\nsplitListOnMatch(): matchFunctor = "+str(matchFunctor))
  matchList = []
  nomatchList = []
  for ele in inputList:
    if matchFunctor(ele): matchList.append(ele)
    else: nomatchList.append(ele)
  return (matchList, nomatchList)


# DECORATOR match functor class that negates the match of a stored functor.
#
class NotMatchFunctor(object):

  # Construct with another functor to negate
  def __init__(self, matchFunctor):
    self.__matchFunctor = matchFunctor

  # Convert to string rep for debugging/etc.
  def __str__(self):
    myStr = "NotMatchFunctor{"+str(self.__matchFunctor)+"}"
    return myStr

  # Negate the matchFunctor
  def __call__(self, item):
    return (self.__matchFunctor(item) == False)


# Apply a functor to transform every element in a list
#
# The object transformFunctor is applied as:
#
#   list_inout[i] = transformFunctor(list_inout[i])
#
# If the elements are small value-type objects, then the assignment is needed.
# However, if the list elements are handled with reference semantics like a
# list [] or a dict {} then really the object is being modified in place and
# the assignment is not needed but it cheap and harmess in that case.
#
# This returns the input list transformed but the return object can be ignored
# because it modifies the input list object's elements in place.
#
def foreachTransform(list_inout, transformFunctor):
  for i in range(len(list_inout)):
    list_inout[i] = transformFunctor(list_inout[i])
  return list_inout


# Remove elements from a list given a list of indexes
#
# This modifies the original list inplace but also returns it.  Therefore, if
# you want to keep the original list, you better create a copy of the base
# list object before passing it in.
#
def removeElementsFromListGivenIndexes(list_inout, indexesToRemoveList_in):
  indexesToRemoveList = copy.copy(indexesToRemoveList_in)
  indexesToRemoveList.sort()
  numRemoved = 0
  for index in indexesToRemoveList:
    del list_inout[index-numRemoved]
    numRemoved += 1
  return list_inout


# Class CsvFileStructure
#
class CsvFileStructure(object):

  def __init__(self, headersList, rowsList):
    self.headersList = headersList
    self.rowsList = rowsList


# Write a CsvFileStructure data to a string
#
def writeCsvFileStructureToStr(csvFileStruct):
  csvFileStr = ", ".join(csvFileStruct.headersList)+"\n"
  for rowFieldsList in csvFileStruct.rowsList:
    csvFileStr += ", ".join(rowFieldsList)+"\n"
  return csvFileStr


########################################
# CDash Specific stuff
########################################


#
# Reporting policy, data, and defaults
#


# Collection of data used to create the final HTML CDash report that is
# updated and queried by various functions.
#
# NOTE: This is put into a class object so that these vars can be updated in
# place when passed to a function.
#
class CDashReportData(object):
  def __init__(self):
    # Gives the final result (assume passing by default)
    self.globalPass = True
    # This is the top of the HTML body
    self.htmlEmailBodyTop = ""
    # This is the bottom of the email body
    self.htmlEmailBodyBottom = ""
    # This var will store the list of data numbers for the summary line
    self.summaryLineDataNumbersList = []
  def reset(self):
    self.globalPass = True
    self.htmlEmailBodyTop = ""
    self.htmlEmailBodyBottom = ""
    self.summaryLineDataNumbersList = []


# Define standard CDash colors
def cdashColorPassed(): return 'green'
def cdashColorFailed(): return 'red'
def cdashColorNotRun(): return 'orange'
def cdashColorMissing(): return 'gray'
# ToDo: Make the above return different colors for a color-blind palette


def getStandardTestsetAcroList():
  return ['twoif', 'twoinr', 'twip', 'twim', 'twif', 'twinr']


# Aggregate info about a test set used for generating the summary and table
# and to determine global pass/fail.
#
# Members:
#
# * testsetAcro: e.g. 'twoif'
# * testsetDescr: e.g. "Tests without issue trackers Failed"
# * testsetTableType: Values: 'nopass', 'pass', 'missing'
# * testsetColor: e.g.  'red', 'green' (whatever is accepted by function
#   colorHtmlText())
# * existanceTriggersGlobalFail: If 'True' and any of tests fall into this
#   test-set category, then it shoulid trigger a global 'False'
class TestsetTypeInfo(object):

  def __init__(self, testsetAcro, testsetDescr, testsetTableType, testsetColor,
      existanceTriggersGlobalFail=True,
    ):
    self.testsetAcro = testsetAcro
    self.testsetDescr = testsetDescr
    self.testsetTableType = testsetTableType
    self.testsetColor = testsetColor
    self.existanceTriggersGlobalFail = existanceTriggersGlobalFail


# Return the TestsetTypeInfo object for the standard types of test sets that get
# their own tables.
#
# testsetArco [in] Acronym for the standard test set (e.g. 'twoif')
#
# testsetColor [in] Gives the color to use for the summary line and the table
# header.  If 'None' is passed in (the default), then a standard color is
# used.  If the empty string is passed in '', then no color will be applied.
#
def getStandardTestsetTypeInfo(testsetAcro, testsetColor=None):
  if testsetAcro == "twoif":
    tsti = TestsetTypeInfo(testsetAcro, "Tests without issue trackers Failed", 'nopass',
      cdashColorFailed())
  elif testsetAcro == "twoinr":
    tsti = TestsetTypeInfo(testsetAcro, "Tests without issue trackers Not Run", 'nopass',
      cdashColorNotRun())
  elif testsetAcro == "twip":
    tsti = TestsetTypeInfo(testsetAcro, "Tests with issue trackers Passed", 'pass',
      cdashColorPassed(), existanceTriggersGlobalFail=False)
  elif testsetAcro == "twim":
    tsti = TestsetTypeInfo(testsetAcro, "Tests with issue trackers Missing", 'missing',
      cdashColorMissing(), existanceTriggersGlobalFail=False)
  elif testsetAcro == "twif":
    tsti = TestsetTypeInfo(testsetAcro, "Tests with issue trackers Failed", 'nopass',
      cdashColorFailed())
  elif testsetAcro == "twinr":
    tsti = TestsetTypeInfo(testsetAcro, "Tests with issue trackers Not Run", 'nopass',
      cdashColorNotRun())
  else:
    raise Exception("Error, testsetAcro = '"+str(testsetAcro)+"' is not supported!")

  if testsetColor != None:
    tsti.testsetColor = testsetColor

  return tsti


# Return the 'status' field from a test dict
#
# Return 'Not Run' if the 'status' field is missing.  (This happens with one
# customer's tests apparently, see SESW-383.)
#
def getTestDictStatusField(testDict):
  return testDict.get('status', 'Not Run')


# Get the Test-set acronym from the fields of a test dict
#
def getTestsetAcroFromTestDict(testDict):
  issueTracker = testDict.get('issue_tracker', None)
  if isTestFailed(testDict) and issueTracker == None:
    return 'twoif'
  if isTestNotRun(testDict) and issueTracker == None:
    return 'twoinr'
  if isTestPassed(testDict) and issueTracker != None:
    return 'twip'
  if isTestMissing(testDict) and issueTracker != None:
    return 'twim'
  if isTestFailed(testDict) and issueTracker != None:
    return 'twif'
  if isTestNotRun(testDict) and issueTracker != None:
    return 'twinr'
  raise Exception(
    "Error, testDict = '"+str(testDict)+"' with fields"+\
    " status = '"+str(testDict.get('status', None))+"' and"+\
    " issue_tracker = '"+str(testDict.get('issue_tracker', None))+"'"+\
    " is not a supported test-set type!")


# Returns True if a test has 'status' 'Passed'
def isTestPassed(testDict):
  return (testDict.get('status', None) == 'Passed')


# Returns True if a test has 'status' 'Failed'
def isTestFailed(testDict):
  return (testDict.get('status', None) == 'Failed')


# Returns True if a test has 'status' 'Not Run'
def isTestNotRun(testDict):
  return (testDict.get('status', None) == 'Not Run')


# Return True if a test is missing
def isTestMissing(testDict):
  status = testDict.get('status', None)
  if status == 'Missing': return True
  if status == 'Missing / Failed': return True
  return False


# Define default test dicts sort order in tables
def getDefaultTestDictsSortKeyList() : return ['testname', 'buildName', 'site']


#
# Implementation functions
#


# Given a CDash query URL PHP page that returns JSON data, return the JSON
# data converged to a Python data-structure.
#
# The returned Python object will be a simple nested set of Python dicts and
# lists.
#
# NOTE: This function can't really be unit tested because it actually gets
# data from CDash.  Therefore, the code below will be structured such that it
# we can avoid getting call it in any automated tests.
#
def extractCDashApiQueryData(cdashApiQueryUrl):
  if sys.version_info < (2,7,5):
    raise Exception("Error: Must be using Python 2.7.5 or newer")
  # NOTE: If we use Python 2.6.6. then the urllib2 function crashes!
  response = urlopen(cdashApiQueryUrl)
  return json.load(response)


# Read a CSV file into a list of dictionaries for each row where the rows of
# the output list are dicts with the column names as keys.
#
# For example, for the CSV file:
#
#   col_0, col_1, col_2
#   val_00, val_01, val_02
#   val_10, val_11, val_12
#
# the returned list of dicts will be:
#
#  [
#    { 'col_0':'val_00', 'col_1':'val_01', 'col_2':'val_02' },
#    { 'col_0':'val_10', 'col_1':'val_11', 'col_2':'val_12' },
#    ]
#
# This function can also allow the user to assert that the included columns
# match a set of required and optional headers.  For example, that above CSV
# file would match:
#
#   requiredColumnHeadersList = [ 'col_0', 'col_1', 'col_2' ]
#
# or:
#
#   requiredColumnHeadersList = [ 'col_0', 'col_1' ]
#   optionalColumnHeadersList = [ 'col_2', 'col_3', ]
#
# The requiredColumnHeadersList and optionalColumnHeadersList argument lists
# are optional.
#
# Also, the columns can be appear in any order as long as they match all of
# the required headers and don't contain any headers not in the list of
# expected headers.
#
def readCsvFileIntoListOfDicts(csvFileName, requiredColumnHeadersList=[],
    optionalColumnHeadersList=[],
  ):
  listOfDicts = []
  with open(csvFileName, 'r') as csvFile:
    csvReader = csv.reader(csvFile)
    columnHeadersList = getColumnHeadersFromCsvFileReader(csvFileName, csvReader)
    assertExpectedColumnHeadersFromCsvFile(csvFileName, requiredColumnHeadersList,
      optionalColumnHeadersList, columnHeadersList)
    # Read the rows of the CSV file into dicts
    dataRow = 0
    for lineList in csvReader:
      if not lineList: continue # Ignore blank line
      stripWhiltespaceFromStrList(lineList)
      assertExpectedNumColsFromCsvFile(csvFileName, dataRow, lineList,
        columnHeadersList)
      # Read the row entries into a new dict
      rowDict = {}
      for j in range(len(columnHeadersList)):
        rowDict.update( { columnHeadersList[j] : lineList[j] } )
      listOfDicts.append(rowDict)
      # Update for next row
      dataRow += 1
  # Return the constructed object
  return listOfDicts


def getColumnHeadersFromCsvFileReader(csvFileName, csvReader):
  try:
    columnHeadersList = csvReaderNext(csvReader)
    stripWhiltespaceFromStrList(columnHeadersList)
    return columnHeadersList
  except StopIteration:
    raise Exception(
      "Error, CSV file '"+csvFileName+"' is empty which is not allowed!"
      )


def csvReaderNext(csvReader):
  if sys.version_info < (3,):
    return csvReader.next()
  else:
    return next(csvReader)


def assertExpectedColumnHeadersFromCsvFile(csvFileName, requiredColumnHeadersList,
    optionalColumnHeadersList, columnHeadersList,
  ):

  if not requiredColumnHeadersList and not optionalColumnHeadersList:
    return  # No expected column headers to assert against!

  requiredAndOptionalHeadersSet = set(requiredColumnHeadersList)
  requiredAndOptionalHeadersSet.update(optionalColumnHeadersList)
  columnHeadersSet = set(columnHeadersList)

  # Assert that each column header is expected
  for colHeader in columnHeadersList:
    if not colHeader in requiredAndOptionalHeadersSet:
      raise Exception(
        "Error, for CSV file '"+csvFileName+"' the"+\
        " column header '"+str(colHeader)+"' is not in the set"+\
        " of required column headers '"+str(requiredColumnHeadersList)+"'"+\
        " or optional column headers '"+str(optionalColumnHeadersList)+"'!"+\
        ""
        )

  # Assert that all of the required headers are present
  for requiredHeader in requiredColumnHeadersList:
    if not requiredHeader in columnHeadersSet:
      raise Exception(
        "Error, for CSV file '"+csvFileName+"' the"+\
        " required header '"+str(requiredHeader)+"' is missing from the"+\
        " set of included column headers '"+str(columnHeadersList)+"'!"+\
        ""
        )


def assertExpectedNumColsFromCsvFile(csvFileName, dataRow, lineList,
    columnHeadersList,
  ):
  if len(lineList) != len(columnHeadersList):
    raise Exception(
      "Error, for CSV file '"+csvFileName+"' the data row"+\
      " "+str(dataRow)+" "+str(lineList)+" has"+\
      " "+str(len(lineList))+" entries which does not macth"+\
      " the number of column headers "+str(len(columnHeadersList))+"!")


def stripWhiltespaceFromStrList(strListInOut):
  for i in range(len(strListInOut)): strListInOut[i] = strListInOut[i].strip()


g_expectedBuildsCsvFileHeadersRequired = \
  ('group', 'site', 'buildname')


def getExpectedBuildsListOfDictsfromCsvFile(expectedBuildsFileName):
  return readCsvFileIntoListOfDicts(expectedBuildsFileName,
    g_expectedBuildsCsvFileHeadersRequired)


def getExpectedBuildsListOfDictsFromCsvFileArg(expectedBuildsFileArg):
  expectedBuildsLOD = []
  if expectedBuildsFileArg:
    expectedBuildsFilenameList = expectedBuildsFileArg.split(",")
    for expectedBuildsFilename in expectedBuildsFilenameList:
      expectedBuildsLOD.extend(
        getExpectedBuildsListOfDictsfromCsvFile(expectedBuildsFilename))
  return expectedBuildsLOD


# Write list of builds from a builds LOD to a CSV file structure meant to
# match the expected builds CSV file.
#
def expectedBuildsListOfDictsToCsvFileStructure(buildsLOD):
  csvFileHeadersList = copy.deepcopy(g_expectedBuildsCsvFileHeadersRequired)
  csvFileRowsList = []
  for buildDict in buildsLOD:
    csvFileRow = (
      buildDict['group'],
      buildDict['site'],
      buildDict['buildname'],
      )
    csvFileRowsList.append(csvFileRow)
  return CsvFileStructure(csvFileHeadersList, csvFileRowsList)


# Write list of builds from a builds LOD to a CSV file meant to match the
# expected builds CSV file.
#
def writeExpectedBuildsListOfDictsToCsvFile(buildsLOD, csvFileName):
  csvFileStruct = expectedBuildsListOfDictsToCsvFileStructure(buildsLOD)
  with open(csvFileName, 'w') as csvFile:
    csvFile.write(writeCsvFileStructureToStr(csvFileStruct))


g_testsWithIssueTrackersCsvFileHeadersRequired = \
  ('site', 'buildName', 'testname', 'issue_tracker_url', 'issue_tracker')


def getTestsWtihIssueTrackersListFromCsvFile(testsWithIssueTrackersFile):
  return readCsvFileIntoListOfDicts(testsWithIssueTrackersFile,
    g_testsWithIssueTrackersCsvFileHeadersRequired)


# Write list of tests from a Tests LOD to a CSV file structure meant to match
# tests with issue trackers CSV file.
#
def writeTestsListOfDictsToCsvFileStructure(testsLOD,
    issueTrackerUrl="", issueTracker="",
  ):
  csvFileHeadersList = copy.deepcopy(g_testsWithIssueTrackersCsvFileHeadersRequired)
  csvFileRowsList = []
  for testDict in testsLOD:
    csvFileRow = (
      testDict['site'],
      testDict['buildName'],
      testDict['testname'],
      issueTrackerUrl,  # issue_tracker_url
      issueTracker,  # issue_tracker
      )
    csvFileRowsList.append(csvFileRow)
  return CsvFileStructure(csvFileHeadersList, csvFileRowsList)


# Write list of tests from a Tests LOD to a CSV file meant to match tests with
# issue trackers CSV file.
#
def writeTestsListOfDictsToCsvFile(testsLOD, csvFileName):
  csvFileStruct = writeTestsListOfDictsToCsvFileStructure(testsLOD)
  with open(csvFileName, 'w') as csvFile:
    csvFile.write(writeCsvFileStructureToStr(csvFileStruct))


# Pretty print a nested Python data-structure to a file
#
# ToDo: Reimplement this to create a better looking set of indented that that
# involves less right-drift and the expense of more vertical space.
#
def pprintPythonDataToFile(pythonData, filePath):
  with open(filePath,'w') as fileObj:
    pp = pprint.PrettyPrinter(stream=fileObj, indent=2)
    pp.pprint(pythonData)


# Get data off CDash and cache it or read from previously cached data
#
# If useCachedCDashData == True, then the file cdashQueryDataCacheFile must
# exist and will be used to get the data instead of calling CDash
#
# If alwaysUseCacheFileIfExists==True and the file cdashQueryDataCacheFile
# already exists, then the file cdashQueryDataCacheFile will be used to get
# the dta instead of calling CDash.
#
# Otherwise, CDash will be called at cdashQueryUrl to get the data and then
# the data will be written to the the file cdashQueryDataCacheFile if
# cdashQueryDataCacheFile != None.
#
# This function can be used to get data off of CDash using any page on CDash
# including cdash/api/v1/index.php, cdash/api/v1/queryTests.php and anything
# other PHP page that returns a JSON data structure (which is all of the
# cdash/api/v1/XXX.php pages).
#
def getAndCacheCDashQueryDataOrReadFromCache(
  cdashQueryUrl,
  cdashQueryDataCacheFile,  # File name
  useCachedCDashData,  # If 'True', then cdasyQueryDataCacheFile must be non-null
  alwaysUseCacheFileIfExists = False,
  verbose = False,
  extractCDashApiQueryData_in=extractCDashApiQueryData,
  ):
  if (
      alwaysUseCacheFileIfExists \
      and cdashQueryDataCacheFile \
      and os.path.exists(cdashQueryDataCacheFile) \
    ):
    if verbose:
      print("  Since the file exists, using cached data from file:\n"+\
        "    "+cdashQueryDataCacheFile )
    with open(cdashQueryDataCacheFile, 'r') as cacheFile:
      cdashQueryData=eval(cacheFile.read())
  elif useCachedCDashData:
    if verbose:
      print("  Using cached data from file:\n    "+cdashQueryUrl )
    with open(cdashQueryDataCacheFile, 'r') as cacheFile:
      cdashQueryData=eval(cacheFile.read())
  else:
    if verbose:
      print("  Downloading CDash data from:\n    "+cdashQueryUrl )
    cdashQueryData = extractCDashApiQueryData_in(cdashQueryUrl)
    if cdashQueryDataCacheFile:
      if verbose:
        print("  Caching data downloaded from CDash to file:\n    "+\
          cdashQueryDataCacheFile)
      pprintPythonDataToFile(cdashQueryData, cdashQueryDataCacheFile)
  return cdashQueryData


def normalizeUrlStrings(*args):
  return [urlquote(x) for x in args]


# Construct full cdash/api/v1/index.php query URL to pull data down given the
# pieces
def getCDashIndexQueryUrl(cdashUrl, projectName, date, filterFields):
  # for legacy reasons, this function assumes we normalized projectName
  projectName, = normalizeUrlStrings(projectName,)
  if date: dateArg = "&date="+date
  else: dateArg = ""
  return cdashUrl+"/api/v1/index.php?project="+projectName+dateArg \
      + "&"+filterFields


# Construct full cdash/index.php browser URL given the pieces
def getCDashIndexBrowserUrl(cdashUrl, projectName, date, filterFields):
  # for legacy reasons, this function assumes we normalized projectName
  projectName, = normalizeUrlStrings(projectName,)
  if date: dateArg = "&date="+date
  else: dateArg = ""
  return cdashUrl+"/index.php?project="+projectName+dateArg \
      + "&"+filterFields


# Construct full cdash/api/v1/queryTests.php query URL given the pieces
def getCDashQueryTestsQueryUrl(cdashUrl, projectName, date, filterFields):
  # for legacy reasons, this function assumes we normalized projectName
  projectName, = normalizeUrlStrings(projectName,)
  if date: dateArg = "&date="+date
  else: dateArg = ""
  cdashTestUrl = cdashUrl+"/api/v1/queryTests.php?project="+projectName+dateArg+"&"+filterFields
  return cdashTestUrl


# Construct full cdash/queryTests.php browser URL given the pieces
def getCDashQueryTestsBrowserUrl(cdashUrl, projectName, date, filterFields):
  # for legacy reasons, this function assumes we normalized projectName
  projectName, = normalizeUrlStrings(projectName,)
  if date: dateArg = "&date="+date
  else: dateArg = ""
  return cdashUrl+"/queryTests.php?project="+projectName+dateArg+"&"+filterFields


# Copy a key/value pair from one dict to another if it eixsts
def copyKeyDictIfExists(sourceDict_in, keyName_in, dict_inout):
  value = sourceDict_in.get(keyName_in, None)
  if value:
    dict_inout.update( { keyName_in : value } )


# Extend the set of fields for a CDash index.phpb build dict
#
# buildDict_in [in]: The build dict gotten from cdash/index.php.  This will be
# modified in place.
#
# Returns the modified build dict.
#
# Change this to get all of the fields and add the 'group' field as well.
#
def extendCDashIndexBuildDict(buildDict_in, groupName):
  buildDict = buildDict_in
  buildDict[u'group'] = groupName
  return buildDict


# Given the full Python JSON data-structure returned from the page
# cdash/api/v1/index.php query from extractCDashApiQueryData(), return a
# flattened-out data-structure that is easier to manipulate.
#
# This function takes in the JSON data-structure (as a nested set of Python
# dicts and listed) directly returned from a query gotten from the page
# cdash/api/v1/index.php with some filters.
#
# The input full CDash index.php JSON data-structure has the following
# structure and fields of interest:
#
#  fullCDashIndexBuildsJson =
#  {
#    'all_buildgroups': [ {'id':1,'name:"Nightly"}, ...],
#    'buildgroups': [
#      {
#        'name':"???",   # group name, e.g. Nightly
#        'builds":[
#          {
#            'site':"???"
#            'buildname':"???",
#            'update': {'errors':???, ...},
#            'configure':{'error': ???, ...},
#            'compilation':{'error': ???, ...},
#            'test': {'fail':???, 'notrun':???, 'pass':???, ...},
#            ...
#            },
#            ...
#          ]
#        },
#        ...
#      ...
#      ]
#      },
#      ...
#    }
#
# This function gets the data from *all* of the collapsed builds and returns
# the flatten-out list of dicts for each build with the 'group' field added in
# as:
#
#   [
#     {
#       'group':"???",
#       'site':"???",
#       'buildname':"???",
#       'update': {'errors':???, ...},
#       'configure':{'error': ???, ...},
#       'compilation':{'error': ???, ...},
#       'test': {'fail':???, 'notrun':???, 'pass':???, ...},
#       ...
#       },
#     ...
#     ]
#
# This collects *all* of the builds from all of the build groups provided by
# that data-structure, not just the 'Nighlty' build group.  Therefore, if you
# want to only consider one set of build groups, you need to add that to the
# CDash query URL (e.g. group='Nighlty').
#
def flattenCDashIndexBuildsToListOfDicts(fullCDashIndexBuildsJson):
  summaryCDashIndexBuilds = []
  for buildgroup in fullCDashIndexBuildsJson["buildgroups"]:
    groupName = buildgroup["name"]
    for build in buildgroup["builds"]:
      summaryBuild = extendCDashIndexBuildDict(build, groupName)
      summaryCDashIndexBuilds.append(summaryBuild)
  return summaryCDashIndexBuilds


# Given the full JSON data-structure returned from the page
# cdash/api/v1/queryTests.php query from extractCDashApiQueryData(), return a
# flattened-out data-structure that is easier to manipulate.
#
# This function takes in the JSON data-structure (as a nested set of Python
# dicts and listed) directly returned from a query gotten from the page
# cdash/api/v1/queryTests.php with some filters.
#
# The input full CDash queryTests.php JSON data-structure has the following
# structure and fields of interest:
#
#  fullCDashQueryTestsJson =
#  {
#    'version':???,
#    'feed_enabled':???,
#    ...
#    'builds': [
#      {
#        'buildName': 'Trilinos-atdm-mutrino-intel-opt-openmp-HSW',
#        'buildSummaryLink': 'buildSummary.php?buildid=4109735',
#        'buildstarttime': '2018-10-29T05:54:03 UTC',
#        'details': 'Completed (Failed)\n',
#        'nprocs': 4,
#        'prettyProcTime': '40s 400ms',
#        'prettyTime': '10s 100ms',
#        'procTime': 40.4,
#        'site': 'mutrino',
#        'siteLink': 'viewSite.php?siteid=223',
#        'status': 'Failed',
#        'statusclass': 'error',
#        'testDetailsLink': 'testDetails.php?test=57925465&build=4109735',
#        'testname': 'Anasazi_Epetra_BKS_norestart_test_MPI_4',
#        'time': 10.1
#        },
#      ...
#      ],
#    ...
#    }
#
# This function gets the data from *all* of the tests and returns the
# flatten-out list of dicts with some additional fields for each test of the
# form:
#
#   [
#     {
#       'buildName': 'Trilinos-atdm-mutrino-intel-opt-openmp-HSW',
#       'buildSummaryLink': 'buildSummary.php?buildid=4109735',
#       'buildstarttime': '2018-10-29T05:54:03 UTC',
#       'details': 'Completed (Failed)\n',
#       'nprocs': 4,
#       'prettyProcTime': '40s 400ms',
#       'prettyTime': '10s 100ms',
#       'procTime': 40.4,
#       'site': 'mutrino',
#       'siteLink': 'viewSite.php?siteid=223',
#       'status': 'Failed',
#       'statusclass': 'error',
#       'testDetailsLink': 'testDetails.php?test=57925465&build=4109735',
#       'testname': 'Anasazi_Epetra_BKS_norestart_test_MPI_4',
#       'time': 10.1,
#       },
#     ...
#     ]
#
# NOTE: This does a shallow copy so any modifications to the returned list and
# dicts will modify the original data-structure fullCDashQueryTestsJson.  If
# that is a problem, then make sure and do a deep copy before passing in
# fullCDashQueryTestsJson.
#
# This collects *all* of the tests from all of the "build" list provided by
# the CDash JSON data-structure.  Therefore, if you want to only consider one
# set of build groups, you need to add that to the CDash query URL
# (e.g. buildName='<build-name>').
#
def flattenCDashQueryTestsToListOfDicts(fullCDashQueryTestsJson):
  testsListOfDicts = []
  for testDict in fullCDashQueryTestsJson['builds']:
    testsListOfDicts.append(testDict)
  return testsListOfDicts


# Create a lookup dict for a list of dicts
#
# listOfDicts [in/out]: List of dict objects that have keys that one will want
# to lookup the dict based on their values.  May have 100% duplicate elements
# removed from the list.
#
# listOfKeys [in]: List of the names of keys in these dicts that are used to
# build a search dict data-structure which is returned from this function.
#
# removeExactDuplicateElements [in]: If True, then dict elements that are 100%
# duplicates and have the exact same key/value pairs will be removed from
# listOfDicts. (default False)
#
# checkDictsAreSame_in [in]: Allows specialization of the check for exact dict
# matches and reporting the differences.  The default value is the function
# checkDictsAreSame().  Any Python object that has the __call__() operator
# function defined that takes those same arguments and returns the same
# outputs as the function checkDictsAreSame() can be passed in.
#
# If listOfDicts has any elements that are 100% complete duplicates with the
# same exact key/value pairs, then the later elements will be removed from the
# list.  But if just the key/value pairs listed in listOfKeys are duplicated
# but one or more of the other key/value pairs is different, then then an
# exception is thrown.
#
# NOTE: This is an implementation function that is used in the class
# SearchableListOfDicts.  Please use that class instead of this raw function.
#
def createLookupDictForListOfDicts(listOfDicts, listOfKeys,
    removeExactDuplicateElements=False, checkDictsAreSame_in=checkDictsAreSame,
  ):
  # Build the lookup dict data-structure. Also, optionally mark any 100%
  # duplicate elements if asked to remove 100% duplicate elements.
  lookupDict = {} ; idx = 0 ; numRemoved = 0 ; duplicateIndexesToRemoveFromList = []
  for dictEle in listOfDicts:
    # Create the structure of recursive dicts for the keys in order
    currentLookupDictRef = lookupDict
    lastLookupDictRef = None
    for key in listOfKeys:
      keyValue = dictEle[key]
      lastLookupDictRef = currentLookupDictRef
      nextLookupDictRef = currentLookupDictRef.setdefault(keyValue, {})
      currentLookupDictRef = nextLookupDictRef
    addEle = True
    # Check to see if this dict has already been added
    if currentLookupDictRef:
      lookedUpDict = currentLookupDictRef.get('dict', None)
      lookedUpIdx = currentLookupDictRef.get('idx', None)
      (hasSameKeyValuePairs, dictDiffErrorMsg) = checkDictsAreSame_in(
        dictEle, "listOfDicts["+str(idx)+"]",
        lookedUpDict, "listOfDicts["+str(lookedUpIdx)+"]" )
      if hasSameKeyValuePairs and removeExactDuplicateElements:
        # This is a 100% duplicate element to one previously added.
        # Therefore, mark this duplicate element to be removed.
        duplicateIndexesToRemoveFromList.append(idx)
        addEle = False
      else:
        raiseDuplicateDictEleException(idx, dictEle, listOfKeys, lookedUpIdx,
          lookedUpDict, dictDiffErrorMsg)
    # Need to go back and reset the dict on the last dict in the
    # data-structure so that modifications to the dicts that are looked up
    # will modify the original list.
    if addEle:
      currentLookupDictRef.update({'dict':dictEle, 'idx':idx-numRemoved})
    else:
      numRemoved += 1
    idx += 1
  # Remove 100% duplicate elements marked above
  removeElementsFromListGivenIndexes(listOfDicts, duplicateIndexesToRemoveFromList)
  return  lookupDict


def raiseDuplicateDictEleException(idx, dictEle, listOfKeys,
    lookedUpIdx, lookedUpDict, dictDiffErrorMsg,
  ):
  raise Exception(
    "Error, The element\n\n"+\
    "    listOfDicts["+str(idx)+"] =\n\n"+\
    "      "+sorted_dict_str(dictEle)+"\n\n"+\
    "  has duplicate values for the list of keys\n\n"+\
    "    "+str(listOfKeys)+"\n\n"+\
    "  with the element already added\n\n"+\
    "    listOfDicts["+str(lookedUpIdx)+"] =\n\n"+\
    "      "+sorted_dict_str(lookedUpDict)+"\n\n"+\
    "  and differs by at least the key/value pair\n\n"+\
    "    "+str(dictDiffErrorMsg))


# Lookup a dict (and optionally also its index location) in a list of dicts
# given a lookup dict returned from createLookupDictForListOfDicts() where the
# key/value pairs match
#
# lookupDict [in]: A dict created by createLookupDictForListOfDicts() given
# the same listOfKeys used in that function.
#
# listOfKeys [in]: List of keys that was used used to create lookupDict.
#
# listOfValues [in]: A list of values for the given list of keys in
# listOfKeys.
#
# alsoReturnIdx [in]: If True, then the index of the located dict in the
# original listOfDicts will be returned as well.  (default False)
#
# If the matching dict is found, then it will be returned as:
#
#   matchingDict = lookupDictGivenLookupDict(...)
#
# If alsoReturnIdx==True, then also the index will be returned as:
#
#   (matchingDict, idx) = lookupDictGivenLookupDict(...)
#
# If the matching dict is not found, then None will be returned or the tuple
# (None, None) if alsoReturnIdx==True.
#
# NOTE: This is an implementation function that is used in the class
# SearchableListOfDicts.  Please use that class instead of this raw function.
#
def lookupDictGivenLookupDict(lookupDict, listOfKeys, listOfValues,
    alsoReturnIdx=False,
  ):
  #print("\nlookupDict = "+str(lookupDict))
  #print("\nlistOfKeys = "+str(listOfKeys))
  #print("\ndictToFind = "+str(dictToFind))
  if len(listOfKeys) != len(listOfValues):
    raise Exception("Error, len(listOfKeys)="+str(len(listOfKeys))+\
    " != len(listOfValues)="+str(len(listOfValues))+" where"+\
    " listOfKeys="+str(listOfKeys)+\
    " and listOfValues="+str(listOfValues)+"!")
  currentSubLookupDict = lookupDict
  idx = 0
  for idx in range(len(listOfValues)):
    key = listOfKeys[idx]
    #print("\nkey = '"+key+"'")
    keyValueToFind = listOfValues[idx]
    #print("keyValueToFind = '"+str(keyValueToFind)+"'")
    #print("currentSubLookupDict = "+str(currentSubLookupDict))
    keyValueLookedUp = currentSubLookupDict.get(keyValueToFind, None)
    #print("keyValueLookedUp = "+str(keyValueLookedUp))
    if not keyValueLookedUp:
      if alsoReturnIdx: return (None, None)
      return None
    currentSubLookupDict = keyValueLookedUp
  if keyValueLookedUp:
    if alsoReturnIdx:
      return (keyValueLookedUp.get('dict'), keyValueLookedUp.get('idx'))
    return keyValueLookedUp.get('dict')
  return None


# Class that encapsulates a list of dicts and an efficient lookup of a dict
# given a list key/value pairs to match
#
# Once created, this object acts like a list of dicts in most cases but also
# contains functions to search for speicfic dicts given a set of key/value
# pairs.
#
# Any modifications to the dicts looked up with this object will edit the
# dicts in the underlying list of dicts.  This therefore makes this class act
# as a type of multi-key lookup dict using the member function
# lookupDictGivenKeyValuesList(['keyval0', 'keyval1', ...]).  This provides a
# handy way to access and edit the underlying dicts that require
# multi-key/value pairs to find them.
#
# NOTE: The key values for the list of keys given in listOfKeys must be
# unique!  If it is not, then an exception will be thrown.
#
class SearchableListOfDicts(object):

  # Constructor
  #
  # listOfDicts [stored, may be modified]: List of dicts that a search
  # data-structure will be created for.
  #
  # listOfKeys [stored, will not be modified]: List of the names of keys in
  # the dicts of listOfDicts that a search data-structure will be created for
  # and defines the set of key/value pairs used to look up up dicts in
  # listOfDicts.
  #
  # removeExactDuplicateElements [in]: If True, then exact duplicate dicts
  # based on the key/value pairs in listOfKeys will be removed from
  # listOfDicts (which is modified in place, not a copy). (default False)
  #
  # keyMapList [in]: Optional list of key names in the input
  # keyValueDictToFind to pull out and used the match the key/value pairs in
  # the listOfKeys.  This allows a mapping from int input dict key names to
  # the output dict key names. (default None)
  #
  # checkDictsAreSame_in [in]: Allows specialization of the check for exact
  # dict matches and reporting the differences.  The default value is the
  # function checkDictsAreSame().  Any Python object that has the __call__()
  # operator function defined that takes those same arguments and returns the
  # same outputs as the function checkDictsAreSame() can be passed in.
  #
  def __init__(self, listOfDicts, listOfKeys,
      removeExactDuplicateElements=False, keyMapList=None,
      checkDictsAreSame_in=checkDictsAreSame,
    ):
    if keyMapList:
      if len(listOfKeys) != len(keyMapList):
        raise Exception("Error, listOfKeys="+str(listOfKeys)+\
          " keyMapList="+str(listOfKeys)+" have different lengths!" )
    self.__listOfDicts = listOfDicts
    self.__listOfKeys = listOfKeys
    self.__keyMapList = keyMapList
    self.__checkDictsAreSame = checkDictsAreSame_in
    self.__lookupDict = createLookupDictForListOfDicts(
      self.__listOfDicts, self.__listOfKeys,
      removeExactDuplicateElements=removeExactDuplicateElements,
      checkDictsAreSame_in=checkDictsAreSame_in)

  # Convert to string rep
  def __str__(self):
    myStr = "SearchableListOfDicts{listOfDicts="+str(self.__listOfDicts)+\
      ", listOfKeys="+str(self.__listOfKeys)+", lookupDict="+str(self.__lookupDict)+"}"
    return myStr

  # Return listOfDicts passed into Constructor
  def getListOfDicts(self):
    return self.__listOfDicts

  # Return listOfKeys passed to Constructor
  def getListOfKeys(self):
    return self.__listOfKeys

  # Return keyMapList passed to Constructor
  def getKeyMapList(self):
    return self.__keyMapList

  # Lookup a dict given a dict with same key/value pairs for keys listed in
  # listOfKeys.
  def lookupDictGivenKeyValueDict(self, keyValueDictToFind, alsoReturnIdx=False):
    if self.__keyMapList:
      keyListToUse = self.__keyMapList
    else:
      keyListToUse = self.__listOfKeys
    keyValuesListToFind = []
    for idx in range(len(keyListToUse)):
      keyValuesListToFind.append(keyValueDictToFind.get(keyListToUse[idx]))
    lookupRtn = self.lookupDictGivenKeyValuesList(keyValuesListToFind, alsoReturnIdx)
    return lookupRtn

  # Lookup a dict given a flat list of values for the keys
  #
  # Must be in same order self.getListOfKeys().
  #
  def lookupDictGivenKeyValuesList(self, keyValuesListToFind, alsoReturnIdx=False):
    lookupRtn = lookupDictGivenLookupDict(self.__lookupDict, self.__listOfKeys,
      keyValuesListToFind, alsoReturnIdx)
    return lookupRtn

  # Functions to allow this to act like a list
  def __len__(self):
    return len(self.__listOfDicts)
  def __getitem__(self, index_in):
    return self.__listOfDicts[index_in]


# Create a SearchableListOfDicts object for a list of builds dicts that allows
# lookups of builds given the keys "group" => "site" => "buildname" :
# build_dict.
def createSearchableListOfBuilds(buildsListOfDicts):
  return SearchableListOfDicts(buildsListOfDicts, ['group', 'site', 'buildname'])


# Create a SearchableListOfDicts object for a list of tests with issue
# trackers that allows lookups of tests given the keys "site" => "buildName"
# => "testname" : test_dict.
def createSearchableListOfTests( testsListOfDicts,
    removeExactDuplicateElements=False,
    checkDictsAreSame_in=checkDictsAreSame,
  ):
  return SearchableListOfDicts(testsListOfDicts, ['site', 'buildName', 'testname'],
    removeExactDuplicateElements=removeExactDuplicateElements,
    checkDictsAreSame_in=checkDictsAreSame_in )


# Create a SearchableListOfDicts object for a list of build dicts allows
# lookups that match the 'site' and 'buildname' fields but uses input for the
# search that are test dicts that have the fields 'site' and 'buildName'.
def createTestToBuildSearchableListOfDicts(buildsLOD,
    removeExactDuplicateElements=False,
  ):
  return SearchableListOfDicts( buildsLOD, ('site', 'buildname'),
    removeExactDuplicateElements=removeExactDuplicateElements,
    keyMapList=('site', 'buildName') )
  # NOTE: The extra keyMapList is needed because CDash used the key name
  # 'buildname' for the build name returned form the cdash/index.php page
  # while it gave the build name the key name 'buildName' for the data
  # returned from cdash/queryTests.php.


# Match functor that returns true if the input dict has key/values that
# matches one of the dicts in the input SearchableListOfDicts.
class MatchDictKeysValuesFunctor(object):

  # Construct with a SearchableListOfDicts object
  def __init__(self, searchableListOfDict):
    self.__searchableListOfDict = searchableListOfDict

  # Convert to string rep for debugging/etc.
  def __str__(self):
    myStr = "MatchDictKeysValuesFunctor{"+str(self.__searchableListOfDict)+"}"
    return myStr

  # Return 'true' if the key/value pairs in dict_in match the key/value pairs
  # in one of the dicts in the searchableListOfDict object.
  def __call__(self, dict_in):
    matchingDict = self.__searchableListOfDict.lookupDictGivenKeyValueDict(dict_in)
    if matchingDict:
      return True
    return False


# Transform functor that adds issue tracker info and URL to an existing test
# dict.
#
# This functor looks up the test based on 'site', 'buildName', and 'testname'
# keys to find the entry in the list of known issues with issue trackers and
# then it copies the issue issue tracker fields to the input/output test dict.
class AddIssueTrackerInfoToTestDictFunctor(object):

  # Construct with a SearchableListOfDicts object that has issue tracker info.
  # This object testsWithIssueTrackersSLOD must have been constructed using
  # the function createSearchableListOfTests() so it will allow lookups based
  # on the 'site', 'buildName', and 'testname' keys.
  def __init__(self, testsWithIssueTrackersSLOD, addEmptyOnNoMatch=True):
    self.__testsWithIssueTrackersSLOD = testsWithIssueTrackersSLOD
    self.__addEmptyOnNoMatch = addEmptyOnNoMatch

  # Lookup the issue tracker info and add it as new key/value pairs to
  # testDict_inout.
  def __call__(self, testDict_inout):
    # Look up the entry for the test tracker info based on the 'site',
    # 'buildName', and 'testname' key/value pairs in testDict_inout.
    matchingDict = \
      self.__testsWithIssueTrackersSLOD.lookupDictGivenKeyValueDict(testDict_inout)
    if matchingDict:
      issue_tracker = matchingDict['issue_tracker']
      issue_tracker_url = matchingDict['issue_tracker_url']
    else:
      if self.__addEmptyOnNoMatch:
        issue_tracker = ""
        issue_tracker_url = ""
      else:
        raise Exception(
         "Error, testDict_inout="+str(testDict_inout)+\
          " does not have an assigned issue tracker!")
    testDict_inout[u'issue_tracker'] = issue_tracker
    testDict_inout[u'issue_tracker_url'] = issue_tracker_url
    return testDict_inout


# Split a list of Test Dicts based on if they match the set of expected builds
# or not.
#
# testsWithIssueTrackersLOD [in]: List of dicts of tests with issue trackers.
#   Here, only the fields 'site', 'buildName', and 'testname' are significant.
#
# expectedBuildsLOD [in]: List of dicts of expected builds with fields 'site'
#   and 'buildname'.  Here, the key/value pairs 'site' and 'buildname' must be
#   unique.  The 'group' field is ignored (because cdash/queryTests.php does
#   not give the 'group' of each test).
#
# Returns the tuple (testsWithIssueTrackersMatchingExpectedBuildsLOD,
# testsWithIssueTrackersNotMatchingExpectedBuildsLOD).
#
def splitTestsOnMatchExpectedBuilds( testsWithIssueTrackersLOD,
    testToExpectedBuildsSLOD,
  ):
  return splitListOnMatch(testsWithIssueTrackersLOD,
    MatchDictKeysValuesFunctor(testToExpectedBuildsSLOD) )


# Check if the list of tests with issue trackers matches the expected builds.
#
# testsWithIssueTrackersLOD [in]: List of dicts of tests with issue trackers.
#   Here, only the fields 'site', 'buildName', and 'testname' are significant.
#
# expectedBuildsLOD [in]: List of dicts of expected builds with fields 'site'
#   and 'buildname'.  Here, the key/value pairs 'site' and 'buildname' must be
#   unique.  The 'group' field is ignored (because cdash/queryTests.php does
#   not give the 'group' of each test).
#
# This returns a tuple (matches, errMsg).  If all of the tests match, then
# 'matches' will be 'True' and errMsg=="".  If one or more of the tests don't
# match the expected builds then 'matches' will be 'False' and 'errMsg' will
# give a message about which tests are missing.
#
def doTestsWithIssueTrackersMatchExpectedBuilds( testsWithIssueTrackersLOD,
    testToExpectedBuildsSLOD,
  ):
  # Get list of tests matching and non-matching expected builds (ignoring
  # matching tests list)
  (_, nomatchingTestsLOD) = \
     splitTestsOnMatchExpectedBuilds(testsWithIssueTrackersLOD,
        testToExpectedBuildsSLOD)
  # Gather up all of the tests that don't match one of the expected builds
  nonmatchingTestsWithIssueTrackersLOD = []
  for testDict in nomatchingTestsLOD:
    nonmatchingTestsWithIssueTrackersLOD.append(
      {'site':testDict['site'], 'buildName':testDict['buildName'],
       'testname':testDict['testname']} )
  # If all tests matched, return True
  if len(nonmatchingTestsWithIssueTrackersLOD) == 0:
    return (True, "")
  # One or more tests did not match so build an error message and return False
  errMsg = \
    "Error: The following tests with issue trackers did not match 'site' and"+\
    " 'buildName' in one of the expected builds:\n"
  for testDict in nonmatchingTestsWithIssueTrackersLOD:
    errMsg += \
      "  {'site'='"+testDict['site']+"'"+\
      ", 'buildName'="+testDict['buildName']+"'"+\
      ", 'testname'="+testDict['testname']+"'}\n"
  return (False, errMsg)


# Extract just the date from the testDict['buildstartdate'] field
def dateFromBuildStartTime(buildStartTime):
  return buildStartTime.split('T')[0]


# Sort list of test history dicts and get statistics
#
# Inputs:
#
#   testHistoryLOD [in]: List of test dicts for the same test.  Neither this
#   list nor its elements are modified in this call.  The base list object is
#   shallow copied before it is sorted and returned.
#
#   currentTestDate [in]: The current testing day (as a string "YYYY-MM-DD").
#   This is needed to define a frame of reference for interpreting if the test
#   is currently 'Passed', 'Failed', 'Not Run', or is 'Missing' (i.e. does not
#   have any test results for current testing date).
#
#   testingDayStartTimeUtc [in]: The CDash project testing day start time
#   "<hh>:<mm>" in UTC.  For example, if the CDash project testing day start
#   time is 6 PM MDT (18:00 MDT), then the testing day start time is "02:00"
#   (UTC) (which is the next calendar day).
#
#   daysOfHistory [in]: Number of days of history that were requested.
#
# Note that len(testHistoryLOD) may be less than daysOfHistory which is
# allowed and handled in function.  Any days in that range missing contribute
# to testHistoryStats['missing_last_x_days'].
#
# Also note that this function will remove any duplicate tests (which seem to
# occur sometimes due to a defect in CDash).
#
# Returns:
#
#   (sortedTestHistoryLOD, testHistoryStats, testStatus)
#
# where:
#
#   sortedTestHistoryLOD: The sorted list of test dicts with most recent dict
#   at the top. New list object with references to the same test dict
#   elements.  Therefore, if the list elements themselves are modified after
#   the function returns, then the elements in the original list testHistoryLOD
#   will be modified as well.
#
#   testHistoryStats: Dict that gives statistics for the test with fields:
#     - 'pass_last_x_days': Number of times test 'Passed'
#     - 'nopass_last_x_days': Number of times the not 'Passed'
#     - 'missing_last_x_days': Number of days there was no test data
#     - 'consec_pass_days': Number of times the test consecutively passed
#     - 'consec_nopass_days': Number of times the test consecutively did not pass
#     - 'consec_missing_days': Number of days test is missing
#     - 'previous_nopass_date': Before current date, the previous nopass date in UTC
#
#   testStatus: The status of the test for the current testing day with values:
#     - 'Passed': Most recent test 'Passed' had date matching curentTestDate
#     - 'Failed': Most recent test 'Failed' had date matching curentTestDate
#     - 'Not Run': Most recent test 'Not Run' had date matching curentTestDate
#     - 'Missing': Most recent test has date before matching curentTestDate
#
def sortTestHistoryGetStatistics(testHistoryLOD,
    currentTestDate, testingDayStartTimeUtc,
    daysOfHistory,
  ):

  # Helper functions
  def incr(testDict, key): testDict[key] = testDict[key] + 1
  def decr(testDict, key): testDict[key] = testDict[key] - 1

  # Initialize outputs assuming no history (i.e. test is missing for
  # alldaysOfHistory of history)
  sortedTestHistoryLOD = []
  testHistoryStats = {
    'pass_last_x_days': 0,
    'nopass_last_x_days': 0,
    'missing_last_x_days': daysOfHistory,
    'consec_pass_days': 0,
    'consec_nopass_days': 0,
    'consec_missing_days': 0,
    'previous_nopass_date': 'None'
    }
  testStatus = "Missing"

  # Return if there is no test history
  if len(testHistoryLOD) == 0:
    testHistoryStats['consec_missing_days'] = daysOfHistory
    return (sortedTestHistoryLOD, testHistoryStats, testStatus)

  # Sort the test history by the buildstarttime (most current date at top)
  sortedTestHistoryLOD = copy.copy(testHistoryLOD)
  sortedTestHistoryLOD.sort(reverse=True, key=DictSortFunctor(['buildstarttime']))

  # Remove duplicate tests from list of dicts
  sortedTestHistoryLOD = getUniqueSortedTestsHistoryListOfDicts(sortedTestHistoryLOD)

  # Get testing day/time helper object
  testingDayTimeObj = CBTD.CDashProjectTestingDay(currentTestDate, testingDayStartTimeUtc)
  currentTestDateDT = testingDayTimeObj.getCurrentTestingDayDateDT()
  #print("currentTestDateDT = "+str(currentTestDateDT))

  # Top (most recent) test history data
  topTestDict = sortedTestHistoryLOD[0]
  #print("topTestDict['buildstarttime'] = "+topTestDict['buildstarttime'])

  # Get the CDash testing date of the most recent test
  topTestDictTestingDayDT = testingDayTimeObj.getTestingDayDateFromBuildStartTimeDT(
    topTestDict['buildstarttime'] )
  #print("topTestDictTestingDayDT ="+str(topTestDictTestingDayDT))

  # testStatus (for this test based on history)
  if topTestDictTestingDayDT == currentTestDateDT:
    testStatus = getTestDictStatusField(topTestDict)
  else:
    testStatus = "Missing"
  #print("testStatus = '"+testStatus+"'")

  # testHistoryStats

  # Set up for counting num of consecutive pass, nopass, or missing

  if testStatus == "Missing":
    # The test is missing so see how many consecutive days that it is missing
    testHistoryStats['consec_missing_days'] = \
      (currentTestDateDT - topTestDictTestingDayDT).days
    # There are no initial consecutive passing or nopassing days
    initialTestStatusHasChanged = True
  else:
    # Count number of consecutive days that test is either passing or
    # nopasssing
    initialTestStatusHasChanged = False

  if testStatus == 'Passed': previousTestStatusPassed = True
  else: previousTestStatusPassed = False

  previousNopassDate = None

  # Loop over test history for each of the kth entries and update quantities
  for pastTestDict_k in sortedTestHistoryLOD:
    pastTestStatus_k = getTestDictStatusField(pastTestDict_k)
    pastTestDateUtc_k = testingDayTimeObj.getTestingDayDateFromBuildStartTimeStr(
      pastTestDict_k['buildstarttime'])
    # Count the initial consecutive streaks for passing and nonpassing
    if (
       (pastTestStatus_k=='Passed') == previousTestStatusPassed \
       and not initialTestStatusHasChanged \
      ):
      # The initial consecutive streak for passing or nonpassing continues!
      if pastTestStatus_k == 'Passed':
        incr(testHistoryStats, 'consec_pass_days')
      else:
        incr(testHistoryStats, 'consec_nopass_days')
    else:
      # The initial consecutive streak has been broken
      initialTestStatusHasChanged = True
    # Count total pass/nopass/missing tests
    decr(testHistoryStats, 'missing_last_x_days') # Test not missing this day!
    if pastTestStatus_k == 'Passed':
      incr(testHistoryStats, 'pass_last_x_days')
    else:
      incr(testHistoryStats, 'nopass_last_x_days')
    # Find most recent previous nopass test date
    if (
        previousNopassDate == None \
        and pastTestDateUtc_k != currentTestDate \
        and pastTestStatus_k != 'Passed' \
      ):
      previousNopassDate = pastTestDateUtc_k
      testHistoryStats['previous_nopass_date'] = previousNopassDate

  # Return the computed stuff
  return (sortedTestHistoryLOD, testHistoryStats, testStatus)


# Get a new list with unique entries from an input sorted list of test dicts.
#
# The returned list is new and does not modify any of the entries in the input
# sorted inputSortedTestHistoryLOD object.
#
def getUniqueSortedTestsHistoryListOfDicts(inputSortedTestHistoryLOD):

  if len(inputSortedTestHistoryLOD) == 0:
    return inputSortedTestHistoryLOD

  uniqueSortedTestHistoryLOD = []

  lastUniqueTestDict = inputSortedTestHistoryLOD[0]
  uniqueSortedTestHistoryLOD.append(lastUniqueTestDict)

  idx = 1
  while idx < len(inputSortedTestHistoryLOD):
    candidateTestDict = inputSortedTestHistoryLOD[idx]

    if not checkCDashTestDictsAreSame(candidateTestDict, "a", lastUniqueTestDict, "b")[0]:
      uniqueSortedTestHistoryLOD.append(candidateTestDict)
      lastUniqueTestDict = candidateTestDict
    # Else, this is duplicate test entry so skip
    idx += 1

  return uniqueSortedTestHistoryLOD


# Extract testid and buildid from 'testDetailsLink' CDash test dict
# field.
def extractTestIdAndBuildIdFromTestDetailsLink(testDetailsLink):
  testDetailsLinkList = testDetailsLink.split('?')
  if (len(testDetailsLinkList) > 1):
    # Older CDash implementations have ?testid=<testid>&buildid=<buildid>
    phpArgsList = testDetailsLinkList[1].split('&')
    testidArgList = phpArgsList[0].split("=")
    buildidArgList = phpArgsList[1].split("=")
    testId = testidArgList[1]
    buildId = buildidArgList[1]
  else:
    # Newer CDash implementations have
    testDetailsLinkList = testDetailsLink.split('/')
    testId = testDetailsLinkList[1]
    buildId = ""
  return (testId, buildId)


# Check if two test dicts returned from CDash are the same, accounting for
# possible CDash defects allowing duplicate tests except for different test
# IDs and small changes in 'time' (strange defects in CDash).
#
# Has the same calling conventions and return value as the function
# checkDictsAreSame().
#
# Returns tuple (hasSameKeyValuePairs, errMsg).  If
# hasSameKeyValuePairs==True, then errMsg==None.  Otherwise, if
# hasSameKeyValuePairs==False, then errMsg gives a string that explains how
# they are different.
#
# This improves on a simple check dict_1 == dict_2 in that shows exactly why
# the dicts are different for a single key/value pair.
#
def checkCDashTestDictsAreSame(testDict_1, testDict_1_name,
    testDict_2, testDict_2_name,
  ):
  # Check the easy case where they are exactly the same
  if testDict_1 == testDict_2:
    return (True, None)
  # Check to see if 'testDetailsLink' is there in both and then check contents
  sameBuildIdDifferentTestIds = False
  if (
      ('testDetailsLink' in testDict_1.keys()) \
      and \
      ('testDetailsLink' in testDict_2.keys()) \
    ):
    (test1d_1, buildid_1) = \
      extractTestIdAndBuildIdFromTestDetailsLink(testDict_1['testDetailsLink'])
    (test1d_2, buildid_2) = \
      extractTestIdAndBuildIdFromTestDetailsLink(testDict_2['testDetailsLink'])
    if (buildid_1 == buildid_2) and (test1d_1 != test1d_2):
      # This is the special case that we are writing this function for!
      sameBuildIdDifferentTestIds = True
  # Set up copy to allow dropping out fields for comparison
  testDict_1_copy = copy.deepcopy(testDict_1)
  testDict_2_copy = copy.deepcopy(testDict_2)
  # If buildIds are the same but the testIds are different, then check the
  # rest of the key/value pairs to determine if they are the same:
  if sameBuildIdDifferentTestIds:
    testDict_1_copy.pop('testDetailsLink', None)
    testDict_2_copy.pop('testDetailsLink', None)
  # Don't require the same test times
  testDict_1_copy.pop('time', None)
  testDict_2_copy.pop('time', None)
  testDict_1_copy.pop('prettyTime', None)
  testDict_2_copy.pop('prettyTime', None)
  testDict_1_copy.pop('procTime', None)
  testDict_2_copy.pop('procTime', None)
  testDict_1_copy.pop('prettyProcTime', None)
  testDict_2_copy.pop('prettyProcTime', None)
  # Don't require identical matching output
  testDict_1_copy.pop('matchingoutput', None)
  testDict_2_copy.pop('matchingoutput', None)
  # Compare what ever fields are left that may be different and just use the
  # standard comparison that will give a good error message for differences.
  return checkDictsAreSame(testDict_1_copy, testDict_1_name,
    testDict_2_copy, testDict_2_name )


# Get the test history CDash cache filename
#
# Note: this takes care of things like having '/' in the test name
#
def getTestHistoryCacheFileName(date, site, buildName, testname, daysOfHistory):
  testHistoryFileName = \
    date+"-"+site+"-"+buildName+"-"+testname+"-HIST-"+str(daysOfHistory)+".json"
  return testHistoryFileName.replace('/', '_')


# Transform functor that computes and add detailed test history to an existing
# test dict so that it can be printed in the table
# createCDashTestHtmlTableStr().
#
# ToDo: Document the fields set by this functor
#
class AddTestHistoryToTestDictFunctor(object):

  # Constructor which takes additional data needed to get the test history and
  # other stuff.
  #
  # By default, this will always read the data from the cache file if that file
  # already exists.
  #
  def __init__(self, cdashUrl, projectName, date, testingDayStartTimeUtc, daysOfHistory,
    testCacheDir, useCachedCDashData=True, alwaysUseCacheFileIfExists=True,
    verbose=False, printDetails=False, requireMatchTestTopTestHistory=True,
    extractCDashApiQueryData_in=extractCDashApiQueryData, # For unit testing
    ):
    self.__cdashUrl = cdashUrl
    self.__projectName = projectName
    self.__date = date
    self.__testingDayStartTimeUtc = testingDayStartTimeUtc
    self.__daysOfHistory = daysOfHistory
    self.__testCacheDir = testCacheDir
    self.__useCachedCDashData = useCachedCDashData
    self.__alwaysUseCacheFileIfExists = alwaysUseCacheFileIfExists
    self.__verbose = verbose
    self.__printDetails = printDetails
    self.__requireMatchTestTopTestHistory = requireMatchTestTopTestHistory
    self.__extractCDashApiQueryData_in = extractCDashApiQueryData_in

  # Get test history off CDash and add test history info and URL to info we
  # find out from that test history
  #
  def __call__(self, testDict):

    #pp = pprint.PrettyPrinter(indent=2)

    #print("\ntestDict:\n")
    #pp.pprint(testDict)

    # Get short names for data inside of this functor
    cdashUrl = self.__cdashUrl
    projectName = self.__projectName
    testDayDate = validateAndConvertYYYYMMDD(self.__date)
    daysOfHistory = self.__daysOfHistory

    # Get basic info about the test from the from the testDict
    site = testDict["site"]
    buildName = testDict["buildName"]
    testname = testDict["testname"]

    # Determine if this test has data from CDash or if it does not
    if testDict.get('buildstarttime', None):
      testAlreadyHasCDashData = True
    else:
      testAlreadyHasCDashData = False

    # Get the date range for CDash queries
    dateRangeBeginDT = testDayDate - datetime.timedelta(days=(daysOfHistory-1))
    dateRangeBeginDateStr = CBTD.getDateStrFromDateTime(dateRangeBeginDT)
    dateRangeEndDateStr = self.__date
    beginEndUrlFields = "begin="+dateRangeBeginDateStr+"&end="+dateRangeEndDateStr

    # normalize names for query
    projectName_url, buildName_url, testname_url, site_url = normalizeUrlStrings(
      projectName, buildName, testname, site)

    # Define queryTests.php query filters for test history
    testHistoryQueryFilters = \
      beginEndUrlFields+"&"+\
      "filtercombine=and&filtercombine=&filtercount=3&showfilters=1&filtercombine=and"+\
      "&field1=buildname&compare1=61&value1="+buildName_url+\
      "&field2=testname&compare2=61&value2="+testname_url+\
      "&field3=site&compare3=61&value3="+site_url

    # URL used to get the history of the test in JSON form
    testHistoryQueryUrl = \
      getCDashQueryTestsQueryUrl(cdashUrl, projectName_url, None, testHistoryQueryFilters)

    # URL to embed in email to show the history of the test to humans
    testHistoryBrowserUrl = \
      getCDashQueryTestsBrowserUrl(cdashUrl, projectName_url, None, testHistoryQueryFilters)

    # URL for to the build summary on index.php page
    buildHistoryEmailUrl = getCDashIndexBrowserUrl(
      cdashUrl, projectName_url, None,
      beginEndUrlFields+"&"+\
      "filtercombine=and&filtercombine=&filtercount=2&showfilters=1&filtercombine=and"+\
      "&field1=buildname&compare1=61&value1="+buildName_url+\
      "&field2=site&compare2=61&value2="+site_url
      )
    # ToDo: Replace this with the the URL to just this one build the index.php
    # page.  To do that, get the build stamp from the list of builds on CDash
    # and then create a URL link for this one build given 'site', 'buildName',
    # and 'buildStamp'.  (NOTE: We can't use 'buildstarttime' without
    # replacing ':' with '%' or the URL will not work with CDash.)

    # Set the names of the cached files so we can check if they exists and
    # write them out otherwise
    testHistoryCacheFileFullName = \
      getTestHistoryCacheFileName(self.__date,site,buildName,testname,daysOfHistory)
    # Possibly compress the file name if it is too long
    testHistoryCacheFilePath = \
     self.__testCacheDir+"/"+\
      getCompressedFileNameIfTooLong(testHistoryCacheFileFullName,self.__date+"-","json")

    if self.__verbose:
      gettingTestHistoryMsg = \
        "Getting "+str(daysOfHistory)+" days of history for "+testname+\
        " in the build "+buildName+" on "+site
      if os.path.exists(testHistoryCacheFilePath):
        gettingTestHistoryMsg += " from cache file"
      else:
        gettingTestHistoryMsg += " from CDash"
      print(gettingTestHistoryMsg)

    # Get the test history off of CDash (or from reading the cache file)
    testHistoryLOD = downloadTestsOffCDashQueryTestsAndFlatten(
      testHistoryQueryUrl, testHistoryCacheFilePath,
      useCachedCDashData=self.__useCachedCDashData,
      alwaysUseCacheFileIfExists=self.__alwaysUseCacheFileIfExists,
      verbose=self.__printDetails,
      extractCDashApiQueryData_in=self.__extractCDashApiQueryData_in
      )

    # Sort and get test history stats and update core testDict fields

    (testHistoryLOD, testHistoryStats, testStatus) = sortTestHistoryGetStatistics(
      testHistoryLOD, self.__date, self.__testingDayStartTimeUtc, daysOfHistory)

    # Assert and update the status

    #print("\ntestStatus = "+str(testStatus))
    #print("\ntestHistoryLOD[0] = "+str(testHistoryLOD[0]))

    if testStatus == "Missing":
      testDict = setTestDictAsMissing(testDict)
    elif testStatus == "Passed":
      testDict.update(testHistoryLOD[0])
      testDict['status_color'] = cdashColorPassed()
    else:
      # If we get here, there should be at least one test dict in
      # testHistoryLOD and this should be a Failed or Not Run test
      # testHistoryLOD[0] should be an exact duplicate of testDict.  The below
      # check confirms that to make sure that CDash is giving us consistent
      # data.
      if self.__requireMatchTestTopTestHistory:
        if testDict.get('status', None) != testStatus:
          raise Exception(
            "Error, test testDict['status'] = '"+str(testDict.get('status',None))+"'"+\
            " != "+\
            "top test history testStatus = '"+testStatus+"'"+\
            " where:\n\n"+\
            "   testDict = "+sorted_dict_str(testDict)+"\n\n"+\
            "   top test history dict = "+sorted_dict_str(testHistoryLOD[0])+"\n\n" )
        if testDict.get('buildstarttime', None) != testHistoryLOD[0]['buildstarttime']:
          raise Exception(
            "Error, testDict['buildstarttime'] = '"+\
            str(testDict.get('buildstarttime',None))+"'"+\
            " != "+\
            "top test history 'buildstarttime' = "+\
            "'"+testHistoryLOD[0]['buildstarttime']+"'"+\
            " where:\n\n"+\
            "   testDict = "+sorted_dict_str(testDict)+"\n\n"+\
            "   top test history dict = "+sorted_dict_str(testHistoryLOD[0])+"\n\n" )
      if testDict.get('status', None) == None and testStatus == "Failed":
        # This is a test missing in the outer list of nonpassing tests but is
        # shown to be failing for the current testing day when looking at the
        # test history.  This can happen when the outer CDash query filters
        # out random system failures (see documentation for option
        # --require-test-history-match-nonpassing-tests).
        testDict = setTestDictAsMissing(testDict)
        testDict.update(testHistoryLOD[0])       # Overwrites 'status' = "Failed"
        testDict['status'] = "Missing / Failed"  # Show this special status!
      elif testStatus == "Failed":
        testDict['status_color'] = cdashColorFailed()
      elif testStatus == "Not Run":
        testDict['status_color'] = cdashColorNotRun()

    # ToDo: Lookup the matching build info so that we can get the buildstamp
    # in order to build a good link to the build on CDash?

    # Get the link to the test details if it exists
    testDetailsLink = testDict.get('testDetailsLink', None)
    if testDetailsLink:
      fullTestDetailsLink = cdashUrl+"/"+testDetailsLink
    else:
      fullTestDetailsLink = None

    # Assign all of the new test dict fields that need to be added
    testDict["site_url"] = ""
    testDict['buildName_url'] = buildHistoryEmailUrl # ToDo: Change to one build
    if fullTestDetailsLink:
      testDict['testname_url'] = fullTestDetailsLink
      testDict['status_url'] = fullTestDetailsLink
    testDict['test_history_num_days'] = daysOfHistory
    testDict['test_history_query_url'] = testHistoryQueryUrl
    testDict['test_history_browser_url'] = testHistoryBrowserUrl
    testDict['test_history_list'] = testHistoryLOD
    testDict.update(testHistoryStats)
    testDict['pass_last_x_days_color'] = cdashColorPassed()
    testDict['pass_last_x_days_url'] = testHistoryBrowserUrl
    testDict['nopass_last_x_days_color'] = cdashColorFailed()
    testDict['nopass_last_x_days_url'] = testHistoryBrowserUrl
    testDict['missing_last_x_days_color'] = cdashColorMissing()
    testDict['missing_last_x_days_url'] = testHistoryBrowserUrl
    testDict['consec_pass_days_color'] = cdashColorPassed()
    testDict['consec_pass_days_url'] = testHistoryBrowserUrl
    testDict['consec_nopass_days_color'] = cdashColorFailed()
    testDict['consec_nopass_days_url'] = testHistoryBrowserUrl
    testDict['consec_missing_days_color'] = cdashColorMissing()
    testDict['consec_missing_days_url'] = testHistoryBrowserUrl

    if testDict.get('status', None) == None:
      print("\ntestStatus = "+testStatus)
      print("\ntestDict:")
      pp = pprint.PrettyPrinter(indent=2)
      pp.pprint(testDict)
      raise Exception("Error, testDict['status']==None for testDict="+str(testDict))

    # Return the updated test dict with the new fields
    return testDict


def setTestDictAsMissing(testDict):
  testDict['status'] = "Missing"
  testDict['status_color'] = cdashColorMissing()
  testDict['details'] = "Missing"
  return testDict


# Transform functor that sets the 'cdash_testing_day' field in a test dict.
#
class AddCDashTestingDayFunctor(object):

  def __init__(self, cdash_testing_day):
    self.cdash_testing_day = cdash_testing_day

  def __call__(self, testDict):
    testDict[u('cdash_testing_day')] = u(self.cdash_testing_day)
    return testDict


# Gather up a list of the missing builds
#
# Inputs:
#
#   buildLookupDict [in]: Lookup dict of build summary dicts gotten off CDash
#
#   expectedBuildsList [in]: List of expected builds dict with fields 'group',
#   'site', and 'buildname'.
#
# Returns an array of dicts of missing expected builds with list elements:
#
#    {'group':"???", 'site':"???", 'buildname':"???", 'status':"???", ...}
#
# where the '...' will be the rest of the fields for builds that exist on CDash
# but don't have full results.
#
# The field 'status' will either be given either:
#
#   "Missing ALL"
#
# or will be:
#
#   "Missing [update], [configure], [build], [tests]"
#
def getMissingExpectedBuildsList(buildsSearchableListOfDicts, expectedBuildsList):
  missingExpectedBuildsList = []
  for expectedBuildDict in expectedBuildsList:
    #print("\nexpectedBuildDict = "+str(expectedBuildDict))
    buildSummaryDict = \
      buildsSearchableListOfDicts.lookupDictGivenKeyValueDict(expectedBuildDict)
    #print("\nbuildSummaryDict = "+str(buildSummaryDict))
    if not buildSummaryDict:
      # No part of the expected build is found!
      missingExpectedBuildDict = copy.deepcopy(expectedBuildDict)
      missingExpectedBuildDict.update({'status':"Missing ALL"})
      #print("missingExpectedBuildDict = "+str(missingExpectedBuildDict))
      missingExpectedBuildsList.append(missingExpectedBuildDict)
    else:
      # This build has some build results so see if all of the results are
      # there or not
      missingPartsList = []
      #if not buildSummaryDict.get('update', None):
      #  missingPartsList.append("update")
      if not buildSummaryDict.get('configure', None):
        missingPartsList.append("configure")
      if not buildSummaryDict.get('compilation', None):
        missingPartsList.append("build")
      if not buildSummaryDict.get('test', None):
        compilationDict = buildSummaryDict.get('compilation', None)
        if compilationDict and compilationDict['time'] == '0s':
          missingPartsList.append("build") # See NOTE below for explanation!
        missingPartsList.append("tests")
      # See if any parts are missing and report if there are as a missing
      # build
      if len(missingPartsList) > 0:
        missingBuildStr = "Missing "+", ".join(missingPartsList)
        missingExpectedBuildDict = copy.deepcopy(expectedBuildDict)
        missingExpectedBuildDict.update({'status':missingBuildStr})
        #print("missingExpectedBuildDict = "+str(missingExpectedBuildDict))
        missingExpectedBuildsList.append(missingExpectedBuildDict)
      else:
        # All parts of the expected build exists so don't report this as an
        # expected build.
        None
  # Return the list of missing expected builds and status
  return missingExpectedBuildsList

  # NOTE: Above uses a heuristic for reporting missing build/compilation data
  # that if the 'compilation' 'time' is '0s' and there is no test data, then
  # we will assume that no build data was submitted to the dashboard.
  # Otherwise, we can't detect when build results are missing because if just
  # 'configure' data is updated to CDash, then CDash will automatically create
  # a 'compilation' subdict and will set the initial time to '0s'.  So we
  # can't directly tell if build has missing build/compilation data or if it
  # just has configure data.  (I consider this to be a defect in CDash and we
  # may ask Kitware to fix this but for now we can just use this heuristic.)
  # But we don't want to assume that a time of '0s' means that no
  # build/compilation data was submitted to CDash since a rebuild with no
  # targets getting built can take under 0s.  However, it is unlikely that the
  # build/compilation time is '0s' and just happens to be missing test
  # results.  In that case, it is likely that the job that was doing the build
  # died and never submitted build results in the first place.  Therefore, we
  # will never report missing build results unless there are also missing test
  # results.  And if there are missing test results, then the build will be
  # listed and missing anyway.  And I don't see a lot of harm in failing to
  # report missing build results if there are test results and all of the
  # tests pass.  The only gap is this logic, therefore, are expected builds
  # that have zero tests defined and submitted 0 tests to CDash (so the 'test'
  # subdict exists).  So this is not perfect but I think this is the best we
  # can do until Kitware fixes CDash to not create a 'compilation' subdict
  # unless build data is actually submitted to CDash.

  # NOTE: Above, we are skipping checking for missing 'update' results for now
  # because we want to allow some builds to not have update results.  In the
  # future, it may be better to allow an expected build to say that it does
  # not have expected 'update' results.


# Download set of builds from CDash builds and return flattened list of dicts
#
# The cdash/api/v1/index.php query selecting the set of builds is provided by
# cdashIndexBuildsQueryUrl.
#
# If cdashIndexBuildsQueryCacheFile != None, then the raw JSON data-structure
# downloaded from CDash will be written to the file
# cdashIndexBuildsQueryCacheFile or read from that file if
# useCachedCDashData==True.
#
# If alwaysUseCacheFileIfExists==True, then if the file
# cdashIndexBuildsQueryCacheFile already exists, it will always be read to get
# data instead of communicating with CDash even if useCachedCDashData==False.
#
# The list of builds pulled off of CDash is flattended and extracted using the
# function flattenCDashIndexBuildsToListOfDicts().
#
# NOTE: The optional argument extractCDashApiQueryData_in is used in unit
# testing to avoid calling CDash.
#
def downloadBuildsOffCDashAndFlatten(
    cdashIndexBuildsQueryUrl,
    fullCDashIndexBuildsJsonCacheFile=None,
    useCachedCDashData=False,
    alwaysUseCacheFileIfExists = False,
    verbose=True,
    extractCDashApiQueryData_in=extractCDashApiQueryData,
  ):
  # Get the query data
  fullCDashIndexBuildsJson = getAndCacheCDashQueryDataOrReadFromCache(
    cdashIndexBuildsQueryUrl, fullCDashIndexBuildsJsonCacheFile, useCachedCDashData,
    alwaysUseCacheFileIfExists, verbose=verbose,
    extractCDashApiQueryData_in=extractCDashApiQueryData_in )
  # Get trimmed down set of builds
  buildsListOfDicts = \
    flattenCDashIndexBuildsToListOfDicts(fullCDashIndexBuildsJson)
  return buildsListOfDicts


# Download set of tests from cdash/api/v1/ctest/queryTests.php and return
# flattened list of dicts
#
# cdashQueryTestsUrl [in]: String URL for cdash/api/v1/ctest/queryTests.php
# with filters.
#
# If verbose==True, the the CDash query URL will be printed to STDOUT.
# Otherwise, this function is silent and will not return any output to STDOUT.
#
# If fullCDashQueryTestsJsonCacheFile != None, then the raw JSON
# data-structure will be written to that file.
#
# If useCachedCDashData==True, then data will not be pulled off of CDash and
# instead the list of builds will be read from the file cdashQueryCacheFile
# which must already exist from a prior call to this function (mostly for
# debugging and unit testing purposes).
#
# If alwaysUseCacheFileIfExists==True, then if the file
# cdashIndexBuildsQueryCacheFile already exists, it will always be read to get
# data instead of communicating with CDash even if useCachedCDashData==False.
#
# The list of tests pulled off CDash is flattended and returned by the
# function flattenCDashQueryTestsToListOfDicts().
#
# NOTE: The optional argument extractCDashApiQueryData_in is used in unit
# testing to avoid calling CDash.
#
def downloadTestsOffCDashQueryTestsAndFlatten(
    cdashQueryTestsUrl,
    fullCDashQueryTestsJsonCacheFile=None,
    useCachedCDashData=False,
    alwaysUseCacheFileIfExists = False,
    verbose=True,
    extractCDashApiQueryData_in=extractCDashApiQueryData,
  ):
  # Get the query data
  fullCDashQueryTestsJson = getAndCacheCDashQueryDataOrReadFromCache(
    cdashQueryTestsUrl, fullCDashQueryTestsJsonCacheFile, useCachedCDashData,
    alwaysUseCacheFileIfExists, verbose=verbose,
    extractCDashApiQueryData_in=extractCDashApiQueryData_in )
  # Get flattened set of tests
  testsListOfDicts = \
    flattenCDashQueryTestsToListOfDicts(fullCDashQueryTestsJson)
  return testsListOfDicts


# Returns True if a build has configure failures
def buildHasConfigureFailures(buildDict):
  configureDict = buildDict.get('configure', None)
  if configureDict and configureDict['error'] > 0:
    return True
  return False


# Returns True if a build has compilation/build failures
def buildHasBuildFailures(buildDict):
  compilationDict = buildDict.get('compilation', None)
  if compilationDict and compilationDict['error'] > 0:
    return True
  return False


# Functor class to sort a row of dicts by multiple columns of string data.
class DictSortFunctor(object):
  def __init__(self, sortKeyList):
    self.sortKeyList = sortKeyList
  def __call__(self, dict_in):
    sortKeyStr=""
    for key in self.sortKeyList:
      keyData = dict_in.get(key)
      if sortKeyStr:
        sortKeyStr += "-"+str(keyData)
      else:
        sortKeyStr = keyData
    return sortKeyStr


# Sort and limit a list of dicts
#
# Arguments:
#
# listOfDicts [in]: List of dicts that will be sorted according to keys.
#
# sortKeyList [in]: List of dict keys that define the sort order for the data
# in the list.  The default is None which means that no sort is performed.
#
# limitRowsToDisplay [in]: The max number of rows to display.  The default is
# None which will result in no limit to the number of rows displayed.  The top
# limitRowsToDisplay items will be displayed after the list is sorted.
#
def sortAndLimitListOfDicts(listOfDicts, sortKeyList = None,
   limitRowsToDisplay = None\
  ):
  # Sort the list
  if sortKeyList:
    listOfDictsOrdered = copy.copy(listOfDicts)  # Shallow copy
    listOfDictsOrdered.sort(key=DictSortFunctor(sortKeyList))
  else:
    listOfDictsOrdered = listOfDicts  # No sort being done
  # Limit rows
  if limitRowsToDisplay == None:
    listOfDictsLimited = listOfDictsOrdered
  else:
    listOfDictsLimited = listOfDictsOrdered[0:limitRowsToDisplay]
  # Return the final sorted limited list
  return listOfDictsLimited


# Create a final summary line of global passfail
#
# cdashReportData [in]: Report data of type CDashReportData
#
# buildsetName [in]: The name of the set of builds to report on
#
# date [in]: Date string in format "YYYY-MM-DD"
#
def getOverallCDashReportSummaryLine(cdashReportData, buildsetName, date):
  if cdashReportData.globalPass:
    summaryLine = "PASSED"
  else:
    summaryLine = "FAILED"

  if cdashReportData.summaryLineDataNumbersList:
    summaryLine += " (" + ", ".join(cdashReportData.summaryLineDataNumbersList) + ")"

  summaryLine += ": "+buildsetName+" on "+date

  return summaryLine


################################################################################
# HTML Support Code
################################################################################


def getFullCDashHtmlReportPageStr(cdashReportData, pageTitle="", pageStyle="",
    detailsBlockSummary=None,
  ):

  htmlPage = \
    "<html>\n\n"

  if pageStyle:
    htmlPage += \
      "<head>\n"+\
      pageStyle+\
      "</head>\n\n"

  htmlPage += \
    "<body>\n\n"

  if pageTitle:
    htmlPage += \
      "<h2>"+pageTitle+"</h2>\n\n"

  htmlPage += \
    cdashReportData.htmlEmailBodyTop+\
    "\n"

  if detailsBlockSummary:
    htmlPage += \
      "<details>\n\n"+\
      "<summary><b>"+detailsBlockSummary+":</b> (click to expand)</b></summary>\n\n"
  htmlPage += \
    cdashReportData.htmlEmailBodyBottom+\
    "\n"

  if detailsBlockSummary:
    htmlPage += \
      "</details>\n\n"

  htmlPage += \
    "</body>\n\n"+\
    "</html>\n"

  return htmlPage


def getDefaultHtmlPageStyleStr():
  return \
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
    "</style>\n"


# Class to store dict key and table header
class TableColumnData(object):

  # Class data
  validColAlignList=["left","right","center","justify","char"]

  # Constructor
  def __init__(self, colHeader, dictKey, colAlign="left"):
    self.colHeader = colHeader
    self.dictKey = dictKey
    if not colAlign in self.validColAlignList:
      raise Exception(
        "Error, colAlign="+colAlign+" not valid.  Please choose from"+\
        " the list ['" + "', '".join(validColAlignList) + "']!" )
    self.colAlign = colAlign


# Color HTML text supported color
def colorHtmlText(htmlText, color_in):
  if color_in == None or color_in == "":
    return htmlText
  elif color_in == "red":
    None # Okay!
  elif color_in == "green":
    None # Okay!
  elif color_in == "gray":
    None # Okay!
  elif color_in == "orange":
    None # Okay!
  else:
    raise Exception("Error, color='"+color_in+"' is invalid."+\
      "  Only 'red', 'green', 'gray' and 'orange' are supported!")
  return("<font color=\""+color_in+"\">"+htmlText+"</font>")


# Add soft word breaks for '_' chars and at other places to allow word wrap
def addHtmlSoftWordBreaks(text_in):
  text_out = text_in.replace('_', '_&shy;')
  return text_out


# Create an HTML table string from a list of dicts and column headers
#
# Arguments:
#
# tableTitle [in]: String for the name of the table included at the top of the
# table.
#
# colDataList [in]: List of TableColumnData objects where
#   colDataList[j].dictKey gives the name of the key for that column of data,
#   colDataList[j].colHeader is the text name for the column header and
#   colDataList[j].colAlign gives the HTML alignment.  That columns in the
#   table will listed in the order given in this list.
#
# rowDataList [in]: List of dicts that provide the data from the table.  The
#   dict in each row must have the keys specified by colData[j].dictKey.  In
#   addition, if (key_url=rowDataList[i].get(colData[j].dictKey+"_url",#
#   None))!=None, then the table entry will be an HTML link <a
#   href="dataRowList[i].get(key_url)">dataRowList[i].get(key)</a>.
#
# htmlStyle [in]: The HTML style data (between <style></style>.  If None is
# passed in then a default style is provided internally.  NOTE: The default
# table style uses CSS formatting for boarders but also sets the <table>
# 'boarder' property since some email clients like Gmail ignore the CSS style
# sections.  To not set a style at all, pass in the empty string "" (not
# None).
#
# htmlTableStyle [in]: The style for the HTML table used in <table
#   style=htmlTableStyle>.  If set to None, then a default style is used.  To
#   not set a style, pass in the empty string "" (not None).
#
# This will also put in soft work breaks for chars like '_' to allow for
# compressing the produced tables.
#
def createHtmlTableStr(tableTitle, colDataList, rowDataList,
    htmlStyle=None, htmlTableStyle=None \
  ):

  htmlStr = ""

  # style options for the table
  defaultHtmlStyle=\
    "<style>table, th, td {\n"+\
    "  padding: 5px;\n"+\
    "  border: 1px solid black;\n"+\
    "  border-collapse: collapse;\n"+\
    "}\n"+\
    "tr:nth-child(even) {background-color: #eee;}\n"+\
    "tr:nth-child(odd) {background-color: #fff;}\n"+\
    "</style>"
  if htmlStyle == "": htmlStyleUsed = ""
  elif htmlStyle != None: htmlStyleUsed = htmlStyle
  else: htmlStyleUsed = defaultHtmlStyle
  if htmlStyleUsed:
    htmlStr+=htmlStyleUsed+"\n"

  # Table title and <table style=...>
  htmlStr+="<h3>"+tableTitle+"</h3>\n"
  if htmlTableStyle != None: htmlTableStyleUsed = htmlTableStyle
  else: htmlTableStyleUsed = "style=\"width:100%\" boarder=\"1\""
  htmlStr+="<table "+htmlTableStyleUsed+">\n\n"

  # Column headings:
  htmlStr+="<tr>\n"
  for colData in colDataList:
    htmlStr+="<th>"+colData.colHeader+"</th>\n"
  htmlStr+="</tr>\n\n"

  # Rows for the table
  row_i = 0
  for rowData in rowDataList:
    htmlStr+="<tr>\n"
    col_j = 0
    for colData in colDataList:
      dictKey = colData.dictKey
      # Get the raw entry for this column
      entry = rowData.get(dictKey, None)
      if entry == None:
        raise Exception(
          "Error, column "+str(col_j)+" dict key='"+colData.dictKey+"'"+\
          " row "+str(row_i)+" entry is 'None' which is not allowed!\n\n"+\
          "Row dict = "+str(rowData))
      # Add soft word breaks to allow line breaks for table compression
      entry = addHtmlSoftWordBreaks(str(entry).strip())
      # Add color if defined for this field
      entryColor = rowData.get(dictKey+"_color", None)
      if entryColor:
        entry = colorHtmlText(entry, entryColor)
      # See if the _url key also exists
      entry_url = rowData.get(dictKey+"_url", None)
      # Set the text for this row/column entry with or without the hyperlink
      if entry_url:
        entryStr = "<a href=\""+entry_url+"\">"+str(entry)+"</a>"
      else:
        entryStr = entry
      # Set the row entry in the HTML table
      htmlStr+=\
        "<td align=\""+colData.colAlign+"\">"+entryStr+"</td>\n"
      col_j += 1
    htmlStr+="</tr>\n\n"
    row_i += 1

  # End of table
  htmlStr+="</table>\n\n"  # Use two newlines makes for good formatting!
  return(htmlStr)


# Get string for table title for CDash data to display
#
# Arguments:
#
# Arguments:
#
# dataTitle [in]: Name of the data category.
#
# dataCountAcronym [in]: Acronym for the type of data being displayed
# (e.g. 'twoi' for "Tests With Out issue trackers").  This is printed in the
# table title in the form dataCoutAcronym=len(rowDataList).
#
# numItems [in]: The number of items of data
#
def getCDashDataSummaryHtmlTableTitleStr(dataTitle, dataCountAcronym, numItems,
    limitRowsToDisplay=None,
  ):
  tableTitle = dataTitle
  if limitRowsToDisplay:
    tableTitle += " (limited to "+str(limitRowsToDisplay)+")"
  tableTitle += ": "+dataCountAcronym+"="+str(numItems)
  return tableTitle


# Create an html table string for CDash summary data.
#
# Arguments:
#
# dataTitle [in]: Name of the data that we be included in the table title.
#
# dataCountAcronym [in]: Acronym for the type of data being displayed
# (e.g. 'twoi' for "Tests With Out issue trackers").  This is printed in the
# table title in the form dataCoutAcronym=len(rowDataList).
#
# colDataList [in]: List of TableColumnData objects where
#   colDataList[j].dictKey gives the name of the key for that column of data,
#   colDataList[j].colHeader is the text name for the column header, and
#   colDataList[j].colAlign gives the HTML alignment.
#   The columns in the table will listed in the order given in this list.
#
# rowDataList [in]: List of dicts that provide the data from the table.  The
#   dict in each row must have the keys specified by colData[j].dictKey.
#
# sortKeyList [in]: List of dict keys that define the sort order for the data
# in the list.  The default is None which means that no sort is performed.
#
# limitRowsToDisplay [in]: The max number of rows to display.  The default is
# None which will result in no limit to the number of rows displayed.  The top
# limitRowsToDisplay items will be displayed after the list is sorted.
#
# htmlStyle [in]: The HTML style data (between <style></style>.  If None is
# passed in then a default style is provided internally (see
# createHtmlTableStr().
#
# htmlTableStyle [in]: The style for the HTML table used in <table
#   style=htmlTableStyle>.  The default is None in which case a default is
#   picked by createHtmlTableStr(().
#
# NOTE: If len(rowDataList) == 0, then the empty string "" is returned.
#
def createCDashDataSummaryHtmlTableStr( dataTitle, dataCountAcronym,
    colDataList, rowDataList, sortKeyList=None, limitRowsToDisplay=None,
    htmlStyle=None, htmlTableStyle=None, titleColor=None,
  ):
  # If no rows, don't create a table
  if len(rowDataList) == 0:
    return ""
  # Sort the list and limit the list
  rowDataListDisplayed = sortAndLimitListOfDicts(
    rowDataList, sortKeyList, limitRowsToDisplay)
  # Table title
  tableTitle = colorHtmlText(
    getCDashDataSummaryHtmlTableTitleStr(
      dataTitle, dataCountAcronym, len(rowDataList), limitRowsToDisplay ),
    titleColor )
  # Create and return the table
  return createHtmlTableStr( tableTitle,
    colDataList, rowDataListDisplayed, htmlStyle, htmlTableStyle )


# Create a tests HTML table string
#
# testsetTypeInfo [in]: Information about the testset of type TestsetTypeInfo
#
# testTypeCountNum [in]: Number of total items for the test type, before
# limiting (e.g. 25)
#
# testsLOD [in]: List of dicts of the test data typically first first
# downloaded from CDash.  Each dict in this list must also have been operated
# on by the functors AddIssueTrackerInfoToTestDictFunctor and
# AddTestHistoryToTestDictFunctor in order to have all of the data needed to
# print in this table.
#
# daysOfHistory [in]: Number of days of test history being displayed.  This is
# needed for one of the table column headers.  (ToDo: Remove this and get this
# from the data).
#
# limitRowsToDisplay [in]: Limit of the number of rows to display.  If this
# limited then this argument is needed in order to print "(limited it ???)" in
# the table title.  Should be 'None' if this listing is not limited. (default
# None)
#
# htmlStyle [in]: HTML style for the entire table (see createHtmlTableStr())
# (default None)
#
# htmlTableStyle [in]: Style inside of <table ... > (see createHtmlTableStr())
# (default None)
#
def createCDashTestHtmlTableStr(testsetTypeInfo, testTypeCountNum, testsLOD,
    limitRowsToDisplay=None, htmlStyle=None, htmlTableStyle=None,
  ):
  # Return empty string if no tests
  if len(testsLOD) == 0:
    return ""
  # Table title
  tableTitle = colorHtmlText(
    getCDashDataSummaryHtmlTableTitleStr(
      testsetTypeInfo.testsetDescr, testsetTypeInfo.testsetAcro,
      testTypeCountNum, limitRowsToDisplay ),
    testsetTypeInfo.testsetColor )
  # Consecutive nopass/pass/missing column
  consecCol = getCDashTestHtmlTableConsecColData(testsetTypeInfo.testsetTableType)
  # Get daysOfHistory out of the data
  daysOfHistory = testsLOD[0]['test_history_num_days']
  # Create column headers
  tcd = TableColumnData
  testsColDataList = [
    tcd("Site", "site"),
    tcd("Build Name", "buildName"),
    tcd("Test Name", "testname"),
    tcd("Status", "status"),
    tcd("Details", "details"),
    consecCol,
    tcd("Non-pass Last "+str(daysOfHistory)+" Days", 'nopass_last_x_days', "right"),
    tcd("Pass Last "+str(daysOfHistory)+" Days", 'pass_last_x_days', "right"),
    tcd("Issue Tracker", "issue_tracker", "right"),
    ]
  # Return the HTML table
  return createHtmlTableStr( tableTitle,
    testsColDataList, testsLOD,
    htmlStyle=htmlStyle, htmlTableStyle=htmlTableStyle )


def getCDashTestHtmlTableConsecColData(testsetTableType):
  tcd = TableColumnData
  if testsetTableType == 'nopass':
    consecCol = tcd("Consec&shy;utive Non-pass Days", 'consec_nopass_days', 'right')
  elif testsetTableType == 'pass':
    consecCol = tcd("Consec&shy;utive Pass Days", 'consec_pass_days', 'right')
  elif testsetTableType == 'missing':
    consecCol = tcd("Consec&shy;utive Missing Days", 'consec_missing_days', 'right')
  else:
    raise Exception("Error, invalid testsetTableType="+str(testsetTableType))
  return consecCol


# Class to generate the data for an HTML report for all test-sets for a given
# issue tracker.
#
# This class can be reused for multiple issue issue trackers.
#
class IssueTrackerTestsStatusReporter(object):

  # Constructor
  #
  # cdashReportData [persisting]: Data used to create the final report (of type
  # CDashReportData).
  #
  def __init__(self, verbose=True):
    self.cdashReportData = CDashReportData()
    self.testsetsReporter = TestsetsReporter(self.cdashReportData, htmlStyle="",
      verbose=verbose )
    self.issueTracker = None
    self.cdashTestingDay = None


  # Generate a report about the status of all of the tests for one issue
  # tracker
  #
  # Inputs:
  #
  #   testsLOD [in]: The list of test dicts for one issue tracker.  (The isuse
  #   tracker info will be extracted from the test dicts field 'issue_tracker'
  #   and all of the test dicts must have the same value for 'issue_tracker'
  #   this will throw.)
  #
  # Return 'True' if all of the tests match an (internally defined) passing
  # criteria such that the issue tracker could be closed.  Otherwise, returns
  # 'False' which means that the tests have not yet mt the passing criteria.
  # If len(testsLOD) == 0, then 'True' will be returned (which assumes that
  # there are no tests remaining related to the issue tracker).
  #
  # Postconditions:
  #
  # * A report about the status of the tests is returned in the function
  #   self.getIssueTrackerTestsStatusReport().
  #
  def reportIssueTrackerTestsStatus(self, testsLOD):
    self.cdashReportData.reset()
    if (len(testsLOD) == 0):
      self.issueTracker = None
      return True
    (self.issueTracker, _) = getIssueTrackerFieldsAndAssertAllSame(testsLOD)
    self.cdashTestingDay = testsLOD[0]['cdash_testing_day']
    self.testsetsReporter.reportTestsets(testsLOD)
    # Return the final status
    return False  # ToDo: Add logic to verify if issue can be closed!


  # Generate a pass/fail report HTML string for the last call to
  # reportIssueTrackerTestsStatus()
  def getIssueTrackerTestsStatusReport(self):
    if self.issueTracker == None:
      return None
    testsSummaryTitle = \
      "Test results for issue "+self.issueTracker+" as of "+self.cdashTestingDay
    return self.testsetsReporter.getTestsHtmlReportStr(testsSummaryTitle,
      detailsBlockSummary="Detailed test results")


# Get the issue tracker fields from a list of test dicts and assert they are
# all the same.
#
# Formal Parameters:
#
#   testsLOD [in]: List of test dicts the must have the 'issue_tracker' and
#   the 'issuer_tracker_url' fields and they must all be identical to each
#   other.
#
# Returns:
#
#   (issue_tracker, issue_tracker_url)
#
# Throws:
#
#   IssueTrackerFieldError: If the 'issue_tracker' or the 'issuer_tracker_url'
#   fields are missing from a test dict or if they don't match each other.
#
def getIssueTrackerFieldsAndAssertAllSame(testsLOD):
  if len(testsLOD) == 0:
    return None
  issue_tracker = None
  issue_tracker_url = None
  i = 0
  for testDict in testsLOD:
    (issue_tracker_i, issue_tracker_url_i) = \
      getIssueTrackerFieldsFromTestDict(testDict, i)
    if issue_tracker == None and issue_tracker_url == None:
      issue_tracker = issue_tracker_i
      issue_tracker_url = issue_tracker_url_i
    else:
      assertSameIssueTracker(issue_tracker_i, issue_tracker, testDict, i)
      assertSameIssueTrackerUrl(issue_tracker_url_i, issue_tracker_url,
        testDict, i)
    i += 1
  return (issue_tracker, issue_tracker_url)


def getIssueTrackerFieldsFromTestDict(testDict, idx):
    issue_tracker = testDict.get('issue_tracker', None)
    issue_tracker_url = testDict.get('issue_tracker_url', None)
    if issue_tracker == None:
      raise IssueTrackerFieldError(
        "Error, the test dict "+sorted_dict_str(testDict)+" at index "+str(idx)+\
        " is missing the 'issue_tracker' field!" )
    if issue_tracker_url == None:
      raise IssueTrackerFieldError(
        "Error, the test dict "+sorted_dict_str(testDict)+" at index "+str(idx)+\
        " is missing the 'issue_tracker_url' field!" )
    return (issue_tracker, issue_tracker_url)


def assertSameIssueTracker(issue_tracker, issue_tracker_expected, testDict, idx):
  if issue_tracker != issue_tracker_expected:
    raise IssueTrackerFieldError(
      "Error, the test dict "+sorted_dict_str(testDict)+" at index "+str(idx)+\
      " has a different 'issue_tracker' field '"+str(issue_tracker)+"' than the"+\
      " expected value of '"+str(issue_tracker_expected)+"'!" )


def assertSameIssueTrackerUrl(issue_tracker_url, issue_tracker_url_expected,
    testDict, idx,
  ):
  if issue_tracker_url != issue_tracker_url_expected:
    raise IssueTrackerFieldError(
      "Error, the test dict "+sorted_dict_str(testDict)+" at index "+str(idx)+\
      " has a different 'issue_tracker_url' field '"+str(issue_tracker_url)+"' than the"+\
      " expected value of '"+str(issue_tracker_url_expected)+"'!" )


class IssueTrackerFieldError(Exception):
  pass


# Class to report a single build-set.
#
# NOTE: The reason this is a class is that the cdashReportData and
# addTestHistoryStrategy objects are set once and are used for multiple calls
# to reportSingleBuildset().
#
class SingleBuildsetReporter(object):

  # Constructor
  #
  def __init__(self, cdashReportData,
      htmlStyle=None, htmlTableStyle=None,
      verbose=True,
    ):
    self.cdashReportData = cdashReportData
    self.htmlStyle = htmlStyle
    self.htmlTableStyle = htmlTableStyle
    self.verbose = verbose
    self.groupSiteBuildNameSortOrder = ['group', 'site', 'buildname']

  # Report on a given build-set and write info to self.cdashReportData
  #
  # Input arguments:
  #
  #   buildsetDescr [in] Description for the set of builds
  #
  #   buildsetAcro [in] Short acronym for the build set
  #
  #   buildsetLOD [in] List of builds in this build set
  #
  #   buildsetGlobalPass [in] If set to False and len(buildsetLOD) > 0, then
  #   global pass is set to False (see below).
  #
  #   buildsetColor [in] Color used for the text (see cdashColorXXX() values).
  #
  #   buildsetColDataList [in] List of TableColumnData entries for each of the
  #   columns to include in the table.  Default is 'None' which gives a
  #   default set of entries.
  #
  #   verbose [in] If set to True then some more verbose info is printed to
  #   STDOUT.  If False, then nothing is printed to STDOUT (which is useful
  #   for unit testing). The default is True.
  #
  # On output, self.cdashReportData data will be updated with the summary and
  # table of this given build-set.  In particular, the following
  # cdashReportData fields will be written to:
  #
  #   cdashReportData.summaryLineDataNumbersList: List will be appended with
  #   lines for each build-set in order called.
  #
  #   cdashReportData.htmlEmailBodyTop: The name of the table 'buildsetDescr',
  #   the acronym 'buildsetAcro' and the size will be written on one line
  #   ending with ``<br>\n``.
  #
  #   cdashReportData.htmlEmailBodyBottom: Summary HTML table (with title)
  #   will be written, along with formatting.
  #
  #   cdashReportData.globalPass: Set to False if buildsetGlobalPass==True and
  #   len(buildsetLOD) > 0.
  #
  def reportSingleBuildset(self, buildsetDescr, buildsetAcro, buildsetLOD,
      buildsetGlobalPass, buildsetColor, buildsetColDataList=None, verbose=True,
    ):

    buildsetNum = len(buildsetLOD)

    buildsetSummaryStr = \
      getCDashDataSummaryHtmlTableTitleStr(buildsetDescr,  buildsetAcro,
        buildsetNum)

    if self.verbose:
      print(buildsetSummaryStr)

    if buildsetNum > 0:

      if not buildsetGlobalPass:
        self.cdashReportData.globalPass = False

      self.cdashReportData.summaryLineDataNumbersList.append(
        buildsetAcro+"="+str(buildsetNum))

      self.cdashReportData.htmlEmailBodyTop += \
        colorHtmlText(buildsetSummaryStr,buildsetColor)+"<br>\n"

      if not buildsetColDataList:
        tcd = TableColumnData
        buildsetColDataList = [
          tcd("Group", 'group'),
          tcd("Site", 'site'),
          tcd("Build Name", 'buildname'),
          ]

      self.cdashReportData.htmlEmailBodyBottom += \
        createCDashDataSummaryHtmlTableStr(
          buildsetDescr,  buildsetAcro, buildsetColDataList, buildsetLOD,
          sortKeyList=self.groupSiteBuildNameSortOrder,
          titleColor=buildsetColor)


# Class to optionally get test history and then analyze and report a single
# test-set.
#
# NOTE: The reason this is a class is that the cdashReportData and
# addTestHistoryStrategy objects are set once and are used for multiple calls
# to reportSingleTestset().
#
class SingleTestsetReporter(object):

  # Constructor
  #
  # cdashReportData [persisting]: Data used to create the final report (of type
  # CDashReportData).
  #
  # testDictsSortKeyList [persisting]: The sort order array tests dicts (input
  # to DictSortFunctor()). (Default getDefaultTestDictsSortKeyList()).
  #
  # addTestHistoryStrategy [persisting]: Strategy object that can set the test
  # history on a list of dicts.  Must have member function
  # getTestHistory(testsLOD).  (Default 'None')
  #
  def __init__(self, cdashReportData,
      testDictsSortKeyList=getDefaultTestDictsSortKeyList(),
      addTestHistoryStrategy=None,
      htmlStyle=None, htmlTableStyle=None,
      verbose=True,
    ):
    self.cdashReportData = cdashReportData
    self.testDictsSortKeyList = testDictsSortKeyList
    self.addTestHistoryStrategy = addTestHistoryStrategy
    self.htmlStyle = htmlStyle
    self.htmlTableStyle = htmlTableStyle
    self.verbose = verbose

  # Report on a given test-set and write info to self.cdashReportData
  #
  # On output, self.cdashReportData data will be updated with the summary and
  # table of this given testsetTypeInfo data.  In particular, the following
  # cdashReportData fields will be written to:
  #
  #   cdashReportData.summaryLineDataNumbersList: List will be appended with
  #   the entry ``testsetTypeInfo.testsetAcro+"="+testsetTotalSize``.
  #
  #   cdashReportData.htmlEmailBodyTop: The name of the table from
  #   testsetTypeInfo.testsetDescr, the acronym testsetTypeInfo.testsetAcro and the
  #   size will be written on one line ending with ``<br>\n``.
  #
  #   cdashReportData.htmlEmailBodyBottom: Summary HTML table (with title)
  #   will be written, along with formatting.
  #
  def reportSingleTestset(self, testsetTypeInfo, testsetTotalSize, testsetLOD,
      sortTests=True,
      limitTableRows=None,   # Change to 'int' > 0 to limit table rows
      getTestHistory=False,
    ):

    testsetSummaryStr = \
      getCDashDataSummaryHtmlTableTitleStr(testsetTypeInfo.testsetDescr,
        testsetTypeInfo.testsetAcro, testsetTotalSize)

    if self.verbose:
      print("")
      print(testsetSummaryStr)

    if testsetTotalSize > 0:

      if testsetTypeInfo.existanceTriggersGlobalFail:
        self.cdashReportData.globalPass = False

      self.cdashReportData.summaryLineDataNumbersList.append(
        testsetTypeInfo.testsetAcro+"="+str(testsetTotalSize))

      self.cdashReportData.htmlEmailBodyTop += \
        colorHtmlText(testsetSummaryStr, testsetTypeInfo.testsetColor)+"<br>\n"

      if sortTests or limitTableRows:
        testsetSortedLimitedLOD = sortAndLimitListOfDicts(
          testsetLOD, self.testDictsSortKeyList, limitTableRows )
      else:
        testsetSortedLimitedLOD = testsetLOD

      if getTestHistory and self.addTestHistoryStrategy:
        self.addTestHistoryStrategy.getTestHistory(testsetSortedLimitedLOD)

      self.cdashReportData.htmlEmailBodyBottom += createCDashTestHtmlTableStr(
        testsetTypeInfo, testsetTotalSize, testsetSortedLimitedLOD,
        limitRowsToDisplay=limitTableRows,
        htmlStyle=self.htmlStyle, htmlTableStyle=self.htmlTableStyle )


# Class to generate the data for an HTML report for all test-sets represented
# in a list of test dicts.
#
class TestsetsReporter(object):

  # Constructor
  #
  # cdashReportData [persisting]: Data used to create the final report (of type
  # CDashReportData).
  #
  def __init__(self, cdashReportData, testsetAcroList=getStandardTestsetAcroList(),
      htmlStyle=None, htmlTableStyle=None, verbose=True,
    ):
    self.cdashReportData = cdashReportData
    self.testsetAcroList = testsetAcroList
    self.singleTestsetReporter = SingleTestsetReporter(cdashReportData,
      htmlStyle=htmlStyle, htmlTableStyle=htmlTableStyle, verbose=verbose)


  # Separate out and report on all of the test-sets in the input list of test
  # dicts.
  #
  def reportTestsets(self, testsLOD):
    testDictsByTestsetAcro = binTestDictsByTestsetAcro(testsLOD)
    for testsetAcro in self.testsetAcroList:
      testsetLOD = testDictsByTestsetAcro.get(testsetAcro, None)
      if testsetLOD:
        testsetTypeInfo = getStandardTestsetTypeInfo(testsetAcro)
        self.singleTestsetReporter.reportSingleTestset(testsetTypeInfo,
          len(testsetLOD), testsetLOD)
  # ToDo: Modify the above to assert that there are no unexpected test-set
  # types!


  # Generate a pass/fail report for all of the testset analyzed by
  # reportTestsets().
  def getTestsHtmlReportStr(self, testsSummaryTitle, detailsBlockSummary=None):
    cdashReportData = self.cdashReportData
    cdashReportData.htmlEmailBodyTop = \
      "<p>\n" + cdashReportData.htmlEmailBodyTop + "</p>\n"
    return getFullCDashHtmlReportPageStr(
      cdashReportData, pageTitle=testsSummaryTitle, pageStyle="",
      detailsBlockSummary=detailsBlockSummary)
    # NOTE: Above, cdashReportData is a shallow copy of some fields like lists
    # and other objects but the string fields are deep copied.  Therefore,
    # this function does not change the state of 'self' at all.!


# Bin a list of test dicts by issue tracker
#
# Input arguments:
#
#   testsLOD [in]: Tests list of dicts
#
# Returns (testDictsByIssueTracker, testsWithoutIssueTrackersLOD)
#
#  testDictsByIssueTracker [out]: A dict where the keys are the issue tracker
#  strings (e.g. '#1234') and the values are the sublist test dicts that have
#  that issue tracker.
#
#  testsWithoutIssueTrackersLOD [out]: The remaining list of test dicts that
#  don't have an issue tracker.
#
# For example, the input list of dicts:
#
#  testsLOD = [
#    { 'testname':'test1' ... 'issue_tracker':'#1234' },
#    { 'testname':'test2' ... },
#    { 'testname':'test3' ... 'issue_tracker':'#1235' },
#    { 'testname':'test4' ... 'issue_tracker':'#1234' },
#    { 'testname':'test5' ... },
#    { 'testname':'test6' ... 'issue_tracker':'#1235' },
#    { 'testname':'test7' ... 'issue_tracker':'#1236' },
#
# would yield:
#
#   testDictsByIssueTracker = {
#     '#1234' : [
#       { 'testname':'test1' ... 'issue_tracker':'#1234' },
#       { 'testname':'test4' ... 'issue_tracker':'#1234' },
#       ],
#     '#1235' : [
#       { 'testname':'test3' ... 'issue_tracker':'#1235' },
#       { 'testname':'test6' ... 'issue_tracker':'#1235' },
#       ],
#     '#1236' : [
#       { 'testname':'test7' ... 'issue_tracker':'#1236' },
#       ],
#     ]
#
#   testsWithoutIssueTrackersLOD = [
#    { 'testname':'test2' ... },
#    { 'testname':'test5' ... },
#    ]
#
def binTestDictsByIssueTracker(testsLOD):
  testDictsByIssueTracker = {}
  testsWithoutIssueTrackersLOD = []
  for testDict in testsLOD:
    issueTracker = testDict.get('issue_tracker', None)
    if issueTracker:
      issueTrackerBinTestsLOD = testDictsByIssueTracker.setdefault(issueTracker, [])
      issueTrackerBinTestsLOD.append(testDict)
    else:
      testsWithoutIssueTrackersLOD.append(testDict)
  return (testDictsByIssueTracker, testsWithoutIssueTrackersLOD)


# Bin a list of test dicts based on their test-set acronym
#
# Input arguments:
#
#   testsLOD [in]: Tests list of dicts
#
# Returns testDictsByTestsetAcro
#
#  testDictsByTestsetAcro [out]: A dict where the keys are the standard
#  test-set acronyms ('twoif', 'twoinr', 'twip', etc) and the value for each
#  is the list of test dicts that below to that test-set category.
#
# For example, the input list of dicts:
#
#   testsLOD = [
#     { 'testname':'test1', 'status':'Failed' ... 'issue_tracker':'#1234' },
#     { 'testname':'test2', 'status':'Failed' ... },
#     { 'testname':'test3', 'status':'Passed' ... 'issue_tracker':'#1235' },
#     { 'testname':'test4', 'status':'Failed' ... 'issue_tracker':'#1234' },
#     { 'testname':'test5', 'status':'Not Run' ... },
#     { 'testname':'test6', 'status':'Missing' ... 'issue_tracker':'#1235' },
#     { 'testname':'test7', 'status':'Not Run' ... 'issue_tracker':'#1236' },
#
# would yield:
#
#   testDictsByTestsetAcro = {
#     'twoif' : [
#       { 'testname':'test2', 'status':'Failed' ... },
#       ],
#     'twoinr' : [
#       { 'testname':'test5', 'status':'Not Run' ... },
#       ],
#     'twip' : [
#       { 'testname':'test3', 'status':'Passed' ... 'issue_tracker':'#1235' },
#       ],
#     'twim' : [
#       { 'testname':'test6', 'status':'Missing' ... 'issue_tracker':'#1235' },
#       ],
#     'twif' : [
#       { 'testname':'test1', 'status':'Failed' ... 'issue_tracker':'#1234' },
#       { 'testname':'test4', 'status':'Failed' ... 'issue_tracker':'#1234' },
#       ],
#     'twinr' : [
#     { 'testname':'test7', 'status':'Not Run' ... 'issue_tracker':'#1236' },
#       ],
#
def binTestDictsByTestsetAcro(testsLOD):
  testDictsByTestsetAcro = {}
  for testDict in testsLOD:
    testsetAcron = getTestsetAcroFromTestDict(testDict)
    testsetAcronBinTestsLOD = testDictsByTestsetAcro.setdefault(testsetAcron, [])
    testsetAcronBinTestsLOD.append(testDict)
  return testDictsByTestsetAcro


#########################################
# HTML Email stuff
#########################################


#
# Create an HTML MIME Email
#

import smtplib

from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.utils import formatdate
from email.charset import Charset, QP
import email.message


# Create MINE formatted email object (but don't send it)
#
def createHtmlMimeEmail(fromAddress, toAddress, subject, textBody, htmlBody,
    doRemoveSoftHyphens=False \
  ):

  # Create message container - the correct MIME type is multipart/alternative.
  msg = MIMEMultipart('alternative')
  msg['From'] = fromAddress
  msg['To'] = toAddress
  msg['Date'] = formatdate(localtime=True)
  msg['Subject'] = subject
  msg['Content-Type'] = "text/html; charset=utf-8"

  # RFC 821 states that the text line size including <CRLF> is 1000 characters.
  # This text line size is causing spaces to be inserted into the email body and
  # break long hyperlinks. Below, we use quoted-printable
  # "Content-Transfer-Encoding" to ensure integrity of the data being sent
  # through the gateway when html text lines are longer than 1000 characters.
  msg['Content-Transfer-Encoding'] = "quoted-printable"

  # Remove hidden soft hyphens from htmlBody so triagers can easily copy paste
  # long build and test names into CDash queries.
  if doRemoveSoftHyphens:
    htmlBody = htmlBody.replace('&shy;', '')

  # Record the MIME types of both parts - text/plain and text/html.
  part1 = MIMEText(textBody.encode('utf-8'), 'plain', 'utf-8')
  part2 = MIMEText(htmlBody.encode('utf-8'), 'html', 'utf-8')

  # Attach parts into message container.  According to RFC 2046, the last part
  # of a multipart message, in this case the HTML message, is best and
  # preferred.
  msg.attach(part1)
  msg.attach(part2)

  return msg


# Send a MIME formatted email
#
def sendMineEmail(mimeEmail):
  # Send the message via local SMTP server.
  s = smtplib.SMTP('localhost')
  # sendmail function takes 3 arguments: sender's address, recipient's address
  # and message to send - here it is sent as one string.
  s.sendmail(mimeEmail['From'], mimeEmail['To'], mimeEmail.as_string())
  s.quit()
