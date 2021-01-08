#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
from decimal import Decimal

from FindTribitsCiSupportDir import *
import GeneralScriptSupport as GSS
import CDashQueryAnalyzeReport as CDQAR


#
# Helper functions
#


# Round a float to n decimal places
#
def roundNum(numIn, numDecPlaces):
  return float(round(Decimal(numIn), numDecPlaces))


# Read a CSV file of build stats into a dict of lists for just the fields we
# want.
#
# Returns a dict of lists where each key is the column/field name and the
# value is an array of data for that field.
#
def readBuildStatsCsvFileIntoDictOfLists(buildStatusCsvFileName,
    computeStdScaledFields=True, normalizeFileName=True,
  ):
  buildStatsDOL = readCsvFileIntoDictOfLists(buildStatusCsvFileName,
    getStdBuildStatsColsAndTypesList() )
  if computeStdScaledFields:
    addStdScaledBuildStatsFields(buildStatsDOL)
  if normalizeFileName:
    normalizeFileNameField(buildStatsDOL)
  return buildStatsDOL


# Standard set of build stats fields we want to read in
#
def getStdBuildStatsColsAndTypesList():
  return [
    ColNameAndType('max_resident_size_Kb', 'float'),
    ColNameAndType('elapsed_real_time_sec', 'float'),
    ColNameAndType('FileName', 'string'),
    ColNameAndType('FileSize', 'float'),
    ]
# NOTE: Above, we use type 'float' instead of 'int' for fields that are ints
# because we want to allow a very large size.


# Read in a CSV file as a dict of lists.
#
def readCsvFileIntoDictOfLists(csvFileName, colNameAndTypeList):
  dictOfLists = {}
  with open(csvFileName, 'r') as csvFile:
    csvReader = csv.reader(csvFile)
    # Get the list of col headers and the index to the col headers we want 
    columnHeadersList = \
      CDQAR.getColumnHeadersFromCsvFileReader(csvFileName, csvReader)
    colNameTypeIdxList = \
      getColNameTypeIdxListGivenColNameAndTypeList(csvFileName, columnHeadersList,
        colNameAndTypeList)
    # Initial empty lists for each column to hold the data
    for colNameTypeIdx in colNameTypeIdxList:
      dictOfLists.update( { colNameTypeIdx.colName() : [] } )
    # Fill the columns of data
    dataRowIdx = 0
    for lineList in csvReader:
      if not lineList: continue # Ignore blank line
      CDQAR.stripWhiltespaceFromStrList(lineList)
      assertNumExpectedCsvFileLineEntries(csvFileName, columnHeadersList,
        dataRowIdx, lineList)
      # Read the row entries
      for colNameTypeIdx in colNameTypeIdxList:
        dictOfLists[colNameTypeIdx.colName()].append(
          colNameTypeIdx.convertFromStr(lineList[colNameTypeIdx.getIdx()]) )
      # Update for next row
      dataRowIdx += 1
  # Return completed dict of lists
  return dictOfLists


def assertNumExpectedCsvFileLineEntries(csvFileName, columnHeadersList,
    dataRowIdx, csvLineList,
  ):
  if len(columnHeadersList) != len(csvLineList):
    raise Exception(
      "Error, the CSV file '"+csvFileName+"' has "+str(len(columnHeadersList))+\
      " column headers but data row "+str(dataRowIdx)+" only has "+\
       str(len(csvLineList))+" entries!" )


def getColNameTypeIdxListGivenColNameAndTypeList(csvFileName, columnHeadersList,
    colNameAndTypesToGetList,
  ):
  colNameTypeIdxList = []
  for colNameAndTypeToGet in colNameAndTypesToGetList:
    colIdx = GSS.findInSequence(columnHeadersList, colNameAndTypeToGet.colName())
    if colIdx != -1:
      colNameTypeIdxList.append(ColNameTypeIdx(colNameAndTypeToGet, colIdx))
    else:
      raise Exception(
        "Error, the CSV file column header '"+colNameAndTypeToGet.colName()+"'"+\
        " does not exist in the list of column headers "+str(columnHeadersList)+\
        " from the CSV file '"+csvFileName+"'!")
  return colNameTypeIdxList


class ColNameAndType(object):
  def __init__(self, colName, colType):
    self.__colName = colName
    self.__colType = colType
    self.assertType()
  def colName(self):
    return self.__colName
  def colType(self):
    return self.__colType
  def __repr__(self):
    myStr = "ColNameAndType{"+self.__colName+","+str(self.__colType)+"}"
    return myStr
  def convertFromStr(self, strIn):
    if self.__colType == "string":
      return strIn
    elif self.__colType == "int":
      return int(strIn)
    elif self.__colType == "float":
      return float(strIn)
  def assertType(self):
    supportedTypes = [ "string", "int", "float" ]
    if -1 == GSS.findInSequence(supportedTypes, self.__colType):
      raise Exception(
        "Error, type '"+str(self.__colType)+"' is not supported!  Supported types"+\
        " include "+str(supportedTypes)+"!")
  def __eq__(self, other):
    return((self.__colName,self.__colType)==(other.__colName,other.__colType))
  def __ne__(self, other):
    return((self.__colName,self.__colType)!=(other.__colName,other.__colType))


class ColNameTypeIdx(object):
  def __init__(self, colNameAndType, colIdx):
    self.__colNameAndType = colNameAndType
    self.__colIdx = colIdx
  def colName(self):
    return self.__colNameAndType.colName()
  def getIdx(self):
    return self.__colIdx
  def convertFromStr(self, strIn):
    return self.__colNameAndType.convertFromStr(strIn)
  def __repr__(self):
    myStr = "ColNameTypeIdx{"+str(self.__colNameAndType)+","+str(self.__colIdx)+"}"
    return myStr
  def __eq__(self, other):
    return ((self.__colNameAndType,self.__colIdx)==(other.__colNameAndType,other.__colIdx))
  def __ne__(self, other):
    return ((self.__colNameAndType,self.__colIdx)!=(other.__colNameAndType,other.__colIdx))


# Add standard scaled fields to read-in build stats dict of lists
#
def addStdScaledBuildStatsFields(buildStatsDOL):
  addNewFieldByScalingExistingField(buildStatsDOL, 'max_resident_size_Kb',
    1.0/1024, 2, 'max_resident_size_mb')
  addNewFieldByScalingExistingField(buildStatsDOL, 'FileSize',
    1.0/(1024*1024), 2, 'file_size_mb')


# Scale an existing field to create a new field
#
def addNewFieldByScalingExistingField(dictOfLists, existingFieldName,
    scaleFactor, decimalPlaces, newFieldName,
  ):
  existingFieldDataList = dictOfLists[existingFieldName]
  newFieldDataList = []
  for entry in existingFieldDataList:
    newEntry = roundNum(scaleFactor*entry, decimalPlaces)
    newFieldDataList.append(newEntry)
  dictOfLists.update( {newFieldName : newFieldDataList} )


# Normize the FileName field to remove beginning './'
#
def normalizeFileNameField(dictOfLists):
  fileNameList = dictOfLists.get('FileName')
  i = 0
  while i < len(fileNameList):
    fileName = fileNameList[i]
    if (len(fileName) > 1) and (fileName[0:2] == './'):
      fileNameList[i] = fileName[2:]
    i += 1


# Bin the build stats dict of lists by subdirs under a given list of dirs for
# the 'FileName' field.
#
def binBuildStatsDictOfListsBySubdirUnderDirs(
    buildStatsDOL,
    binBySubdirsUnderDirsList,
  ):
  binnedBuildStatsDOL_dict = {}
  numTotalRows = len(buildStatsDOL.get('FileName'))
  row_idx = 0
  while row_idx < numTotalRows:
    fileName = buildStatsDOL.get('FileName')[row_idx]
    for baseDir in binBySubdirsUnderDirsList:
      if fileName.startswith(baseDir+"/"):
        subdirUnderBaseDir = getSubdirUnderBaseDir(baseDir, fileName)
        fileNameBuildStatsDOL = \
          binnedBuildStatsDOL_dict.setdefault(subdirUnderBaseDir, {})
        addRowFromDictOfListsToDictOfLists(buildStatsDOL, row_idx,
          fileNameBuildStatsDOL)
    row_idx += 1
  #
  return BuildStatsBinnedBySubdirs(buildStatsDOL, binnedBuildStatsDOL_dict)


class BuildStatsBinnedBySubdirs(object):
  def __init__(self, fullBuildStatsDOL, binnedBuildStatsDOL_dict):
    self.fullBuildStatsDOL = fullBuildStatsDOL
    self.binnedBuildStatsDOL_dict = binnedBuildStatsDOL_dict


def getSubdirUnderBaseDir(baseDir, fileName):
  lenBaseDir = len(baseDir)
  positionOfDirCharAfterSubDir = fileName.find("/", lenBaseDir+2)
  subdirUnderBaseDir = fileName[lenBaseDir+1 : positionOfDirCharAfterSubDir]
  return subdirUnderBaseDir


def addRowFromDictOfListsToDictOfLists(inputDOL, row_idx, inoutDOL):
  keysList = inputDOL.keys()
  for key in keysList:
    keyValList = inoutDOL.setdefault(key, [])
    keyValList.append(inputDOL.get(key)[row_idx])


# Compute summary info about a sinlgle build stat from a dict of list of build
# stats
#
def computeBuildStatusSummaryForOneField(buildStatsDOL, fieldName, decimalPlaces):
  buildStatList = buildStatsDOL[fieldName]
  fileNameList = buildStatsDOL['FileName']
  # Set easy fields
  buildStatSummary = BuildStatSummary(fieldName)
  buildStatSummary.numValues = len(buildStatList)
  buildStatSummary.sumValue = roundNum(sum(buildStatList), decimalPlaces)
  # Compute max and the corresponding filename
  maxValue = 0
  maxFileName = ""
  for i in range(buildStatSummary.numValues):
    buildStat = buildStatList[i]
    fileName = fileNameList[i]
    if buildStat > maxValue:
      maxValue = buildStat
      maxFileName = fileName
  buildStatSummary.maxValue = maxValue
  buildStatSummary.maxFileName = maxFileName
  # Return
  return buildStatSummary


class BuildStatSummary(object):
  def __init__(self, fieldName):
    self.fieldName = fieldName
    self.numValues = None
    self.sumValue = None
    self.maxValue = None
    self.maxFileName = None
  def __str__(self):
    return "{"+\
      "fieldName="+self.fieldName+", " + \
      "numValues="+str(self.numValues)+", " + \
      "sumValue="+str(self.sumValue)+", " + \
      "maxValue="+str(self.maxValue)+", " + \
      "maxFileName="+self.maxFileName + \
      "}"


# Compute and return a list of standard build stats summaries from a dict of
# lists of build stats.
#
def computeStdBuildStatsSummariesSingleDOL(buildStatsDOL):
  buildStatsSummariesList = []
  buildStatsSummariesList.append(
    computeBuildStatusSummaryForOneField(buildStatsDOL, 'max_resident_size_mb', 2) )
  buildStatsSummariesList.append(
    computeBuildStatusSummaryForOneField(buildStatsDOL, 'elapsed_real_time_sec', 2) )
  buildStatsSummariesList.append(
    computeBuildStatusSummaryForOneField(buildStatsDOL, 'file_size_mb', 2) )
  return buildStatsSummariesList


# Compute and return the lists of standard build stats summaries for the full
# project as well as those binned by subdirs.
#
def computeStdBuildStatsSummaries(buildStatsBinnedBySubdirs):
  fullBuildStatsSummariesList = \
    computeStdBuildStatsSummariesSingleDOL(buildStatsBinnedBySubdirs.fullBuildStatsDOL)
  binnedBuildStatsSummariesList_dict = {}
  for subdir in buildStatsBinnedBySubdirs.binnedBuildStatsDOL_dict.keys():
    binnedBuildStatsSummariesList_dict.update(
      {
        subdir :
        computeStdBuildStatsSummariesSingleDOL(
          buildStatsBinnedBySubdirs.binnedBuildStatsDOL_dict.get(subdir) )
        }
      )
  return BuildStatsSummariesBinnedBySubdirs(fullBuildStatsSummariesList,
    binnedBuildStatsSummariesList_dict)


class BuildStatsSummariesBinnedBySubdirs(object):
  def __init__(self, fullBuildStatsSummariesList, binnedBuildStatsSummariesList_dict):
    self.fullBuildStatsSummariesList = fullBuildStatsSummariesList
    self.binnedBuildStatsSummariesList_dict = binnedBuildStatsSummariesList_dict


# Create an ASCII text report block for a list of build stats summaries for a
# single list of stats.
#
def createAsciiReportOfBuildStatsSummariesSingleSet(buildStatsSummariesList,
    buildStatsSetName,
  ):
  asciiReportStr = ""
  for buildStatsSummary in buildStatsSummariesList:
    asciiReportStr += createAsciiReportOfOneBuildStatsSummary(buildStatsSummary,
      buildStatsSetName)
  return asciiReportStr


def createAsciiReportOfOneBuildStatsSummary(buildStatsSummary, buildStatsSetName):
  # Shorter names for below
  bss = buildStatsSummary
  bssn = buildStatsSetName
  # Create and return the report str
  asciiReportStr = \
    bssn+": sum("+bss.fieldName+") = "+str(bss.sumValue)+\
      " ("+str(bss.numValues)+" entries)\n"+\
    bssn+": max("+bss.fieldName+") = "+str(bss.maxValue)+" ("+bss.maxFileName+")\n"
  return asciiReportStr



# Create an ASCII text report block for a list of build stats summaries for a
# single list of stats.
#
def createAsciiReportOfBuildStatsSummaries(buildStatsSummariesBinnedBySubdirs):
  asciiReportStr = ""
  asciiReportStr += createAsciiReportOfBuildStatsSummariesSingleSet(
    buildStatsSummariesBinnedBySubdirs.fullBuildStatsSummariesList,
    "Full Project")
  binnedBuildStatsSummariesList_dict = \
    buildStatsSummariesBinnedBySubdirs.binnedBuildStatsSummariesList_dict
  for subdir in sorted(binnedBuildStatsSummariesList_dict.keys()):
    asciiReportStr += "\n"
    asciiReportStr += createAsciiReportOfBuildStatsSummariesSingleSet(
      binnedBuildStatsSummariesList_dict.get(subdir), subdir )
  return asciiReportStr


#
# Help message
#

usageHelp = r"""summarize_build_stats.py --build-stats-csv-file=<csv-file>
  --bin-by-subdirs-under-dirs=<basedir0>,<basedir1>,...

Summarize gathered build stats from the the build stats CSV file and print the
report as ASCII text to STDOUT.  This prints a report like:

Full Project: sum(max_resident_size_mb) = ??? (??? entries)
Full Project: max(max_resident_size_mb) = ??? (<file-name>)
Full Project: max(elapsed_real_time_sec) = ??? (<file-name>)
Full Project: sum(elapsed_real_time_sec) = ??? (??? entries)
Full Project: sum(file_size_mb) = ??? (??? entries)
Full Project: max(file_size_mb) = ??? (<file-name>)

<subdir0>: sum(max_resident_size_mb) = ??? (??? entries)
<subdir0>: max(max_resident_size_mb) = ??? (<file-name>)
<subdir0>: max(elapsed_real_time_sec) = ??? (<file-name>)
<subdir0>: sum(elapsed_real_time_sec) = ??? (??? entries)
<subdir0>: sum(file_size_mb) = ??? (??? entries)
<subdir0>: max(file_size_mb) = ??? (<file-name>)

...
"""


#
# Helper functions for main()
#


def injectCmndLineOptionsInParser(clp, gitoliteRootDefault=""):
  
  clp.add_option(
    "--build-stats-csv-file", dest="buildStatsCsvFile", type="string", default="",
    help="The build status CSV file created by build wappers and gathered up." )
  
  clp.add_option(
    "--bin-by-subdirs-under-dirs", dest="binBySubdirsUnderDirsStr", type="string",
    default="",
    help="List of base dirs to group results by subdir under."+\
      " Format '<basedir0>,<basedir1>,..." )


def getCmndLineOptions():
  from optparse import OptionParser
  clp = OptionParser(usage=usageHelp)
  injectCmndLineOptionsInParser(clp)
  (options, args) = clp.parse_args()
  if options.buildStatsCsvFile == "":
    raise Exception(
      "Error, input argument --build-stats-csv-file must be set!")
  if not os.path.exists(options.buildStatsCsvFile):
    raise Exception(
      "Error, file '"+options.buildStatsCsvFile+"' does not exist!")
  return options


#
#  Main()
# 

if __name__ == '__main__':

  inOptions = getCmndLineOptions()

  buildStatsDOL = readBuildStatsCsvFileIntoDictOfLists(inOptions.buildStatsCsvFile)
  addStdScaledBuildStatsFields(buildStatsDOL)
  buildStatsBinnedBySubdirs = binBuildStatsDictOfListsBySubdirUnderDirs(
    buildStatsDOL, inOptions.binBySubdirsUnderDirsStr.split(',') )
  buildStatsSummariesBinnedBySubdirs = computeStdBuildStatsSummaries(
    buildStatsBinnedBySubdirs )
  buildStatsAsciiReport = createAsciiReportOfBuildStatsSummaries(
    buildStatsSummariesBinnedBySubdirs )

  print(buildStatsAsciiReport)
