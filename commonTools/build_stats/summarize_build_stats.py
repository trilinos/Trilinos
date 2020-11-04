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


# Read a CSV file of build stats into a dict of lists for just the fields we
# want.
#
# Returns a dict of lists where each key is the column/field name and the
# value is an array of data for that field.
#
def readBuildStatsCsvFileIntoDictOfLists(buildStatusCsvFileName):
  return readCsvFileIntoDictOfLists(buildStatusCsvFileName,
    getStdBuildStatsColsAndTypesList() )


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
    newEntry = round(Decimal(scaleFactor*entry), decimalPlaces) 
    newFieldDataList.append(newEntry)
  dictOfLists.update( {newFieldName : newFieldDataList} )


# Compute summary info about a sinlgle build stat from a dict of list of build
# stats
#
def computeBuildStatusSummaryForOneField(buildStatsDOL, fieldName, decimalPlaces):
  buildStatList = buildStatsDOL[fieldName]
  fileNameList = buildStatsDOL['FileName']
  # Set easy fields
  buildStatSummary = BuildStatSummary(fieldName)
  buildStatSummary.numValues = len(buildStatList)
  buildStatSummary.sumValue = round(Decimal(sum(buildStatList)), decimalPlaces)
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


# Compute and return a list of standard build stats summaries from a dict of
# lists of build stats.
#
def computeStdBuildStatsSummaries(buildStatsDOL):
  buildStatsSummariesList = []
  buildStatsSummariesList.append(
    computeBuildStatusSummaryForOneField(buildStatsDOL, 'max_resident_size_mb', 2) )
  buildStatsSummariesList.append(
    computeBuildStatusSummaryForOneField(buildStatsDOL, 'elapsed_real_time_sec', 2) )
  buildStatsSummariesList.append(
    computeBuildStatusSummaryForOneField(buildStatsDOL, 'file_size_mb', 2) )
  return buildStatsSummariesList


# Create an ASCII text report block for a list of build stats summaries
#
def createAsciiReportOfBuildStatsSummaries(buildStatsSummariesList,
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

#
# Help message
#

usageHelp = r"""summarize_build_stats.py --build-stats-csv-file=<csv-file>

Summarize gathered build stats from the the build stats CSV file and print the
report as ASCII text to STDOUT.  This prints a report like:

Full Project: sum(max_resident_size_size_mb) = ??? (??? entries)
Full Project: max(max_resident_size_size_mb) = ??? (<file-name>)
Full Project: max(elapsed_real_time_sec) = ??? (<file-name>)
Full Project: sum(elapsed_real_time_sec) = ??? (??? entries)
Full Project: sum(file_size_mb) = ??? (??? entries)
Full Project: max(file_size_mb) = ??? (<file-name>)
"""


#
# Helper functions for main()
#


def injectCmndLineOptionsInParser(clp, gitoliteRootDefault=""):
  
  clp.add_option(
    "--build-stats-csv-file", dest="buildStatsCsvFile", type="string", default="",
    help="The build status CSV file created by build wappers and gathered up." )


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
  buildStatsSummariesList = computeStdBuildStatsSummaries(buildStatsDOL)
  buildStatsAsciiReport = createAsciiReportOfBuildStatsSummaries(
    buildStatsSummariesList, "Full Project")

  print(buildStatsAsciiReport)
