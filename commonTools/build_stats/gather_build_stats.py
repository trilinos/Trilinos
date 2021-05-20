#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv

from FindTribitsCiSupportDir import *
import GeneralScriptSupport as GSS
import CDashQueryAnalyzeReport as CDQAR

from BuildStatsData import *


#
# Helper functions
#


# Robustly read all CSV build stats *.timing files under a base dir and return
# as a list of dicts (LOD).
#
# This robustly deals with *.timing files and discards any *.timing files that
# have any problems.
#
def readAllValidTimingFiles(baseDir, printErrMsg=True, printStats=False):
  listOfAllTimingFiles = getListOfAllTimingFiles(baseDir)
  if printStats:
    print("Number of *.timing files found = "+str(len(listOfAllTimingFiles)))
  allValidTimingFilesLOD = []
  for timingFile in listOfAllTimingFiles:
    timingFileFullPath = baseDir+"/"+timingFile
    (buildStatsTimingDict, errMsg) = \
      readBuildStatsTimingFileIntoDict(timingFileFullPath)
    if errMsg != "" and printErrMsg:
      print(errMsg)
    if not buildStatsTimingDict == None:
      allValidTimingFilesLOD.append(buildStatsTimingDict)
  if printStats:
    print("Number of valid *.timing files found = "+str(len(allValidTimingFilesLOD)))
  return allValidTimingFilesLOD


# Robustly read a CSV build stats timing file created by magic_wrapper.py and
# return as dict.
#
# Returns tuple:
#
#  (timingBuildStatsDict, errorMsg)
#
# If the timing build stats file 'buildStatsTimingFile' exists and has valid
# data, then 'timingBuildStatsDict' will be a simple dict with the contents of
# the CSV file.  Otherwise, 'timingBuildStatsDict' will be 'None' and 'errMsg'
# will contain the error message.
#
# The provides for a very robust reading of these timing build stats files in
# case there are problems with the running of the magic_wrapper.py tool.
#
def readBuildStatsTimingFileIntoDict(buildStatsTimingFile):

  # Output data initialization
  buildStatsTimingDict = None
  errMsg = ""

  (listOfDicts, errMsg) = robustReadCsvFileIntoListOfDicts(buildStatsTimingFile)

  if errMsg == "" and not len(listOfDicts) == 1:
    errMsg = buildStatsTimingFile+": ERROR: Contains "+\
      str(len(listOfDicts))+" != 1 data rows!"

  if listOfDicts != None and errMsg == "":
     # No errors found to this point, so grab the first row as the build stats dict
     buildStatsTimingDict = listOfDicts[0]

  if buildStatsTimingDict != None and errMsg == "":
    errMsgBody = checkBuildStatsTimingDictHasError(buildStatsTimingDict)
    if errMsgBody != "":
      errMsg = buildStatsTimingFile+": "+errMsgBody

  if buildStatsTimingDict != None and errMsg == "":
    normalizeFileNameFieldInDict(buildStatsTimingDict)

  if errMsg != "":
    buildStatsTimingDict = None

  return (buildStatsTimingDict, errMsg)


# Call readCsvFileIntoListOfDicts() but make robust to basic read errors.
#
# Returns:
#
#   (listOfDicts, errMsg)
#
# Returns a valid list of dicts listOfDicts!=None unless some error occurs.
# If some error occured, then errMsg will be sets to a string descrdibing what
# the problem was.
#
# No exception should ever be thrown from this function.  This is useful to
# use in cases where the existance or basic structure of a CSV file may be
# broken and we want to ignore or gracefully deal with invalid files.
#
def robustReadCsvFileIntoListOfDicts(csvFile):
   listOfDicts = None
   errMsg = ""
   if os.path.exists(csvFile):
     try:
       listOfDicts = CDQAR.readCsvFileIntoListOfDicts(csvFile)
     except Exception as exceptObj:
       if str(exceptObj).find("is empty which is not allowed") != -1:
         errMsg = csvFile+": ERROR: File is empty!"
       else:
         errMsg = csvFile+": ERROR: "+str(exceptObj)
       # NOTE: The above check is tied pretty tighlty to the implementation of
       # readCsvFileIntoListOfDicts() in looking for a specific substring in
       # the error message but it will still capture any other error as well
       # and report it through errMsg.
   else:
     errMsg = csvFile+": ERROR: File does not exist!"
   return (listOfDicts, errMsg)
# ToDo: Move the above function to CsvFileUtils.py!


# Assert that a build stats timing dict contains the required fields and has
# valid data in those required field.
#
# Returns:
#
#   errMsg
#
# Returns errMsg=="" if there is no error.  Otherwise, errMsg describes the
# nature of the error.
#
def checkBuildStatsTimingDictHasError(buildStatsTimingDict):
  errMsg = ""
  for stdBuildStatColAndType in getStdBuildStatsColsAndTypesList():
    requiredFieldName = stdBuildStatColAndType.colName()
    requiredFieldType = stdBuildStatColAndType.colType()
    strVal = buildStatsTimingDict.get(requiredFieldName, None)
    if strVal == None:
      errMsg = "ERROR: The required field '"+requiredFieldName+"' is missing!"
      break
    try:
      convertedVal = stdBuildStatColAndType.convertFromStr(strVal)
    except Exception as exceptObj:
      errMsg = "ERROR: For field '"+requiredFieldName+"' the string value '"+strVal+"'"+\
        " could not be converted to the expected type '"+requiredFieldType+"'!"
  return errMsg


# Normalize the 'FileName' field value
#
def normalizeFileNameFieldInDict(aDict):
  fileName = aDict.get('FileName')
  if fileName.startswith("./"):
    aDict.update({'FileName':fileName[2:]})


# Get list of all *.timing files under baseDir and return paths relative to
# baseDir.
#
# This does not read the contents of any of the timing files, it just returns
# a list of all of them.
#
def getListOfAllTimingFiles(baseDir):
  listOfAllTimingFiles = []
  for root, subdirs, files in os.walk(baseDir):
    if root == baseDir:  relRoot = ""
    else:                relRoot = root.replace(baseDir+"/","")
    for aFile in files:
      if aFile.endswith(".timing"):
        aFileRelPath = os.path.join(relRoot, aFile)
        listOfAllTimingFiles.append(aFileRelPath)
  return listOfAllTimingFiles


# Fill in dict of lists for combined info from a list of dicts
#
# The output dict of lists will have the superset of keys from all of the
# input dicts in the listOfDicts and any non-existent values will be given the
# empty string "" instead of `None`.
#
def getDictOfListsFromListOfDicts(listOfDicts, printStats=False):
  numTotalRows = len(listOfDicts)
  supersetOfFieldNamesList = getSupersetOfFieldNamesList(listOfDicts)
  if printStats:
    print(
      "Combined build-stats keys sorted:\n"+\
      "  "+str(supersetOfFieldNamesList) )
  dictOfLists = {}
  # Create dict of lists with all empty values
  for keyName in supersetOfFieldNamesList:
    dictOfLists.update( { keyName : [""] * numTotalRows } )
  # Fill in the values from the dicts in the list
  for i in range(numTotalRows):
    aDict = listOfDicts[i]
    for key, value in aDict.items():
      dictOfLists.get(key)[i] = value
  # Return the completed data-structure
  return dictOfLists


# Get superset of all of the field names for a list of dicts
#
def getSupersetOfFieldNamesList(listOfDicts):
  supersetOfFieldNames = set()
  for aDict in listOfDicts:
    supersetOfFieldNames = supersetOfFieldNames.union(aDict.keys())
  return sorted(list(supersetOfFieldNames))


# Write a dict of lists to a CSV File
#
# Note, this writes the column names (keys) in sorted order.
#
def writeDictOfListsToCsvFile(dictOfLists, csvFile):
  keysList = sorted(dictOfLists.keys())
  if len(keysList) > 0:
    numTotalRows = len(dictOfLists.get(keysList[0]))  # All lists are same length
  else:
    numTotalRows = 0
  numTotalKeys = len(keysList)
  with open(csvFile, "w") as csvFileHandle:
    csvWriter = csv.writer(csvFileHandle, delimiter=",", lineterminator="\n")
    csvWriter.writerow(keysList)
    for i in range(numTotalRows):
      rowList = []
      for aKey in keysList:
        rowList.append(dictOfLists.get(aKey)[i])
      csvWriter.writerow(rowList)
# ToDo: Move the above function to CsvFileUtils.py!


#
# Helper functions for main()
#


#
# Help message
#


def getRequiredFieldsAndTypesDocStr():
  docStr = ""
  for stdBuildStatColAndType in getStdBuildStatsColsAndTypesList():
    requiredFieldName = stdBuildStatColAndType.colName()
    requiredFieldType = stdBuildStatColAndType.colType()
    docStr += "  "+requiredFieldName+" : "+requiredFieldType+"\n"
  return docStr


usageHelp = r"""

Gather up build stats from *.timing CSV files under the given base directory
created by the magic_wrapper.py tool as a byproduct of building the various
targets in a project.

This will discard the data from any *.timing file that does not have valid
values for the required minimum column headers/fields with types:

"""+getRequiredFieldsAndTypesDocStr()+r"""

or if any other problems are found with a *.timing file.

The column headers in all of the *.timing files are combined into one superset
in the generated 'buildStatsCsvFile' file and the columns are listed in sorted
order.  (The values for any fields missing in a *.timing file are given the null
string ''.)
"""


def injectCmndLineOptionsInParser(clp):

  clp.add_argument(
    "--base-dir", "-d", dest="baseDir", default="",
    help="Base directory for project build dir containing the *.timing files."+\
      " [default is current working directory]" )

  clp.add_argument(
    "--verbose", "-v", dest="verbose", action="store_true", default=False,
    help="Provide verbose output." )

  clp.add_argument("buildStatsCsvFile", nargs='?', default="build_stats.csv",
    help="The build status CSV file to created on output."+\
      "  [default is 'build_stats.csv' in current working directory]" )


def getCmndLineOptions():
  from argparse import ArgumentParser, RawDescriptionHelpFormatter
  clp = ArgumentParser(description=usageHelp,
    formatter_class=RawDescriptionHelpFormatter)
  injectCmndLineOptionsInParser(clp)
  options = clp.parse_args(sys.argv[1:])
  if options.baseDir == "":
    options.baseDir = os.getcwd()
  elif not os.path.exists(options.baseDir):
    print("Error, the base dir '"+options.baseDir+"' does not exist!")
  return options


#
#  Main()
#

if __name__ == '__main__':
  inOptions = getCmndLineOptions()
  if inOptions.verbose:
    print("Reading all *.timing files from under '"+inOptions.baseDir+"' ...")
  allValidTimingFilesListOfDicts = readAllValidTimingFiles(inOptions.baseDir,
    printStats=inOptions.verbose)
  allValidTimingFilesDictOfLists = \
    getDictOfListsFromListOfDicts(allValidTimingFilesListOfDicts,
      printStats=inOptions.verbose)
  writeDictOfListsToCsvFile(allValidTimingFilesDictOfLists,
    inOptions.buildStatsCsvFile)
  if inOptions.verbose:
    print("Wrote file '"+inOptions.buildStatsCsvFile+"'")
