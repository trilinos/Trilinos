# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

import os
import sys
import copy
import shutil
import unittest
import pprint

thisScriptsDir = os.path.dirname(os.path.abspath(__file__))
g_testBaseDir = thisScriptsDir
sys.path = [thisScriptsDir+"/.."] + sys.path
import gather_build_stats as GBS
import FindTribitsCiSupportDir
import GeneralScriptSupport as GSS

g_pp = pprint.PrettyPrinter(indent=2)

# Shared test data

g_listOfDicts = [
  {'field1':'11', 'field2':'12', 'field4':'14'},
  {'field1':'21', 'field2':'22', 'field3':'23', 'field5':"25"},
  ]


#############################################################################
#
# Test gather_build_stats.readAllValidTimingFiles()
#
#############################################################################


class test_readAllValidTimingFiles(unittest.TestCase):

  def test_1(self):
    baseDir = g_testBaseDir+"/dummy_build_dir"
    allValidTimingFiles = GBS.readAllValidTimingFiles(baseDir, printErrMsg=False)
    allValidTimingFiles_expected = [
      {'FileName': 'target4.o',
       'FileSize': '260000',
       'elapsed_real_time_sec': '1.9',
       'max_resident_size_Kb': '2000'},
      {'FileName': 'packages/pkga/src/target2.lib',
       'FileSize': '870000',
       'cpu_sec_user_mode': '1.38',
       'elapsed_real_time_sec': '1.5',
       'max_resident_size_Kb': '180000'},
      {'FileName': 'some/base/dir/target1.o',
       'FileSize': '3300000',
       'elapsed_real_time_sec': '3.5',
       'max_resident_size_Kb': '240000',
       'num_filesystem_outputs': '20368',
       'num_involuntary_context_switch': '46'}]
    # NOTE: The bad timign file 'some/base/target3.timing' was gracefully
    # skipped!
    allValidTimingFiles.sort(key=lambda item: item.get('FileName')) # Avoid system-dependent behavior
    allValidTimingFiles_expected.sort(key=lambda item: item.get('FileName'))
    self.assertEqual(allValidTimingFiles, allValidTimingFiles_expected)


#############################################################################
#
# Test gather_build_stats.readBuildStatsTimingFileIntoDict()
#
#############################################################################


def readBuildStatsTimingFileIntoDictTest(testObj, buildStatsTimingFile,
    numKeys_expected, buildStatsTimingDict_expected, errMsg_expected,
  ):
  (buildStatsTimingDict, errMsg) = GBS.readBuildStatsTimingFileIntoDict(
    buildStatsTimingFile)
  testObj.assertEqual(errMsg, errMsg_expected)
  if numKeys_expected > 0:
    testObj.assertEqual(len(buildStatsTimingDict.keys()), numKeys_expected)
  testObj.assertEqual(buildStatsTimingDict, buildStatsTimingDict_expected)


class test_readBuildStatsTimingFileIntoDict(unittest.TestCase):

  def test_correct(self):
    buildStatsTimingFile = \
      g_testBaseDir+"/dummy_build_dir/some/base/dir/target1.timing"
    numKeys_expected = 6
    buildStatsTimingDict_expected = {
      'FileName': 'some/base/dir/target1.o',
      'FileSize': '3300000',
      'elapsed_real_time_sec': '3.5',
      'max_resident_size_Kb': '240000',
      'num_filesystem_outputs': '20368',
      'num_involuntary_context_switch': '46',
      }
    errMsg_expected = ""
    readBuildStatsTimingFileIntoDictTest(self, buildStatsTimingFile,
      numKeys_expected, buildStatsTimingDict_expected, errMsg_expected)

  def test_missing_fail(self):
    buildStatsTimingFile = \
      g_testBaseDir+"/file_does_not_exist.timing"
    numKeys_expected = 0
    buildStatsTimingDict_expected = None
    errMsg_expected = buildStatsTimingFile+": ERROR: File does not exist!"
    readBuildStatsTimingFileIntoDictTest(self, buildStatsTimingFile,
      numKeys_expected, buildStatsTimingDict_expected, errMsg_expected)

  def test_two_data_rows_fail(self):
    buildStatsTimingFile = \
      g_testBaseDir+"/bad_timing_build_stats_files/target1.timing.two_data_rows"
    numKeys_expected = 0
    buildStatsTimingDict_expected = None
    errMsg_expected = buildStatsTimingFile+": ERROR: Contains 2 != 1 data rows!"
    readBuildStatsTimingFileIntoDictTest(self, buildStatsTimingFile,
      numKeys_expected, buildStatsTimingDict_expected, errMsg_expected)

  def test_empty_fail(self):
    buildStatsTimingFile = \
      g_testBaseDir+"/bad_timing_build_stats_files/target1.timing.empty"
    numKeys_expected = 0
    buildStatsTimingDict_expected = None
    errMsg_expected = buildStatsTimingFile+": ERROR: File is empty!"
    readBuildStatsTimingFileIntoDictTest(self, buildStatsTimingFile,
      numKeys_expected, buildStatsTimingDict_expected, errMsg_expected)

  def test_junk_fail(self):
    buildStatsTimingFile = \
      g_testBaseDir+"/bad_timing_build_stats_files/target1.timing.junk"
    numKeys_expected = 0
    buildStatsTimingDict_expected = None
    errMsg_expected = buildStatsTimingFile+": ERROR: Error, for CSV file"+\
     " '"+buildStatsTimingFile+"' the data row 0 ['for this garbage'] has 1 entries"+\
     " which does not macth the number of column headers 3!"
    readBuildStatsTimingFileIntoDictTest(self, buildStatsTimingFile,
      numKeys_expected, buildStatsTimingDict_expected, errMsg_expected)
    # NOTE: The above test is very much tied to the implementation of
    # readCsvFileIntoListOfDicts() for the error message it puts out.  That is
    # very

  def test_missing_col_filename_fail(self):
    buildStatsTimingFile = \
      g_testBaseDir+"/bad_timing_build_stats_files/target1.timing.missing_col_filename"
    numKeys_expected = 0
    buildStatsTimingDict_expected = None
    errMsg_expected = \
      buildStatsTimingFile+": ERROR: The required field 'FileName' is missing!"
    readBuildStatsTimingFileIntoDictTest(self, buildStatsTimingFile,
      numKeys_expected, buildStatsTimingDict_expected, errMsg_expected)

  def test_bad_type_filesize_fail(self):
    buildStatsTimingFile = \
      g_testBaseDir+"/bad_timing_build_stats_files/target1.timing.bad_type_filesize"
    numKeys_expected = 0
    buildStatsTimingDict_expected = None
    errMsg_expected = \
      buildStatsTimingFile+": ERROR: For field 'FileSize' the string value"+\
       " 'bad size type' could not be converted to the expected type 'float'!"
    readBuildStatsTimingFileIntoDictTest(self, buildStatsTimingFile,
      numKeys_expected, buildStatsTimingDict_expected, errMsg_expected)


#############################################################################
#
# Test gather_build_stats.writeDictOfListsToCsvFile()
#
#############################################################################

class test_writeDictOfListsToCsvFile(unittest.TestCase):

  def test_1(self):
    dictOfLists = GBS.getDictOfListsFromListOfDicts(g_listOfDicts)
    csvFile = "test_writeDictOfListsToCsvFile_build_stats.csv"
    csvFileText_expected = \
      "field1,field2,field3,field4,field5\n11,12,,14,\n21,22,23,,25\n"
    GBS.writeDictOfListsToCsvFile(dictOfLists, csvFile)
    with open(csvFile, 'r') as csvFileHandle:
      csvFileText = csvFileHandle.read()
    self.assertEqual(csvFileText, csvFileText_expected)


#############################################################################
#
# Test gather_build_stats.getListOfAllTimingFiles()
#
#############################################################################


class test_getListOfAllTimingFiles(unittest.TestCase):

  def test_1(self):
    baseDir = g_testBaseDir+"/dummy_build_dir"
    listOfAllTimingFiles = GBS.getListOfAllTimingFiles(baseDir)
    listOfAllTimingFiles_expected = [
      'packages/pkga/src/target2.timing',
      'some/base/dir/target1.timing',
      'some/base/target3.timing',
      'target4.timing',
      ]
    listOfAllTimingFiles.sort() # Avoid system-dependent behavior
    listOfAllTimingFiles_expected.sort()
    self.assertEqual(listOfAllTimingFiles, listOfAllTimingFiles_expected)


#############################################################################
#
# Test gather_build_stats.getDictOfListsFromListOfDicts()
#
#############################################################################

class test_getDictOfListsFromListOfDicts(unittest.TestCase):

  def test_1(self):
    dictOfLists = GBS.getDictOfListsFromListOfDicts(g_listOfDicts)
    dictOfLists_expected = {
      'field1': ['11', '21'],
      'field2': ['12', '22'],
      'field3': ['', '23'],
      'field4': ['14', ''],
      'field5': ['', '25'],
      }
    self.assertEqual(dictOfLists, dictOfLists_expected)


#############################################################################
#
# Test gather_build_stats.getSupersetOfFieldNamesList()
#
#############################################################################

class test_getSupersetOfFieldNamesList(unittest.TestCase):

  def test_1(self):
    supersetOfFieldNamesList = GBS.getSupersetOfFieldNamesList(g_listOfDicts)
    supersetOfFieldNamesList_expected = \
      ['field1', 'field2', 'field3', 'field4', 'field5']
    supersetOfFieldNamesList.sort() # Make system independent
    supersetOfFieldNamesList_expected.sort()
    self.assertEqual(supersetOfFieldNamesList, supersetOfFieldNamesList_expected)



#############################################################################
#
# Test gather_build_stats.py
#
#############################################################################


csvFileText_expected = \
  "FileName,FileSize,cpu_sec_user_mode,elapsed_real_time_sec,max_resident_size_Kb,num_filesystem_outputs,num_involuntary_context_switch\n"+\
  "target4.o,260000,,1.9,2000,,\n"+\
  "some/base/dir/target1.o,3300000,,3.5,240000,20368,46\n"+\
  "packages/pkga/src/target2.lib,870000,1.38,1.5,180000,,\n"


def gather_build_stats_py_expected_output(csvFile):
  return \
    "Reading all *.timing files from under '"+g_testBaseDir+"/dummy_build_dir' ...\n"+\
    "Number of *.timing files found = 4\n"+\
    g_testBaseDir+"/dummy_build_dir/some/base/target3.timing: ERROR: Contains 0 != 1 data rows!\n"+\
    "Number of valid *.timing files found = 3\n"+\
    "Combined build-stats keys sorted:\n"+\
    "  ['FileName', 'FileSize', 'cpu_sec_user_mode', 'elapsed_real_time_sec', 'max_resident_size_Kb', 'num_filesystem_outputs', 'num_involuntary_context_switch']\n"+\
    "Wrote file '"+csvFile+"'\n"


def sortCsvFileTextList(csvFileText):
  csvFileTextList_orig = csvFileText.split('\n')
  csvFileTextList = []
  csvFileTextList.append(csvFileTextList_orig[0]) # Headers
  csvFileTextList.extend(sorted(csvFileTextList_orig[1:])) # Rows
  return csvFileTextList


def test_gather_build_stats_py_body(testObj, csvFile, cmnd, silentStdout=False):
  output = GSS.getCmndOutput(cmnd)
  #print("output:\n"+output)
  with open(csvFile, 'r') as csvFileHandle:
    csvFileText = csvFileHandle.read()
  testObj.assertEqual(
    sortCsvFileTextList(csvFileText),
    sortCsvFileTextList(csvFileText_expected))
  if not silentStdout:
    testObj.assertEqual(
      output.split('\n'),
      gather_build_stats_py_expected_output(csvFile).split('\n'))


class test_gather_build_stats_py(unittest.TestCase):

  def test_help(self):
    cmnd = thisScriptsDir+"/../gather_build_stats.py --help"
    output = GSS.getCmndOutput(cmnd)
    #print("output:\n"+output+"\n")
    self.assertTrue(output.find("Gather up build stats from *.timing CSV files")!=-1)
    self.assertTrue(output.find("max_resident_size_Kb : float")!=-1)
    self.assertTrue(output.find("FileName : string")!=-1)
    self.assertTrue(output.find("The column headers in all of the *.timing files are combined")!=-1)

  def test_default_out_file(self):
    csvFile = "build_stats.csv"
    cmnd = thisScriptsDir+"/../gather_build_stats.py"+\
      " -d "+g_testBaseDir+"/dummy_build_dir"
    test_gather_build_stats_py_body(self, csvFile, cmnd, silentStdout=True)

  def test_default_out_file_verbose(self):
    csvFile = "build_stats.csv"
    cmnd = thisScriptsDir+"/../gather_build_stats.py -v"+\
      " -d "+g_testBaseDir+"/dummy_build_dir"
    test_gather_build_stats_py_body(self, csvFile, cmnd)

  def test_explicit_out_file_verbose(self):
    csvFile = "test_gather_build_stats_py_build_stats.csv"
    cmnd = thisScriptsDir+"/../gather_build_stats.py -v"+\
      " -d "+g_testBaseDir+"/dummy_build_dir "+csvFile
    test_gather_build_stats_py_body(self, csvFile, cmnd)


#
# Run the unit tests!
#

if __name__ == '__main__':

  unittest.main()
