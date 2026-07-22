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
from decimal import Decimal

thisScriptsDir = os.path.dirname(os.path.abspath(__file__))
g_testBaseDir = thisScriptsDir
sys.path = [thisScriptsDir+"/.."] + sys.path
import summarize_build_stats as SBS
import FindTribitsCiSupportDir
import GeneralScriptSupport as GSS

g_pp = pprint.PrettyPrinter(indent=2)


# Get a copied dict of lists of build stats read from input file
#
def getBuildStatsForTests(computeStdScaledFields=True):
  global g_buildStatsDOL
  if not g_buildStatsDOL:
    g_buildStatsDOL = SBS.readBuildStatsCsvFileIntoDictOfLists(
      g_testBaseDir+"/build_stats.big.small.csv",
      computeStdScaledFields=computeStdScaledFields,
      )
  return copy.deepcopy(g_buildStatsDOL)

g_buildStatsDOL = None

# Note: It is structured this way above because we don't want the unit test
# for the function readBuildStatsCsvFileIntoDictOfLists() to be unable to run if a defect is added to it.


#############################################################################
#
# Test summarize_build_stats.readBuildStatsCsvFileIntoDictOfLists()
#
#############################################################################


class test_readBuildStatsCsvFileIntoDictOfLists(unittest.TestCase):

  def test_build_stats_big_little(self):
    buildStatsDOL = SBS.readBuildStatsCsvFileIntoDictOfLists(
      g_testBaseDir+"/build_stats.big.small.csv", computeStdScaledFields=False )
    numCols_expected = 4
    numRows_expected = 21
    self.assertEqual(len(buildStatsDOL.keys()), numCols_expected)
    self.assertEqual(len(buildStatsDOL['max_resident_size_Kb']),numRows_expected)
    self.assertEqual(len(buildStatsDOL['elapsed_real_time_sec']), numRows_expected)
    self.assertEqual(len(buildStatsDOL['FileName']), numRows_expected)
    self.assertEqual(len(buildStatsDOL['FileSize']), numRows_expected)
    self.assertEqual(buildStatsDOL['max_resident_size_Kb'][0], 240000)
    self.assertEqual(buildStatsDOL['max_resident_size_Kb'][11], 730000)
    self.assertEqual(buildStatsDOL['max_resident_size_Kb'][20], 77000)
    self.assertEqual(buildStatsDOL['elapsed_real_time_sec'][0], 3.5)
    self.assertEqual(buildStatsDOL['elapsed_real_time_sec'][11], 48.2)
    self.assertEqual(buildStatsDOL['elapsed_real_time_sec'][20], 0.4)
    self.assertEqual(buildStatsDOL['FileName'][0],
      "commonTools/gtest/CMakeFiles/gtest.dir/gtest/gtest-all.cc.o")
    self.assertEqual(buildStatsDOL['FileName'][11],
      "packages/rol/adapters/epetra/test/sol/CMakeFiles/ROL_adapters_epetra_test_sol_EpetraSROMSampleGenerator.dir/test_02.cpp.o")
    self.assertEqual(buildStatsDOL['FileName'][20],
      "packages/adelus/test/vector_random/Adelus_vector_random.exe")
    self.assertEqual(buildStatsDOL['FileSize'][0], 3300000)
    self.assertEqual(buildStatsDOL['FileSize'][11], 17000000)
    self.assertEqual(buildStatsDOL['FileSize'][20], 5200000)

# NOTE: Above tests also indireclty tests
# summarize_build_stats.readCsvFileIntoDictOfLists()!


#############################################################################
#
# Test summarize_build_stats.readCsvFileIntoDictOfLists()
#
#############################################################################


def cnat(colName, colType="string"):
  return SBS.ColNameAndType(colName, colType)


class test_readCsvFileIntoDictOfLists(unittest.TestCase):

  def test_invalid_col_header(self):
    threwExcept = True
    try:
      buildStatsDOL = SBS.readCsvFileIntoDictOfLists(
        g_testBaseDir+"/build_stats.big.small.csv",
        [
          cnat('missing_header', 'float'),
          ]
        )
      threwExcept = False
    except Exception as exceptObj:
      errMsg = str(exceptObj)
      subStrIdx = errMsg.find(
        "Error, the CSV file column header 'missing_header' does not exist")
      self.assertEqual(subStrIdx, 0)
      # ToDo: Do a better match with fail msg above
    if not threwExcept:
      self.assertFalse("ERROR: Did not thown an excpetion")


  def test_row_missing_elements(self):
    threwExcept = True
#    buildStatsDOL = SBS.readCsvFileIntoDictOfLists(
#      g_testBaseDir+"/build_stats.incomplete_row.csv",
#      [
#        cnat('max_resident_size_Kb', 'float'),
#        cnat('elapsed_real_time_sec', 'float'),
#        ]
#      )
    try:
      buildStatsDOL = SBS.readCsvFileIntoDictOfLists(
        g_testBaseDir+"/build_stats.incomplete_row.csv",
        [
          cnat('max_resident_size_Kb', 'float'),
          cnat('elapsed_real_time_sec', 'float'),
          ]
        )
      threwExcept = False
    except Exception as exceptObj:
      errMsg = str(exceptObj)
      errMsgSubStrExpected = \
        "build_stats.incomplete_row.csv' has 10 column headers but data row 2"+\
        " only has 5 entries"
      subStrIdx = errMsg.find(errMsgSubStrExpected)
      self.assertNotEqual(subStrIdx, -1)
      # ToDo: Do a better match with fail msg above
    if not threwExcept:
      self.assertFalse("ERROR: Did not thown an excpetion")


#############################################################################
#
# Test summarize_build_stats.getColNameTypeIdxListGivenColNameAndTypeList()
#
#############################################################################


def cnti(colName, colType, colIdx):
  return SBS.ColNameTypeIdx(SBS.ColNameAndType(colName, colType), colIdx)


class test_getColNameTypeIdxListGivenColNameAndTypeList(unittest.TestCase):


  def test_subset(self):
    colNameAndIdxList = SBS.getColNameTypeIdxListGivenColNameAndTypeList(
      "dummyFileName",
      [ "aaa", "bbb", "ccc", "ddd", "eee", "fff", "ggg", "hhh" ],
      [ cnat("bbb",'int'), cnat("ccc",'float'), cnat("fff",'string') ] )
    colNameAndIdxList_expected = [
       cnti("bbb",'int',1) , cnti("ccc",'float',2), cnti("fff",'string',5) ]
    self.assertEqual(colNameAndIdxList, colNameAndIdxList_expected)


  def test_missing_col_header(self):
    threwExcept = True
    try:
      colNameAndIdxList = SBS.getColNameTypeIdxListGivenColNameAndTypeList(
        "dummyFileName",
        [ "aaa", "bbb", "ccc", "ddd", "eee", "fff", "ggg", "hhh" ],
        [ cnat("bbb"), cnat("ccc"), cnat("ffg", "int") ] )
      threwExcept = False
    except Exception as errMsg:
      self.assertEqual( str(errMsg),
        "Error, the CSV file column header 'ffg' does not exist in the list"+\
        " of column headers ['aaa', 'bbb', 'ccc', 'ddd', 'eee', 'fff', 'ggg', 'hhh']"+\
        " from the CSV file 'dummyFileName'!" )
    if not threwExcept:
      self.assertFalse("ERROR: Did not thown an excpetion")



#############################################################################
#
# Test summarize_build_stats.ColNameTypeIdx
#
#############################################################################

class test_ColNameTypeIdx(unittest.TestCase):

  def test_float(self):
    colNameTypeIdx = SBS.ColNameTypeIdx(SBS.ColNameAndType("name", "float"), 5)
    asStr = str(colNameTypeIdx)
    self.assertEqual(asStr, "ColNameTypeIdx{ColNameAndType{name,float},5}")
    self.assertEqual(colNameTypeIdx.colName(), "name")
    self.assertEqual(colNameTypeIdx.getIdx(), 5)
    self.assertEqual(colNameTypeIdx.convertFromStr("10.5"), 10.5)

  def test_int(self):
    colNameTypeIdx = SBS.ColNameTypeIdx(SBS.ColNameAndType("name", "int"), 4)
    asStr = str(colNameTypeIdx)
    self.assertEqual(asStr, "ColNameTypeIdx{ColNameAndType{name,int},4}")
    self.assertEqual(colNameTypeIdx.colName(), "name")
    self.assertEqual(colNameTypeIdx.getIdx(), 4)
    self.assertEqual(colNameTypeIdx.convertFromStr("12"), 12)

  def test_string(self):
    colNameTypeIdx = SBS.ColNameTypeIdx(SBS.ColNameAndType("name", "string"), 3)
    asStr = str(colNameTypeIdx)
    self.assertEqual(asStr, "ColNameTypeIdx{ColNameAndType{name,string},3}")
    self.assertEqual(colNameTypeIdx.colName(), "name")
    self.assertEqual(colNameTypeIdx.getIdx(), 3)
    self.assertEqual(colNameTypeIdx.convertFromStr("some str"), "some str")

  def test_invalid_type(self):
    threwExcept = True
    try:
      colNameTypeIdx = SBS.ColNameTypeIdx(SBS.ColNameAndType("name", "invalid"), 2)
      threwExcept = False
    except Exception as errMsg:
      self.assertEqual( str(errMsg),
        "Error, type 'invalid' is not supported!  Supported types include"+\
        " ['string', 'int', 'float']!" )
    if not threwExcept:
      self.assertFalse("ERROR: Did not thown an excpetion")


#############################################################################
#
# Test summarize_build_stats.addStdScaledBuildStatsFields()
#
#############################################################################

class test_addStdScaledBuildStatsFields(unittest.TestCase):

  def test_read_in_and_create_new_fields(self):
    buildStatsDOL = SBS.readBuildStatsCsvFileIntoDictOfLists(
      g_testBaseDir+"/build_stats.big.small.csv" )
    SBS.addStdScaledBuildStatsFields(buildStatsDOL)
    self.assertEqual(len(buildStatsDOL), 6)
    self.assertEqual(len(buildStatsDOL['max_resident_size_mb']), 21)
    self.assertEqual(len(buildStatsDOL['file_size_mb']), 21)
    self.assertEqual(buildStatsDOL['max_resident_size_mb'][0], 234.38)
    self.assertEqual(buildStatsDOL['max_resident_size_mb'][11], 712.89)
    self.assertEqual(buildStatsDOL['max_resident_size_mb'][20], 75.20)
    self.assertEqual(buildStatsDOL['file_size_mb'][0], 3.15)
    self.assertEqual(buildStatsDOL['file_size_mb'][11], 16.21)
    self.assertEqual(buildStatsDOL['file_size_mb'][20], 4.96)


#############################################################################
#
# Test summarize_build_stats.addNewFieldByScalingExistingField()
#
#############################################################################

class test_addNewFieldByScalingExistingField(unittest.TestCase):

  def test_add_field_1(self):
    dictOfLists = { 'field_1' : [ 1.1, 2.2, 3.3 ] }
    SBS.addNewFieldByScalingExistingField(dictOfLists, 'field_1', 0.1, 2,
      'scaled_field')
    self.assertEqual(len(dictOfLists.keys()), 2)
    self.assertEqual(dictOfLists['field_1'], [ 1.1, 2.2, 3.3 ])
    self.assertEqual(dictOfLists['scaled_field'], [ 0.11, 0.22, 0.33 ])


#############################################################################
#
# Test summarize_build_stats.binBuildStatsDictOfListsBySubdirUnderDirs()
#
#############################################################################


dummy1BuildStatsDOL = {
  'field1' : [ "00", "01", "02", "03", "04" ],
  'FileName' : [
    "basedir/pkg0/some_file0",
    "basedir/pkg1/some_file1",
    "basedir/pkg0/some_file2",
    "basedir/pkg1/some_file3",
    "basedir/pkg2/some_file4",
    ],
  'field2' : [ "10", "11", "12", "13", "14" ],
  }

binnedDummy1BuildStatsDOL_dict = {
  "pkg0": {
    'field1' : [ "00", "02" ],
    'FileName' : [
      "basedir/pkg0/some_file0",
      "basedir/pkg0/some_file2",
      ],
    'field2' : [ "10", "12" ],
    },
  "pkg1" : {
    'field1' : [ "01", "03" ],
    'FileName' : [
      "basedir/pkg1/some_file1",
      "basedir/pkg1/some_file3",
       ],
    'field2' : [ "11", "13" ],
    },
  "pkg2" : {
    'field1' : [ "04" ],
    'FileName' : [
      "basedir/pkg2/some_file4",
      ],
    'field2' : [ "14" ],
     },
  }

dummy2BuildStatsDOL = {
  'field1' : [ "00", "01", "02", "03", "04" ],
  'FileName' : [
    "dir2/pkg0/some_file0",
    "basedir/pkg1/some_file1",
    "dir2/pkg0/some_file2",
    "basedir/pkg1/some_file3",
    "basedir/pkg2/some_file4",
    ],
  'field2' : [ "10", "11", "12", "13", "14" ],
  }

binnedDummy2BuildStatsDOL_dict = {
  "pkg0": {
    'field1' : [ "00", "02" ],
    'FileName' : [
      "dir2/pkg0/some_file0",
      "dir2/pkg0/some_file2",
      ],
    'field2' : [ "10", "12" ],
    },
  "pkg1" : {
    'field1' : [ "01", "03" ],
    'FileName' : [
      "basedir/pkg1/some_file1",
      "basedir/pkg1/some_file3",
       ],
    'field2' : [ "11", "13" ],
    },
  "pkg2" : {
    'field1' : [ "04" ],
    'FileName' : [
      "basedir/pkg2/some_file4",
      ],
    'field2' : [ "14" ],
     },
  }


class test_binBuildStatsDictOfListsBySubdirUnderDirs(unittest.TestCase):

  def test_1(self):
    buildStatsDOL = dummy1BuildStatsDOL
    buildStatsBinnedBySubdirs = SBS.binBuildStatsDictOfListsBySubdirUnderDirs(
      buildStatsDOL, [ "basedir" ] )
    self.assertEqual(buildStatsBinnedBySubdirs.fullBuildStatsDOL, dummy1BuildStatsDOL)
    binnedBuildStatsDOL_dict = buildStatsBinnedBySubdirs.binnedBuildStatsDOL_dict
    self.assertEqual(len(binnedBuildStatsDOL_dict.keys()), 3)
    self.assertEqual(binnedBuildStatsDOL_dict, binnedDummy1BuildStatsDOL_dict)

  def test_2(self):
    buildStatsDOL = dummy2BuildStatsDOL
    buildStatsBinnedBySubdirs = SBS.binBuildStatsDictOfListsBySubdirUnderDirs(
      buildStatsDOL, [ "basedir", "dir2" ] )
    self.assertEqual(buildStatsBinnedBySubdirs.fullBuildStatsDOL, dummy2BuildStatsDOL)
    binnedBuildStatsDOL_dict = buildStatsBinnedBySubdirs.binnedBuildStatsDOL_dict
    self.assertEqual(len(binnedBuildStatsDOL_dict.keys()), 3)
    self.assertEqual(binnedBuildStatsDOL_dict, binnedDummy2BuildStatsDOL_dict)


#############################################################################
#
# Test summarize_build_stats.computeBuildStatusSummaryForOneField()
#
#############################################################################

class test_computeBuildStatusSummaryForOneField(unittest.TestCase):

  def test_field_1(self):
    buildStatsDOL = getBuildStatsForTests()
    buildStatSummary = \
      SBS.computeBuildStatusSummaryForOneField(buildStatsDOL, 'max_resident_size_mb', 2)
    self.assertEqual(buildStatSummary.fieldName, 'max_resident_size_mb')
    self.assertEqual(buildStatSummary.numValues, 21)
    self.assertEqual(buildStatSummary.sumValue, 10023.45)
    self.assertEqual(buildStatSummary.maxValue, 2400000/1024.0)
    self.assertEqual(buildStatSummary.maxFileName,
      'packages/tpetra/classic/NodeAPI/CMakeFiles/tpetraclassicnodeapi.dir/Kokkos_DefaultNode.cpp.o' )


#############################################################################
#
# Test summarize_build_stats.computeStdBuildStatsSummariesSingleDOL()
#
#############################################################################


class test_computeStdBuildStatsSummariesSingleDOL(unittest.TestCase):

  def test_big_small(self):
    buildStatsDOL = getBuildStatsForTests()
    bssl = SBS.computeStdBuildStatsSummariesSingleDOL(buildStatsDOL)
    self.assertEqual(len(bssl), 3)
    self.assertEqual(bssl[0].fieldName, 'max_resident_size_mb')
    self.assertEqual(bssl[0].numValues, 21)
    self.assertEqual(bssl[0].sumValue, 10023.45)
    self.assertEqual(bssl[0].maxValue, 2400000/1024.0)
    self.assertEqual(bssl[0].maxFileName,
      'packages/tpetra/classic/NodeAPI/CMakeFiles/tpetraclassicnodeapi.dir/Kokkos_DefaultNode.cpp.o' )
    self.assertEqual(bssl[1].fieldName, 'elapsed_real_time_sec')
    self.assertEqual(bssl[1].numValues, 21)
    self.assertEqual(bssl[1].sumValue, 157.9)
    self.assertEqual(bssl[1].maxValue, 48.2)
    self.assertEqual(bssl[1].maxFileName,
     'packages/rol/adapters/epetra/test/sol/CMakeFiles/ROL_adapters_epetra_test_sol_EpetraSROMSampleGenerator.dir/test_02.cpp.o' )
    self.assertEqual(bssl[2].fieldName, 'file_size_mb')
    self.assertEqual(bssl[2].numValues, 21)
    self.assertEqual(bssl[2].sumValue, 157.19)
    self.assertEqual(bssl[2].maxValue, SBS.roundNum(45000000/(1024.0*1024.0),2))
    self.assertEqual(bssl[2].maxFileName,
      'packages/panzer/adapters-stk/example/CurlLaplacianExample/CMakeFiles/PanzerAdaptersSTK_CurlLaplacianExample.dir/main.cpp.o' )
  # NOTE: Above is a white-box test and we want to validate the order as that
  # is also the order these stats will be displayed.


#############################################################################
#
# Test summarize_build_stats.computeStdBuildStatsSummaries()
#
#############################################################################

class test_computeStdBuildStatsSummaries(unittest.TestCase):

  def test_big_small(self):
    buildStatsDOL = getBuildStatsForTests()
    buildStatsBinnedBySubdirs = SBS.binBuildStatsDictOfListsBySubdirUnderDirs(
      buildStatsDOL, [ "commonTools", "packages" ] )
    #print("\nbuildStatsBinnedBySubdirs.fullBuildStatsDOL:")
    #g_pp.pprint(buildStatsBinnedBySubdirs.fullBuildStatsDOL)
    #print("buildStatsBinnedBySubdirs.binnedBuildStatsDOL_dict:")
    #g_pp.pprint(buildStatsBinnedBySubdirs.binnedBuildStatsDOL_dict)
    buildStatsSummariesBinnedBySubdirs = SBS.computeStdBuildStatsSummaries(
      buildStatsBinnedBySubdirs )
    # Full project build stats
    bssl = buildStatsSummariesBinnedBySubdirs.fullBuildStatsSummariesList
    self.assertEqual(len(bssl), 3)
    self.assertEqual(bssl[0].fieldName, 'max_resident_size_mb')
    self.assertEqual(bssl[0].numValues, 21)
    self.assertEqual(bssl[0].sumValue, 10023.45)
    self.assertEqual(bssl[0].maxValue, 2400000/1024.0)
    self.assertEqual(bssl[0].maxFileName,
      'packages/tpetra/classic/NodeAPI/CMakeFiles/tpetraclassicnodeapi.dir/Kokkos_DefaultNode.cpp.o' )
    self.assertEqual(bssl[1].fieldName, 'elapsed_real_time_sec')
    self.assertEqual(bssl[1].numValues, 21)
    self.assertEqual(bssl[1].sumValue, 157.9)
    self.assertEqual(bssl[1].maxValue, 48.2)
    self.assertEqual(bssl[1].maxFileName,
     'packages/rol/adapters/epetra/test/sol/CMakeFiles/ROL_adapters_epetra_test_sol_EpetraSROMSampleGenerator.dir/test_02.cpp.o' )
    self.assertEqual(bssl[2].fieldName, 'file_size_mb')
    self.assertEqual(bssl[2].numValues, 21)
    self.assertEqual(bssl[2].sumValue, 157.19)
    self.assertEqual(bssl[2].maxValue, SBS.roundNum(45000000/(1024.0*1024.0),2))
    self.assertEqual(bssl[2].maxFileName,
      'packages/panzer/adapters-stk/example/CurlLaplacianExample/CMakeFiles/PanzerAdaptersSTK_CurlLaplacianExample.dir/main.cpp.o' )
    # Verify number of build stats summaries binned by subdirs
    self.assertEqual(
      len(buildStatsSummariesBinnedBySubdirs.binnedBuildStatsSummariesList_dict.keys()),
      6)
    self.assertEqual(
      sorted(buildStatsSummariesBinnedBySubdirs.binnedBuildStatsSummariesList_dict.keys()),
      ['adelus', 'gtest', 'panzer', 'rol', 'thyra', 'tpetra'])
    # Build stats for 'adelus'
    bssl = buildStatsSummariesBinnedBySubdirs.binnedBuildStatsSummariesList_dict['adelus']
    #print("\nbssl[0]:"); g_pp.pprint(str(bssl[0]))
    self.assertEqual(len(bssl), 3)
    self.assertEqual(bssl[0].fieldName, 'max_resident_size_mb')
    self.assertEqual(bssl[0].numValues, 4)
    self.assertEqual(bssl[0].sumValue, SBS.roundNum((680000+35000+64000+77000)/1024.0,2))
    self.assertEqual(bssl[0].maxValue, SBS.roundNum(680000/1024.0,2))
    self.assertEqual(bssl[0].maxFileName, 'packages/adelus/src/CMakeFiles/zadelus.dir/Adelus_pcomm.cpp.o')
    #print("\nbssl[1]:"); g_pp.pprint(str(bssl[1]))
    self.assertEqual(bssl[1].fieldName, 'elapsed_real_time_sec')
    self.assertEqual(bssl[1].numValues, 4)
    self.assertEqual(bssl[1].sumValue, SBS.roundNum(0.5+0.2+0.3+0.4,2))
    self.assertEqual(bssl[1].maxValue, 0.5)
    self.assertEqual(bssl[1].maxFileName, 'packages/adelus/src/CMakeFiles/zadelus.dir/Adelus_pcomm.cpp.o')
    self.assertEqual(bssl[2].fieldName, 'file_size_mb')
    self.assertEqual(bssl[2].numValues, 4)
    self.assertEqual(bssl[2].sumValue, 5.73)
    self.assertEqual(bssl[2].maxValue, SBS.roundNum(5200000/(1024.0*1024.0),2))
    self.assertEqual(bssl[2].maxFileName, 'packages/adelus/test/vector_random/Adelus_vector_random.exe')
    # Build stats for 'gtest'
    bssl = buildStatsSummariesBinnedBySubdirs.binnedBuildStatsSummariesList_dict['gtest']
    #print("\nbssl[0]:"); g_pp.pprint(str(bssl[0]))
    self.assertEqual(len(bssl), 3)
    self.assertEqual(bssl[0].fieldName, 'max_resident_size_mb')
    self.assertEqual(bssl[0].numValues, 1)
    self.assertEqual(bssl[0].sumValue, SBS.roundNum(240000/1024.0,2))
    self.assertEqual(bssl[0].maxValue, SBS.roundNum(240000/1024.0,2))
    self.assertEqual(bssl[0].maxFileName, 'commonTools/gtest/CMakeFiles/gtest.dir/gtest/gtest-all.cc.o')
    self.assertEqual(bssl[1].fieldName, 'elapsed_real_time_sec')
    self.assertEqual(bssl[1].numValues, 1)
    self.assertEqual(bssl[1].sumValue, 3.5)
    self.assertEqual(bssl[1].maxValue, 3.5)
    self.assertEqual(bssl[1].maxFileName, 'commonTools/gtest/CMakeFiles/gtest.dir/gtest/gtest-all.cc.o')
    self.assertEqual(bssl[2].fieldName, 'file_size_mb')
    self.assertEqual(bssl[2].numValues, 1)
    self.assertEqual(bssl[2].sumValue, SBS.roundNum(3300000/(1024.0*1024.0),2))
    self.assertEqual(bssl[2].maxValue, SBS.roundNum(3300000/(1024.0*1024.0),2))
    self.assertEqual(bssl[2].maxFileName,  'commonTools/gtest/CMakeFiles/gtest.dir/gtest/gtest-all.cc.o')
    # Build stats for 'panzer' (don't bother checking values, above is good enough)
    bssl = buildStatsSummariesBinnedBySubdirs.binnedBuildStatsSummariesList_dict['panzer']
    #print("\nbssl[0]:"); g_pp.pprint(str(bssl[0]))
    self.assertEqual(len(bssl), 3)
    self.assertEqual(bssl[0].fieldName, 'max_resident_size_mb')
    self.assertEqual(bssl[0].numValues, 4)
    self.assertEqual(bssl[1].fieldName, 'elapsed_real_time_sec')
    self.assertEqual(bssl[1].numValues, 4)
    self.assertEqual(bssl[2].fieldName, 'file_size_mb')
    self.assertEqual(bssl[2].numValues, 4)
    # Build stats for 'rol'  (don't bother checking values, above is good enough)
    bssl = buildStatsSummariesBinnedBySubdirs.binnedBuildStatsSummariesList_dict['rol']
    #print("\nbssl[0]:"); g_pp.pprint(str(bssl[0]))
    self.assertEqual(len(bssl), 3)
    self.assertEqual(bssl[0].fieldName, 'max_resident_size_mb')
    self.assertEqual(bssl[0].numValues, 4)
    self.assertEqual(bssl[1].fieldName, 'elapsed_real_time_sec')
    self.assertEqual(bssl[1].numValues, 4)
    self.assertEqual(bssl[2].fieldName, 'file_size_mb')
    self.assertEqual(bssl[2].numValues, 4)
    # Build stats for 'thyra' (don't bother checking values, above is good enough)
    bssl = buildStatsSummariesBinnedBySubdirs.binnedBuildStatsSummariesList_dict['thyra']
    #print("\nbssl[0]:"); g_pp.pprint(str(bssl[0]))
    self.assertEqual(len(bssl), 3)
    self.assertEqual(bssl[0].fieldName, 'max_resident_size_mb')
    self.assertEqual(bssl[0].numValues, 4)
    self.assertEqual(bssl[1].fieldName, 'elapsed_real_time_sec')
    self.assertEqual(bssl[1].numValues, 4)
    self.assertEqual(bssl[2].fieldName, 'file_size_mb')
    self.assertEqual(bssl[2].numValues, 4)
    # Build stats for 'tpetra' (don't bother checking values, above is good enough)
    bssl = buildStatsSummariesBinnedBySubdirs.binnedBuildStatsSummariesList_dict['tpetra']
    #print("\nbssl[0]:"); g_pp.pprint(str(bssl[0]))
    self.assertEqual(len(bssl), 3)
    self.assertEqual(bssl[0].fieldName, 'max_resident_size_mb')
    self.assertEqual(bssl[0].numValues, 4)
    self.assertEqual(bssl[1].fieldName, 'elapsed_real_time_sec')
    self.assertEqual(bssl[1].numValues, 4)
    self.assertEqual(bssl[2].fieldName, 'file_size_mb')
    self.assertEqual(bssl[2].numValues, 4)


#############################################################################
#
# Test summarize_build_stats.createAsciiReportOfBuildStatsSummariesSingleSet()
#
#############################################################################


class test_createAsciiReportOfBuildStatsSummariesSingleSet(unittest.TestCase):

  def test_big_small(self):
    buildStatsDOL = getBuildStatsForTests()
    buildStatsSummariesList = SBS.computeStdBuildStatsSummariesSingleDOL(buildStatsDOL)
    buildStatsAsciiReport = SBS.createAsciiReportOfBuildStatsSummariesSingleSet(
      buildStatsSummariesList, "Full Project")
    self.assertEqual(buildStatsAsciiReport,
      "Full Project: sum(max_resident_size_mb) = 10023.45 (21 entries)\n"+\
      "Full Project: max(max_resident_size_mb) = 2343.75 (packages/tpetra/classic/NodeAPI/CMakeFiles/tpetraclassicnodeapi.dir/Kokkos_DefaultNode.cpp.o)\n"+\
      "Full Project: sum(elapsed_real_time_sec) = 157.9 (21 entries)\n"+\
      "Full Project: max(elapsed_real_time_sec) = 48.2 (packages/rol/adapters/epetra/test/sol/CMakeFiles/ROL_adapters_epetra_test_sol_EpetraSROMSampleGenerator.dir/test_02.cpp.o)\n"+\
      "Full Project: sum(file_size_mb) = 157.19 (21 entries)\n"+\
      "Full Project: max(file_size_mb) = 42.92 (packages/panzer/adapters-stk/example/CurlLaplacianExample/CMakeFiles/PanzerAdaptersSTK_CurlLaplacianExample.dir/main.cpp.o)\n" )


#############################################################################
#
# Test summarize_build_stats.createAsciiReportOfBuildStatsSummariesSingleSet()
#
#############################################################################


big_small_summary_full_project_ascii = \
r"""Full Project: sum(max_resident_size_mb) = 10023.45 (21 entries)
Full Project: max(max_resident_size_mb) = 2343.75 (packages/tpetra/classic/NodeAPI/CMakeFiles/tpetraclassicnodeapi.dir/Kokkos_DefaultNode.cpp.o)
Full Project: sum(elapsed_real_time_sec) = 157.9 (21 entries)
Full Project: max(elapsed_real_time_sec) = 48.2 (packages/rol/adapters/epetra/test/sol/CMakeFiles/ROL_adapters_epetra_test_sol_EpetraSROMSampleGenerator.dir/test_02.cpp.o)
Full Project: sum(file_size_mb) = 157.19 (21 entries)
Full Project: max(file_size_mb) = 42.92 (packages/panzer/adapters-stk/example/CurlLaplacianExample/CMakeFiles/PanzerAdaptersSTK_CurlLaplacianExample.dir/main.cpp.o)
"""

big_small_summary_ascii = \
big_small_summary_full_project_ascii + \
"\n" + \
r"""adelus: sum(max_resident_size_mb) = 835.94 (4 entries)
adelus: max(max_resident_size_mb) = 664.06 (packages/adelus/src/CMakeFiles/zadelus.dir/Adelus_pcomm.cpp.o)
adelus: sum(elapsed_real_time_sec) = 1.4 (4 entries)
adelus: max(elapsed_real_time_sec) = 0.5 (packages/adelus/src/CMakeFiles/zadelus.dir/Adelus_pcomm.cpp.o)
adelus: sum(file_size_mb) = 5.73 (4 entries)
adelus: max(file_size_mb) = 4.96 (packages/adelus/test/vector_random/Adelus_vector_random.exe)

gtest: sum(max_resident_size_mb) = 234.38 (1 entries)
gtest: max(max_resident_size_mb) = 234.38 (commonTools/gtest/CMakeFiles/gtest.dir/gtest/gtest-all.cc.o)
gtest: sum(elapsed_real_time_sec) = 3.5 (1 entries)
gtest: max(elapsed_real_time_sec) = 3.5 (commonTools/gtest/CMakeFiles/gtest.dir/gtest/gtest-all.cc.o)
gtest: sum(file_size_mb) = 3.15 (1 entries)
gtest: max(file_size_mb) = 3.15 (commonTools/gtest/CMakeFiles/gtest.dir/gtest/gtest-all.cc.o)

panzer: sum(max_resident_size_mb) = 3828.13 (4 entries)
panzer: max(max_resident_size_mb) = 1660.16 (packages/panzer/adapters-stk/example/CurlLaplacianExample/CMakeFiles/PanzerAdaptersSTK_CurlLaplacianExample.dir/main.cpp.o)
panzer: sum(elapsed_real_time_sec) = 68.5 (4 entries)
panzer: max(elapsed_real_time_sec) = 37.9 (packages/panzer/adapters-stk/example/CurlLaplacianExample/CMakeFiles/PanzerAdaptersSTK_CurlLaplacianExample.dir/main.cpp.o)
panzer: sum(file_size_mb) = 91.65 (4 entries)
panzer: max(file_size_mb) = 42.92 (packages/panzer/adapters-stk/example/CurlLaplacianExample/CMakeFiles/PanzerAdaptersSTK_CurlLaplacianExample.dir/main.cpp.o)

rol: sum(max_resident_size_mb) = 1982.42 (4 entries)
rol: max(max_resident_size_mb) = 712.89 (packages/rol/adapters/epetra/test/sol/CMakeFiles/ROL_adapters_epetra_test_sol_EpetraSROMSampleGenerator.dir/test_02.cpp.o)
rol: sum(elapsed_real_time_sec) = 75.7 (4 entries)
rol: max(elapsed_real_time_sec) = 48.2 (packages/rol/adapters/epetra/test/sol/CMakeFiles/ROL_adapters_epetra_test_sol_EpetraSROMSampleGenerator.dir/test_02.cpp.o)
rol: sum(file_size_mb) = 51.88 (4 entries)
rol: max(file_size_mb) = 18.12 (packages/rol/adapters/belos/test/vector/CMakeFiles/ROL_adapters_belos_test_vector_BelosInterface.dir/test_01.cpp.o)

thyra: sum(max_resident_size_mb) = 732.42 (4 entries)
thyra: max(max_resident_size_mb) = 195.31 (packages/thyra/adapters/epetra/example/CMakeFiles/ThyraEpetraAdapters_sillyPowerMethod_epetra.dir/sillyPowerMethod_epetra.cpp.o)
thyra: sum(elapsed_real_time_sec) = 7.4 (4 entries)
thyra: max(elapsed_real_time_sec) = 2.2 (packages/thyra/adapters/epetra/example/CMakeFiles/ThyraEpetraAdapters_sillyPowerMethod_epetra.dir/sillyPowerMethod_epetra.cpp.o)
thyra: sum(file_size_mb) = 4.5 (4 entries)
thyra: max(file_size_mb) = 1.43 (packages/thyra/adapters/epetra/example/CMakeFiles/ThyraEpetraAdapters_sillyPowerMethod_epetra.dir/sillyPowerMethod_epetra.cpp.o)

tpetra: sum(max_resident_size_mb) = 2410.16 (4 entries)
tpetra: max(max_resident_size_mb) = 2343.75 (packages/tpetra/classic/NodeAPI/CMakeFiles/tpetraclassicnodeapi.dir/Kokkos_DefaultNode.cpp.o)
tpetra: sum(elapsed_real_time_sec) = 1.4 (4 entries)
tpetra: max(elapsed_real_time_sec) = 1.0 (packages/tpetra/classic/NodeAPI/CMakeFiles/tpetraclassicnodeapi.dir/Kokkos_DefaultNode.cpp.o)
tpetra: sum(file_size_mb) = 0.28 (4 entries)
tpetra: max(file_size_mb) = 0.16 (packages/tpetra/classic/NodeAPI/CMakeFiles/tpetraclassicnodeapi.dir/Kokkos_DefaultNode.cpp.o)
"""


empty_summary_ascii = \
r"""No build statistics to summarize!"""


class test_createAsciiReportOfBuildStatsSummaries(unittest.TestCase):

  def test_big_small(self):
    buildStatsDOL = getBuildStatsForTests()
    buildStatsBinnedBySubdirs = SBS.binBuildStatsDictOfListsBySubdirUnderDirs(
      buildStatsDOL, [ "commonTools", "packages" ] )
    buildStatsSummariesBinnedBySubdirs = SBS.computeStdBuildStatsSummaries(
      buildStatsBinnedBySubdirs )
    buildStatsAsciiReport = SBS.createAsciiReportOfBuildStatsSummaries(
      buildStatsSummariesBinnedBySubdirs )
    self.assertEqual(buildStatsAsciiReport,
      big_small_summary_ascii )


#############################################################################
#
# Test summarize_build_stats.py
#
#############################################################################

class test_summarize_build_stats_py(unittest.TestCase):

  def test_big_small_full_project(self):
    cmnd = thisScriptsDir+"/../summarize_build_stats.py"+\
      " "+g_testBaseDir+"/build_stats.big.small.csv"
    output = GSS.getCmndOutput(cmnd)
    self.assertEqual(GSS.s(output), GSS.s(big_small_summary_full_project_ascii+"\n"))

  def test_big_small_by_subdir(self):
    cmnd = thisScriptsDir+"/../summarize_build_stats.py"+\
      " --bin-by-subdirs-under-dirs=commonTools,packages"+\
      " "+g_testBaseDir+"/build_stats.big.small.csv"
    output = GSS.getCmndOutput(cmnd)
    self.assertEqual(GSS.s(output), GSS.s(big_small_summary_ascii+"\n"))

  def test_empty_build_stats(self):
    cmnd = thisScriptsDir+"/../summarize_build_stats.py"+\
      " --bin-by-subdirs-under-dirs=commonTools,packages"+\
      " "+g_testBaseDir+"/build_stats.empty.csv"
    output = GSS.getCmndOutput(cmnd)
    self.assertEqual(GSS.s(output), GSS.s(empty_summary_ascii+"\n"))


#
# Run the unit tests!
#

if __name__ == '__main__':

  unittest.main()
