// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include <catch2/catch_test_macros.hpp>

TEST_CASE_METHOD(CatalystTestFixture, "CatalystCGNSVariableComponents",
                 "[catalyst partitioned data set collection]")
{

  runCatalystMultiBlockMeshTest("aero_blunt_wedge_test3.cgns");
  checkTestOutputFileExists("iossDatabaseCatalystMesh_time_0.vtpc");
  checkTestOutputFileExists("iossDatabaseCatalystMesh_time_1.vtpc");
  checkTestOutputFileExists("iossDatabaseCatalystMesh_time_2.vtpc");
  checkTestOutputFileExists("iossDatabaseCatalystMesh_time_3.vtpc");
  checkTestOutputFileExists("iossDatabaseCatalystMesh_time_4.vtpc");

  CatalystTestFixture::VarAndCompCountVec cellVariables  = {{"density", 1},     {"pressure", 1},
                                                            {"temperature", 1}, {"velocity", 3},
                                                            {"object_id", 1},   {"cell_ids", 1}};
  CatalystTestFixture::VarAndCompCountVec pointVariables = {{"cell_node_ids", 1}};
  CatalystTestFixture::VarAndCompCountVec globalVariables;

  checkMeshOutputVariables("iossDatabaseCatalystMesh_time_4.vtpc", cellVariables, pointVariables,
                           globalVariables, "/IOSS/structured_blocks");
}

TEST_CASE_METHOD(CatalystTestFixture, "CatalystExodusVariableComponents",
                 "[catalyst partitioned data set collection]")
{

  runCatalystMultiBlockMeshTest("block_crush_1.ex2");
  checkTestOutputFileExists("iossDatabaseCatalystMesh_time_0.vtpc");
  checkTestOutputFileExists("iossDatabaseCatalystMesh_time_10.vtpc");

  CatalystTestFixture::VarAndCompCountVec cellVariables = {
      {"vonmises", 1}, {"stress", 6}, {"object_id", 1}, {"ids", 1}};
  CatalystTestFixture::VarAndCompCountVec pointVariables = {
      {"cetan", 1}, {"cftan", 1}, {"vel", 3}, {"displ", 3}, {"ids", 1}};
  CatalystTestFixture::VarAndCompCountVec globalVariables = {{"momentum", 3}, {"kineticenergy", 1}};

  checkMeshOutputVariables("iossDatabaseCatalystMesh_time_0.vtpc", cellVariables, pointVariables,
                           globalVariables, "/IOSS/element_blocks");
}

TEST_CASE_METHOD(CatalystTestFixture, "CatalystCGNSPartitionedDataSetCollectionStructure",
                 "[catalyst partitioned data set collection]")
{

  runCatalystMultiBlockMeshTest("aero_blunt_wedge_test3.cgns");
  StringVec partitions    = {"blk-1"};
  StringVec searchQueries = {"/IOSS/structured_blocks/blk-1"};
  checkPartitionedDataSetCollectionStructure("iossDatabaseCatalystMesh_time_4.vtpc", partitions, 16,
                                             searchQueries);
}

TEST_CASE_METHOD(CatalystTestFixture, "CatalystExodusPartitionedDataSetCollectionStructure",
                 "[catalyst partitioned data set collection]")
{

  runCatalystMultiBlockMeshTest("block_crush_1.ex2");
  StringVec partitions    = {"block_1", "block_2", "rigidbodyoutputpart"};
  StringVec searchQueries = {"/IOSS/element_blocks/block_1", "/IOSS/element_blocks/block_2",
                             "/IOSS/element_blocks/rigidbodyoutputpart"};
  checkPartitionedDataSetCollectionStructure("iossDatabaseCatalystMesh_time_10.vtpc", partitions,
                                             132, searchQueries);
}
