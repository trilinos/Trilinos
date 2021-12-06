// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "mpi.h"

#include "UnitTestIotmTextMeshFixture.h"

namespace
{
TEST_F(TestTextMesh, singleHex)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  setup_text_mesh(meshDesc);

  verify_num_elements(1);
  verify_single_element(1u, "HEX_8", EntityIdVector{1,2,3,4,5,6,7,8});
}

TEST_F(TestTextMesh, singleHexWithCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  EntityIdVector nodeIds = EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };

  setup_text_mesh(meshDesc, coordinates);

  verify_num_elements(1);
  verify_single_element(1u, "HEX_8", nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestTextMesh, singleHexWithCoordinates_separatedNodeIds)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,11,12,13,14";
  EntityIdVector nodeIds = EntityIdVector{1, 2, 3, 4, 11, 12, 13, 14};
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };

  setup_text_mesh(meshDesc, coordinates);

  verify_num_elements(1);
  verify_single_element(1u, "HEX_8", nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestTextMesh, twoHexesSerial)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12";
  setup_text_mesh(meshDesc);

  verify_num_elements(2);
  verify_single_element(1u, "HEX_8", EntityIdVector{1,2,3,4,5,6,7,8});
  verify_single_element(2u, "HEX_8", EntityIdVector{5,6,7,8,9,10,11,12});
}

TEST_F(TestTextMesh, twoTet10Serial)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,TET_10,1,2,3,4,5,6,7,8,9,10\n"
                         "0,2,TET_10,2,11,3,4,12,13,6,9,14,10";
//                                       1       2      3        4          5          6         7           8
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 0.5,1,0, 0.5,0.5,1, 0.5,0,0, 0.75,0.5,0, 0.25,0.5,0, 0.25,0.25,0.5,
//                                       9              10            11         12           13        14
                                      0.75,0.25,0.5, 0.5,0.75,0.5, 1.5,0.5,0, 1.25,0.25,0, 1,0.75,0, 1,0.5,0.5 };

  setup_text_mesh(meshDesc, coordinates);

  verify_num_elements(2);
  verify_single_element(1u, "TET_10", EntityIdVector{1,2,3,4,5,6,7,8,9,10});
  verify_single_element(2u, "TET_10", EntityIdVector{2,11,3,4,12,13,6,9,14,10});
}

TEST_F(TestTextMesh, twoHexDisconnectedWithCoordinatesAndParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_1\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,block_2";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                      0,0,1, 1,0,1, 1,1,1, 0,1,1, 0,0,2, 1,0,2, 1,1,2, 0,1,2 };

  setup_text_mesh(meshDesc, coordinates);

  verify_num_elements(2);
  verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  verify_single_element(2u, "HEX_8", EntityIdVector{9, 10, 11, 12, 13, 14, 15, 16});
  verify_coordinates(EntityIdVector{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}, coordinates);
  verify_part_membership({{"block_1", {1u}}, {"block_2", {2u}}});
}

TEST_F(TestTextMesh, twoHexDisconnectedWithDefaultParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16";
  setup_text_mesh(meshDesc);

  verify_num_elements(2);
  verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  verify_single_element(2u, "HEX_8", EntityIdVector{9, 10, 11, 12, 13, 14, 15, 16});
  verify_part_membership({{"block_" + get_topology_name("HEX_8"), {1u, 2u}}});
}

TEST_F(TestTextMesh, threeTriShellsWithCoordinatesAndParts)
{
  //      4-----5             //
  //      |\  2 |\            //
  //      |  \  |  \          //
  //      | 1  \| 3  \        //
  //      1-----2-----3       //

  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,SHELL_TRI_3,1,2,4,block_1\n"
                         "0,2,SHELL_TRI_3,2,5,4,block_2\n"
                         "0,3,SHELL_TRI_3,2,3,5,block_2";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 2,0,0,
    0,1,0, 1,1,0
  };

  setup_text_mesh(meshDesc, coordinates);

  verify_num_elements(3);
  verify_single_element(1u, "SHELL_TRI_3", EntityIdVector{1, 2, 4});
  verify_single_element(2u, "SHELL_TRI_3", EntityIdVector{2, 5, 4});
  verify_single_element(3u, "SHELL_TRI_3", EntityIdVector{2, 3, 5});
  verify_coordinates(EntityIdVector{1, 2, 3, 4, 5}, coordinates);
  verify_part_membership({{"block_1", {1u}}, {"block_2", {2u, 3u}}});
}

TEST_F(TestTextMesh, threeTriShellsWithDefaultParts)
{
  //      4-----5             //
  //      |\  2 |\            //
  //      |  \  |  \          //
  //      | 1  \| 3  \        //
  //      1-----2-----3       //

  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,SHELL_TRI_3,1,2,4\n"
                         "0,2,SHELL_TRI_3,2,5,4\n"
                         "0,3,SHELL_TRI_3,2,3,5";
  setup_text_mesh(meshDesc);

  verify_num_elements(3);
  verify_single_element(1u, "SHELL_TRI_3", EntityIdVector{1, 2, 4});
  verify_single_element(2u, "SHELL_TRI_3", EntityIdVector{2, 5, 4});
  verify_single_element(3u, "SHELL_TRI_3", EntityIdVector{2, 3, 5});
  verify_part_membership({{"block_" + get_topology_name("SHELL_TRI_3"), {1u, 2u, 3u}}});
}


TEST_F(TestTextMesh, partIds_oneDefaultPartOneElem)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"block_" + get_topology_name("HEX_8"), 1u}});
}

TEST_F(TestTextMesh, partIds_oneDefaultPartTwoElems)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"block_" + get_topology_name("HEX_8"), 1u}});
}

TEST_F(TestTextMesh, partIds_twoDefaultParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16\n"
                         "0,3,SHELL_QUAD_4,17,18,19,20";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"block_" + get_topology_name("HEX_8"), 1u}, {"block_" + get_topology_name("SHELL_QUAD_4"), 2u}});
}

TEST_F(TestTextMesh, partIds_oneDefaultPartOneUserSpecifiedPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,my_cool_part";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"block_" + get_topology_name("HEX_8"), 1u}, {"my_cool_part", 2u}});
}

TEST_F(TestTextMesh, partIds_samePartTwice)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,my_cool_part\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,my_cool_part";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"my_cool_part", 1u}});
}

TEST_F(TestTextMesh, partIds_orderingIsByLine)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,2,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partOne\n"
                         "0,1,HEX_8,9,10,11,12,13,14,15,16,partTwo";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"partOne", 1u}, {"partTwo", 2u}});
}

TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_onePart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_101";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"block_101", 101u}});
}

TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_twoParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_101\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,block_201";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"block_101", 101u}, {"block_201", 201u}});
}

TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_withDefaultPartFirst)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,block_101";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"block_" + get_topology_name("HEX_8"), 1u}, {"block_101", 101u}});
}

TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_withDefaultPartSecond)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_101\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"block_101", 101u}, {"block_" + get_topology_name("HEX_8"), 1u}});
}

TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_withDefaultPartIdCollision)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,block_1";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"block_" + get_topology_name("HEX_8"), 2u}, {"block_1", 1u}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_onePart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,partThree,3";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"partThree", 3u}});
  verify_num_elements(1);
  verify_single_element(1u, "HEX_8", EntityIdVector{1,2,3,4,5,6,7,8});
  verify_part_membership({{"partThree", {1u}}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_twoParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,partFive,5";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"partThree", 3u}, {"partFive", 5u}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_twoPartsSameId)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partFour,4\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,partFour,4";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"partFour", 4u}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_samePartDifferentIds)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partFour,4\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,partFour,5";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_withDefaultPartIdCollision)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,partOne,1";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"block_" + get_topology_name("HEX_8"), 2u}, {"partOne", 1u}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_withDefaultPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"partThree", 3u}, {"block_" + get_topology_name("HEX_8"), 1u}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_withExodusPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,block_4";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"partThree", 3u}, {"block_4", 4u}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_collidesWithExodusPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,block_3";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_collidesWithPreviousSpec)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThreeA,3\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,partThreeB,3";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_forExodusPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_2,3\n";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, partIds_shortPartNamesAreValid)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,a\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,b,3";
  setup_text_mesh(meshDesc);

  verify_part_ids({{"a", 1u}, {"b", 3u}});
}

TEST_F(TestTextMesh, partIds_integerPartNamesAreInvalid)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}


TEST_F(TestTextMesh, twoHexesParallel)
{
  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "1,2,HEX_8,5,6,7,8,9,10,11,12";
  setup_text_mesh(meshDesc);

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, "HEX_8", EntityIdVector{1,2,3,4,5,6,7,8});
    verify_shared_nodes(EntityIdVector{5,6,7,8}, 1);
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, "HEX_8", EntityIdVector{5,6,7,8,9,10,11,12});
    verify_shared_nodes(EntityIdVector{5,6,7,8}, 0);
  }
}

TEST_F(TestTextMesh, twoQuadShellsWithCoordinatesParallel)
{
  //      4-----5-----6
  //      |     |     |
  //      |  1  |  2  |
  //      |     |     |
  //      1-----2-----3

  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,5,4\n"
                         "1,2,SHELL_QUAD_4,2,3,6,5";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 2,0,0,
    0,1,0, 1,1,0, 2,1,0
  };

  setup_text_mesh(meshDesc, coordinates);

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, "SHELL_QUAD_4", EntityIdVector{1, 2, 5, 4});
    verify_coordinates(EntityIdVector{1, 2, 4, 5}, {0,0,0, 1,0,0, 0,1,0, 1,1,0});
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, "SHELL_QUAD_4", EntityIdVector{2, 3, 6, 5});
    verify_coordinates(EntityIdVector{2, 3, 5, 6}, {1,0,0, 2,0,0, 1,1,0, 2,1,0});
  }
}

TEST_F(TestTextMesh, threeTriShellsWithCoordinatesParallel)
{
  //      4-----5             //
  //      |\  2 |\            //
  //      |  \  |  \          //
  //      | 1  \| 3  \        //
  //      1-----2-----3       //

  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc = "0,1,SHELL_TRI_3,1,2,4\n"
                         "0,2,SHELL_TRI_3,2,5,4\n"
                         "1,3,SHELL_TRI_3,2,3,5";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 2,0,0,
    0,1,0, 1,1,0
  };

  setup_text_mesh(meshDesc, coordinates);

  if (rank == 0) {
    verify_num_elements(2);
    verify_single_element(1u, "SHELL_TRI_3", EntityIdVector{1, 2, 4});
    verify_single_element(2u, "SHELL_TRI_3", EntityIdVector{2, 5, 4});
    verify_coordinates(EntityIdVector{1, 2, 4, 5}, {0,0,0, 1,0,0, 0,1,0, 1,1,0});
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(3u, "SHELL_TRI_3", EntityIdVector{2, 3, 5});
    verify_coordinates(EntityIdVector{2, 3, 5}, {1,0,0, 2,0,0, 1,1,0});
  }
}

TEST_F(TestTextMesh, singleHexWithSpaces)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0, 1, HEX_8, 1, 2, 3, 4, 5, 6, 7, 8";
  setup_text_mesh(meshDesc);

  verify_num_elements(1);
  verify_single_element(1u, "HEX_8", EntityIdVector{1,2,3,4,5,6,7,8});
}

TEST_F(TestTextMesh, singleHexWithCoordinatesAndSpaces)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0, 1, HEX_8, 1, 2, 3, 4, 5, 6, 7, 8";
  EntityIdVector nodeIds = EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };

  setup_text_mesh(meshDesc, coordinates);

  verify_num_elements(1);
  verify_single_element(1u, "HEX_8", nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestTextMesh, singleHexWithLowerCase)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,Hex_8,1,2,3,4,5,6,7,8";
  setup_text_mesh(meshDesc);

  verify_num_elements(1);
  verify_single_element(1u, "HEX_8", EntityIdVector{1,2,3,4,5,6,7,8});
}

TEST_F(TestTextMesh, singleHexWithCoordinatesAndLowerCase)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,Hex_8,1,2,3,4,5,6,7,8";
  EntityIdVector nodeIds = EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };

  setup_text_mesh(meshDesc, coordinates);

  verify_num_elements(1);
  verify_single_element(1u, "HEX_8", nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestTextMesh, tooFewNodes)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooFewNodesWithCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1 };
  EXPECT_THROW(setup_text_mesh(meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, tooFewCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1 };
  EXPECT_THROW(setup_text_mesh(meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, tooManyNodes)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9,10";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooManyNodesWithCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9,10";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1, 2,0,0, 2,1,0 };
  EXPECT_THROW(setup_text_mesh(meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, tooManyCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1, 52 };
  EXPECT_THROW(setup_text_mesh(meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleData_empty)
{
  std::string meshDesc = "";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleData_startsWithString)
{
  std::string meshDesc = "hi";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleData_noGlobalId)
{
  std::string meshDesc = "0 ";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleData_noTopology)
{
  std::string meshDesc = "0,1,";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleData_noNodes)
{
  std::string meshDesc = "0,1,HEX_8";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleDataWithCoordinates)
{
  std::string meshDesc = "0,1,";
  std::vector<double> coordinates;
  EXPECT_THROW(setup_text_mesh(meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, invalidTopology)
{
  std::string meshDesc = "0,1,invalid,1";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, invalidTopologyWithCoordinates)
{
  std::string meshDesc = "0,1,invalid,1";
  std::vector<double> coordinates = { 0,0,0 };
  EXPECT_THROW(setup_text_mesh(meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, mixedSpatialDim)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,QUAD_4_2D,5,6,7,8";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, mixedSpatialDimWithCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,QUAD_4_2D,5,6,7,8";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };
  EXPECT_THROW(setup_text_mesh(meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, spatialDimInconsistentWithMetaData)
{
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, spatialDimInconsistentWithMetaDataWithCoordinates)
{
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
  std::vector<double> coordinates = { 0,0, 1,0, 1,1, 0,1 };
  EXPECT_THROW(setup_text_mesh(meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, endingWithNewlineIsOk)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n";
  EXPECT_NO_THROW(setup_text_mesh(meshDesc));
}

TEST_F(TestTextMesh, stringAfterPartNameIsError)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1,bogus\n";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, particleHex)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,PARTICLE,1\n"
                         "0,2,HEX_8,2,3,4,5,6,7,8,9";
  EXPECT_NO_THROW(setup_text_mesh(meshDesc));

  verify_num_elements(2);
  verify_single_element(1u, "PARTICLE", EntityIdVector{1});
  verify_single_element(2u, "HEX_8", EntityIdVector{2,3,4,5,6,7,8,9});
}

TEST_F(TestTextMesh, particleHexWithCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,PARTICLE,1\n"
                         "0,2,HEX_8,2,3,4,5,6,7,8,9";
  std::vector<double> coordinates = { 2,0,0, 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };
  EXPECT_NO_THROW(setup_text_mesh(meshDesc, coordinates));

  verify_num_elements(2);
  verify_single_element(1u, "PARTICLE", EntityIdVector{1});
  verify_single_element(2u, "HEX_8", EntityIdVector{2,3,4,5,6,7,8,9});
  verify_coordinates(EntityIdVector{1,2,3,4,5,6,7,8,9}, coordinates);
}

TEST_F(TestTextMesh2d, singleQuad)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
  setup_text_mesh(meshDesc);

  verify_num_elements(1);
  verify_single_element(1u, "QUAD_4_2D", EntityIdVector{1,2,3,4});
}

TEST_F(TestTextMesh2d, twoSprings)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,SPRING_2,1,2\n"
                         "0,2,SPRING_2,2,3";
  setup_text_mesh(meshDesc);

  verify_num_elements(2);
  verify_single_element(1u, "SPRING_2", EntityIdVector{1,2});
  verify_single_element(2u, "SPRING_2", EntityIdVector{2,3});
}

TEST_F(TestTextMesh2d, threeQuadsWithCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4\n"
                         "0,2,QUAD_4_2D,2,3,5,6\n"
                         "0,3,QUAD_4_2D,5,7,8,6";
  std::vector<double> coordinates = { 0,0, 1,0, 1,1, 0,1, 2,0, 2,1, 3,0, 3,1 };

  setup_text_mesh(meshDesc, coordinates);

  verify_num_elements(3);
  verify_single_element(1u, "QUAD_4_2D", EntityIdVector{1,2,3,4});
  verify_single_element(2u, "QUAD_4_2D", EntityIdVector{2,3,5,6});
  verify_single_element(3u, "QUAD_4_2D", EntityIdVector{5,7,8,6});
  verify_coordinates(EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8}, coordinates);
}

TEST_F(TestTextMesh2d, twoQuadsWithCoordinatesParallel)
{
  //      4-----5-----6
  //      |     |     |
  //      |  1  |  2  |
  //      |     |     |
  //      1-----2-----3

  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,5,4\n"
                         "1,2,QUAD_4_2D,2,3,6,5";
  std::vector<double> coordinates = {
    0,0, 1,0, 2,0,
    0,1, 1,1, 2,1
  };

  setup_text_mesh(meshDesc, coordinates);

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, "QUAD_4_2D", EntityIdVector{1, 2, 5, 4});
    verify_shared_nodes(EntityIdVector{2,5}, 1);
    verify_coordinates(EntityIdVector{1, 2, 4, 5}, {0,0, 1,0, 0,1, 1,1});
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, "QUAD_4_2D", EntityIdVector{2, 3, 6, 5});
    verify_shared_nodes(EntityIdVector{2,5}, 0);
    verify_coordinates(EntityIdVector{2, 3, 5, 6}, {1,0, 2,0, 1,1, 2,1});
  }
}

TEST_F(TestTextMesh2d, twoQuadOneShellParallel)
{
  if (get_parallel_size() != 3) return;
  int rank = get_parallel_rank();

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4\n"
                         "1,2,QUAD_4_2D,3,4,5,6\n"
                         "2,3,SHELL_LINE_2,3,4";
  setup_text_mesh(meshDesc);

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, "QUAD_4_2D", EntityIdVector{1,2,3,4});
    verify_shared_nodes(EntityIdVector{3,4}, 1);
    verify_shared_nodes(EntityIdVector{3,4}, 2);
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, "QUAD_4_2D", EntityIdVector{3,4,5,6});
    verify_shared_nodes(EntityIdVector{3,4}, 0);
    verify_shared_nodes(EntityIdVector{3,4}, 2);
  }
  else if (rank == 2) {
    verify_num_elements(1);
    verify_single_element(3u, "SHELL_LINE_2", EntityIdVector{3,4});
    verify_shared_nodes(EntityIdVector{3,4}, 0);
    verify_shared_nodes(EntityIdVector{3,4}, 1);
  }
}


TEST_F(TestTextMesh1d, oneDimensionNotSupported)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,LINE_2_1D,1,2";
  EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
}


}  // namespace
