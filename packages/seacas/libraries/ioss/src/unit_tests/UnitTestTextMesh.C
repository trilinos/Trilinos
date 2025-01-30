// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <string>
#include <vector>

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"

#include "UnitTestIotmTextMeshFixture.h"

namespace {
  std::string add_coords_to_connectivity(const std::string         &textMeshConnectivityDesc,
                                         const std::vector<double> &coordVec)
  {
    std::stringstream coords;
    coords << "|coordinates:";

    for (double coord : coordVec) {
      coords << coord << ",";
    }

    std::string fullTextMeshDesc = textMeshConnectivityDesc + coords.str();
    return fullTextMeshDesc;
  }

  TEST_F(TestTextMesh, singleHex)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  }

  TEST_F(TestTextMesh, singleHexWithCoordinates)
  {
    if (get_parallel_size() != 1)
      return;

    std::string         meshDesc        = "0,1,HEX_8,1,2,3,4,5,6,7,8"
                                          "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1";
    std::vector<double> goldCoordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                           0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    EntityIdVector      nodeIds         = EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};

    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_single_element(1u, "HEX_8", nodeIds);
    verify_coordinates(nodeIds, goldCoordinates);
  }

  TEST_F(TestTextMesh, singleHexWithCoordinates_separatedNodeIds)
  {
    if (get_parallel_size() != 1)
      return;

    std::string         meshDesc        = "0,1,HEX_8,1,2,3,4,11,12,13,14";
    EntityIdVector      nodeIds         = EntityIdVector{1, 2, 3, 4, 11, 12, 13, 14};
    std::vector<double> goldCoordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                           0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};

    setup_text_mesh(add_coords_to_connectivity(meshDesc, goldCoordinates));

    verify_num_elements(1);
    verify_single_element(1u, "HEX_8", nodeIds);
    verify_coordinates(nodeIds, goldCoordinates);
  }

  TEST_F(TestTextMesh, twoHexesSerial)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12";
    setup_text_mesh(meshDesc);

    verify_num_elements(2);
    verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
    verify_single_element(2u, "HEX_8", EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});
  }

  TEST_F(TestTextMesh, twoTet10Serial)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,TET_10,1,2,3,4,5,6,7,8,9,10\n"
                           "0,2,TET_10,2,11,3,4,12,13,6,9,14,10";
    //                                       1       2      3        4          5          6 7 8
    std::vector<double> goldCoordinates = {0, 0, 0, 1, 0, 0, 0.5, 1, 0, 0.5, 0.5, 1, 0.5, 0, 0,
                                           0.75, 0.5, 0, 0.25, 0.5, 0, 0.25, 0.25, 0.5,
                                           //                                       9 10 11 12 13 14
                                           0.75, 0.25, 0.5, 0.5, 0.75, 0.5, 1.5, 0.5, 0, 1.25, 0.25,
                                           0, 1, 0.75, 0, 1, 0.5, 0.5};

    setup_text_mesh(add_coords_to_connectivity(meshDesc, goldCoordinates));

    verify_num_elements(2);
    verify_single_element(1u, "TET_10", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    verify_single_element(2u, "TET_10", EntityIdVector{2, 11, 3, 4, 12, 13, 6, 9, 14, 10});
  }

  TEST_F(TestTextMesh, twoHexDisconnectedWithCoordinatesAndParts)
  {
    if (get_parallel_size() != 1)
      return;

    std::string         meshDesc        = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_1\n"
                                          "0,2,HEX_8,9,10,11,12,13,14,15,16,block_2";
    std::vector<double> goldCoordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1,
                                           0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1,
                                           1, 0, 1, 1, 0, 0, 2, 1, 0, 2, 1, 1, 2, 0, 1, 2};

    setup_text_mesh(add_coords_to_connectivity(meshDesc, goldCoordinates));

    verify_num_elements(2);
    verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
    verify_single_element(2u, "HEX_8", EntityIdVector{9, 10, 11, 12, 13, 14, 15, 16});
    verify_coordinates(EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},
                       goldCoordinates);
    verify_part_membership({{"block_1", {1u}}, {"block_2", {2u}}});
  }

  TEST_F(TestTextMesh, twoHexDisconnectedWithDefaultParts)
  {
    if (get_parallel_size() != 1)
      return;

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

    if (get_parallel_size() != 1)
      return;

    std::string         meshDesc        = "0,1,SHELL_TRI_3,1,2,4,block_1\n"
                                          "0,2,SHELL_TRI_3,2,5,4,block_2\n"
                                          "0,3,SHELL_TRI_3,2,3,5,block_2";
    std::vector<double> goldCoordinates = {0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1, 0};

    setup_text_mesh(add_coords_to_connectivity(meshDesc, goldCoordinates));

    verify_num_elements(3);
    verify_single_element(1u, "SHELL_TRI_3", EntityIdVector{1, 2, 4});
    verify_single_element(2u, "SHELL_TRI_3", EntityIdVector{2, 5, 4});
    verify_single_element(3u, "SHELL_TRI_3", EntityIdVector{2, 3, 5});
    verify_coordinates(EntityIdVector{1, 2, 3, 4, 5}, goldCoordinates);
    verify_part_membership({{"block_1", {1u}}, {"block_2", {2u, 3u}}});
  }

  TEST_F(TestTextMesh, threeTriShellsWithDefaultParts)
  {
    //      4-----5             //
    //      |\  2 |\            //
    //      |  \  |  \          //
    //      | 1  \| 3  \        //
    //      1-----2-----3       //

    if (get_parallel_size() != 1)
      return;

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
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"block_" + get_topology_name("HEX_8"), 1u}});
  }

  TEST_F(TestTextMesh, partIds_oneDefaultPartTwoElems)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"block_" + get_topology_name("HEX_8"), 1u}});
  }

  TEST_F(TestTextMesh, partIds_twoDefaultParts)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16\n"
                           "0,3,SHELL_QUAD_4,17,18,19,20";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"block_" + get_topology_name("HEX_8"), 1u},
                     {"block_" + get_topology_name("SHELL_QUAD_4"), 2u}});
  }

  TEST_F(TestTextMesh, partIds_oneDefaultPartOneUserSpecifiedPart)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,my_cool_part";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"block_" + get_topology_name("HEX_8"), 1u}, {"my_cool_part", 2u}});
  }

  TEST_F(TestTextMesh, partIds_samePartTwice)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,my_cool_part\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,my_cool_part";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"my_cool_part", 1u}});
  }

  TEST_F(TestTextMesh, partIds_orderingIsByLine)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,2,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partOne\n"
                           "0,1,HEX_8,9,10,11,12,13,14,15,16,partTwo";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"partOne", 1u}, {"partTwo", 2u}});
  }

  TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_onePart)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_101";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"block_101", 101u}});
  }

  TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_twoParts)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_101\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,block_201";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"block_101", 101u}, {"block_201", 201u}});
  }

  TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_withDefaultPartFirst)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,block_101";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"block_" + get_topology_name("HEX_8"), 1u}, {"block_101", 101u}});
  }

  TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_withDefaultPartSecond)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_101\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"block_101", 101u}, {"block_" + get_topology_name("HEX_8"), 1u}});
  }

  TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_withDefaultPartIdCollision)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,block_1";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"block_" + get_topology_name("HEX_8"), 2u}, {"block_1", 1u}});
  }

  TEST_F(TestTextMesh, partIds_userSpecifiedPartId_onePart)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,partThree,3";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"partThree", 3u}});
    verify_num_elements(1);
    verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
    verify_part_membership({{"partThree", {1u}}});
  }

  TEST_F(TestTextMesh, partIds_userSpecifiedPartId_twoParts)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,partFive,5";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"partThree", 3u}, {"partFive", 5u}});
  }

  TEST_F(TestTextMesh, partIds_userSpecifiedPartId_twoPartsSameId)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partFour,4\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,partFour,4";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"partFour", 4u}});
  }

  TEST_F(TestTextMesh, partIds_userSpecifiedPartId_samePartDifferentIds)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partFour,4\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,partFour,5";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, partIds_userSpecifiedPartId_withDefaultPartIdCollision)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,partOne,1";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"block_" + get_topology_name("HEX_8"), 2u}, {"partOne", 1u}});
  }

  TEST_F(TestTextMesh, partIds_userSpecifiedPartId_withDefaultPart)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"partThree", 3u}, {"block_" + get_topology_name("HEX_8"), 1u}});
  }

  TEST_F(TestTextMesh, partIds_userSpecifiedPartId_withExodusPart)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,block_4";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"partThree", 3u}, {"block_4", 4u}});
  }

  TEST_F(TestTextMesh, partIds_userSpecifiedPartId_collidesWithExodusPart)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,block_3";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, partIds_userSpecifiedPartId_collidesWithPreviousSpec)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThreeA,3\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,partThreeB,3";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, partIds_userSpecifiedPartId_forExodusPart)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_2,3\n";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, partIds_shortPartNamesAreValid)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,a\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,b,3";
    setup_text_mesh(meshDesc);

    verify_part_ids({{"a", 1u}, {"b", 3u}});
  }

  TEST_F(TestTextMesh, partIds_integerPartNamesAreInvalid)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, twoHexesParallel)
  {
    if (get_parallel_size() != 2)
      return;
    int rank = get_parallel_rank();

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "1,2,HEX_8,5,6,7,8,9,10,11,12";
    setup_text_mesh(meshDesc);

    if (rank == 0) {
      verify_num_elements(1);
      verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
      verify_shared_nodes(EntityIdVector{5, 6, 7, 8}, 1);
    }
    else if (rank == 1) {
      verify_num_elements(1);
      verify_single_element(2u, "HEX_8", EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});
      verify_shared_nodes(EntityIdVector{5, 6, 7, 8}, 0);
    }
  }

  TEST_F(TestTextMesh, twoQuadShellsWithCoordinatesParallel)
  {
    //      4-----5-----6
    //      |     |     |
    //      |  1  |  2  |
    //      |     |     |
    //      1-----2-----3

    if (get_parallel_size() != 2)
      return;
    int rank = get_parallel_rank();

    std::string         meshDesc        = "0,1,SHELL_QUAD_4,1,2,5,4\n"
                                          "1,2,SHELL_QUAD_4,2,3,6,5";
    std::vector<double> goldCoordinates = {0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1, 0, 2, 1, 0};

    setup_text_mesh(add_coords_to_connectivity(meshDesc, goldCoordinates));

    if (rank == 0) {
      verify_num_elements(1);
      verify_single_element(1u, "SHELL_QUAD_4", EntityIdVector{1, 2, 5, 4});
      verify_coordinates(EntityIdVector{1, 2, 4, 5}, {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0});
    }
    else if (rank == 1) {
      verify_num_elements(1);
      verify_single_element(2u, "SHELL_QUAD_4", EntityIdVector{2, 3, 6, 5});
      verify_coordinates(EntityIdVector{2, 3, 5, 6}, {1, 0, 0, 2, 0, 0, 1, 1, 0, 2, 1, 0});
    }
  }

  TEST_F(TestTextMesh, threeTriShellsWithCoordinatesParallel)
  {
    //      4-----5             //
    //      |\  2 |\            //
    //      |  \  |  \          //
    //      | 1  \| 3  \        //
    //      1-----2-----3       //

#ifndef NDEBUG
    GTEST_SKIP();
#endif

    if (get_parallel_size() != 2)
      return;
    int rank = get_parallel_rank();

    std::string         meshDesc        = "0,1,SHELL_TRI_3,1,2,4\n"
                                          "0,2,SHELL_TRI_3,2,5,4\n"
                                          "1,3,SHELL_TRI_3,2,3,5";
    std::vector<double> goldCoordinates = {0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1, 0};

    setup_text_mesh(add_coords_to_connectivity(meshDesc, goldCoordinates));

    if (rank == 0) {
      verify_num_elements(2);
      verify_single_element(1u, "SHELL_TRI_3", EntityIdVector{1, 2, 4});
      verify_single_element(2u, "SHELL_TRI_3", EntityIdVector{2, 5, 4});
      verify_coordinates(EntityIdVector{1, 2, 4, 5}, {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0});
    }
    else if (rank == 1) {
      verify_num_elements(1);
      verify_single_element(3u, "SHELL_TRI_3", EntityIdVector{2, 3, 5});
      // parallel consistency checks in IOSS fails because field access in
      // parallel needs to test against same GroupingEntity type
      verify_single_element(3u, "SHELL_TRI_3", EntityIdVector{2, 3, 5});
      verify_coordinates(EntityIdVector{2, 2, 3, 5}, {1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 1, 0});
    }
  }

  TEST_F(TestTextMesh, singleHexWithSpaces)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0, 1, HEX_8, 1, 2, 3, 4, 5, 6, 7, 8";
    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  }

  TEST_F(TestTextMesh, singleHexWithCoordinatesAndSpaces)
  {
    if (get_parallel_size() != 1)
      return;

    std::string         meshDesc        = "0, 1, HEX_8, 1, 2, 3, 4, 5, 6, 7, 8";
    EntityIdVector      nodeIds         = EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<double> goldCoordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                           0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};

    setup_text_mesh(add_coords_to_connectivity(meshDesc, goldCoordinates));

    verify_num_elements(1);
    verify_single_element(1u, "HEX_8", nodeIds);
    verify_coordinates(nodeIds, goldCoordinates);
  }

  TEST_F(TestTextMesh, singleHexWithLowerCase)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,Hex_8,1,2,3,4,5,6,7,8";
    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  }

  TEST_F(TestTextMesh, singleHexWithCoordinatesAndLowerCase)
  {
    if (get_parallel_size() != 1)
      return;

    std::string         meshDesc        = "0,1,Hex_8,1,2,3,4,5,6,7,8";
    EntityIdVector      nodeIds         = EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<double> goldCoordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                           0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};

    setup_text_mesh(add_coords_to_connectivity(meshDesc, goldCoordinates));

    verify_num_elements(1);
    verify_single_element(1u, "HEX_8", nodeIds);
    verify_coordinates(nodeIds, goldCoordinates);
  }

  TEST_F(TestTextMesh, tooFewNodes)
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, tooFewNodesWithCoordinates)
  {
    std::string         meshDesc    = "0,1,HEX_8,1,2,3,4,5,6,7";
    std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1,
                                       0, 0, 0, 1, 1, 0, 1, 1, 1, 1};
    EXPECT_THROW(setup_text_mesh(add_coords_to_connectivity(meshDesc, coordinates)),
                 std::logic_error);
  }

  TEST_F(TestTextMesh, tooFewCoordinates)
  {
    std::string         meshDesc    = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                       0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1};
    EXPECT_THROW(setup_text_mesh(add_coords_to_connectivity(meshDesc, coordinates)),
                 std::logic_error);
  }

  TEST_F(TestTextMesh, tooManyNodes)
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9,10";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, tooManyNodesWithCoordinates)
  {
    std::string         meshDesc    = "0,1,HEX_8,1,2,3,4,5,6,7,8,9,10";
    std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1,
                                       1, 0, 1, 1, 1, 1, 0, 1, 1, 2, 0, 0, 2, 1, 0};
    EXPECT_THROW(setup_text_mesh(add_coords_to_connectivity(meshDesc, coordinates)),
                 std::logic_error);
  }

  TEST_F(TestTextMesh, tooManyCoordinates)
  {
    std::string         meshDesc    = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0,
                                       0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 52};
    EXPECT_THROW(setup_text_mesh(add_coords_to_connectivity(meshDesc, coordinates)),
                 std::logic_error);
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
    std::string         meshDesc = "0,1,";
    std::vector<double> coordinates;
    EXPECT_THROW(setup_text_mesh(add_coords_to_connectivity(meshDesc, coordinates)),
                 std::logic_error);
  }

  TEST_F(TestTextMesh, invalidTopology)
  {
    std::string meshDesc = "0,1,invalid,1";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, invalidTopologyWithCoordinates)
  {
    std::string         meshDesc    = "0,1,invalid,1";
    std::vector<double> coordinates = {0, 0, 0};
    EXPECT_THROW(setup_text_mesh(add_coords_to_connectivity(meshDesc, coordinates)),
                 std::logic_error);
  }

  TEST_F(TestTextMesh, mixedSpatialDim)
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,QUAD_4_2D,5,6,7,8";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, mixedSpatialDimWithCoordinates)
  {
    std::string         meshDesc    = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                                      "0,2,QUAD_4_2D,5,6,7,8";
    std::vector<double> coordinates = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                       0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    EXPECT_THROW(setup_text_mesh(add_coords_to_connectivity(meshDesc, coordinates)),
                 std::logic_error);
  }

  TEST_F(TestTextMesh, spatialDimInconsistentWithMetaData)
  {
    std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, spatialDimInconsistentWithMetaDataWithCoordinates)
  {
    std::string         meshDesc    = "0,1,QUAD_4_2D,1,2,3,4";
    std::vector<double> coordinates = {0, 0, 1, 0, 1, 1, 0, 1};
    EXPECT_THROW(setup_text_mesh(add_coords_to_connectivity(meshDesc, coordinates)),
                 std::logic_error);
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

  TEST_F(TestTextMesh, incompatibleParsedAndEnforcedDimensions)
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|dimension:2";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, multipleDimensionSpecification)
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|dimension:3|dimension:2";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, invalidDimensionOptionSyntax)
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|dimension:3:2";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, particleWithCoordinates)
  {
    if (get_parallel_size() != 1)
      return;

    std::string         meshDesc        = "0,1,PARTICLE,1|coordinates:2,0,0";
    std::vector<double> goldCoordinates = {2, 0, 0};
    EntityIdVector      nodeIds         = EntityIdVector{1};
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_elements(1);
    verify_single_element(1u, "PARTICLE", EntityIdVector{1});
    verify_coordinates(nodeIds, goldCoordinates);
  }

  TEST_F(TestTextMesh, particleWithCoordinatesAlreadyProvided)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,PARTICLE,1|coordinates:2,0,0|coordinates:3,0,0";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, particleHex)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,PARTICLE,1\n"
                           "0,2,HEX_8,2,3,4,5,6,7,8,9";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_elements(2);
    verify_single_element(1u, "PARTICLE", EntityIdVector{1});
    verify_single_element(2u, "HEX_8", EntityIdVector{2, 3, 4, 5, 6, 7, 8, 9});
  }

  TEST_F(TestTextMesh, particleHexWithCoordinates)
  {
    if (get_parallel_size() != 1)
      return;

    std::string         meshDesc        = "0,1,PARTICLE,1\n"
                                          "0,2,HEX_8,2,3,4,5,6,7,8,9";
    std::vector<double> goldCoordinates = {2, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1,
                                           0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    EXPECT_NO_THROW(setup_text_mesh(add_coords_to_connectivity(meshDesc, goldCoordinates)));

    verify_num_elements(2);
    verify_single_element(1u, "PARTICLE", EntityIdVector{1});
    verify_single_element(2u, "HEX_8", EntityIdVector{2, 3, 4, 5, 6, 7, 8, 9});
    verify_coordinates(EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8, 9}, goldCoordinates);
  }

  TEST_F(TestTextMesh, hexWithSideset_standardName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=surface_1; data=1,3";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1, SideVector{{1, 3}});
  }

  TEST_F(TestTextMesh, hexWithSideset_nonStandardName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=my_ss; data=1,2";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_sidesets(1);
    verify_single_sideset("my_ss", 1, SideVector{{1, 2}});
  }

  TEST_F(TestTextMesh, hexWithSideset_noName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:data=1,5";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1, SideVector{{1, 5}});
  }

  TEST_F(TestTextMesh, hexWithSideset_noData)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_sidesets(0);
  }

  TEST_F(TestTextMesh, hexWithSideset_noDataMixedCase)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sIdESeT:";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_sidesets(0);
  }

  TEST_F(TestTextMesh, hexWithSideset_onlyName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=my_ss";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_sidesets(1);
    verify_single_sideset("my_ss", 1, SideVector{});
  }

  TEST_F(TestTextMesh, hexWithSideset_invalidPairs)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=my_ss; data=1,1, 2";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithSideset_invalidName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=block_1; data=1,1";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithSideset_nameCollisionWithElementBlock)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,my_part|sideset:name=my_part; data=1,1";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithSideset_nameCollisionWithSideset)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:my_part,1,1|sideset:name=my_part; data=1,2";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithSideset_nameCollisionWithSideset_MixedCase)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:my_part,1,1|sideset:name=MY_PART; data=1,2";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithSideset_standardNameWithInvalidId)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=surface_0; data=1,1";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithSideset_nonStandardNamesPreserveIdOrder)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=first_ss; "
                           "data=1,1|sideset:name=second_ss; data=1,2";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_sidesets(2);
    verify_single_sideset("first_ss", 1, SideVector{{1, 1}});
    verify_single_sideset("second_ss", 2, SideVector{{1, 2}});
  }

  TEST_F(TestTextMesh, hexWithSideset_nonStandardNameDoesNotAffectAssignedId)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=ss_5; data=1,1";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_sidesets(1);
    verify_single_sideset("ss_5", 1, SideVector{{1, 1}});
  }

  TEST_F(TestTextMesh, hexWithSideset_standardNameHasPrecedenceForIdAssignment)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=my_ss; data=1,1|sideset:name=surface_1; data=1,2";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_sidesets(2);
    verify_single_sideset("surface_1", 1, SideVector{{1, 2}});
    verify_single_sideset("my_ss", 2, SideVector{{1, 1}});
  }

  TEST_F(TestTextMesh, hexWithSideset_standardNamesPreserveIds)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=surface_5; "
                           "data=1,5|sideset:name=surface_1; data=1,1";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_sidesets(2);
    verify_single_sideset("surface_1", 1, SideVector{{1, 1}});
    verify_single_sideset("surface_5", 5, SideVector{{1, 5}});
  }

  TEST_F(TestTextMesh, hexWithSideset_idAssignmentOrderIsStandardThenEmptyThenNonStandardName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=my_ss;data=1,5|sideset:name="
                           "surface_1;data=1,1|sideset:data=1,3";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_sidesets(3);
    verify_single_sideset("surface_1", 1, SideVector{{1, 1}});
    verify_single_sideset("surface_2", 2, SideVector{{1, 3}});
    verify_single_sideset("my_ss", 3, SideVector{{1, 5}});
  }

  TEST_F(TestTextMesh, hexWithSideset_referenceToInvalidElement)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=surface_1; data=2,1";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithSideset_referenceToInvalidOrdinal)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=surface_1; data=2,-1";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, twoHexesWithSpanningSideset)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12|sideset:name=surface_1; data=1,1, 2,1";
    setup_text_mesh(meshDesc);

    verify_num_elements(2);
    verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
    verify_single_element(2u, "HEX_8", EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1, SideVector{{1, 1}, {2, 1}});
  }

  TEST_F(TestTextMesh, twoHexesWithSpanningSidesetAndStandardName_splitByBlock)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
        "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2|sideset:name=surface_1; data=1,1, 2,1; split=block";
    setup_text_mesh(meshDesc);

    verify_num_elements(2);
    verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
    verify_single_element(2u, "HEX_8", EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});

    std::vector<std::string> sidesetSubsets;
    std::string              prefix("surface_");
    sidesetSubsets.push_back(prefix + std::string("block_1") + "_" + get_topology_name("QUAD_4") +
                             "_1");
    sidesetSubsets.push_back(prefix + std::string("block_2") + "_" + get_topology_name("QUAD_4") +
                             "_1");

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1, sidesetSubsets, SideVector{{1, 1}, {2, 1}});
  }

  TEST_F(TestTextMesh, twoHexesWithSpanningSidesetAndStandardName_splitByTopology)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2|sideset:name=surface_1; data=1,1, "
                           "2,1; split=topology";
    setup_text_mesh(meshDesc);

    verify_num_elements(2);
    verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
    verify_single_element(2u, "HEX_8", EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});

    std::vector<std::string> sidesetSubsets;
    sidesetSubsets.push_back("surface_" + get_topology_name("HEX_8") + "_" +
                             get_topology_name("QUAD_4") + "_1");

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1, sidesetSubsets, SideVector{{1, 1}, {2, 1}});
  }

  TEST_F(TestTextMesh, twoHexesWithSpanningSidesetAndNonStandardName_splitByBlock)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
        "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2|sideset:name=my_ss; data=1,1, 2,1; split=block";
    setup_text_mesh(meshDesc);

    verify_num_elements(2);
    verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
    verify_single_element(2u, "HEX_8", EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});

    std::vector<std::string> sidesetSubsets;
    std::string              prefix("my_ss_");
    sidesetSubsets.push_back(prefix + std::string("block_1") + "_" + get_topology_name("QUAD_4"));
    sidesetSubsets.push_back(prefix + std::string("block_2") + "_" + get_topology_name("QUAD_4"));

    verify_num_sidesets(1);
    verify_single_sideset("my_ss", 1, sidesetSubsets, SideVector{{1, 1}, {2, 1}});
  }

  TEST_F(TestTextMesh, twoHexesWithSpanningSidesetAndNonStandardName_splitByTopology)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
        "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2|sideset:name=my_ss; data=1,1, 2,1; split=topology";
    setup_text_mesh(meshDesc);

    verify_num_elements(2);
    verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
    verify_single_element(2u, "HEX_8", EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});

    std::vector<std::string> sidesetSubsets;
    sidesetSubsets.push_back("my_ss_" + get_topology_name("HEX_8") + "_" +
                             get_topology_name("QUAD_4"));

    verify_num_sidesets(1);
    verify_single_sideset("my_ss", 1, sidesetSubsets, SideVector{{1, 1}, {2, 1}});
  }

  TEST_F(TestTextMesh, twoHexesWithSpanningSideset_invalidSplitType)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2|sideset:name=surface_1; data=1,1, "
                           "2,1; split=blahblah";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, pyramidWithStandardName_splitByBlock)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,PYRAMID_5,1,2,3,4,5,block_1|sideset:name=surface_1; data=1,1, 1,2, "
                           "1,3, 1,4, 1,5; split=block";
    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_single_element(1u, "PYRAMID_5", EntityIdVector{1, 2, 3, 4, 5});

    std::vector<std::string> sidesetSubsets;
    std::string              prefix("surface_");
    sidesetSubsets.push_back(prefix + std::string("block_1") + "_" + get_topology_name("QUAD_4") +
                             "_1");
    sidesetSubsets.push_back(prefix + std::string("block_1") + "_" + get_topology_name("TRI_3") +
                             "_1");

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1, sidesetSubsets,
                          SideVector{{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}});
  }

  TEST_F(TestTextMesh, pyramidWithStandardName_splitByTopology)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,PYRAMID_5,1,2,3,4,5,block_1|sideset:name=surface_1; data=1,1, 1,2, "
                           "1,3, 1,4, 1,5; split=topology";
    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_single_element(1u, "PYRAMID_5", EntityIdVector{1, 2, 3, 4, 5});

    std::vector<std::string> sidesetSubsets;
    std::string              prefix("surface_");
    sidesetSubsets.push_back(prefix + get_topology_name("PYRAMID_5") + "_" +
                             get_topology_name("QUAD_4") + "_1");
    sidesetSubsets.push_back(prefix + get_topology_name("PYRAMID_5") + "_" +
                             get_topology_name("TRI_3") + "_1");

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1, sidesetSubsets,
                          SideVector{{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}});
  }

  TEST_F(TestTextMesh, twoHexesWithSpanningSidesetAndPartiallyOverlayedSideset)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12|sideset:name=surface_1; data=1,1, "
                           "2,1|sideset:name=surface_2; data=2,1";
    setup_text_mesh(meshDesc);

    verify_num_elements(2);
    verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
    verify_single_element(2u, "HEX_8", EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});

    verify_num_sidesets(2);
    verify_single_sideset("surface_1", 1, SideVector{{1, 1}, {2, 1}});
    verify_single_sideset("surface_2", 2, SideVector{{2, 1}});
  }

  TEST_F(TestTextMesh, twoHexesWithSpanningSideset_Parallel)
  {
    if (get_parallel_size() != 2)
      return;
    int rank = get_parallel_rank();

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "1,2,HEX_8,5,6,7,8,9,10,11,12|sideset:name=surface_1; data=1,1, 2,1";
    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_num_sidesets(1);

    if (rank == 0) {
      verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
      verify_single_sideset("surface_1", 1, SideVector{{1, 1}});
    }
    else if (rank == 1) {
      verify_single_element(2u, "HEX_8", EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});
      verify_single_sideset("surface_1", 1, SideVector{{2, 1}});
    }
  }

  TEST_F(TestTextMesh, hexWithNodeset_standardName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=nodelist_1; data=1,2,3";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_nodesets(1);
    verify_single_nodeset("nodelist_1", 1, EntityIdVector{1, 2, 3});
  }

  TEST_F(TestTextMesh, hexWithNodeset_nonStandardName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=my_ns; data=1,2";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_nodesets(1);
    verify_single_nodeset("my_ns", 1, EntityIdVector{1, 2});
  }

  TEST_F(TestTextMesh, hexWithNodeset_noName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:data=1,5";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_nodesets(1);
    verify_single_nodeset("nodelist_1", 1, EntityIdVector{1, 5});
  }

  TEST_F(TestTextMesh, hexWithNodeset_noData)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_nodesets(0);
  }

  TEST_F(TestTextMesh, hexWithNodeset_noDataMixedCase)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nOdESeT:";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_nodesets(0);
  }

  TEST_F(TestTextMesh, hexWithNodeset_onlyName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=my_ns";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_nodesets(1);
    verify_single_nodeset("my_ns", 1, EntityIdVector{});
  }

  TEST_F(TestTextMesh, hexWithNodeset_invalidName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=block_1; data=2,3";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithNodeset_invalidName2)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=surface_1; data=3,4";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithNodeset_nameCollisionWithElementBlock)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,my_part|nodeset:name=my_part;data=3,6";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithNodeset_nameCollisionWithSideset)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:my_part,1,2|nodeset:name=my_part;data=8";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithNodeset_nameCollisionWithNodeset_MixedCase)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=my_part;data=1,1|nodeset:name=MY_PART;data=3,4";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithNodeset_standardNameWithInvalidId)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=nodelist_0; data=1,1";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithNodeset_nonStandardNamesPreserveIdOrder)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=first_ss;data=1,2|nodeset:name="
                           "second_ss;data=3,4,5";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_nodesets(2);
    verify_single_nodeset("first_ss", 1, EntityIdVector{1, 2});
    verify_single_nodeset("second_ss", 2, EntityIdVector{3, 4, 5});
  }

  TEST_F(TestTextMesh, hexWithNodeset_nonStandardNameDoesNotAffectAssignedId)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=ns_5; data=7,8";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_nodesets(1);
    verify_single_nodeset("ns_5", 1, EntityIdVector{7, 8});
  }

  TEST_F(TestTextMesh, hexWithNodeset_standardNameHasPrecedenceForIdAssignment)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=my_ns; data=1,2|nodeset:name=nodelist_1; data=3,4";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_nodesets(2);
    verify_single_nodeset("nodelist_1", 1, EntityIdVector{3, 4});
    verify_single_nodeset("my_ns", 2, EntityIdVector{1, 2});
  }

  TEST_F(TestTextMesh, hexWithNodeset_standardNamesPreserveIds)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=nodelist_5; "
                           "data=1,5|nodeset:name=nodelist_1; data=3,4";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_nodesets(2);
    verify_single_nodeset("nodelist_1", 1, EntityIdVector{3, 4});
    verify_single_nodeset("nodelist_5", 5, EntityIdVector{1, 5});
  }

  TEST_F(TestTextMesh, hexWithNodeset_idAssignmentOrderIsStandardThenEmptyThenNonStandardName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=my_ns;data=5,6,7,8|nodeset:name="
                           "nodelist_1;data=1,2|nodeset:data=3,4";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_nodesets(3);
    verify_single_nodeset("nodelist_1", 1, EntityIdVector{1, 2});
    verify_single_nodeset("nodelist_2", 2, EntityIdVector{3, 4});
    verify_single_nodeset("my_ns", 3, EntityIdVector{5, 6, 7, 8});
  }

  TEST_F(TestTextMesh, hexWithNodeset_referenceToInvalidNode)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=nodelist_1; data=9";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, twoHexesWithSpanningNodeset)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
        "0,2,HEX_8,5,6,7,8,9,10,11,12|nodeset:name=nodelist_1; data=1,2,5,6,9,10";
    setup_text_mesh(meshDesc);

    verify_num_elements(2);
    verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
    verify_single_element(2u, "HEX_8", EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});

    verify_num_nodesets(1);
    verify_single_nodeset("nodelist_1", 1, EntityIdVector{1, 2, 5, 6, 9, 10});
  }

  TEST_F(TestTextMesh, twoHexesWithSpanningNodeset_Parallel)
  {
    if (get_parallel_size() != 2)
      return;
    int rank = get_parallel_rank();

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
        "1,2,HEX_8,5,6,7,8,9,10,11,12|nodeset:name=nodelist_1; data=1,2,5,6,9,10";
    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_num_nodesets(1);

    if (rank == 0) {
      verify_single_element(1u, "HEX_8", EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
      verify_single_nodeset("nodelist_1", 1, EntityIdVector{1, 2, 5, 6});
    }
    else if (rank == 1) {
      verify_single_element(2u, "HEX_8", EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});
      verify_single_nodeset("nodelist_1", 1, EntityIdVector{5, 6, 9, 10});
    }
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_standardName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1|assembly:name=assembly_1; type=block; member=block_1";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_assemblies(1);
    verify_single_assembly("assembly_1", 1, {"block_1"});
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_nonStandardName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1|assembly:name=my_assembly; type=block; member=block_1";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_assemblies(1);
    verify_single_assembly("my_assembly", 1, {"block_1"});
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_noName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1|assembly:type=block; member=block_1";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_assemblies(1);
    verify_single_assembly("assembly_1", 1, {"block_1"});
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_noData)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|assembly:";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_assemblies(0);
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_noDataMixedCase)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|aSsEmBlY:";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_assemblies(0);
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_mustSpecifyTypeWithName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1|assembly:name=my_assembly;";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_invalidName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|assembly:name=block_1;";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_invalidName2)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|assembly:name=surface_1;";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_invalidName3)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|assembly:name=nodelist_1;";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_nameCollisionWithElementBlock)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,my_part|assembly:name=my_part; type=block";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_nameCollisionWithSideset)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8|sideset:name=my_part;data=1,2|assembly:name=my_part; type=block";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_nameCollisionWithNodeset)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=my_part;data=1,2|assembly:name=my_part; type=block";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_nameCollisionWithNodeset_MixedCase)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8|nodeset:name=my_part;data=1,2|assembly:name=MY_PART; type=block";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_standardNameWithInvalidId)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|assembly:name=assembly_0; type=block";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_nonStandardNamesPreserveIdOrder)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|assembly:name=first; "
                           "type=block|assembly:name=second; type=block";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_assemblies(2);
    verify_single_assembly("first", 1, {});
    verify_single_assembly("second", 2, {});
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_nonStandardNameDoesNotAffectAssignedId)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1|assembly:name=ass_5; type=block; member=block_1";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_assemblies(1);
    verify_single_assembly("ass_5", 1, {"block_1"});
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_standardNameHasPrecedenceForIdAssignment)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|assembly:name=my_assembly; "
                           "type=block|assembly:name=assembly_1; type=block";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_assemblies(2);
    verify_single_assembly("assembly_1", 1, {});
    verify_single_assembly("my_assembly", 2, {});
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_standardNamesPreserveIds)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8|assembly:name=assembly_5; "
                           "type=block|assembly:name=assembly_1; type=block";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_assemblies(2);
    verify_single_assembly("assembly_1", 1, {});
    verify_single_assembly("assembly_5", 5, {});
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_idAssignmentOrderIsStandardThenEmptyThenNonStandardName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8"
                           "|assembly:name=my_assembly; type=block|assembly:name=assembly_1; "
                           "type=block|assembly:type=block";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_assemblies(3);
    verify_single_assembly("assembly_1", 1, {});
    verify_single_assembly("assembly_2", 2, {});
    verify_single_assembly("my_assembly", 3, {});
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_referenceToInvalidBlock)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1|assembly:name=assembly_1; type=block; member=block_2";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_deeplyNestedHomogeneous)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8, 1, 2, 3, 4, 5, 6, 7, 8,block_1\n"
        "0,2,HEX_8, 9,10,11,12,13,14,15,16,block_2\n"
        "0,3,HEX_8,17,18,19,20,21,22,23,24,block_3\n"
        "0,4,HEX_8,25,26,27,28,29,30,31,32,block_4"
        "|assembly:name=assembly_9000; type=assembly; member=assembly_9001,assembly_9002"
        "|assembly:name=assembly_9001; type=assembly; member=assembly_9003"
        "|assembly:name=assembly_9002; type=block   ; member=block_3"
        "|assembly:name=assembly_9003; type=assembly; member=assembly_9004,assembly_9005"
        "|assembly:name=assembly_9004; type=block   ; member=block_1,block_2"
        "|assembly:name=assembly_9005; type=block   ; member=block_4";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_assemblies(6);
    verify_single_assembly("assembly_9000", 9000, {"block_1", "block_2", "block_3", "block_4"});
    verify_single_assembly("assembly_9001", 9001, {"block_1", "block_2", "block_4"});
    verify_single_assembly("assembly_9002", 9002, {"block_3"});
    verify_single_assembly("assembly_9003", 9003, {"block_1", "block_2", "block_4"});
    verify_single_assembly("assembly_9004", 9004, {"block_1", "block_2"});
    verify_single_assembly("assembly_9005", 9005, {"block_4"});
  }

  TEST_F(TestTextMesh, hexWithBlockAssembly_deeplyNestedNonHomogeneous)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1"
        "|sideset:name=surface_1; data=1,3"
        "|sideset:name=surface_2; data=1,5"
        "|nodeset:name=nodelist_1; data=1,2,3,4"
        "|assembly:name=assembly_9000; type=assembly; member=assembly_9001,assembly_9002"
        "|assembly:name=assembly_9001; type=assembly; member=assembly_9003"
        "|assembly:name=assembly_9002; type=nodeset ; member=nodelist_1"
        "|assembly:name=assembly_9003; type=assembly; member=assembly_9004,assembly_9005"
        "|assembly:name=assembly_9004; type=sideset ; member=surface_1,surface_2"
        "|assembly:name=assembly_9005; type=block   ; member=block_1";
    EXPECT_NO_THROW(setup_text_mesh(meshDesc));

    verify_num_assemblies(6);
    verify_single_assembly("assembly_9000", 9000,
                           {"block_1", "nodelist_1", "surface_1", "surface_2"});
    verify_single_assembly("assembly_9001", 9001, {"block_1", "surface_1", "surface_2"});
    verify_single_assembly("assembly_9002", 9002, {"nodelist_1"});
    verify_single_assembly("assembly_9003", 9003, {"block_1", "surface_1", "surface_2"});
    verify_single_assembly("assembly_9004", 9004, {"surface_1", "surface_2"});
    verify_single_assembly("assembly_9005", 9005, {"block_1"});
  }

  TEST_F(TestTextMesh2d, singleQuad)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_single_element(1u, "QUAD_4_2D", EntityIdVector{1, 2, 3, 4});
  }

  TEST_F(TestTextMesh2d, twoSprings)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,SPRING_2,1,2\n"
                           "0,2,SPRING_2,2,3";
    setup_text_mesh(meshDesc);

    verify_num_elements(2);
    verify_single_element(1u, "SPRING_2", EntityIdVector{1, 2});
    verify_single_element(2u, "SPRING_2", EntityIdVector{2, 3});
  }

  TEST_F(TestTextMesh2d, threeQuadsWithCoordinates)
  {
    if (get_parallel_size() != 1)
      return;

    std::string         meshDesc        = "0,1,QUAD_4_2D,1,2,3,4\n"
                                          "0,2,QUAD_4_2D,2,3,5,6\n"
                                          "0,3,QUAD_4_2D,5,7,8,6";
    std::vector<double> goldCoordinates = {0, 0, 1, 0, 1, 1, 0, 1, 2, 0, 2, 1, 3, 0, 3, 1};

    setup_text_mesh(add_coords_to_connectivity(meshDesc, goldCoordinates));

    verify_num_elements(3);
    verify_single_element(1u, "QUAD_4_2D", EntityIdVector{1, 2, 3, 4});
    verify_single_element(2u, "QUAD_4_2D", EntityIdVector{2, 3, 5, 6});
    verify_single_element(3u, "QUAD_4_2D", EntityIdVector{5, 7, 8, 6});
    verify_coordinates(EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8}, goldCoordinates);
  }

  TEST_F(TestTextMesh2d, twoQuadsWithCoordinatesParallel)
  {
    //      4-----5-----6
    //      |     |     |
    //      |  1  |  2  |
    //      |     |     |
    //      1-----2-----3

    if (get_parallel_size() != 2)
      return;
    int rank = get_parallel_rank();

    std::string         meshDesc        = "0,1,QUAD_4_2D,1,2,5,4\n"
                                          "1,2,QUAD_4_2D,2,3,6,5";
    std::vector<double> goldCoordinates = {0, 0, 1, 0, 2, 0, 0, 1, 1, 1, 2, 1};

    setup_text_mesh(add_coords_to_connectivity(meshDesc, goldCoordinates));

    if (rank == 0) {
      verify_num_elements(1);
      verify_single_element(1u, "QUAD_4_2D", EntityIdVector{1, 2, 5, 4});
      verify_shared_nodes(EntityIdVector{2, 5}, 1);
      verify_coordinates(EntityIdVector{1, 2, 4, 5}, {0, 0, 1, 0, 0, 1, 1, 1});
    }
    else if (rank == 1) {
      verify_num_elements(1);
      verify_single_element(2u, "QUAD_4_2D", EntityIdVector{2, 3, 6, 5});
      verify_shared_nodes(EntityIdVector{2, 5}, 0);
      verify_coordinates(EntityIdVector{2, 3, 5, 6}, {1, 0, 2, 0, 1, 1, 2, 1});
    }
  }

  TEST_F(TestTextMesh2d, twoQuadOneShellParallel)
  {
    if (get_parallel_size() != 3)
      return;
    int rank = get_parallel_rank();

    std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4\n"
                           "1,2,QUAD_4_2D,3,4,5,6\n"
                           "2,3,SHELL_LINE_2,3,4";
    setup_text_mesh(meshDesc);

    if (rank == 0) {
      verify_num_elements(1);
      verify_single_element(1u, "QUAD_4_2D", EntityIdVector{1, 2, 3, 4});
      verify_shared_nodes(EntityIdVector{3, 4}, 1);
      verify_shared_nodes(EntityIdVector{3, 4}, 2);
    }
    else if (rank == 1) {
      verify_num_elements(1);
      verify_single_element(2u, "QUAD_4_2D", EntityIdVector{3, 4, 5, 6});
      verify_shared_nodes(EntityIdVector{3, 4}, 0);
      verify_shared_nodes(EntityIdVector{3, 4}, 2);
    }
    else if (rank == 2) {
      verify_num_elements(1);
      verify_single_element(3u, "SHELL_LINE_2", EntityIdVector{3, 4});
      verify_shared_nodes(EntityIdVector{3, 4}, 0);
      verify_shared_nodes(EntityIdVector{3, 4}, 1);
    }
  }

  TEST_F(TestTextMesh2d, singleQuadWithSideset_noName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4|sideset:data=1,1";
    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_single_element(1u, "QUAD_4_2D", EntityIdVector{1, 2, 3, 4});

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1, SideVector{{1, 1}});
  }

  TEST_F(TestTextMesh2d, singleQuadWithSideset_standardName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4|sideset:name=surface_2;data=1,1";
    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_single_element(1u, "QUAD_4_2D", EntityIdVector{1, 2, 3, 4});

    verify_num_sidesets(1);
    verify_single_sideset("surface_2", 2, SideVector{{1, 1}});
  }

  TEST_F(TestTextMesh2d, singleQuadWithSideset_nonStandardName)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4|sideset:name=my_ss;data=1,1";
    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_single_element(1u, "QUAD_4_2D", EntityIdVector{1, 2, 3, 4});

    verify_num_sidesets(1);
    verify_single_sideset("my_ss", 1, SideVector{{1, 1}});
  }

  TEST_F(TestTextMesh1d, oneDimensionNotSupported)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,LINE_2_1D,1,2";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMeshGraph, singleHex)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    setup_text_mesh_graph(meshDesc);

    Adjacency goldHexNeighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, -1});
    verify_side_adjacency({goldHexNeighbors});
  }

  TEST_F(TestTextMeshGraph, singleHexConnectedShell)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,SHELL_QUAD_4,5,6,7,8";
    setup_text_mesh_graph(meshDesc);

    Adjacency goldHexNeighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 1});
    Adjacency goldShellNeighbors(1u, Adjacency::SimpleNeighborVector{-1, 0});
    verify_side_adjacency({goldHexNeighbors, goldShellNeighbors});
  }

  TEST_F(TestTextMeshGraph, hexShellShell)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,SHELL_QUAD_4,5,6,7,8\n"
                           "0,3,SHELL_QUAD_4,5,6,7,8";
    ;
    setup_text_mesh_graph(meshDesc);

    Adjacency goldHexNeighbors(
        0u, Adjacency::NeighborVector{{0, -1}, {1, -1}, {2, -1}, {3, -1}, {4, -1}, {5, 1}, {5, 2}});
    Adjacency goldShell1Neighbors(1u, Adjacency::SimpleNeighborVector{-1, 0});
    Adjacency goldShell2Neighbors(2u, Adjacency::SimpleNeighborVector{-1, 0});
    verify_side_adjacency({goldHexNeighbors, goldShell1Neighbors, goldShell2Neighbors});
  }

  TEST_F(TestTextMeshGraph, hexShellShellHex)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,SHELL_QUAD_4,5,6,7,8\n"
                           "0,3,SHELL_QUAD_4,5,6,7,8\n"
                           "0,4,HEX_8,5,6,7,8,9,10,11,12\n";
    setup_text_mesh_graph(meshDesc);

    Adjacency goldHex1Neighbors(
        0u, Adjacency::NeighborVector{{0, -1}, {1, -1}, {2, -1}, {3, -1}, {4, -1}, {5, 1}, {5, 2}});
    Adjacency goldShell1Neighbors(1u, Adjacency::SimpleNeighborVector{3, 0});
    Adjacency goldShell2Neighbors(2u, Adjacency::SimpleNeighborVector{3, 0});
    Adjacency goldHex2Neighbors(
        3u, Adjacency::NeighborVector{{0, -1}, {1, -1}, {2, -1}, {3, -1}, {4, 1}, {4, 2}, {5, -1}});

    verify_side_adjacency(
        {goldHex1Neighbors, goldShell1Neighbors, goldShell2Neighbors, goldHex2Neighbors});
  }

  TEST_F(TestTextMeshGraph, hexShellShellHex_splitCoincidentShells)
  {
    if (get_parallel_size() != 2)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,SHELL_QUAD_4,5,6,7,8\n"
                           "1,3,SHELL_QUAD_4,5,6,7,8\n"
                           "1,4,HEX_8,5,6,7,8,9,10,11,12\n";
    EXPECT_THROW(setup_text_mesh_graph(meshDesc, {}, get_parallel_rank()), std::logic_error);
  }

  TEST_F(TestTextMeshGraph, hexShellShellHex_parallelWithLocalGraph)
  {
    if (get_parallel_size() != 3)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "1,2,SHELL_QUAD_4,5,6,7,8\n"
                           "1,3,SHELL_QUAD_4,5,6,7,8\n"
                           "2,4,HEX_8,5,6,7,8,9,10,11,12\n";
    setup_text_mesh_graph(meshDesc, {}, get_parallel_rank());

    Adjacency goldHex1Neighbors(
        0u, Adjacency::NeighborVector{{0, -1}, {1, -1}, {2, -1}, {3, -1}, {4, -1}, {5, 1}, {5, 2}});
    Adjacency goldShell1Neighbors(1u, Adjacency::SimpleNeighborVector{3, 0});
    Adjacency goldShell2Neighbors(2u, Adjacency::SimpleNeighborVector{3, 0});
    Adjacency goldHex2Neighbors(
        3u, Adjacency::NeighborVector{{0, -1}, {1, -1}, {2, -1}, {3, -1}, {4, 1}, {4, 2}, {5, -1}});

    verify_side_adjacency(
        {goldHex1Neighbors, goldShell1Neighbors, goldShell2Neighbors, goldHex2Neighbors});
  }

  TEST_F(TestTextMeshGraph, singleHexConnectedShell_highOrder)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20\n"
                           "0,2,SHELL_QUAD_8,5,6,7,8,17,18,19,20";
    setup_text_mesh_graph(meshDesc);

    Adjacency goldHexNeighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 1});
    Adjacency goldShellNeighbors(1u, Adjacency::SimpleNeighborVector{-1, 0});
    verify_side_adjacency({goldHexNeighbors, goldShellNeighbors});
  }

  TEST_F(TestTextMeshGraph, singleHexConnectedShell_highOrderWithRotation)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_20,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20\n"
                           "0,2,SHELL_QUAD_8,6,7,8,5,18,19,20,17";
    setup_text_mesh_graph(meshDesc);

    Adjacency goldHexNeighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 1});
    Adjacency goldShellNeighbors(1u, Adjacency::SimpleNeighborVector{-1, 0});
    verify_side_adjacency({goldHexNeighbors, goldShellNeighbors});
  }

  TEST_F(TestTextMeshGraph, twoHexConnected)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12";
    setup_text_mesh_graph(meshDesc);

    Adjacency goldHex1Neighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 1});
    Adjacency goldHex2Neighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 0, -1});
    verify_side_adjacency({goldHex1Neighbors, goldHex2Neighbors});
  }

  TEST_F(TestTextMeshGraph, twoHexConnected_parallelWithLocalGraph)
  {
    if (get_parallel_size() != 2)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "1,2,HEX_8,5,6,7,8,9,10,11,12";
    setup_text_mesh_graph(meshDesc, {}, get_parallel_rank());

    Adjacency goldHex1Neighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 1});
    Adjacency goldHex2Neighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 0, -1});
    verify_side_adjacency({goldHex1Neighbors, goldHex2Neighbors});
  }

  TEST_F(TestTextMeshGraph, twoHexConnected_parallelGlobal)
  {
    if (get_parallel_size() != 2)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "1,2,HEX_8,5,6,7,8,9,10,11,12";
    setup_text_mesh_graph(meshDesc, {});

    Adjacency goldHex1Neighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 1});
    Adjacency goldHex2Neighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 0, -1});
    verify_side_adjacency({goldHex1Neighbors, goldHex2Neighbors});
  }

  TEST_F(TestTextMeshGraph, threeHexConnected)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                           "0,3,HEX_8,9,10,11,12,13,14,15,16";
    setup_text_mesh_graph(meshDesc);

    Adjacency goldHex1Neighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 1});
    Adjacency goldHex2Neighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 0, 2});
    Adjacency goldHex3Neighbors(2u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 1, -1});
    verify_side_adjacency({goldHex1Neighbors, goldHex2Neighbors, goldHex3Neighbors});
  }

  TEST_F(TestTextMeshGraph, threeHexConnected_parallelWithLocalGraph)
  {
    if (get_parallel_size() != 3)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "1,2,HEX_8,5,6,7,8,9,10,11,12\n"
                           "2,3,HEX_8,9,10,11,12,13,14,15,16";
    setup_text_mesh_graph(meshDesc, {}, get_parallel_rank());

    if (get_parallel_rank() == 0) {
      Adjacency goldHex1Neighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 1});
      Adjacency goldHex2Neighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 0, -1});

      verify_side_adjacency({goldHex1Neighbors, goldHex2Neighbors});
    }
    else if (get_parallel_rank() == 1) {
      Adjacency goldHex1Neighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 1});
      Adjacency goldHex2Neighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 0, 2});
      Adjacency goldHex3Neighbors(2u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 1, -1});

      verify_side_adjacency({goldHex1Neighbors, goldHex2Neighbors, goldHex3Neighbors});
    }
    else if (get_parallel_rank() == 2) {
      Adjacency goldHex2Neighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 2});
      Adjacency goldHex3Neighbors(2u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 1, -1});

      verify_side_adjacency({goldHex2Neighbors, goldHex3Neighbors});
    }
  }

  TEST_F(TestTextMeshGraph, threeHexConnected_parallelGlobal)
  {
    if (get_parallel_size() != 3)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12\n"
                           "0,3,HEX_8,9,10,11,12,13,14,15,16";
    setup_text_mesh_graph(meshDesc, {});

    Adjacency goldHex1Neighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 1});
    Adjacency goldHex2Neighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 0, 2});
    Adjacency goldHex3Neighbors(2u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 1, -1});
    verify_side_adjacency({goldHex1Neighbors, goldHex2Neighbors, goldHex3Neighbors});
  }

  TEST_F(TestTextMeshGraph, twoHexConnected_separateBlocks)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
    setup_text_mesh_graph(meshDesc, {"block_1"});

    Adjacency goldHex1Neighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, -1});
    verify_side_adjacency({goldHex1Neighbors});
  }

  TEST_F(TestTextMeshGraph, twoHexConnected_separateBlocksWithInvalidAdjacency)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
    setup_text_mesh_graph(meshDesc, {"block_1"});

    Adjacency goldHex2Neighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, -1});
    EXPECT_THROW(verify_side_adjacency({goldHex2Neighbors}), std::logic_error);
  }

  TEST_F(TestTextMeshGraph, twoHexDisconnected)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_1\n"
                           "0,2,HEX_8,9,10,11,12,13,14,15,16,block_2";
    setup_text_mesh_graph(meshDesc);

    Adjacency goldHex1Neighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, -1});
    Adjacency goldHex2Neighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, -1});
    verify_side_adjacency({goldHex1Neighbors, goldHex2Neighbors});
  }

  TEST_F(TestTextMeshGraph, twoHexConnectedAndShellDisconnected)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                           "0,2,HEX_8,5, 6, 7, 8, 9,10,11,12\n"
                           "0,3,SHELL_QUAD_4,17,18,19,20";
    setup_text_mesh_graph(meshDesc);

    Adjacency goldHex1Neighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 1});
    Adjacency goldHex2Neighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 0, -1});
    Adjacency goldShellNeighbors(2u, Adjacency::SimpleNeighborVector{-1, -1});
    verify_side_adjacency({goldHex1Neighbors, goldHex2Neighbors, goldShellNeighbors});
  }

  TEST_F(TestTextMeshGraph, twoHexConnectedAndShellExternallyConnected)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                           "0,2,HEX_8,5, 6, 7, 8, 9,10,11,12\n"
                           "0,3,SHELL_QUAD_4,9,10,11,12";
    setup_text_mesh_graph(meshDesc);

    Adjacency goldHex1Neighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 1});
    Adjacency goldHex2Neighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 0, 2});
    Adjacency goldShellNeighbors(2u, Adjacency::SimpleNeighborVector{-1, 1});
    verify_side_adjacency({goldHex1Neighbors, goldHex2Neighbors, goldShellNeighbors});
  }

  TEST_F(TestTextMeshGraph, twoHexConnectedAndShellExternallyConnected_shellInSeparateBlock)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_1\n"
                           "0,2,HEX_8,5, 6, 7, 8, 9,10,11,12,block_1\n"
                           "0,3,SHELL_QUAD_4,9,10,11,12,block_2";
    setup_text_mesh_graph(meshDesc, {"block_2"});

    Adjacency goldShellNeighbors(2u, Adjacency::SimpleNeighborVector{-1, -1});
    verify_side_adjacency({goldShellNeighbors});
  }

  TEST_F(TestTextMeshGraph, twoHexConnectedAndShellInternallyConnected)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                           "0,2,HEX_8,5, 6, 7, 8, 9,10,11,12\n"
                           "0,3,SHELL_QUAD_4,5,6,7,8";
    setup_text_mesh_graph(meshDesc);

    Adjacency goldHex1Neighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 2});
    Adjacency goldHex2Neighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 2, -1});
    Adjacency goldShellNeighbors(2u, Adjacency::SimpleNeighborVector{1, 0});
    verify_side_adjacency({goldHex1Neighbors, goldHex2Neighbors, goldShellNeighbors});
  }

  TEST_F(TestTextMeshGraph, twoHexConnectedAndShellInternallyConnected_reorderShell)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                           "0,3,HEX_8,5, 6, 7, 8, 9,10,11,12\n"
                           "0,2,SHELL_QUAD_4,5,6,7,8";
    setup_text_mesh_graph(meshDesc);

    Adjacency goldHex1Neighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, 1});
    Adjacency goldHex2Neighbors(2u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 1, -1});
    Adjacency goldShellNeighbors(1u, Adjacency::SimpleNeighborVector{2, 0});
    verify_side_adjacency({goldHex1Neighbors, goldHex2Neighbors, goldShellNeighbors});
  }

  TEST_F(TestTextMeshGraph, twoHexConnectedAndShellInternallyConnected_shellInSeperateBlock)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_1\n"
                           "0,2,HEX_8,5, 6, 7, 8, 9,10,11,12,block_1\n"
                           "0,3,SHELL_QUAD_4,5,6,7,8,block_2";
    setup_text_mesh_graph(meshDesc, {"block_2"});

    Adjacency goldShellNeighbors(2u, Adjacency::SimpleNeighborVector{-1, -1});
    verify_side_adjacency({goldShellNeighbors});
  }

  TEST_F(TestTextMeshGraph, connectedPyramidAndHexAndTet)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,PYRAMID_5,1,2,3,4,5\n"
                           "0,2,HEX_8,1,4,3,2,6,9,8,7\n"
                           "0,3,TET_4,2,3,5,10";
    setup_text_mesh_graph(meshDesc);

    Adjacency goldPyramidNeighbors(0u, Adjacency::SimpleNeighborVector{-1, 2, -1, -1, 1});
    Adjacency goldHexNeighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 0, -1});
    Adjacency goldTetNeighbors(2u, Adjacency::SimpleNeighborVector{-1, -1, -1, 0});
    verify_side_adjacency({goldPyramidNeighbors, goldHexNeighbors, goldTetNeighbors});
  }

  TEST_F(TestTextMeshGraph, connectedPyramidAndHexAndTet_parallel)
  {
    if (get_parallel_size() != 3)
      return;

    std::string meshDesc = "0,1,PYRAMID_5,1,2,3,4,5\n"
                           "1,2,HEX_8,1,4,3,2,6,9,8,7\n"
                           "2,3,TET_4,2,3,5,10";
    setup_text_mesh_graph(meshDesc, {}, get_parallel_rank());

    Adjacency goldPyramidNeighbors(0u, Adjacency::SimpleNeighborVector{-1, 2, -1, -1, 1});
    Adjacency goldHexNeighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 0, -1});
    Adjacency goldTetNeighbors(2u, Adjacency::SimpleNeighborVector{-1, -1, -1, 0});
    verify_side_adjacency({goldPyramidNeighbors, goldHexNeighbors, goldTetNeighbors});
  }

  TEST_F(TestTextMeshGraph, connectedPyramidAndHexAndTet_graphPyramidAndHexBlocks)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,PYRAMID_5,1,2,3,4,5,block_1\n"
                           "0,2,HEX_8,1,4,3,2,6,9,8,7,block_2\n"
                           "0,3,TET_4,2,3,5,10,block_3";
    setup_text_mesh_graph(meshDesc, {"block_1", "block_2"});

    Adjacency goldPyramidNeighbors(0u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 1});
    Adjacency goldHexNeighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, 0, -1});
    verify_side_adjacency({goldPyramidNeighbors, goldHexNeighbors});
  }

  TEST_F(TestTextMeshGraph, connectedPyramidAndHexAndTet_graphHexAndTetBlocks)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,PYRAMID_5,1,2,3,4,5,block_1\n"
                           "0,2,HEX_8,1,4,3,2,6,9,8,7,block_2\n"
                           "0,3,TET_4,2,3,5,10,block_3";
    setup_text_mesh_graph(meshDesc, {"block_2", "block_3"});

    Adjacency goldHexNeighbors(1u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1, -1, -1});
    Adjacency goldTetNeighbors(2u, Adjacency::SimpleNeighborVector{-1, -1, -1, -1});
    verify_side_adjacency({goldHexNeighbors, goldTetNeighbors});
  }

  TEST_F(TestTextMeshSkin, singleHex_skin)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1|sideset:name=skinned; skin=block_1";
    setup_text_mesh(meshDesc);

    verify_num_sidesets(1);
    verify_single_sideset("skinned", 1, SideVector{{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}});
  }

  TEST_F(TestTextMeshSkin, singleHex_invalidComboSkinAndData)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1|sideset:name=skinned; data=1,2; skin=block_1";
    EXPECT_THROW(setup_text_mesh(meshDesc), std::logic_error);
  }

  TEST_F(TestTextMeshSkin, pyramidWithStandardName_skinAndSplitByBlock)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,PYRAMID_5,1,2,3,4,5,block_1|sideset:name=surface_1; skin=block_1; split=block";
    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_single_element(1u, "PYRAMID_5", EntityIdVector{1, 2, 3, 4, 5});

    std::vector<std::string> sidesetSubsets;
    std::string              prefix("surface_");
    sidesetSubsets.push_back(prefix + std::string("block_1") + "_" + get_topology_name("QUAD_4") +
                             "_1");
    sidesetSubsets.push_back(prefix + std::string("block_1") + "_" + get_topology_name("TRI_3") +
                             "_1");

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1, sidesetSubsets,
                          SideVector{{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}});
  }

  TEST_F(TestTextMeshSkin, pyramidWithStandardName_skinAndSplitByTopology)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,PYRAMID_5,1,2,3,4,5,block_1|sideset:name=surface_1; skin=block_1; split=topology";
    setup_text_mesh(meshDesc);

    verify_num_elements(1);
    verify_single_element(1u, "PYRAMID_5", EntityIdVector{1, 2, 3, 4, 5});

    std::vector<std::string> sidesetSubsets;
    std::string              prefix("surface_");
    sidesetSubsets.push_back(prefix + get_topology_name("PYRAMID_5") + "_" +
                             get_topology_name("QUAD_4") + "_1");
    sidesetSubsets.push_back(prefix + get_topology_name("PYRAMID_5") + "_" +
                             get_topology_name("TRI_3") + "_1");

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1, sidesetSubsets,
                          SideVector{{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}});
  }

  TEST_F(TestTextMeshSkin, twoHexConnected_separateBlocks_skinOneBlock)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
        "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2|sideset:name=surface_1; skin=block_1";
    setup_text_mesh(meshDesc);

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1, SideVector{{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}});
  }

  TEST_F(TestTextMeshSkin, connectedPyramidAndHexAndTet_skinAll)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,PYRAMID_5,1,2,3,4,5,block_1\n"
                           "0,2,HEX_8,1,4,3,2,6,9,8,7,block_2\n"
                           "0,3,TET_4,2,3,5,10,block_3|sideset:name=surface_1; skin=all";
    setup_text_mesh(meshDesc);

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1,
                          SideVector{{1, 1},
                                     {1, 3},
                                     {1, 4},
                                     {2, 1},
                                     {2, 2},
                                     {2, 3},
                                     {2, 4},
                                     {2, 6},
                                     {3, 1},
                                     {3, 2},
                                     {3, 3}});
  }

  TEST_F(TestTextMeshSkin, singleHexAndConnectedShell_skinAll)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "0,2,SHELL_QUAD_4,5,6,7,8,block_2|sideset:name=surface_1; skin=all";
    setup_text_mesh(meshDesc);

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1,
                          SideVector{{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 1}});
  }

  TEST_F(TestTextMeshSkin, twoHexConnectedTwoBlocks_skinIntoTwoSideSets)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2|sideset:name=surface_1; "
                           "skin=block_1|sideset:name=surface_2; skin=block_2";
    setup_text_mesh(meshDesc);

    verify_num_sidesets(2);
    verify_single_sideset("surface_1", 1,
                          SideVector{{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}});
    verify_single_sideset("surface_2", 2,
                          SideVector{{2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}, {2, 6}});
  }

  TEST_F(TestTextMeshSkin, threeHexConnectedThreeBlocks_skinTwoEndBlocks)
  {
    if (get_parallel_size() != 1)
      return;

    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
        "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2\n"
        "0,3,HEX_8,9,10,11,12,13,14,15,16,block_3|sideset:name=surface_1; skin=block_1,block_3";
    setup_text_mesh(meshDesc);

    verify_num_sidesets(1);
    verify_single_sideset("surface_1", 1,
                          SideVector{{1, 1},
                                     {1, 2},
                                     {1, 3},
                                     {1, 4},
                                     {1, 5},
                                     {1, 6},
                                     {3, 1},
                                     {3, 2},
                                     {3, 3},
                                     {3, 4},
                                     {3, 5},
                                     {3, 6}});
  }

} // namespace
