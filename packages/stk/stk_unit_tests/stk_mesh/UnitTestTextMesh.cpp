#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>  // for Field
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>  // for MetaData, put_field, etc
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/TextMeshFixture.hpp>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "mpi.h"

namespace
{
class TestTextMesh : public stk::unit_test_util::TextMeshFixture
{
protected:
  TestTextMesh() : TextMeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }
};

TEST_F(TestTextMesh, singleHex)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
}

TEST_F(TestTextMesh, singleHexWithCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::mesh::EntityIdVector nodeIds = stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestTextMesh, singleHexWithCoordinates_separatedNodeIds)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,11,12,13,14";
  stk::mesh::EntityIdVector nodeIds = stk::mesh::EntityIdVector{1, 2, 3, 4, 11, 12, 13, 14};
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestTextMesh, twoHexesSerial)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{5,6,7,8,9,10,11,12});
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

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::TET_10, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8,9,10});
  verify_single_element(2u, stk::topology::TET_10, stk::mesh::EntityIdVector{2,11,3,4,12,13,6,9,14,10});
}

TEST_F(TestTextMesh, twoHexDisconnectedWithCoordinatesAndParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_1\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,block_2";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                      0,0,1, 1,0,1, 1,1,1, 0,1,1, 0,0,2, 1,0,2, 1,1,2, 0,1,2 };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{9, 10, 11, 12, 13, 14, 15, 16});
  verify_coordinates(stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}, coordinates);
  verify_part_membership({{"block_1", {1u}}, {"block_2", {2u}}});
}

TEST_F(TestTextMesh, twoHexDisconnectedWithDefaultParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{9, 10, 11, 12, 13, 14, 15, 16});
  verify_part_membership({{"block_HEXAHEDRON_8", {1u, 2u}}});
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

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(3);
  verify_single_element(1u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{1, 2, 4});
  verify_single_element(2u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 5, 4});
  verify_single_element(3u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 3, 5});
  verify_coordinates(stk::mesh::EntityIdVector{1, 2, 3, 4, 5}, coordinates);
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
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(3);
  verify_single_element(1u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{1, 2, 4});
  verify_single_element(2u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 5, 4});
  verify_single_element(3u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 3, 5});
  verify_part_membership({{"block_SHELL_TRIANGLE_3", {1u, 2u, 3u}}});
}


TEST_F(TestTextMesh, partIds_oneDefaultPartOneElem)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"block_HEXAHEDRON_8", 1u}});
}

TEST_F(TestTextMesh, partIds_oneDefaultPartTwoElems)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"block_HEXAHEDRON_8", 1u}});
}

TEST_F(TestTextMesh, partIds_twoDefaultParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16\n"
                         "0,3,SHELL_QUAD_4,17,18,19,20";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"block_HEXAHEDRON_8", 1u}, {"block_SHELL_QUADRILATERAL_4", 2u}});
}

TEST_F(TestTextMesh, partIds_oneDefaultPartOneUserSpecifiedPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,my_cool_part";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"block_HEXAHEDRON_8", 1u}, {"my_cool_part", 2u}});
}

TEST_F(TestTextMesh, partIds_samePartTwice)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,my_cool_part\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,my_cool_part";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"my_cool_part", 1u}});
}

TEST_F(TestTextMesh, partIds_orderingIsByLine)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,2,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partOne\n"
                         "0,1,HEX_8,9,10,11,12,13,14,15,16,partTwo";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"partOne", 1u}, {"partTwo", 2u}});
}

TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_onePart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_101";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"block_101", 101u}});
}

TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_twoParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_101\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,block_201";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"block_101", 101u}, {"block_201", 201u}});
}

TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_withDefaultPartFirst)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,block_101";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"block_HEXAHEDRON_8", 1u}, {"block_101", 101u}});
}

TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_withDefaultPartSecond)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_101\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"block_101", 101u}, {"block_HEXAHEDRON_8", 1u}});
}

TEST_F(TestTextMesh, partIds_respectExodusNamingConvention_withDefaultPartIdCollision)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,block_1";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"block_HEXAHEDRON_8", 2u}, {"block_1", 1u}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_onePart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,partThree,3";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"partThree", 3u}});
  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
  verify_part_membership({{"partThree", {1u}}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_twoParts)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,partFive,5";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"partThree", 3u}, {"partFive", 5u}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_twoPartsSameId)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partFour,4\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,partFour,4";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"partFour", 4u}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_samePartDifferentIds)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partFour,4\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,partFour,5";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_withDefaultPartIdCollision)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,partOne,1";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"block_HEXAHEDRON_8", 2u}, {"partOne", 1u}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_withDefaultPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"partThree", 3u}, {"block_HEXAHEDRON_8", 1u}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_withExodusPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,block_4";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"partThree", 3u}, {"block_4", 4u}});
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_collidesWithExodusPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThree,3\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,block_3";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_collidesWithPreviousSpec)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,partThreeA,3\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,partThreeB,3";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, partIds_userSpecifiedPartId_forExodusPart)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_2,3\n";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, partIds_shortPartNamesAreValid)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,a\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16,b,3";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_part_ids({{"a", 1u}, {"b", 3u}});
}

TEST_F(TestTextMesh, partIds_integerPartNamesAreInvalid)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}


TEST_F(TestTextMesh, twoHexesParallel)
{
  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "1,2,HEX_8,5,6,7,8,9,10,11,12";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
    verify_shared_nodes(stk::mesh::EntityIdVector{5,6,7,8}, 1);
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{5,6,7,8,9,10,11,12});
    verify_shared_nodes(stk::mesh::EntityIdVector{5,6,7,8}, 0);
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

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, stk::topology::SHELL_QUAD_4, stk::mesh::EntityIdVector{1, 2, 5, 4});
    verify_coordinates(stk::mesh::EntityIdVector{1, 2, 4, 5}, {0,0,0, 1,0,0, 0,1,0, 1,1,0});
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, stk::topology::SHELL_QUAD_4, stk::mesh::EntityIdVector{2, 3, 6, 5});
    verify_coordinates(stk::mesh::EntityIdVector{2, 3, 5, 6}, {1,0,0, 2,0,0, 1,1,0, 2,1,0});
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

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  if (rank == 0) {
    verify_num_elements(2);
    verify_single_element(1u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{1, 2, 4});
    verify_single_element(2u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 5, 4});
    verify_coordinates(stk::mesh::EntityIdVector{1, 2, 4, 5}, {0,0,0, 1,0,0, 0,1,0, 1,1,0});
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(3u, stk::topology::SHELL_TRI_3, stk::mesh::EntityIdVector{2, 3, 5});
    verify_coordinates(stk::mesh::EntityIdVector{2, 3, 5}, {1,0,0, 2,0,0, 1,1,0});
  }
}

TEST_F(TestTextMesh, singleHexWithSpaces)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0, 1, HEX_8, 1, 2, 3, 4, 5, 6, 7, 8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
}

TEST_F(TestTextMesh, singleHexWithCoordinatesAndSpaces)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0, 1, HEX_8, 1, 2, 3, 4, 5, 6, 7, 8";
  stk::mesh::EntityIdVector nodeIds = stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestTextMesh, singleHexWithLowerCase)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,Hex_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
}

TEST_F(TestTextMesh, singleHexWithCoordinatesAndLowerCase)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,Hex_8,1,2,3,4,5,6,7,8";
  stk::mesh::EntityIdVector nodeIds = stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::HEX_8, nodeIds);
  verify_coordinates(nodeIds, coordinates);
}

TEST_F(TestTextMesh, tooFewNodes)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooFewNodesWithCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, tooFewCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, tooManyNodes)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9,10";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooManyNodesWithCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9,10";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1, 2,0,0, 2,1,0 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, tooManyCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1, 52 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleData_empty)
{
  std::string meshDesc = "";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleData_startsWithString)
{
  std::string meshDesc = "hi";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleData_noGlobalId)
{
  std::string meshDesc = "0 ";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleData_noTopology)
{
  std::string meshDesc = "0,1,";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleData_noNodes)
{
  std::string meshDesc = "0,1,HEX_8";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, tooLittleDataWithCoordinates)
{
  std::string meshDesc = "0,1,";
  std::vector<double> coordinates;
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, invalidTopology)
{
  std::string meshDesc = "0,1,invalid,1";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, invalidTopologyWithCoordinates)
{
  std::string meshDesc = "0,1,invalid,1";
  std::vector<double> coordinates = { 0,0,0 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, mixedSpatialDim)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,QUAD_4_2D,5,6,7,8";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, mixedSpatialDimWithCoordinates)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,QUAD_4_2D,5,6,7,8";
  std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, spatialDimInconsistentWithMetaData)
{
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, spatialDimInconsistentWithMetaDataWithCoordinates)
{
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
  std::vector<double> coordinates = { 0,0, 1,0, 1,1, 0,1 };
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates), std::logic_error);
}

TEST_F(TestTextMesh, endingWithNewlineIsOk)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n";
  EXPECT_NO_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc));
}

TEST_F(TestTextMesh, stringAfterPartNameIsError)
{
  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1,bogus\n";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

TEST_F(TestTextMesh, particleHex)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,PARTICLE,1\n"
                         "0,2,HEX_8,2,3,4,5,6,7,8,9";
  EXPECT_NO_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc));

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::PARTICLE, stk::mesh::EntityIdVector{1});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{2,3,4,5,6,7,8,9});
}

TEST_F(TestTextMesh, particleHexWithCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,PARTICLE,1\n"
                         "0,2,HEX_8,2,3,4,5,6,7,8,9";
  std::vector<double> coordinates = { 2,0,0, 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 };
  EXPECT_NO_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates));

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::PARTICLE, stk::mesh::EntityIdVector{1});
  verify_single_element(2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{2,3,4,5,6,7,8,9});
  verify_coordinates(stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8,9}, coordinates);
}

class TestTextMeshAura : public stk::unit_test_util::TextMeshFixture
{
protected:
  TestTextMeshAura() : TextMeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }
};

TEST_F(TestTextMeshAura, twoQuadShellWithCoordinates)
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

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::SHELL_QUAD_4, stk::mesh::EntityIdVector{1, 2, 5, 4});
  verify_single_element(2u, stk::topology::SHELL_QUAD_4, stk::mesh::EntityIdVector{2, 3, 6, 5});
  verify_coordinates(stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6}, coordinates);

  if (rank == 0) verify_shared_nodes(stk::mesh::EntityIdVector{2,5}, 1);
  if (rank == 1) verify_shared_nodes(stk::mesh::EntityIdVector{2,5}, 0);
}

class TestTextMesh2d : public stk::unit_test_util::TextMeshFixture
{
protected:
  TestTextMesh2d() : TextMeshFixture(2)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }
};

TEST_F(TestTextMesh2d, singleQuad)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(1);
  verify_single_element(1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1,2,3,4});
}

TEST_F(TestTextMesh2d, twoSprings)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,SPRING_2,1,2\n"
                         "0,2,SPRING_2,2,3";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  verify_num_elements(2);
  verify_single_element(1u, stk::topology::SPRING_2, stk::mesh::EntityIdVector{1,2});
  verify_single_element(2u, stk::topology::SPRING_2, stk::mesh::EntityIdVector{2,3});
}

TEST_F(TestTextMesh2d, threeQuadsWithCoordinates)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4\n"
                         "0,2,QUAD_4_2D,2,3,5,6\n"
                         "0,3,QUAD_4_2D,5,7,8,6";
  std::vector<double> coordinates = { 0,0, 1,0, 1,1, 0,1, 2,0, 2,1, 3,0, 3,1 };

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  verify_num_elements(3);
  verify_single_element(1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1,2,3,4});
  verify_single_element(2u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{2,3,5,6});
  verify_single_element(3u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{5,7,8,6});
  verify_coordinates(stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8}, coordinates);
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

  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coordinates);

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1, 2, 5, 4});
    verify_shared_nodes(stk::mesh::EntityIdVector{2,5}, 1);
    verify_coordinates(stk::mesh::EntityIdVector{1, 2, 4, 5}, {0,0, 1,0, 0,1, 1,1});
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{2, 3, 6, 5});
    verify_shared_nodes(stk::mesh::EntityIdVector{2,5}, 0);
    verify_coordinates(stk::mesh::EntityIdVector{2, 3, 5, 6}, {1,0, 2,0, 1,1, 2,1});
  }
}

TEST_F(TestTextMesh2d, twoQuadOneShellParallel)
{
  if (get_parallel_size() != 3) return;
  int rank = get_parallel_rank();

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4\n"
                         "1,2,QUAD_4_2D,3,4,5,6\n"
                         "2,3,SHELL_LINE_2,3,4";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  if (rank == 0) {
    verify_num_elements(1);
    verify_single_element(1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1,2,3,4});
    verify_shared_nodes(stk::mesh::EntityIdVector{3,4}, 1);
    verify_shared_nodes(stk::mesh::EntityIdVector{3,4}, 2);
  }
  else if (rank == 1) {
    verify_num_elements(1);
    verify_single_element(2u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{3,4,5,6});
    verify_shared_nodes(stk::mesh::EntityIdVector{3,4}, 0);
    verify_shared_nodes(stk::mesh::EntityIdVector{3,4}, 2);
  }
  else if (rank == 2) {
    verify_num_elements(1);
    verify_single_element(3u, stk::topology::SHELL_LINE_2, stk::mesh::EntityIdVector{3,4});
    verify_shared_nodes(stk::mesh::EntityIdVector{3,4}, 0);
    verify_shared_nodes(stk::mesh::EntityIdVector{3,4}, 1);
  }
}

class TestTextMesh1d : public stk::unit_test_util::TextMeshFixture
{
protected:
  TestTextMesh1d() : TextMeshFixture(1)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }
};

TEST_F(TestTextMesh1d, oneDimensionNotSupported)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,LINE_2_1D,1,2";
  EXPECT_THROW(stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc), std::logic_error);
}

void test_get_mesh_spec(unsigned blockCountToDist, const std::vector<unsigned>& numProcs,
                        const std::vector<std::vector<unsigned>>& expectedDist)
{
  EXPECT_EQ(numProcs.size(), expectedDist.size());

  for(unsigned i = 0; i < numProcs.size(); i++) {
    unsigned procCount = numProcs[i];
    std::vector<unsigned> procs;
    stk::unit_test_util::get_block_proc_distribution(blockCountToDist, procCount, procs);

    EXPECT_EQ(expectedDist[i].size(), procs.size());
    for(unsigned j = 0; j < procs.size(); j++) {
      EXPECT_EQ(expectedDist[i][j], procs[j]) << "i,j: (" << i << ", " << j << ")";
    }
  }
}

TEST(GetMeshSpecTest, TestGetMeshSpecWithMultiProc)
{
  unsigned blockCountToDist = 6;
  std::vector<unsigned> numProcs = {1,2,4,6};
  std::vector<std::vector<unsigned>> expectedDist = { {0,0,0,0,0,0},
                                                      {0,0,0,1,1,1},
                                                      {0,0,1,1,2,3},
                                                      {0,1,2,3,4,5} };
  test_get_mesh_spec(blockCountToDist, numProcs, expectedDist);
}

}  // namespace
