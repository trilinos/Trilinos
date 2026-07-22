#include "gtest/gtest.h"
#include <string>
#include <vector>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_balance/internal/StkBalanceUtils.hpp>

namespace {

class ElementConnectivity : public stk::unit_test_util::MeshFixture
{
protected:
  ElementConnectivity()
    : stk::unit_test_util::MeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }
};

TEST_F(ElementConnectivity, NumSharedNodes_DisconnectedHexes)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,9,10,11,12,13,14,15,16";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);

  EXPECT_EQ(stk::balance::internal::get_num_common_nodes_between_elements(get_bulk(), elem1, elem2), 0);
  EXPECT_EQ(stk::balance::internal::get_num_common_nodes_between_elements(get_bulk(), elem2, elem1), 0);
}

TEST_F(ElementConnectivity, NumSharedNodes_HexesConnectedAtOneNode)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,8,9,10,11,12,13,14,15";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);

  EXPECT_EQ(stk::balance::internal::get_num_common_nodes_between_elements(get_bulk(), elem1, elem2), 1);
  EXPECT_EQ(stk::balance::internal::get_num_common_nodes_between_elements(get_bulk(), elem2, elem1), 1);
}

TEST_F(ElementConnectivity, NumSharedNodes_HexesConnectedAtTwoNodes)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,7,8,9,10,11,12,13,14";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);

  EXPECT_EQ(stk::balance::internal::get_num_common_nodes_between_elements(get_bulk(), elem1, elem2), 2);
  EXPECT_EQ(stk::balance::internal::get_num_common_nodes_between_elements(get_bulk(), elem2, elem1), 2);
}

TEST_F(ElementConnectivity, NumSharedNodes_HexesConnectedAtThreeNodes)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,6,7,8,9,10,11,12,13";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);

  EXPECT_EQ(stk::balance::internal::get_num_common_nodes_between_elements(get_bulk(), elem1, elem2), 3);
  EXPECT_EQ(stk::balance::internal::get_num_common_nodes_between_elements(get_bulk(), elem2, elem1), 3);
}

TEST_F(ElementConnectivity, NumSharedNodes_HexesConnectedAtFourNodes)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
                         "0,2,HEX_8,5,6,7,8,9,10,11,12";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);

  EXPECT_EQ(stk::balance::internal::get_num_common_nodes_between_elements(get_bulk(), elem1, elem2), 4);
  EXPECT_EQ(stk::balance::internal::get_num_common_nodes_between_elements(get_bulk(), elem2, elem1), 4);
}

TEST_F(ElementConnectivity, NumSharedNodes_DegenerateHexAndTet)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,5,5,5\n"
                         "0,2,TET_4,4,5,6,7";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);

  EXPECT_EQ(stk::balance::internal::get_num_common_nodes_between_elements(get_bulk(), elem1, elem2), 2);
  EXPECT_EQ(stk::balance::internal::get_num_common_nodes_between_elements(get_bulk(), elem2, elem1), 2);
}

}
