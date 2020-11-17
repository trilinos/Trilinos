#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_unit_test_utils/getOption.h"
#include "DisconnectBlocksMeshConstruction.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include <algorithm>
#include <fstream>
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <libgen.h>
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_tools/mesh_tools/DetectHingesImpl.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocks.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocksImpl.hpp>
#include <stk_tools/mesh_tools/DisconnectUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/ConstructedMesh.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <string>

class TestDisconnectBlocks2D : public stk::unit_test_util::MeshFixture
{
protected:
  TestDisconnectBlocks2D() : stk::unit_test_util::MeshFixture(2)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }

};

TEST_F(TestDisconnectBlocks2D, disconnect_1block_1quad)
{
  stk::mesh::PartVector blocks = setup_mesh_1block_1quad(get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(4u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_1block_1quad.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_user_block_1block_1quad_empty_part)
{
  stk::mesh::Part & block2 = create_part(get_meta(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::PartVector blocks = setup_mesh_1block_1quad(get_bulk());

  stk::tools::BlockPairVector blocksToDisconnect;
  blocksToDisconnect.emplace_back(blocks[0], &block2);

  stk::tools::DisconnectBlocksOption option(stk::tools::DISCONNECT_LOCAL, stk::tools::SNIP_ALL_HINGES);
  EXPECT_NO_THROW(stk::tools::disconnect_user_blocks(get_bulk(), blocksToDisconnect, option));
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_1quad)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_1quad(get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(4u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_1quad.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_2quad)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad(get_bulk());

  output_mesh(get_bulk(), "disconnect_2block_2quad_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_2quad.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_user_blocks_2block_2quad_in_reverse_block_order)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad(get_bulk());
  stk::tools::BlockPairVector blocksToDisconnect;
  blocksToDisconnect.emplace_back(blocks[1], blocks[0]);

  stk::tools::DisconnectBlocksOption option(stk::tools::DISCONNECT_LOCAL, stk::tools::SNIP_ALL_HINGES);

  EXPECT_NO_THROW(stk::tools::disconnect_user_blocks(get_bulk(), blocksToDisconnect, option));
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_2quad_updateGraph)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad(get_bulk());

  get_bulk().initialize_face_adjacent_element_graph();
  stk::tools::disconnect_all_blocks(get_bulk());

  const stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
  if (get_bulk().is_valid(elem1)) {
    EXPECT_EQ(0u, get_bulk().get_face_adjacent_element_graph().get_num_connected_elems(elem1));
  }
  if (get_bulk().is_valid(elem2)) {
    EXPECT_EQ(0u, get_bulk().get_face_adjacent_element_graph().get_num_connected_elems(elem2));
  }

  create_all_sides(get_bulk(), get_meta().locally_owned_part(), {&get_meta().get_topology_root_part(stk::topology::LINE_2)}, false);

  unsigned localNumSides = stk::mesh::count_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::EDGE_RANK));
  unsigned numSides = stk::get_global_sum(get_bulk().parallel(), localNumSides);
  EXPECT_EQ(8u, numSides);
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_2quad_reversed)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad_reversed(get_bulk());

  output_mesh(get_bulk(), "disconnect_2block_2quad_reversed_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_2quad_reversed.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_corner_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 1);

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp1_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp1.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_corner_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 2);

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp2_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp2.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_corner_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 3);

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp3_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp3.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_corner_decomp4)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 4);

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp4_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp4.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_swappedCorner_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 1);

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp1_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp1.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_swappedCorner_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 2);

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp2_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp2.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_swappedCorner_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 3);

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp3_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp3.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_swappedCorner_decomp4)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 4);

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp4_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp4.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder1_decomp1)
{
  test_mesh_3block_4quad(get_bulk(), 1, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder1_decomp2)
{
  test_mesh_3block_4quad(get_bulk(), 1, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder1_decomp3)
{
  test_mesh_3block_4quad(get_bulk(), 1, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder2_decomp1)
{
  test_mesh_3block_4quad(get_bulk(), 2, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder2_decomp2)
{
  test_mesh_3block_4quad(get_bulk(), 2, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder2_decomp3)
{
  test_mesh_3block_4quad(get_bulk(), 2, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder3_decomp1)
{
  test_mesh_3block_4quad(get_bulk(), 3, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder3_decomp2)
{
  test_mesh_3block_4quad(get_bulk(), 3, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder3_decomp3)
{
  test_mesh_3block_4quad(get_bulk(), 3, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder4_decomp1)
{
  test_mesh_3block_4quad(get_bulk(), 4, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder4_decomp2)
{
  test_mesh_3block_4quad(get_bulk(), 4, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder4_decomp3)
{
  test_mesh_3block_4quad(get_bulk(), 4, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder5_decomp1)
{
  test_mesh_3block_4quad(get_bulk(), 5, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder5_decomp2)
{
  test_mesh_3block_4quad(get_bulk(), 5, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder5_decomp3)
{
  test_mesh_3block_4quad(get_bulk(), 5, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder6_decomp1)
{
  test_mesh_3block_4quad(get_bulk(), 6, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder6_decomp2)
{
  test_mesh_3block_4quad(get_bulk(), 6, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder6_decomp3)
{
  test_mesh_3block_4quad(get_bulk(), 6, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_keepLowerRight_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(get_bulk(), 1);

  output_mesh(get_bulk(), "disconnect_3block_4quad_keepLowerRight_decomp1_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_3block_4quad_keepLowerRight_decomp1.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_keepLowerRight_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(get_bulk(), 2);

  output_mesh(get_bulk(), "disconnect_3block_4quad_keepLowerRight_decomp2_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_3block_4quad_keepLowerRight_decomp2.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_keepLowerRight_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(get_bulk(), 3);

  output_mesh(get_bulk(), "disconnect_3block_4quad_keepLowerRight_decomp3_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_3block_4quad_keepLowerRight_decomp3.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_checkerboard_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(get_bulk(), 1);

  output_mesh(get_bulk(), "disconnect_2block_4quad_checkerboard_decomp1_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_checkerboard_decomp1.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_checkerboard_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(get_bulk(), 2);

  output_mesh(get_bulk(), "disconnect_2block_4quad_checkerboard_decomp2_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_checkerboard_decomp2.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_checkerboard_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(get_bulk(), 3);

  output_mesh(get_bulk(), "disconnect_2block_4quad_checkerboard_decomp3_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_checkerboard_decomp3.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_checkerboard_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(get_bulk(), 1);

  output_mesh(get_bulk(), "disconnect_3block_4quad_checkerboard_decomp1_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(15u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_3block_4quad_checkerboard_decomp1.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_checkerboard_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(get_bulk(), 2);

  output_mesh(get_bulk(), "disconnect_3block_4quad_checkerboard_decomp2_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(15u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_3block_4quad_checkerboard_decomp2.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_checkerboard_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(get_bulk(), 3);

  output_mesh(get_bulk(), "disconnect_3block_4quad_checkerboard_decomp3_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(15u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_3block_4quad_checkerboard_decomp3.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_2quad_diagonal)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad_diagonal(get_bulk());

  output_mesh(get_bulk(), "disconnect_2block_2quad_diagonal_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_2quad_diagonal.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_bowtie)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4quad_bowtie(get_bulk());

    output_mesh(get_bulk(), "disconnect_3block_4quad_bowtie_init.g");

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ( 0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(15u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_3block_4quad_bowtie.g");
  }
}

typedef TestDisconnectBlocks2D TestReconnectList2D;

void test_reconnect_list(const stk::mesh::BulkData& bulk, const BlockConnectionVector& disconnectList, const BlockConnectionVector& expectedReconnectList)
{
  stk::tools::impl::PartPairLess comparator;
  stk::tools::BlockPairVector reconnectList = get_local_reconnect_list(bulk, disconnectList);
  EXPECT_EQ(reconnectList.size(), expectedReconnectList.size());

  stk::tools::BlockPairVector convertedExpectedList = convert_connection_vector_to_pair_vector(bulk, expectedReconnectList);

  for(stk::tools::BlockPair& reconnectPair : reconnectList) {
    EXPECT_TRUE(std::binary_search(convertedExpectedList.begin(), convertedExpectedList.end(), reconnectPair, comparator))
        << " did not find: " << reconnectPair.first->name() << " and " << reconnectPair.second->name() << std::endl;
  }
}

TEST_F(TestReconnectList2D, reconnect_2block_2quad)
{
  stk::mesh::BulkData& bulk = get_bulk();
  setup_mesh_2block_2quad(bulk);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2") };
  BlockConnectionVector expectedReconnectList;

  test_reconnect_list(bulk, blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_3block_4quad_permutation1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  setup_mesh_3block_4quad(bulk, 1, 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_3") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_1", "block_2"), BlockConnection("block_2", "block_3") };

  test_reconnect_list(bulk, blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_3block_4quad_permutation2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  setup_mesh_3block_4quad(bulk, 1, 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_1", "block_3"), BlockConnection("block_2", "block_3") };

  test_reconnect_list(bulk, blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_3block_4quad_permutation3)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  setup_mesh_3block_4quad(bulk, 1, 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_2","block_3") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_1", "block_2"), BlockConnection("block_1", "block_3") };

  test_reconnect_list(bulk, blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_4block_4quad_permutation1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4) { return; }
  setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_1", "block_3"), BlockConnection("block_2", "block_3"), BlockConnection("block_2", "block_4"), BlockConnection("block_1", "block_4"), BlockConnection("block_3", "block_4") };

  test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_4block_4quad_permutation2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4) { return; }
  setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2"), BlockConnection("block_1","block_3") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_2", "block_3"), BlockConnection("block_2", "block_4"), BlockConnection("block_1", "block_4"), BlockConnection("block_3", "block_4") };

  test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_4block_4quad_permutation3)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4) { return; }
  setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2"), BlockConnection("block_3", "block_4") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_2", "block_3"), BlockConnection("block_2", "block_4"), BlockConnection("block_1", "block_4"), BlockConnection("block_1", "block_3") };

  test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_4block_4quad_permutation4)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4) { return; }
  setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2"), BlockConnection("block_2", "block_4"), BlockConnection("block_3", "block_4"), BlockConnection("block_1", "block_3") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_2", "block_3"),  BlockConnection("block_1", "block_4") };

  test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_9block_9quad_permutation1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_9block_9quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2"), BlockConnection("block_1", "block_4") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_1", "block_5"), BlockConnection("block_2", "block_4"), BlockConnection("block_2", "block_5"), BlockConnection("block_5", "block_4") };

  test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_9block_9quad_permutation2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_9block_9quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_5","block_2"), BlockConnection("block_5", "block_4"), BlockConnection("block_5", "block_6"), BlockConnection("block_5", "block_8") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_1", "block_5"), BlockConnection("block_1", "block_2"), BlockConnection("block_1", "block_4"), BlockConnection("block_2", "block_4"),
                                               BlockConnection("block_3", "block_5"), BlockConnection("block_2", "block_6"), BlockConnection("block_2", "block_3"), BlockConnection("block_3", "block_6"),
                                               BlockConnection("block_5", "block_9"), BlockConnection("block_6", "block_8"), BlockConnection("block_8", "block_9"), BlockConnection("block_6", "block_9"),
                                               BlockConnection("block_5", "block_7"), BlockConnection("block_4", "block_8"), BlockConnection("block_4", "block_7"), BlockConnection("block_7", "block_8"),
                                             };

  test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
}


class TestDisconnectBlocks : public stk::unit_test_util::MeshFixture
{
protected:
  TestDisconnectBlocks() : stk::unit_test_util::MeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }
};

typedef TestDisconnectBlocks TestReconnectList;

TEST_F(TestReconnectList, reconnect_8block_8hex_permutation1)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnect { BlockConnection("block_1","block_2"), BlockConnection("block_2","block_4"), BlockConnection("block_2","block_6"),
                                                   BlockConnection("block_5","block_7"), BlockConnection("block_5","block_6"), BlockConnection("block_5","block_1"),
                                                   BlockConnection("block_8","block_6"), BlockConnection("block_8","block_4"), BlockConnection("block_8","block_7"),
                                                   BlockConnection("block_3","block_7"), BlockConnection("block_3","block_4"), BlockConnection("block_3","block_1")};
    BlockConnectionVector expectedReconnectList { BlockConnection("block_2","block_3"), BlockConnection("block_2","block_5"), BlockConnection("block_2","block_8"), BlockConnection("block_2","block_7"),
                                                  BlockConnection("block_3","block_5"), BlockConnection("block_3","block_8"), BlockConnection("block_5","block_8"), BlockConnection("block_1","block_8"),
                                                  BlockConnection("block_1","block_4"), BlockConnection("block_1","block_6"), BlockConnection("block_1","block_7"), BlockConnection("block_3","block_6"),
                                                  BlockConnection("block_4","block_6"), BlockConnection("block_4","block_7"), BlockConnection("block_6","block_7"), BlockConnection("block_4","block_5")
                                                };
    test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
  }
}

TEST_F(TestDisconnectBlocks, disconnect_1block_1hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_1block_1hex(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_1block_1hex.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_1hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_2block_1hex(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_1hex.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_2hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2hex(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_2hex.g");
  }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_2block_2hex_with_internal_sides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2hex_with_internal_sides(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_2hex_with_internal_sides.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_2hex_with_external_sides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2hex_with_external_sides(get_bulk());

    EXPECT_EQ(12u, get_num_surface_nodes(get_bulk(), {"surface_1"}));

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(16u, get_num_surface_nodes(get_bulk(), {"surface_1"}));
    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u,get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_2hex_with_external_sides.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_8block_8hex_with_external_sides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_with_external_sides(get_bulk());

    EXPECT_EQ(26u, get_num_surface_nodes(get_bulk(), {"surface_1"}));

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(56u, get_num_surface_nodes(get_bulk(), {"surface_1"}));
    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(64u,get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_8block_8hex_with_external_sides.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_2cubeOfTet)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2cubeOfTet(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_2cubeOfTet.g");
  }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_2block_2cubeOfTet_with_internal_sides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2cubeOfTet_with_internal_sides(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_2cubeOfTet_with_internal_sides.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_3block_4hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex(get_bulk());

    output_mesh(get_bulk(), "disconnect_3block_4hex_init.g");

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_3block_4hex.g");
  }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_3block_4hex_with_internal_sides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex_with_internal_sides(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_3block_4hex_with_internal_sides.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_3block_4cubeOfTet)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4cubeOfTet(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_3block_4cubeOfTet.g");
  }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_3block_4cubeOfTet_with_internal_sides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4cubeOfTet_with_internal_sides(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_3block_4cubeOfTet_with_internal_sides.g");
  }
}

TEST(DisconnectBlocks, input_mesh)
{
  double startTime = stk::wall_time();

  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

  std::string exodusFileName = stk::unit_test_util::get_option("-f", "");
  if (exodusFileName.empty()) return;

  stk::io::fill_mesh_with_auto_decomp(exodusFileName, bulk);

  double meshReadTime = stk::wall_time();

  stk::tools::disconnect_all_blocks(bulk);

  double disconnectTime = stk::wall_time();

  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  EXPECT_EQ(0u,  get_num_intersecting_nodes(bulk, allBlocksInMesh));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  output_mesh(bulk, "disconnected_" + exodusFileName);
  double meshWriteTime = stk::wall_time();

  if (bulk.parallel_rank() == 0) {
    std::cout << " Mesh read time = " << (meshReadTime - startTime) << " s" << std::endl;
    std::cout << "Disconnect time = " << (disconnectTime - meshReadTime) << " s" << std::endl;
    std::cout << "Mesh write time = " << (meshWriteTime - disconnectTime) << " s" << std::endl;
  }
}

typedef TestDisconnectBlocks2D TestReconnectBlocks2D;
typedef TestDisconnectBlocks   TestReconnectBlocks;
typedef TestDisconnectBlocks2D TestDisconnectFullAlgorithm2D;
typedef TestDisconnectBlocks2D TestDisconnectFullAlgorithmPartial2D;
typedef TestDisconnectBlocks   TestDisconnectFullAlgorithm;
typedef TestDisconnectBlocks   TestDisconnectFullAlgorithmPartial;

void test_one_block_reconnect(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  ThrowRequire(allBlocksInMesh.size() == 1u);
  stk::tools::impl::LinkInfo info;
  info.preserveOrphans = true;

  std::vector<stk::tools::BlockPair> blockPairs;

  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::disconnect_block_pairs(bulk, blockPairs, info);
  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::reconnect_block_pairs(bulk, blockPairs, info);
  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));

  output_mesh(bulk);
}

void test_two_block_reconnect(stk::mesh::BulkData& bulk, const stk::mesh::PartVector& blocks,
                              const unsigned expectedCommonNodeCount)
{
  ThrowRequire(blocks.size() >= 2u);
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::impl::LinkInfo info;
  info.preserveOrphans = true;

  std::vector<stk::tools::BlockPair> blockPairs;

  stk::tools::impl::insert_block_pair(blocks[0], blocks[1], blockPairs);

  EXPECT_EQ(expectedCommonNodeCount, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::disconnect_block_pairs(bulk, blockPairs, info);
  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::reconnect_block_pairs(bulk, blockPairs, info);
  EXPECT_EQ(expectedCommonNodeCount, get_num_intersecting_nodes(bulk, allBlocksInMesh));

  output_mesh(bulk);
}

TEST_F(TestReconnectBlocks2D, reconnect_1block_1quad)
{
  stk::mesh::PartVector blocks = setup_mesh_1block_1quad(get_bulk());
  test_one_block_reconnect(get_bulk());
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_2quad)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad(get_bulk());
  test_two_block_reconnect(get_bulk(), blocks, 2u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_2quad_diagonal)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad_diagonal(get_bulk());
  test_two_block_reconnect(get_bulk(), blocks, 1u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_2quad_reversed)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad_reversed(get_bulk());
  test_two_block_reconnect(get_bulk(), blocks, 2u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_corner_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 1);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_corner_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 2);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_corner_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 3);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_corner_decomp4)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 4);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_swappedCorner_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 1);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_swappedCorner_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 2);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_swappedCorner_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 3);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_swappedCorner_decomp4)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 4);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_checkerboard_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(get_bulk(), 1);
  test_two_block_reconnect(get_bulk(), blocks, 5u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_checkerboard_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(get_bulk(), 2);
  test_two_block_reconnect(get_bulk(), blocks, 5u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_checkerboard_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(get_bulk(), 3);
  test_two_block_reconnect(get_bulk(), blocks, 5u);
}

void test_reconnect_block_pairs(stk::mesh::BulkData& bulk, stk::mesh::PartVector& blocks,
                                const std::vector<BlockConnection>& reconnectPairs,
                                unsigned expectedGlobalInitialCommonNodes,
                                unsigned expectedGlobalConnectCommonNodes)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::impl::LinkInfo info;
  info.preserveOrphans = false;

  std::vector<stk::tools::BlockPair> disconnectBlockPairs;
  for(unsigned i = 0; i < allBlocksInMesh.size() - 1; i++) {
    for(unsigned j = i+1; j < allBlocksInMesh.size(); j++) {
      stk::tools::impl::insert_block_pair(allBlocksInMesh[i], allBlocksInMesh[j], disconnectBlockPairs);
    }
  }

  std::vector<stk::tools::BlockPair> reconnectBlockPairs;

  for(const BlockConnection& connectPair : reconnectPairs) {
    ThrowRequire(connectPair.block1 != connectPair.block2);

    stk::mesh::Part* block1 = meta.get_part(connectPair.block1);
    stk::mesh::Part* block2 = meta.get_part(connectPair.block2);

    ThrowRequire(block1 != nullptr && block2 != nullptr);

    stk::tools::impl::insert_block_pair(block1, block2, reconnectBlockPairs);
  }

  EXPECT_EQ(expectedGlobalInitialCommonNodes, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::disconnect_block_pairs(bulk, disconnectBlockPairs, info);
  bulk.dump_mesh_per_proc("dump");
  output_mesh(bulk, "test.g");
  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::reconnect_block_pairs(bulk, reconnectBlockPairs, info);
  EXPECT_EQ(expectedGlobalConnectCommonNodes, get_num_intersecting_nodes(bulk, allBlocksInMesh));

  stk::mesh::PartVector connectedBlockPair(2);
  for(const BlockConnection& connectPair : reconnectPairs) {
    stk::mesh::Part* block1 = meta.get_part(connectPair.block1);
    stk::mesh::Part* block2 = meta.get_part(connectPair.block2);

    connectedBlockPair[0] = block1;
    connectedBlockPair[1] = block2;

    EXPECT_EQ(connectPair.numExpectedIntersectingNodes, get_num_intersecting_nodes(bulk, connectedBlockPair)) << " first block: " << block1->name() << " second block: " << block2->name();
  }
  output_mesh(bulk);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder1_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder1_decomp1_v2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2), BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder1_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder1_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 2, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",2), BlockConnection("block_2","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 2, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",2), BlockConnection("block_2","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp2_v2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 2, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 2);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp2_v3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 2, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 2);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 2, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",2), BlockConnection("block_2","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder3_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 3, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_1",2), BlockConnection("block_3","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder3_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 3, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_1",2), BlockConnection("block_3","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder3_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 3, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_1",2), BlockConnection("block_3","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder4_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 4, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_2",2), BlockConnection("block_3","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder4_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 4, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_2",2), BlockConnection("block_3","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder4_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 4, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_2",2), BlockConnection("block_3","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder5_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 5, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2), BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder5_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 5, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2), BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder5_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 5, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2), BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder6_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 6, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_3",2), BlockConnection("block_1","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder6_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 6, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_3",2), BlockConnection("block_1","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder6_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 6, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_3",2), BlockConnection("block_1","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_keepLowerRight_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_keepLowerRight_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(get_bulk(), 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_keepLowerRight_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(get_bulk(), 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_checkerboard_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",3), BlockConnection("block_1","block_3",3)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 5);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_checkerboard_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(get_bulk(), 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",3), BlockConnection("block_1","block_3",3)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 5);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_checkerboard_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(get_bulk(), 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",3), BlockConnection("block_1","block_3",3)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 5);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_bowtie)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_bowtie(get_bulk());
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",1), BlockConnection("block_1","block_3",1)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 1, 1);
}

TEST_F(TestReconnectBlocks2D, reconnect_4block_4quad_blockorder1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_4",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 2);
}

TEST_F(TestReconnectBlocks2D, reconnect_4block_4quad_blockorder1_v2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_4",2), BlockConnection("block_2","block_3",1)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 2);
}

TEST_F(TestReconnectBlocks2D, reconnect_4block_4quad_blockorder1_v3)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_4",2), BlockConnection("block_3","block_4",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_4block_4quad_blockorder1_v4)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_2","block_3",1), BlockConnection("block_3","block_4",2)/*, BlockConnection("block_1","block_4",1)*/};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_4block_4quad_blockorder1_v5)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_1",2), BlockConnection("block_1","block_2",2), BlockConnection("block_2","block_4",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 4);
}

TEST_F(TestReconnectBlocks2D, reconnect_4block_4quad_blockorder2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 2);
  output_mesh(get_bulk(), "initial.g");
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",1)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 2);
}

TEST_F(TestReconnectBlocks, reconnect_1block_1hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_1block_1hex(get_bulk());
    test_one_block_reconnect(get_bulk());
  }
}

TEST_F(TestReconnectBlocks, reconnect_2block_1hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_2block_1hex(get_bulk());
    test_two_block_reconnect(get_bulk(), blocks, 0u);
  }
}

TEST_F(TestReconnectBlocks, reconnect_2block_2hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2hex(get_bulk());
    test_two_block_reconnect(get_bulk(), blocks, 4u);
  }
}

TEST_F(TestReconnectBlocks, reconnect_2block_2cubeOfTet)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2cubeOfTet(get_bulk());
    test_two_block_reconnect(get_bulk(), blocks, 4u);
  }
}

TEST_F(TestReconnectBlocks, reconnect_3block_4hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex(get_bulk());
    BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",4), BlockConnection("block_2","block_1",4)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 8, 6);
  }
}

TEST_F(TestReconnectBlocks, reconnect_3block_4hex_v2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex(get_bulk());
    BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",4)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 8, 4);
  }
}

TEST_F(TestReconnectBlocks, reconnect_3block_4hex_v3)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex(get_bulk());
    BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",4)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 8, 4);
  }
}

TEST_F(TestReconnectBlocks, reconnect_3block_4cubeOfTet)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4cubeOfTet(get_bulk());
    BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",4), BlockConnection("block_2","block_1",4)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 8, 6);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_4hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4hex(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_2","block_1",4), BlockConnection("block_3","block_2",2) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 10, 4);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_4hex_v2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4hex(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_2","block_3",2), BlockConnection("block_1","block_4",2) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 10, 4);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_4hex_v3)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4hex(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_1","block_2",4), BlockConnection("block_2","block_4",4), BlockConnection("block_4","block_3",4),
                                          BlockConnection("block_2","block_3",2), BlockConnection("block_1","block_4",2)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 10, 8);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_4hex_vertical_stack)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 4) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4hex_vertical_stack(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_1","block_2",4), BlockConnection("block_3","block_4",4) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 12, 8);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_4hex_vertical_stack_v2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 4) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4hex_vertical_stack(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_2","block_3",4) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 12, 4);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_8hex_cube)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_8hex_cube(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_2","block_4",6) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 15, 6);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_8hex_cube_v2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_8hex_cube(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_2","block_3",3) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 15, 3);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_8hex_cube_v3)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_8hex_cube(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_2","block_3",3), BlockConnection("block_3","block_4",6) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 15, 6);
  }
}

TEST_F(TestReconnectBlocks, DISABLED_reconnect_8block_8hex_cube)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_1","block_3",4),
                                          BlockConnection("block_3","block_4",4),
                                          BlockConnection("block_4","block_8",4),
                                          BlockConnection("block_8","block_7",4),
                                          BlockConnection("block_7","block_5",4),
                                          BlockConnection("block_5","block_6",4)
                                        };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 19, 16);
  }
}

TEST_F(TestReconnectBlocks, reconnect_2block_2hex_with_external_sides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2hex_with_external_sides(get_bulk());

    EXPECT_EQ(12u, get_num_surface_nodes(get_bulk(), {"surface_1"}));

    stk::tools::impl::LinkInfo info;
    info.preserveOrphans = true;

    stk::tools::disconnect_all_blocks(get_bulk(), info);

    EXPECT_EQ(16u, get_num_surface_nodes(get_bulk(), {"surface_1"}));
    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u,get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_2hex_with_external_sides.g");

    stk::tools::impl::reconnect_block_pairs(get_bulk(), {stk::tools::BlockPair(blocks[0], blocks[1])}, info);

    EXPECT_EQ(12u, get_num_surface_nodes(get_bulk(), {"surface_1"}));
    EXPECT_EQ(4u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(12u,get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk());
  }
}

TEST_F(TestReconnectBlocks, reconnect_8block_8hex_with_external_sides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_with_external_sides(get_bulk());

    EXPECT_EQ(26u, get_num_surface_nodes(get_bulk(), {"surface_1"}));

    stk::tools::impl::LinkInfo info;
    info.preserveOrphans = true;

    stk::tools::disconnect_all_blocks(get_bulk(), info);

    EXPECT_EQ(56u, get_num_surface_nodes(get_bulk(), {"surface_1"}));
    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(64u,get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_8block_8hex_with_external_sides.g");

    stk::tools::BlockPairVector blockPairs { stk::tools::impl::get_block_pair(blocks[0], blocks[1]),
                                             stk::tools::impl::get_block_pair(blocks[0], blocks[2]),
                                             stk::tools::impl::get_block_pair(blocks[0], blocks[3]),
                                             stk::tools::impl::get_block_pair(blocks[1], blocks[2]),
                                             stk::tools::impl::get_block_pair(blocks[1], blocks[3]),
                                             stk::tools::impl::get_block_pair(blocks[2], blocks[3])
                                           };

    stk::tools::impl::reconnect_block_pairs(get_bulk(), blockPairs, info);

    EXPECT_EQ(45u, get_num_surface_nodes(get_bulk(), {"surface_1"}));
    EXPECT_EQ(10u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(50u,get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk());
  }
}

stk::tools::BlockPairVector convert_connection_vector_to_pair_vector(const stk::mesh::BulkData& bulk, BlockConnectionVector& disconnectConnVector)
{
  stk::tools::BlockPairVector pairVector;
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  for(BlockConnection& connection : disconnectConnVector) {
    stk::mesh::Part* part1 = meta.get_part(connection.block1);
    stk::mesh::Part* part2 = meta.get_part(connection.block2);

    pairVector.push_back(stk::tools::impl::get_block_pair(part1, part2));
  }

  return pairVector;
}

void test_populate_blocks_to_reconnect(stk::mesh::BulkData& bulk, const stk::tools::BlockPairVector& expectedReconnectBlockPairs, BlockConnectionVector& disconnectPairs)
{
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::BlockPairVector orderedBlockPairsInMesh;
  stk::tools::BlockPairVector blockPairsToReconnect;
  stk::tools::BlockPairVector blockPairsToDisconnect = convert_connection_vector_to_pair_vector(bulk, disconnectPairs);

  stk::tools::impl::fill_ordered_block_pairs(allBlocksInMesh, orderedBlockPairsInMesh);

  stk::tools::impl::populate_blocks_to_reconnect(bulk, orderedBlockPairsInMesh, blockPairsToDisconnect, blockPairsToReconnect);
  EXPECT_EQ(expectedReconnectBlockPairs.size(), blockPairsToReconnect.size());

  for(const stk::tools::BlockPair& blockPair : expectedReconnectBlockPairs) {
    EXPECT_TRUE(std::binary_search(blockPairsToReconnect.begin(), blockPairsToReconnect.end(), blockPair, stk::tools::impl::PartPairLess()))
        << "first block name: " << blockPair.first->name() << " second block name: " << blockPair.second->name() << std::endl;
  }
}

TEST_F(TestDisconnectBlocks2D, test_populate_blocks_to_reconnect_2block_2quad)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad(bulk);
  BlockConnectionVector disconnectPairs{BlockConnection("block_1","block_2",2)};
  stk::tools::BlockPairVector expectedReconnectBlockPairs;

  test_populate_blocks_to_reconnect(bulk, expectedReconnectBlockPairs, disconnectPairs);
}

TEST_F(TestDisconnectBlocks2D, test_populate_blocks_to_reconnect_3block_4quad)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(bulk, 1, 1);

  stk::mesh::Part* block1 = meta.get_part("block_1");
  stk::mesh::Part* block2 = meta.get_part("block_2");
  stk::mesh::Part* block3 = meta.get_part("block_3");

  stk::tools::BlockPairVector expectedReconnectBlockPairs;
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block1, block3));
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block2, block3));

  BlockConnectionVector disconnectPairs{BlockConnection("block_1","block_2",2)};

  test_populate_blocks_to_reconnect(bulk, expectedReconnectBlockPairs, disconnectPairs);
}

TEST_F(TestDisconnectBlocks2D, test_populate_blocks_to_reconnect_3block_4quad_v2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(bulk, 1, 1);
  BlockConnectionVector disconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1", "block_3", 2)};

  stk::mesh::Part* block2 = meta.get_part("block_2");
  stk::mesh::Part* block3 = meta.get_part("block_3");

  stk::tools::BlockPairVector expectedReconnectBlockPairs;
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block2, block3));

  test_populate_blocks_to_reconnect(bulk, expectedReconnectBlockPairs, disconnectPairs);
}

TEST_F(TestDisconnectBlocks2D, test_populate_blocks_to_reconnect_4block_4quad)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(bulk, 1);

  stk::mesh::Part* block1 = meta.get_part("block_1");
  stk::mesh::Part* block2 = meta.get_part("block_2");
  stk::mesh::Part* block3 = meta.get_part("block_3");
  stk::mesh::Part* block4 = meta.get_part("block_4");

  stk::tools::BlockPairVector expectedReconnectBlockPairs;
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block1, block3));
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block1, block4));
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block2, block3));
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block3, block4));

  BlockConnectionVector disconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_4","block_2",2)};

  test_populate_blocks_to_reconnect(bulk, expectedReconnectBlockPairs, disconnectPairs);
}

void print_block_pair_list(const stk::tools::BlockPairVector& blockPairs)
{
  return;
  std::cout << "Block pairs to reconnect" << std::endl;
  for(const stk::tools::BlockPair& blockPair : blockPairs) {
    std::cout << "\t block pair: {" << blockPair.first->name() << "," << blockPair.second->name() <<"}" << std::endl;
  }
  std::cout << std::endl;
}

void test_block_disconnect_full_algorithm_partial_disconnect(stk::mesh::BulkData& bulk,
                                                             BlockConnectionVector& blockPairConnectionsToDisconnect,
                                                             unsigned expectedGlobalInitialCommonNodes,
                                                             unsigned expectedGlobalConnectCommonNodesAfterHingeSnip,
                                                             unsigned expectedGlobalTotalNumberOfNodes)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::impl::LinkInfo info;
  info.preserveOrphans = true;

  stk::tools::BlockPairVector blockPairsToDisconnect = convert_connection_vector_to_pair_vector(bulk, blockPairConnectionsToDisconnect);
  stk::tools::BlockPairVector blockPairsToReconnect;
  blockPairsToReconnect = get_local_reconnect_list(bulk, blockPairConnectionsToDisconnect);

  EXPECT_EQ(expectedGlobalInitialCommonNodes, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::disconnect_block_pairs(bulk, blockPairsToDisconnect, info);
  output_mesh(bulk, "dis_seq_disconnected_mesh.g");

  stk::mesh::PartVector connectedBlockPair(2);
  for(const stk::tools::BlockPair& connectPair :  blockPairsToDisconnect) {
    connectedBlockPair[0] = connectPair.first;
    connectedBlockPair[1] = connectPair.second;

    EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, connectedBlockPair)) << " block1: " << connectedBlockPair[0]->name() << " block2: " << connectedBlockPair[1]->name();
  }

  stk::tools::impl::reconnect_block_pairs(bulk, blockPairsToReconnect, info);

  output_mesh(bulk, "dis_seq_reconnected_mesh.g");

  stk::tools::impl::snip_all_hinges_between_blocks(bulk);

  EXPECT_EQ(expectedGlobalConnectCommonNodesAfterHingeSnip, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  EXPECT_EQ(expectedGlobalTotalNumberOfNodes, get_num_total_nodes(bulk));

  for(const BlockConnection& connectPair :  blockPairConnectionsToDisconnect) {
    stk::mesh::Part* block1 = meta.get_part(connectPair.block1);
    stk::mesh::Part* block2 = meta.get_part(connectPair.block2);

    connectedBlockPair[0] = block1;
    connectedBlockPair[1] = block2;

    EXPECT_EQ(connectPair.numExpectedIntersectingNodes, get_num_intersecting_nodes(bulk, connectedBlockPair)) << " block1: " << block1->name() << " block2: " << block2->name();
  }
  output_mesh(bulk, "dis_seq_hinge_snipped_mesh.g");
}

void test_block_disconnect_full_algorithm(stk::mesh::BulkData& bulk,
                                          BlockConnectionVector& blockPairConnectionsToDisconnect,
                                          unsigned expectedGlobalInitialCommonNodes,
                                          unsigned expectedGlobalConnectCommonNodesAfterHingeSnip,
                                          unsigned expectedGlobalTotalNumberOfNodes)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::impl::LinkInfo info;
  info.preserveOrphans = true;

  stk::tools::BlockPairVector orderedBlockPairsInMesh;
  stk::tools::impl::fill_ordered_block_pairs(allBlocksInMesh, orderedBlockPairsInMesh);

  stk::tools::BlockPairVector blockPairsToDisconnect = convert_connection_vector_to_pair_vector(bulk, blockPairConnectionsToDisconnect);
  stk::tools::BlockPairVector blockPairsToReconnect;
  stk::tools::impl::populate_blocks_to_reconnect(bulk, orderedBlockPairsInMesh, blockPairsToDisconnect, blockPairsToReconnect);

  EXPECT_EQ(expectedGlobalInitialCommonNodes, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::disconnect_block_pairs(bulk, orderedBlockPairsInMesh, info);
  output_mesh(bulk, "dis_seq_disconnected_mesh.g");
  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::reconnect_block_pairs(bulk, blockPairsToReconnect, info);
//  EXPECT_EQ(expectedGlobalConnectCommonNodesAfterReconnect, get_num_intersecting_nodes(bulk, allBlocksInMesh));

  output_mesh(bulk, "dis_seq_reconnected_mesh.g");

  stk::tools::impl::snip_all_hinges_between_blocks(bulk);

  EXPECT_EQ(expectedGlobalConnectCommonNodesAfterHingeSnip, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  EXPECT_EQ(expectedGlobalTotalNumberOfNodes, get_num_total_nodes(bulk));

  stk::mesh::PartVector connectedBlockPair(2);
  for(const BlockConnection& connectPair :  blockPairConnectionsToDisconnect) {
    stk::mesh::Part* block1 = meta.get_part(connectPair.block1);
    stk::mesh::Part* block2 = meta.get_part(connectPair.block2);

    connectedBlockPair[0] = block1;
    connectedBlockPair[1] = block2;

    EXPECT_EQ(connectPair.numExpectedIntersectingNodes, get_num_intersecting_nodes(bulk, connectedBlockPair)) << " block1: " << block1->name() << " block2: " << block2->name();
  }
  output_mesh(bulk, "dis_seq_hinge_snipped_mesh.g");
}

void test_connection_pairs(stk::mesh::BulkData& bulk, stk::mesh::PartVector& allBlocksInMesh,
                           BlockConnectionVector& blockPairConnectionsToDisconnect,
                           unsigned expectedFinalCommonNodeCount)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector connectedBlockPair(2);
  for(const BlockConnection& connectPair :  blockPairConnectionsToDisconnect) {
    stk::mesh::Part* block1 = meta.get_part(connectPair.block1);
    stk::mesh::Part* block2 = meta.get_part(connectPair.block2);

    connectedBlockPair[0] = block1;
    connectedBlockPair[1] = block2;

    EXPECT_EQ(connectPair.numExpectedIntersectingNodes, get_num_intersecting_nodes(bulk, connectedBlockPair)) << " block1: " << block1->name() << " block2: " << block2->name();
  }

  output_mesh(bulk);

  EXPECT_EQ(expectedFinalCommonNodeCount, get_num_intersecting_nodes(bulk, allBlocksInMesh));
}

void test_user_block_disconnect(stk::mesh::BulkData& bulk,
                               BlockConnectionVector& blockPairConnectionsToDisconnect,
                               unsigned expectedFinalCommonNodeCount,
                               stk::tools::DisconnectBlocksOption disconnectOptions = stk::tools::DisconnectBlocksOption())
{
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::BlockPairVector blockPairsToDisconnect = convert_connection_vector_to_pair_vector(bulk, blockPairConnectionsToDisconnect);

  stk::tools::disconnect_user_blocks(bulk, blockPairsToDisconnect, disconnectOptions);

  test_connection_pairs(bulk, allBlocksInMesh, blockPairConnectionsToDisconnect, expectedFinalCommonNodeCount);
}

stk::tools::BlockNamePairVector convert_connection_vector_to_name_pair_vector(const stk::mesh::BulkData& bulk, const BlockConnectionVector& blockNamePairConnections)
{
  stk::tools::BlockNamePairVector namePairVector;

  for(const BlockConnection& connection : blockNamePairConnections) {
    namePairVector.push_back(std::make_pair(connection.block1, connection.block2));
  }

  return namePairVector;
}

TEST_F(TestDisconnectFullAlgorithm2D, disconnect_full_algorithm_2block_2quad)
{
  stk::mesh::PartVector initialMesh = setup_mesh_2block_2quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{BlockConnection("block_1","block_2",0)};
  test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnect, 2u, 0u, 8u);
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_2block_2quad)
{
  stk::mesh::PartVector initialMesh = setup_mesh_2block_2quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{BlockConnection("block_1","block_2",0)};
  test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnect, 2u, 0u, 8u);
}

TEST_F(TestDisconnectFullAlgorithm2D, disconnect_full_algorithm_3block_4quad_blockOrder1_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 1);
  BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_1","block_3",0)};
  test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 4u, 2u, 12u);
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_3block_4quad_blockOrder1_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 1);
  BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_1","block_3",0)};
  test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 4u, 2u, 12u);
}

TEST_F(TestDisconnectFullAlgorithm2D, disconnect_full_algorithm_3block_4quad_blockOrder1_decomp1_pacman)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 1);
  BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_2","block_3",1)};
  test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 4u, 3u, 10u);
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_3block_4quad_blockOrder1_decomp1_pacman)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 1);
  BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_2","block_3",1)};
  test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 4u, 3u, 10u);
}

TEST_F(TestDisconnectFullAlgorithm2D, disconnect_full_algorithm_4block_4quad_blockOrder1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) <= 3) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_2","block_4",0), BlockConnection("block_3","block_4",0)};
    test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 5u, 3u, 12u);
  }
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_4block_4quad_blockOrder1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) <= 3) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_2","block_4",0), BlockConnection("block_3","block_4",0)};
    test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 5u, 3u, 12u);
  }
}

TEST_F(TestDisconnectFullAlgorithm2D, disconnect_full_algorithm_6block_6quad)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_6block_6quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_2","block_5",1), BlockConnection("block_3","block_6",0),
                                                BlockConnection("block_2","block_6",0), BlockConnection("block_3","block_5",0) };

  test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnect, 8u, 8u, 14u);
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_6block_6quad)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_6block_6quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_2","block_5",1), BlockConnection("block_3","block_6",0),
                                                BlockConnection("block_2","block_6",0), BlockConnection("block_3","block_5",0) };

  test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnect, 8u, 8u, 14u);
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_9block_9quad)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_9block_9quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2"), BlockConnection("block_1", "block_4") };

  test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnect, 12u, 10u, 19u);
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_9block_9quad_v2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_9block_9quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_2","block_3"), BlockConnection("block_5","block_6"),
                                                BlockConnection("block_3","block_5"), BlockConnection("block_2","block_6"), BlockConnection("block_9", "block_6")};

  test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnect, 12u, 11u, 20u);
}

TEST_F(TestDisconnectFullAlgorithm2D, disconnect_full_algorithm_9block_9quad_v2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_9block_9quad(get_bulk());
//  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_2","block_3"), BlockConnection("block_5", "block_6"), BlockConnection("block_9", "block_6")};
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_2","block_3"), BlockConnection("block_5","block_6"),
                                                BlockConnection("block_3","block_5"), BlockConnection("block_2","block_6"), BlockConnection("block_9", "block_6")};

  test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnect, 12u, 11u, 20u);
}

TEST_F(TestDisconnectFullAlgorithm, disconnect_full_algorithm_4block_8hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_2","block_4",0), BlockConnection("block_1","block_2",0)};
    test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 15u, 9u, 36u);
  }
}

TEST_F(TestDisconnectFullAlgorithmPartial, disconnect_full_algorithm_4block_8hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_2","block_4",0), BlockConnection("block_1","block_2",0)};
    test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 15u, 9u, 36u);
  }
}

TEST_F(TestDisconnectFullAlgorithm, disconnect_full_algorithm_8block_8hex_lower_corner)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0)};
    test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 19u, 16u, 34u);
  }
}

TEST_F(TestDisconnectFullAlgorithmPartial, disconnect_full_algorithm_8block_8hex_lower_corner)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0)};
    test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 19u, 16u, 34u);
  }
}

TEST_F(TestDisconnectFullAlgorithm, disconnect_full_algorithm_8block_8hex_opposite_corners)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0),
                                                      BlockConnection("block_3","block_7",0), BlockConnection("block_7","block_8",0), BlockConnection("block_5","block_7",0)};
    test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 19u, 13u, 41u);
  }
}

TEST_F(TestDisconnectFullAlgorithmPartial, disconnect_full_algorithm_8block_8hex_opposite_corners)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0),
                                                      BlockConnection("block_3","block_7",0), BlockConnection("block_7","block_8",0), BlockConnection("block_5","block_7",0)};
    test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 19u, 13u, 41u);
  }
}

TEST_F(TestDisconnectFullAlgorithm, disconnect_full_algorithm_8block_8hex_opposite_corners2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0),
                                                      BlockConnection("block_5","block_7",0), BlockConnection("block_5","block_6",0), BlockConnection("block_5","block_1",0),
                                                      BlockConnection("block_8","block_6",0), BlockConnection("block_8","block_4",0), BlockConnection("block_8","block_7",0),
                                                      BlockConnection("block_3","block_7",0), BlockConnection("block_3","block_4",0), BlockConnection("block_3","block_1",0)};
    test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 19u, 0u, 64u);
  }
}

TEST_F(TestDisconnectFullAlgorithmPartial, disconnect_full_algorithm_8block_8hex_opposite_corners2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0),
                                                      BlockConnection("block_5","block_7",0), BlockConnection("block_5","block_6",0), BlockConnection("block_5","block_1",0),
                                                      BlockConnection("block_8","block_6",0), BlockConnection("block_8","block_4",0), BlockConnection("block_8","block_7",0),
                                                      BlockConnection("block_3","block_7",0), BlockConnection("block_3","block_4",0), BlockConnection("block_3","block_1",0)};
    test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 19u, 0u, 64u);
  }
}

std::string get_basename(const std::string& path)
{
  return path.substr(path.find_last_of("/\\")+1);
}

TEST(TestDisconnectInputFile, input_mesh)
{
  double startTime = stk::wall_time();

  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);

  std::string exodusFileName = stk::unit_test_util::get_option("-exoFile", "");
  if (exodusFileName.empty()) return;

  std::string disconnectBlockFile = stk::unit_test_util::get_option("-blockFile", "");
  if (disconnectBlockFile.empty()) return;

  stk::io::fill_mesh(exodusFileName, bulk);

  std::ifstream infile(disconnectBlockFile);

  std::string block1;
  std::string block2;
  stk::tools::BlockPairVector disconnectBlockVec;

  while(infile >> block1 >> block2) {
    stk::mesh::Part* block1Part = meta.get_part(block1);
    stk::mesh::Part* block2Part = meta.get_part(block2);

    if(block1Part == nullptr || block2Part == nullptr) {
      if(bulk.parallel_rank() == 0) {
        if(block1Part == nullptr) {
          std::cout << "could not find block name: " << block1 << std::endl;
        }
        if(block2Part == nullptr) {
          std::cout << "could not find block name: " << block2 << std::endl;
        }
      }

      continue;
    }

    disconnectBlockVec.push_back(stk::tools::impl::get_block_pair(block1Part, block2Part));
  }
  infile.close();

  double meshReadTime = stk::wall_time();

  std::cout << "Starting disconnect block sequence" << std::endl;
  stk::tools::disconnect_user_blocks(bulk, disconnectBlockVec, stk::tools::DisconnectBlocksOption(stk::tools::DISCONNECT_LOCAL, stk::tools::PRESERVE_INITIAL_HINGES));

  double disconnectTime = stk::wall_time();

  std::string filename = get_basename(exodusFileName);
  output_mesh(bulk, "disconnected_" + filename);
  double meshWriteTime = stk::wall_time();

  if (bulk.parallel_rank() == 0) {
    std::cout << "Mesh read time = " << (meshReadTime - startTime) << " s" << std::endl;
    std::cout << "Mesh write time = " << (meshWriteTime - disconnectTime) << " s" << std::endl;
  }
}

typedef TestDisconnectFullAlgorithm TestDisconnectUserBlocks;
typedef TestDisconnectFullAlgorithm2D TestDisconnectUserBlocks2D;

TEST_F(TestDisconnectUserBlocks, disconnect_user_blocks_8block_8hex_opposite_corners)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());

    BlockConnectionVector blockPairsToDisconnectVector{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0),
                                                       BlockConnection("block_3","block_7",0), BlockConnection("block_7","block_8",0), BlockConnection("block_5","block_7",0)};
    test_user_block_disconnect(get_bulk(), blockPairsToDisconnectVector, 13u);
  }
}

TEST_F(TestDisconnectUserBlocks, disconnect_user_blocks_preserve_snip_option_4block_4hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4hex(get_bulk());

    BlockConnectionVector blockPairsToDisconnectVector{BlockConnection("block_2","block_1",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_3",0)};
    test_user_block_disconnect(get_bulk(), blockPairsToDisconnectVector, 6u, stk::tools::DisconnectBlocksOption(stk::tools::DISCONNECT_GLOBAL, stk::tools::PRESERVE_INITIAL_HINGES));
  }
}

TEST_F(TestDisconnectUserBlocks, disconnect_user_blocks_snip_all_hinges_snip_option_4block_4hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4hex(get_bulk());

    BlockConnectionVector blockPairsToDisconnectVector{BlockConnection("block_2","block_1",0), BlockConnection("block_2","block_4",0)};
    test_user_block_disconnect(get_bulk(), blockPairsToDisconnectVector, 6u, stk::tools::DisconnectBlocksOption(stk::tools::DISCONNECT_GLOBAL, stk::tools::SNIP_ALL_HINGES));
  }
}

TEST_F(TestDisconnectUserBlocks2D, disconnect_user_blocks_2block_3quad_1hinge)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_3block_3quad_1hinge_linear_stack(get_bulk());

    BlockConnectionVector blockPairsToDisconnectVector{BlockConnection("block_1","block_3",0)};
    test_user_block_disconnect(get_bulk(), blockPairsToDisconnectVector, 3u, stk::tools::DisconnectBlocksOption(stk::tools::DISCONNECT_LOCAL, stk::tools::PRESERVE_INITIAL_HINGES));
  }
}

TEST_F(TestDisconnectUserBlocks2D, disconnect_user_blocks_3blocks_4quad_3proc)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(get_bulk(), 1);

    BlockConnectionVector blockPairsToDisconnectVector{BlockConnection("block_1","block_3",0), BlockConnection("block_2","block_3",0)};
    test_user_block_disconnect(get_bulk(), blockPairsToDisconnectVector, 2u, stk::tools::DisconnectBlocksOption(stk::tools::DISCONNECT_LOCAL, stk::tools::PRESERVE_INITIAL_HINGES));
  }
}

TEST_F(TestDisconnectUserBlocks2D, disconnect_user_blocks_3blocks_4quad_custom_ordinal)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_reverse_ordinal(get_bulk());

  BlockConnectionVector blockPairsToDisconnectVector{BlockConnection("vl","lateral",0), BlockConnection("radax","lateral",0)};
  test_user_block_disconnect(get_bulk(), blockPairsToDisconnectVector, 2u, stk::tools::DisconnectBlocksOption(stk::tools::DISCONNECT_LOCAL, stk::tools::PRESERVE_INITIAL_HINGES));

  stk::io::write_mesh("custom_ordinal_unit_test_mesh.g", get_bulk());
}


typedef TestDisconnectBlocks2D TestBlockPairCreation;

TEST_F(TestBlockPairCreation, test_block_pair_creation)
{
  setup_mesh_3block_4quad(get_bulk(), 1, 1);
  stk::mesh::PartVector parts;
  stk::tools::BlockPairVector partPairs;
  parts.push_back(get_meta().get_part("block_1"));
  parts.push_back(get_meta().get_part("block_2"));
  parts.push_back(get_meta().get_part("block_3"));
  stk::tools::impl::fill_ordered_block_pairs(parts, partPairs);

  EXPECT_EQ(3u, partPairs.size());
  EXPECT_EQ("block_1", partPairs[0].first->name());
  EXPECT_EQ("block_2", partPairs[0].second->name());
  EXPECT_EQ("block_1", partPairs[1].first->name());
  EXPECT_EQ("block_3", partPairs[1].second->name());
  EXPECT_EQ("block_2", partPairs[2].first->name());
  EXPECT_EQ("block_3", partPairs[2].second->name());
}

void create_ngs_jtd_sub_mesh(stk::mesh::BulkData& bulk)
{
  stk::unit_test_util::ConstructedMesh data(3);

  data.set_x_coordinates({ 2.69902090331971, 2.66623140978201, 2.75874573101832, 2.69902090331971,
                           2.69659806972163, 2.68053633989566, 2.66870203682596, 2.7137889539581,
                           2.67272988291214, 2.65367227486113, 2.76694241182748, 2.72228375387999,
                           2.68782101822595, 2.69659805697906, 2.71378890566682, 2.73602958156478});

  data.set_y_coordinates({ -2.54798187240236, -2.625, -2.625, -2.54798187240236, -2.66,
                           -2.67627405696173, -2.71269330410536, -2.625, -2.625, -2.66, -2.625,
                           -2.625, -2.625, -2.625, -2.67580141240548, -2.63747011593858 });

  data.set_z_coordinates({ -7.71929967118528, -7.76787460984145, -7.87919008575732,
                           -7.82910729232061, -7.83153012591869, -7.73778423460932,
                           -7.77420348175294, -7.81433924168221, -7.71000549971963,
                           -7.72979606467204, -7.87099340494817, -7.80584444176032,
                           -7.72193754541886, -7.83153013866126, -7.81433928997349, -7.84971773066203 });

  data.set_node_ids( {953299,  953347,  1046307, 1086643, 1086883, 1543476, 1598788, 3021447,
                      3121015, 3292103, 3357047, 3704977, 3705289, 3705425, 4715689, 5185015} );

  data.set_elem_ids( {25203969, 26788853, 27186527, 29884816, 30162566, 30162567,
                      30822787, 30824958, 30824959, 30824960, 30824961, 26778235, 26788861,
                      29528157, 29886393, 24040308, 24040314, 24045937, 24467153, 24467347,
                      24562607, 24562617, 29861241, 30822789} );

  data.set_elem_distribution( {{30824959, 1}, {26788861, 1}, {29886393, 1}, {24045937, 1}, {24467153, 1},
                               {24467347, 1}, {24562607, 1}, {24562617, 1}, {25203969, 1}, {26788853, 1},
                               {29884816, 1}, {30824960, 1}, {30824961, 1}, {24040308, 1}, {24040314, 1},
                               {30822787, 1}, {30824958, 1}, {26778235, 1}, {30822789, 1}, {30162567, 0},
                               {27186527, 1}, {29861241, 1}, {29528157, 0}, {30162566, 0}} );

  stk::unit_test_util::ConstructedElementBlock block1(stk::topology::TET_4, "block_1", 202, { {2, 8, 9, 10},
                                                                                              {10, 8, 9, 6},
                                                                                              {5, 15, 16, 8},
                                                                                              {2, 8, 10, 5},
                                                                                              {14, 16, 3, 8},
                                                                                              {14, 16, 8, 5},
                                                                                              {14, 5, 8, 2},
                                                                                              {7, 5, 8, 15},
                                                                                              {8, 5, 7, 10},
                                                                                              {7, 8, 6, 15},
                                                                                              {6, 8, 7, 10} });
  data.add_elem_block(block1);

  stk::unit_test_util::ConstructedElementBlock block2(stk::topology::TET_4, "block_2", 223, { {15, 12, 8, 6},
                                                                                              {8, 13, 9, 6},
                                                                                              {12, 16, 15, 8},
                                                                                              {8, 12, 13, 6} });
  data.add_elem_block(block2);

  stk::unit_test_util::ConstructedElementBlock block3(stk::topology::TET_4, "block_3", 245, { {2, 8, 4, 1},
                                                                                              {2, 8, 1, 9},
                                                                                              {8, 11, 3, 4},
                                                                                              {8, 12, 11, 4},
                                                                                              {8, 13, 1, 9},
                                                                                              {8, 12, 4, 1},
                                                                                              {8, 12, 1, 13},
                                                                                              {3, 4, 14, 8},
                                                                                              {4, 14, 8, 2} });
  data.add_elem_block(block3);

  data.create_mesh(bulk);
}

TEST(TestNGSDisconnect, jtd_sub_mesh)
{
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  create_ngs_jtd_sub_mesh(bulk);

  stk::mesh::Part* block1 = meta.get_part("block_1");
  stk::mesh::Part* block2 = meta.get_part("block_2");
  stk::mesh::Part* block3 = meta.get_part("block_3");

  stk::tools::DisconnectBlocksOption disconnectOption(stk::tools::DISCONNECT_LOCAL, stk::tools::PRESERVE_INITIAL_HINGES);
  stk::tools::BlockPairVector disconnectPairs{{block1, block3}, {block2, block3}};
  EXPECT_NO_THROW(stk::tools::disconnect_user_blocks(bulk, disconnectPairs, disconnectOption));
}
