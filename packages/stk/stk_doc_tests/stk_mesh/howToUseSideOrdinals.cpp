#include <gtest/gtest.h>

#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace
{

//-BEGIN
TEST(StkMeshHowTo, useSideOrdinalsOnMeshWithMixedSideRanks)
{
  // ============================================================
  // INITIALIZATION
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(communicator).create();

  // Create a SHELL_QUAD_4 with 2 faces and 2 edges (on opposite sides)
  const std::string meshDesc =
         "textmesh:0,1,SHELL_QUAD_4, 1,2,3,4, block_1|sideset:name=surface_1; data=1,1 ,1,2, 1,3, 1,5";

  stk::io::fill_mesh(meshDesc, *bulk);

  stk::mesh::Entity elem1 = bulk->get_entity(stk::topology::ELEM_RANK, 1);
  ASSERT_EQ(2u, bulk->num_faces(elem1));
  ASSERT_EQ(2u, bulk->num_edges(elem1));

  stk::topology elemTopo = bulk->bucket(elem1).topology();
  EXPECT_EQ(stk::topology::SHELL_QUAD_4, elemTopo);

  // Verify mapping
  const stk::mesh::ConnectivityOrdinal* faceOrdinals = bulk->begin_face_ordinals(elem1);
  const stk::mesh::ConnectivityOrdinal* edgeOrdinals = bulk->begin_edge_ordinals(elem1);

  unsigned goldSideOrdinal;
  unsigned rankedOrdinal, goldRankedOrdinal;
  stk::topology::rank_t rank, goldRank;

  auto sidesets = bulk->get_sidesets();
  EXPECT_EQ(1u, sidesets.size());

  // (0, FACE_RANK) -> 0
  goldSideOrdinal = 0;
  goldRankedOrdinal = 0;
  goldRank = stk::topology::FACE_RANK;
  EXPECT_EQ(goldRankedOrdinal, faceOrdinals[0]);
  EXPECT_EQ(goldSideOrdinal, elemTopo.side_ordinal(goldRankedOrdinal, goldRank));
  // Reverse mapping: 0 -> (0, FACE_RANK)
  elemTopo.ranked_side_ordinal(goldSideOrdinal, rankedOrdinal, rank);
  EXPECT_EQ(goldRankedOrdinal, rankedOrdinal);
  EXPECT_EQ(goldRank, rank);
  // Check sideset entry existence
  EXPECT_TRUE(sidesets[0]->contains(elem1, goldSideOrdinal));

  // (1, FACE_RANK) -> 1
  goldSideOrdinal = 1;
  goldRankedOrdinal = 1;
  goldRank = stk::topology::FACE_RANK;
  EXPECT_EQ(goldRankedOrdinal, faceOrdinals[1]);
  EXPECT_EQ(goldSideOrdinal, elemTopo.side_ordinal(goldRankedOrdinal, goldRank));
  // Reverse mapping: 1 -> (1, FACE_RANK)
  elemTopo.ranked_side_ordinal(goldSideOrdinal, rankedOrdinal, rank);
  EXPECT_EQ(goldRankedOrdinal, rankedOrdinal);
  EXPECT_EQ(goldRank, rank);
  // Check sideset entry existence
  EXPECT_TRUE(sidesets[0]->contains(elem1, goldSideOrdinal));

  // (0, EDGE_RANK) -> 2
  goldSideOrdinal = 2;
  goldRankedOrdinal = 0;
  goldRank = stk::topology::EDGE_RANK;
  EXPECT_EQ(goldRankedOrdinal, edgeOrdinals[0]);
  EXPECT_EQ(goldSideOrdinal, elemTopo.side_ordinal(goldRankedOrdinal, goldRank));
  // Reverse mapping: 2 -> (0, EDGE_RANK)
  elemTopo.ranked_side_ordinal(goldSideOrdinal, rankedOrdinal, rank);
  EXPECT_EQ(goldRankedOrdinal, rankedOrdinal);
  EXPECT_EQ(goldRank, rank);
  // Check sideset entry existence
  EXPECT_TRUE(sidesets[0]->contains(elem1, goldSideOrdinal));

  // (2, EDGE_RANK) -> 4
  goldSideOrdinal = 4;
  goldRankedOrdinal = 2;
  goldRank = stk::topology::EDGE_RANK;
  EXPECT_EQ(goldRankedOrdinal, edgeOrdinals[1]);
  EXPECT_EQ(goldSideOrdinal, elemTopo.side_ordinal(goldRankedOrdinal, goldRank));
  // Reverse mapping: 4 -> (2, EDGE_RANK)
  elemTopo.ranked_side_ordinal(goldSideOrdinal, rankedOrdinal, rank);
  EXPECT_EQ(goldRankedOrdinal, rankedOrdinal);
  EXPECT_EQ(goldRank, rank);
  // Check sideset entry existence
  EXPECT_TRUE(sidesets[0]->contains(elem1, goldSideOrdinal));
}
//-END

}
