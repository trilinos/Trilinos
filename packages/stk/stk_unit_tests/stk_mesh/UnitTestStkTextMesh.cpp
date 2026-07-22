#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "mpi.h"

#include "UnitTestTextMeshFixture.hpp"

// STK specific versions of TextMesh
namespace
{
TEST_F(TestTextMeshAura, twoQuadShellWithCoordinates)
{
  //      4-----5-----6
  //      |     |     |
  //      |  1  |  2  |
  //      |     |     |
  //      1-----2-----3

  if (get_parallel_size() != 2) return;
  int rank = get_parallel_rank();

  std::string meshDesc =
      "0,1,SHELL_QUAD_4,1,2,5,4\n"
      "1,2,SHELL_QUAD_4,2,3,6,5"
      "|coordinates:  0,0,0, 1,0,0, 2,0,0, 0,1,0, 1,1,0, 2,1,0";
  std::vector<double> goldCoordinates = {0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1, 0, 2, 1, 0};

  setup_text_mesh(meshDesc);

  verify_num_elements(2);
  verify_single_element(1u, "SHELL_QUAD_4", EntityIdVector{1, 2, 5, 4});
  verify_single_element(2u, "SHELL_QUAD_4", EntityIdVector{2, 3, 6, 5});
  verify_coordinates(EntityIdVector{1, 2, 3, 4, 5, 6}, goldCoordinates);

  if (rank == 0) verify_shared_nodes(EntityIdVector{2,5}, 1);
  if (rank == 1) verify_shared_nodes(EntityIdVector{2,5}, 0);
}

TEST_F(TestTextMeshAura, sideBlockTopologyForConnectedHexes)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,HEX_8,1, 2, 3, 4, 5, 6, 7, 8,block_1\n"
      "0,2,HEX_8,5, 6, 7, 8, 9,10,11,12,block_2\n"
      "|sideset:name=surface_1; data=1,1, 2,1; split=block";

  setup_text_mesh(meshDesc);

  verify_num_elements(2);

  std::vector<std::string> sidesetSubsets;
  std::string prefix("surface_");
  sidesetSubsets.push_back(prefix + std::string("block_1") + "_" + get_topology_name("QUAD_4") + "_1");
  sidesetSubsets.push_back(prefix + std::string("block_2") + "_" + get_topology_name("QUAD_4") + "_1");

  verify_num_sidesets(1);
  verify_single_sideset("surface_1", 1, sidesetSubsets, SideVector{{1, 1}, {2, 1}});

  for (const std::string& subsetName : sidesetSubsets) {
    stk::mesh::Part* subsetPart = get_meta().get_part(subsetName);
    EXPECT_TRUE(nullptr != subsetPart);
    EXPECT_EQ(stk::topology::QUAD_4, subsetPart->topology());
  }
}

TEST_F(TestTextMeshAura, surfaceToBlockMapping_noSplit)
{
  if (get_parallel_size() != 2) return;

  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
      "1,2,HEX_8,2,9,10,3,6,11,12,7,block_2\n"
      "|sideset:name=left_surf;data=1,4";

  setup_text_mesh(meshDesc);

  verify_num_elements(2);
  verify_num_sidesets(1);

  const stk::mesh::Part* block_1 = get_meta().get_part("block_1");

  EXPECT_TRUE(nullptr != block_1);

  {
    std::string sidesetName("left_surf");

    stk::mesh::Part* sidesetPart = get_meta().get_part(sidesetName);
    EXPECT_TRUE(nullptr != sidesetPart);
    std::vector<const stk::mesh::Part*> goldParts{block_1};
    EXPECT_EQ(goldParts, get_meta().get_blocks_touching_surface(sidesetPart));
  }
}

TEST_F(TestTextMeshAura2d, surfaceToBlockMapping_splitByBlock)
{
  if (get_parallel_size() != 1) return;

  std::string meshDesc =
      "0,1,TRI_3_2D,3,1,4,block_1\n"
      "0,2,TRI_3_2D,1,2,4,block_1\n"
      "0,3,TRI_3_2D,2,5,4,block_1\n"
      "0,4,TRI_3_2D,5,7,4,block_2\n"
      "0,5,TRI_3_2D,7,6,4,block_2\n"
      "0,6,TRI_3_2D,6,3,4,block_2\n"
      "|coordinates: 0,0,0.1,0,0,0.1,0.05,0.1,0.1,0.1,0,0.2,0.1,0.2"
      "|dimension:2"
      "|sideset:name=skinned_surf; skin=all; split=block"
      "|sideset:name=internal_surf; data=1,3,3,2,4,3,6,2; split=block"
      "|sideset:name=shared_surf; data=1,1,6,1; split=block"
      "|sideset:name=owned_surf; data=5,1; split=block";

  setup_text_mesh(meshDesc);

  verify_num_elements(6);
  verify_num_sidesets(4);

  const stk::mesh::Part* block_1 = get_meta().get_part("block_1");
  const stk::mesh::Part* block_2 = get_meta().get_part("block_2");

  EXPECT_TRUE(nullptr != block_1);
  EXPECT_TRUE(nullptr != block_2);

  {
    std::string sidesetName("skinned_surf");

    stk::mesh::Part* sidesetPart = get_meta().get_part(sidesetName);
    EXPECT_TRUE(nullptr != sidesetPart);
    std::vector<const stk::mesh::Part*> goldParts{block_1, block_2};
    EXPECT_EQ(goldParts, get_meta().get_blocks_touching_surface(sidesetPart));
  }

  {
    std::string sidesetName("skinned_surf_" + std::string("block_1") + "_" + get_topology_name("LINE_2"));

    stk::mesh::Part* sidesetPart = get_meta().get_part(sidesetName);
    EXPECT_TRUE(nullptr != sidesetPart);
    std::vector<const stk::mesh::Part*> goldParts{block_1};
    EXPECT_EQ(goldParts, get_meta().get_blocks_touching_surface(sidesetPart));
  }

  {
    std::string sidesetName("skinned_surf_" + std::string("block_2") + "_" + get_topology_name("LINE_2"));

    stk::mesh::Part* sidesetPart = get_meta().get_part(sidesetName);
    EXPECT_TRUE(nullptr != sidesetPart);
    std::vector<const stk::mesh::Part*> goldParts{block_2};
    EXPECT_EQ(goldParts, get_meta().get_blocks_touching_surface(sidesetPart));
  }

  {
    std::string sidesetName("internal_surf");

    stk::mesh::Part* sidesetPart = get_meta().get_part(sidesetName);
    EXPECT_TRUE(nullptr != sidesetPart);
    std::vector<const stk::mesh::Part*> goldParts{block_1, block_2};
    EXPECT_EQ(goldParts, get_meta().get_blocks_touching_surface(sidesetPart));
  }

  {
    std::string sidesetName("internal_surf_" + std::string("block_1") + "_" + get_topology_name("LINE_2"));

    stk::mesh::Part* sidesetPart = get_meta().get_part(sidesetName);
    EXPECT_TRUE(nullptr != sidesetPart);
    std::vector<const stk::mesh::Part*> goldParts{block_1};
    EXPECT_EQ(goldParts, get_meta().get_blocks_touching_surface(sidesetPart));
  }

  {
    std::string sidesetName("internal_surf_" + std::string("block_2") + "_" + get_topology_name("LINE_2"));

    stk::mesh::Part* sidesetPart = get_meta().get_part(sidesetName);
    EXPECT_TRUE(nullptr != sidesetPart);
    std::vector<const stk::mesh::Part*> goldParts{block_2};
    EXPECT_EQ(goldParts, get_meta().get_blocks_touching_surface(sidesetPart));
  }

  {
    std::string sidesetName("shared_surf");

    stk::mesh::Part* sidesetPart = get_meta().get_part(sidesetName);
    EXPECT_TRUE(nullptr != sidesetPart);
    std::vector<const stk::mesh::Part*> goldParts{block_1, block_2};
    EXPECT_EQ(goldParts, get_meta().get_blocks_touching_surface(sidesetPart));
  }

  {
    std::string sidesetName("shared_surf_" + std::string("block_1") + "_" + get_topology_name("LINE_2"));

    stk::mesh::Part* sidesetPart = get_meta().get_part(sidesetName);
    EXPECT_TRUE(nullptr != sidesetPart);
    std::vector<const stk::mesh::Part*> goldParts{block_1};
    EXPECT_EQ(goldParts, get_meta().get_blocks_touching_surface(sidesetPart));
  }

  {
    std::string sidesetName("shared_surf_" + std::string("block_2") + "_" + get_topology_name("LINE_2"));

    stk::mesh::Part* sidesetPart = get_meta().get_part(sidesetName);
    EXPECT_TRUE(nullptr != sidesetPart);
    std::vector<const stk::mesh::Part*> goldParts{block_2};
    EXPECT_EQ(goldParts, get_meta().get_blocks_touching_surface(sidesetPart));
  }

  {
    std::string sidesetName("owned_surf");

    stk::mesh::Part* sidesetPart = get_meta().get_part(sidesetName);
    EXPECT_TRUE(nullptr != sidesetPart);
    std::vector<const stk::mesh::Part*> goldParts{block_2};
    EXPECT_EQ(goldParts, get_meta().get_blocks_touching_surface(sidesetPart));
  }

  {
    std::string sidesetName("owned_surf_" + std::string("block_2") + "_" + get_topology_name("LINE_2"));

    stk::mesh::Part* sidesetPart = get_meta().get_part(sidesetName);
    EXPECT_TRUE(nullptr != sidesetPart);
    std::vector<const stk::mesh::Part*> goldParts{block_2};
    EXPECT_EQ(goldParts, get_meta().get_blocks_touching_surface(sidesetPart));
  }
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
}
