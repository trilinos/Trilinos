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

  std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,5,4\n"
                         "1,2,SHELL_QUAD_4,2,3,6,5";
  std::vector<double> coordinates = {
    0,0,0, 1,0,0, 2,0,0,
    0,1,0, 1,1,0, 2,1,0
  };

  setup_text_mesh(meshDesc, coordinates);

  verify_num_elements(2);
  verify_single_element(1u, "SHELL_QUAD_4", EntityIdVector{1, 2, 5, 4});
  verify_single_element(2u, "SHELL_QUAD_4", EntityIdVector{2, 3, 6, 5});
  verify_coordinates(EntityIdVector{1, 2, 3, 4, 5, 6}, coordinates);

  if (rank == 0) verify_shared_nodes(EntityIdVector{2,5}, 1);
  if (rank == 1) verify_shared_nodes(EntityIdVector{2,5}, 0);
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
