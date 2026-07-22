
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Comm.hpp>       // for comm_mesh_counts
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <vector>                       // for vector
#include "UnitTestSkinMeshUseCaseUtils.hpp"  // for get_skin_parts, etc
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Types.hpp"      // for PartVector
#include "stk_unit_test_utils/ElemGraphTestUtils.hpp"

namespace {

using namespace stk::mesh::impl;
using namespace stk::mesh;

TEST(ElementGraph, skin_exposed_boundary)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) <= 2)
  {
    unsigned spatialDim = 3;

    stk::mesh::MeshBuilder builder(comm);
    builder.set_spatial_dimension(spatialDim);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
    stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
    make_2_hex_mesh_with_element1_inactive(*bulkPtr);
    stk::mesh::PartVector skin_parts = get_skin_parts(meta);
    ElemGraphTestUtils::skin_boundary(*bulkPtr, *meta.get_part("active"), skin_parts);
    ElemGraphTestUtils::test_num_faces_per_element(*bulkPtr, {5u, 0u});
    std::vector<size_t> global_mesh_counts;
    stk::mesh::comm_mesh_counts(*bulkPtr, global_mesh_counts);
    EXPECT_EQ(5u, global_mesh_counts[meta.side_rank()]);
  }
}


} // namespace
