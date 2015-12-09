#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/ElemElemGraph.hpp>
#include <stk_mesh/base/ElemElemGraphImpl.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/environment/ReportHandler.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "stk_unit_test_utils/unittestMeshUtils.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>

#include <stk_unit_tests/stk_mesh/SetupKeyholeMesh.hpp>

#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_mesh/fixtures/heterogeneous_mesh.hpp>
#include <stk_mesh/fixtures/degenerate_mesh.hpp>

#include "UnitTestSkinMeshUseCaseUtils.hpp"

namespace {

using namespace stk::mesh::impl;
using namespace stk::mesh;


TEST(ElementGraph, skin_exposed_boundary)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

     if(stk::parallel_machine_size(comm) <= 2)
     {
         unsigned spatialDim = 3;

         stk::mesh::MetaData meta(spatialDim);
         stk::mesh::BulkData bulkData(meta, comm);
         make_2_hex_mesh_with_element1_inactive(bulkData);
         stk::mesh::PartVector skin_parts = get_skin_parts(meta);
         ElemGraphTestUtils::skin_boundary(bulkData, *meta.get_part("active"), skin_parts);
         ElemGraphTestUtils::test_num_faces_per_element(bulkData, {5u, 0u});
         std::vector<size_t> global_mesh_counts;
         stk::mesh::comm_mesh_counts(bulkData, global_mesh_counts);
         EXPECT_EQ(5u, global_mesh_counts[meta.side_rank()]);
     }
}


} // namespace
