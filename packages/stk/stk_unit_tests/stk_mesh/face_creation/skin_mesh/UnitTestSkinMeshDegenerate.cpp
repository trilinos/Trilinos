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
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "stk_unit_test_utils/unittestMeshUtils.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>

#include "SetupKeyholeMesh.hpp"
#include "stk_unit_test_utils/stk_mesh_fixtures/QuadFixture.hpp"  // for QuadFixture
#include "stk_unit_test_utils/stk_mesh_fixtures/degenerate_mesh.hpp"
#include "stk_unit_test_utils/stk_mesh_fixtures/heterogeneous_mesh.hpp"

#include "UnitTestSkinMeshUseCaseUtils.hpp"

namespace {

using namespace stk::mesh::impl;
using namespace stk::mesh;
using stk::unit_test_util::build_mesh;

TEST(ElementGraph, degenerate_mesh)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) <= 2)
  {
    std::string fileName("degenerate.g");
    {
      std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_SELF, stk::mesh::BulkData::NO_AUTO_AURA);
      stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
      stk::mesh::fixtures::VectorFieldType & node_coord =
          meta_data.declare_field<double>(stk::topology::NODE_RANK, "coordinates");
      stk::mesh::put_field_on_mesh(node_coord, meta_data.universal_part(), 3, nullptr);

      stk::mesh::fixtures::degenerate_mesh_meta_data(meta_data, node_coord);
      meta_data.commit();

      stk::mesh::BulkData& bulk_data = *bulkPtr;
      stk::mesh::fixtures::degenerate_mesh_bulk_data(bulk_data, node_coord);
      if(stk::parallel_machine_rank(comm) == 0)
      {
        stk::io::write_mesh(fileName, bulk_data);
      }
    }
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, comm, stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& bulk_data = *bulkPtr;
    stk::mesh::Part &skin = meta_data.declare_part("skin", meta_data.side_rank());
    stk::io::put_io_part_attribute(skin);
    stk::unit_test_util::read_from_serial_file_and_decompose(fileName, bulk_data, "RIB");
    unlink(fileName.c_str());
    EXPECT_NO_FATAL_FAILURE(ElemGraphTestUtils::skin_boundary(bulk_data, meta_data.locally_owned_part(), {&skin}));
    std::vector<size_t> mesh_counts;
    stk::mesh::comm_mesh_counts(bulk_data, mesh_counts);
    EXPECT_EQ(10u, mesh_counts[meta_data.side_rank()]);
  }
}

}
