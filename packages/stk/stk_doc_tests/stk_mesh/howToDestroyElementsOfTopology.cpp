#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_topology/topology.hpp>
#include <stk_io/FillMesh.hpp>
namespace
{
TEST(StkMeshHowTo, DestroyElementsOfTopology)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { GTEST_SKIP(); }

  std::unique_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();
  stk::io::fill_mesh("generated:1x1x4", *bulkPtr);

  EXPECT_EQ(4u, stk::mesh::count_entities(*bulkPtr, stk::topology::ELEM_RANK,  metaData.universal_part()));
  bulkPtr->destroy_elements_of_topology(stk::topology::HEX_8);
  EXPECT_EQ(0u, stk::mesh::count_entities(*bulkPtr, stk::topology::ELEM_RANK,  metaData.universal_part()));
}
}
