#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
namespace
{
TEST(StkMeshHowTo, DestroyElementsOfTopology)
{
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
  bulkPtr->mesh_meta_data().use_simple_fields();
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulkData = *bulkPtr;
  stk::io::fill_mesh("generated:1x1x4", bulkData);
  EXPECT_GT(stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::ELEM_RANK)), 0u);
  bulkData.destroy_elements_of_topology(stk::topology::HEX_8);
  EXPECT_EQ(0u, stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::ELEM_RANK)));
}
}
