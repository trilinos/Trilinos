#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
namespace
{
TEST(StkMeshHowTo, DestroyElementsOfTopology)
{
    stk::mesh::MetaData metaData;
    stk::mesh::BulkData bulkData(metaData, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData);
    EXPECT_GT(stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::ELEM_RANK)), 0u);
    bulkData.destroy_elements_of_topology(stk::topology::HEX_8);
    EXPECT_EQ(0u, stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::ELEM_RANK)));
}
}
