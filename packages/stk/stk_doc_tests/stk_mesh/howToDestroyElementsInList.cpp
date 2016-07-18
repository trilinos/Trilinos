#include <gtest/gtest.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_mesh/base/DestroyElements.hpp>

namespace
{
//BEGIN_DESTROY_ELEMENTS_IN_LIST
TEST(StkMeshHowTo, DestroyElementsInList)
{
    stk::mesh::MetaData metaData;
    stk::mesh::BulkData bulkData(metaData, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData);
    EXPECT_GT(stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::ELEM_RANK)), 0u);
    stk::mesh::EntityVector elementsToDestroy{bulkData.get_entity(stk::topology::ELEMENT_RANK,1)};
    stk::mesh::destroy_elements(bulkData, elementsToDestroy);

    stk::mesh::EntityVector orphanedNodes{
        bulkData.get_entity(stk::topology::NODE_RANK,1),
        bulkData.get_entity(stk::topology::NODE_RANK,2),
        bulkData.get_entity(stk::topology::NODE_RANK,3),
        bulkData.get_entity(stk::topology::NODE_RANK,4)
    };

    for(stk::mesh::Entity node : orphanedNodes)
        EXPECT_FALSE(bulkData.is_valid(node));
}
//END_DESTROY_ELEMENTS_IN_LIST
}
