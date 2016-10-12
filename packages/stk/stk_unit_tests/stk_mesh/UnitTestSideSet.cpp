#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>

class TestSideSet : public stk::unit_test_util::MeshFixture
{};

TEST_F(TestSideSet, creatingSideOfOneElem_eachProcHasOneSide)
{
    setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);
    get_bulk().initialize_face_adjacent_element_graph();

    stk::mesh::EntityVector localElems;
    stk::mesh::get_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEM_RANK), localElems);
    ASSERT_GE(localElems.size(), 1u);

    stk::mesh::StkSideSet sideSet;
    sideSet.push_back(stk::mesh::SideSetEntry(localElems[0], 1));

    EXPECT_EQ(0u, stk::mesh::count_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(get_meta().side_rank())));
    get_bulk().create_side_entities(sideSet, {});
    EXPECT_EQ(1u, stk::mesh::count_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(get_meta().side_rank())));
}
