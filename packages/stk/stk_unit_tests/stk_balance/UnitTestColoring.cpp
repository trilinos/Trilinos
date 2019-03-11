#include <gtest/gtest.h>

#include "stk_balance/balanceUtils.hpp"
#include <stk_balance/balance.hpp>

#include "stk_mesh/base/GetEntities.hpp"

#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_util/diag/StringUtil.hpp"

namespace {

class BasicColoring : public stk::unit_test_util::MeshFixture {};

void test_adjacent_elements_have_different_coloring(const stk::mesh::BulkData& bulk)
{
    stk::mesh::EntityVector nodes;
    const stk::mesh::MetaData& metaData = bulk.mesh_meta_data();
    stk::mesh::get_selected_entities(metaData.locally_owned_part(), bulk.buckets(stk::topology::NODE_RANK), nodes);

    for (stk::mesh::Entity& node : nodes)
    {
        std::set<stk::mesh::Part*> coloringParts;

        unsigned numElems = bulk.num_elements(node);
        unsigned numLocallyOwnedElems = 0;
        const stk::mesh::Entity* elems = bulk.begin_elements(node);
        for (unsigned i = 0; i < numElems; ++i)
        {
            stk::mesh::Entity elem = elems[i];
            if (bulk.bucket(elem).owned()) {
                ++numLocallyOwnedElems;
                stk::mesh::Part* colorPart = stk::balance::get_coloring_part(bulk, elem);
                if (nullptr != colorPart) {
                    coloringParts.insert(colorPart);
                }
            }
        }
        EXPECT_EQ(numLocallyOwnedElems, coloringParts.size());
    }
}

TEST_F(BasicColoring, checkAdjacentElementsHaveDifferentColors)
{
    if (stk::parallel_machine_size(get_comm()) > 3) return;
    setup_mesh("generated:3x3x3", stk::mesh::BulkData::AUTO_AURA);

    stk::balance::BasicColoringSettings coloringSettings;
    bool meshIsColored = stk::balance::colorStkMesh(coloringSettings, *bulkData);
    EXPECT_TRUE(meshIsColored);
    test_adjacent_elements_have_different_coloring(*bulkData);
}

TEST_F(BasicColoring, checkForCorrectColors)
{
    if (stk::parallel_machine_size(get_comm()) > 1) return;
    setup_mesh("generated:3x3x3", stk::mesh::BulkData::AUTO_AURA);

    stk::balance::BasicColoringSettings coloringSettings;
    bool meshIsColored = stk::balance::colorStkMesh(coloringSettings, *bulkData);
    EXPECT_TRUE(meshIsColored);

    stk::mesh::PartVector coloringParts;
    stk::balance::fill_coloring_parts(*metaData, coloringParts);

    size_t goldNumberOfColors = 9;
    EXPECT_EQ(goldNumberOfColors, coloringParts.size());

    for (size_t i = 1; i <= goldNumberOfColors; ++i)
    {
        std::string goldPartName = stk::balance::construct_coloring_part_name(i);
        bool colorFound = false;
        for (stk::mesh::Part* part : coloringParts)
        {
          if (part->name() == goldPartName)
          {
            colorFound = true;
            break;
          }
        }
        EXPECT_TRUE(colorFound);
    }
}
}
