#include <stddef.h>                     // for size_t
#include "gtest/gtest.h"                // for TEST_F
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"  // for ElemElemGraph
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_unit_test_utils/ElemGraphMultipleSharedSidesUtils.hpp"
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>

namespace {

TEST_F(TwoElemTwoSharedSideTester, skin_mesh)
{
     if(bulkData.parallel_size() <= 2)
     {
         {
             stk::mesh::create_exposed_block_boundary_sides(bulkData, activePart, {&activePart, &skinPart});
         }
         //stk::mesh::skin_mesh( bulkData, activePart, {&activePart, &skinPart});
         test_skinned_mesh(bulkData, 4u);
     }
}

TEST_F(TwoElemTwoSharedSideTester, skin_one_hex)
{
     if(bulkData.parallel_size() <= 2)
     {
         remove_element_from_part(bulkData, 2, activePart);
         stk::mesh::Selector activeSelector = activePart;
         {
             stk::mesh::Selector sel = activeSelector;
             stk::mesh::Selector air = !activeSelector;
             stk::mesh::create_exposed_block_boundary_sides(bulkData, sel, {&activePart, &skinPart}, air);
         }
         //stk::mesh::skin_mesh( bulkData, activePart, {&activePart, &skinPart}, &activeSelector);
         test_total_sides_and_sides_per_element(bulkData, 6u, {6u, 2u});
         stk::io::write_mesh("doublyKissingHexes.e", bulkData);
     }
}


TEST_F(TwoElemThreeSharedSideTester, skin_mesh)
{
     if(bulkData.parallel_size() <= 2)
     {
         {
             stk::mesh::create_exposed_block_boundary_sides(bulkData, activePart, {&activePart, &skinPart});
         }
         //stk::mesh::skin_mesh( bulkData, activePart, {&activePart, &skinPart});
         test_skinned_mesh(bulkData, 3u);
     }
}

TEST_F(TwoElemThreeSharedSideTester, skin_one_hex)
{
     if(bulkData.parallel_size() <= 2)
     {
         remove_element_from_part(bulkData, 2, activePart);
         stk::mesh::Selector activeSelector = activePart;
         {
             stk::mesh::Selector sel = activeSelector;
             stk::mesh::Selector air = !activeSelector;
             stk::mesh::create_exposed_block_boundary_sides(bulkData, sel, {&activePart, &skinPart}, air);
         }
         //stk::mesh::skin_mesh( bulkData, activePart, {&activePart, &skinPart}, &activeSelector);
         test_total_sides_and_sides_per_element(bulkData, 6u, {6u, 3u});
         stk::io::write_mesh("triplyKissingHexes.e", bulkData);
     }
}


TEST_F(TwoElemThreeSharedSideNoAuraTester, skin_mesh)
{
     if(bulkData.parallel_size() <= 2)
     {
         {
             stk::mesh::create_exposed_block_boundary_sides(bulkData, activePart, {&activePart, &skinPart});
         }
         //stk::mesh::skin_mesh( bulkData, activePart, {&activePart, &skinPart});
         test_skinned_mesh(bulkData, 3u);
     }
}

TEST_F(TwoElemThreeSharedSideNoAuraTester, skin_one_hex)
{
     if(bulkData.parallel_size() <= 2)
     {
         remove_element_from_part(bulkData, 2, activePart);
         stk::mesh::Selector activeSelector = activePart;
         {
             stk::mesh::Selector sel = activeSelector;
             stk::mesh::Selector air = !activeSelector;
             stk::mesh::create_exposed_block_boundary_sides(bulkData, sel, {&activePart, &skinPart}, air);
         }
         //stk::mesh::skin_mesh( bulkData, activePart, {&activePart, &skinPart}, &activeSelector);
         test_total_sides_and_sides_per_element(bulkData, 6u, {6u, 3u});
     }
}

}

