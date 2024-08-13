#include <stddef.h>                     // for size_t
#include "gtest/gtest.h"                // for TEST_F
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"  // for ElemElemGraph
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_unit_test_utils/ElemGraphMultipleSharedSidesUtils.hpp"
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>

namespace {

class SkinMesh_TwoElemTwoSharedSide : public TwoElemTwoSharedSideTester {};

TEST_F(SkinMesh_TwoElemTwoSharedSide, skin_mesh)
{
  if(bulkData.parallel_size() <= 2)
  {
    stk::mesh::create_exposed_block_boundary_sides(bulkData, activePart, {&activePart, &skinPart});
    test_skinned_mesh(bulkData, 4u);
  }
}

TEST_F(SkinMesh_TwoElemTwoSharedSide, skin_one_hex)
{
  if(bulkData.parallel_size() <= 2)
  {
    remove_element_from_part(bulkData, 2, activePart);
    stk::mesh::Selector activeSelector = activePart;
    stk::mesh::Selector sel = activeSelector;
    stk::mesh::Selector air = !activeSelector;
    stk::mesh::create_exposed_block_boundary_sides(bulkData, sel, {&activePart, &skinPart}, air);
    test_total_sides_and_sides_per_element(bulkData, 6u, {6u, 2u});
  }
}


class SkinMesh_TwoElemThreeSharedSide : public TwoElemThreeSharedSideTester {};

TEST_F(SkinMesh_TwoElemThreeSharedSide, skin_mesh)
{
  if(bulkData.parallel_size() <= 2)
  {
    stk::mesh::create_exposed_block_boundary_sides(bulkData, activePart, {&activePart, &skinPart});
    test_skinned_mesh(bulkData, 3u);
  }
}

TEST_F(SkinMesh_TwoElemThreeSharedSide, skin_one_hex)
{
  if(bulkData.parallel_size() <= 2)
  {
    remove_element_from_part(bulkData, 2, activePart);
    stk::mesh::Selector activeSelector = activePart;
    stk::mesh::Selector sel = activeSelector;
    stk::mesh::Selector air = !activeSelector;
    stk::mesh::create_exposed_block_boundary_sides(bulkData, sel, {&activePart, &skinPart}, air);
    test_total_sides_and_sides_per_element(bulkData, 6u, {6u, 3u});
  }
}

class SkinMesh_TwoElemThreeSharedSideNoAura : public TwoElemThreeSharedSideNoAuraTester {};

TEST_F(SkinMesh_TwoElemThreeSharedSideNoAura, skin_mesh)
{
  if(bulkData.parallel_size() <= 2)
  {
    stk::mesh::create_exposed_block_boundary_sides(bulkData, activePart, {&activePart, &skinPart});
    test_skinned_mesh(bulkData, 3u);
  }
}

TEST_F(SkinMesh_TwoElemThreeSharedSideNoAura, skin_one_hex)
{
  if(bulkData.parallel_size() <= 2)
  {
    remove_element_from_part(bulkData, 2, activePart);
    stk::mesh::Selector activeSelector = activePart;
    stk::mesh::Selector sel = activeSelector;
    stk::mesh::Selector air = !activeSelector;
    stk::mesh::create_exposed_block_boundary_sides(bulkData, sel, {&activePart, &skinPart}, air);
    test_total_sides_and_sides_per_element(bulkData, 6u, {6u, 3u});
  }
}

}

