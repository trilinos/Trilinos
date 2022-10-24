#include <gtest/gtest.h>
#include "stk_unit_test_utils/ElemGraphMultipleSharedSidesUtils.hpp"

namespace {

class ElemGraph_TwoElemTwoSharedSide : public simple_fields::TwoElemTwoSharedSideTester {};

TEST_F(ElemGraph_TwoElemTwoSharedSide, elem_death)
{
  if(bulkData.parallel_size() <= 2)
  {
    test_element_death_with_multiple_shared_sides(bulkData, activePart, skinPart);
  }
}

TEST_F(ElemGraph_TwoElemTwoSharedSide, double_kissing_hexes)
{
  if(bulkData.parallel_size() <= 2)
  {
    test_elems_kissing_n_times(bulkData, activePart, 2);
  }
}

class ElemGraph_TwoElemThreeSharedSide : public simple_fields::TwoElemThreeSharedSideTester {};

TEST_F(ElemGraph_TwoElemThreeSharedSide, triple_kissing_hexes)
{
  if(bulkData.parallel_size() <= 2)
  {
    test_elems_kissing_n_times(bulkData, activePart, 3);
  }
}

class ElemGraph_TwoElemThreeSharedSideNoAura : public simple_fields::TwoElemThreeSharedSideNoAuraTester {};

TEST_F(ElemGraph_TwoElemThreeSharedSideNoAura, triple_kissing_hexes)
{
  if(bulkData.parallel_size() <= 2)
  {
    test_elems_kissing_n_times(bulkData, activePart, 3);
  }
}

class ElemGraph_TwoElem2dTwoSharedSide : public simple_fields::TwoElem2dTwoSharedSideTester {};

TEST_F(ElemGraph_TwoElem2dTwoSharedSide, double_kissing_quads)
{
  if(bulkData.parallel_size() <= 2)
  {
    test_elems_kissing_n_times(bulkData, activePart, 2);
  }
}


}
