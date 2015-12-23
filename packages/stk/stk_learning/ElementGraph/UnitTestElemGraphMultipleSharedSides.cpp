#include <gtest/gtest.h>
#include "stk_unit_test_utils/ElemGraphMultipleSharedSidesUtils.hpp"

namespace {

TEST_F(TwoElemTwoSharedSideTester, elem_death)
{
     if(bulkData.parallel_size() <= 2)
     {
         test_element_death_with_multiple_shared_sides(bulkData, activePart, skinPart);
     }
}

TEST_F(TwoElemTwoSharedSideTester, double_kissing_hexes)
{
     if(bulkData.parallel_size() <= 2)
     {
         test_elems_kissing_n_times(bulkData, activePart, 2);
     }
}


TEST_F(TwoElemThreeSharedSideTester, triple_kissing_hexes)
{
     if(bulkData.parallel_size() <= 2)
     {
         test_elems_kissing_n_times(bulkData, activePart, 3);
     }
}

TEST_F(TwoElemThreeSharedSideNoAuraTester, triple_kissing_hexes)
{
     if(bulkData.parallel_size() <= 2)
     {
         test_elems_kissing_n_times(bulkData, activePart, 3);
     }
}

TEST_F(TwoElem2dTwoSharedSideTester, double_kissing_quads)
{
     if(bulkData.parallel_size() <= 2)
     {
         test_elems_kissing_n_times(bulkData, activePart, 2);
     }
}


}
