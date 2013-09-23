
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>

namespace stk_search_unit {

STKUNIT_UNIT_TEST(CoarseSearch, 2D)
{
  typedef stk::search::ident::IdentProc<int,unsigned> ID;
  typedef stk::search::box::AxisAlignedBoundingBox<ID, float, 2 /*spatial dimension*/> AABox2;

  typedef std::pair<ID,ID> Domain2Range;
  typedef std::vector<Domain2Range>   DomainToRangeVector;
  typedef std::vector<AABox2>          AABoxVector;

  stk::search::FactoryOrder order;
  DomainToRangeVector domain_2_range;

  AABoxVector range, domain;

  float box[4] = {};

  {
    int id = 0;
    //add to the range
    for(int i=0; i<10; ++i) {
      box[0] = 2*i;     // lower x
      box[1] = 0;       // lower y
      box[2] = 2*(i+1); // upper x
      box[3] = 1;       // upper y

      domain.push_back(AABox2(box,ID(id++,0)));
    }


    //add to the domain
    for(int i=0; i<10; ++i) {
      box[0] = 2*i + 1;     // lower x
      box[1] = 0;           // lower y
      box[2] = 2*(i+1) + 1; // upper x
      box[3] = 1;           // upper y

      range.push_back(AABox2(box,ID(id++,0)));
    }
  }

  stk::search::coarse_search( domain_2_range, range, domain, order);

  DomainToRangeVector expected_domain_2_range;
  expected_domain_2_range.push_back(Domain2Range(ID(0,0),ID(10,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(1,0),ID(10,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(1,0),ID(11,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(2,0),ID(11,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(2,0),ID(12,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(3,0),ID(12,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(3,0),ID(13,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(4,0),ID(13,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(4,0),ID(14,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(5,0),ID(14,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(5,0),ID(15,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(6,0),ID(15,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(6,0),ID(16,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(7,0),ID(16,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(7,0),ID(17,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(8,0),ID(17,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(8,0),ID(18,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(9,0),ID(18,0)));
  expected_domain_2_range.push_back(Domain2Range(ID(9,0),ID(19,0)));

  EXPECT_TRUE(std::equal(expected_domain_2_range.begin(), expected_domain_2_range.end(), domain_2_range.begin()));

}

} // namespace stk_search_unit
