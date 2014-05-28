/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <unit_tests/UnitTestStkSearchFixture.hpp>

namespace stk_search_unit {

OverlappingBoxes::OverlappingBoxes() {

  const int END = 25;
  unsigned current_id = 0;

  for (int i = 2; i<END-2; ++i) {
    for (int j = 2; j<END-2; ++j) {
      for (int k = 2; k<END-2; ++k) {
        BoundingVolume box;
        box.key.ident = current_id++;
        box.box[0] = i-1.5; box.box[1] = j-1.5; box.box[2] = k-1.5;
        box.box[3] = i+0.5; box.box[4] = j+0.5; box.box[5] = k+0.5;
        m_domain.push_back(box);
      }
    }
  }

  for (int i = 0; i<END; ++i) {
    for (int j = 0; j<END; ++j) {
      for (int k = 0; k<END; ++k) {
        BoundingVolume box;
        box.key.ident = current_id++;
        box.box[0] = i-1; box.box[1] = j-1; box.box[2] = k-1;
        box.box[3] = i+1; box.box[4] = j+1; box.box[5] = k+1;
        m_range.push_back(box);
      }
    }
  }
}

bool OverlappingBoxes::check_results(const IdentProcRelation & relation) const
{
  bool result = true;
  std::vector<unsigned> count;
  count.resize(m_domain.size());

  for (unsigned i = 0; i<relation.size(); ++i) {
    count[relation[i].first.ident]++;
  }

  for (unsigned i=0; i<count.size(); i++) {
    if (count[i] < EXPECTED_INTERSECTIONS) {
      result = false;
    }
  }

  return result;
}

}
