/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef UNIT_TEST_STK_SEARCH_FIXTURE_HPP
#define UNIT_TEST_STK_SEARCH_FIXTURE_HPP

#include <stk_search/BoundingBox.hpp>
#include <stk_search/diag/IdentProc.hpp>

namespace stk_search_unit {

class OverlappingBoxes {
  public:

    typedef stk_classic::search::ident::IdentProc<uint64_t,unsigned>            IdentProc;
    typedef std::vector<std::pair<IdentProc, IdentProc> >               IdentProcRelation;
    typedef stk_classic::search::box::AxisAlignedBoundingBox<IdentProc,float,3> BoundingVolume;
    typedef std::vector<BoundingVolume>                                 BoundingVolumeVector;

    OverlappingBoxes();

    inline
    const BoundingVolumeVector & domain() const { return m_domain; }

    inline
    const BoundingVolumeVector & range() const { return m_range; }

    bool check_results(const IdentProcRelation & relations) const;

  private:

    enum {EXPECTED_INTERSECTIONS = 64};

    BoundingVolumeVector m_domain;
    BoundingVolumeVector m_range;

};



}




#endif //UNIT_TEST_STK_SEARCH_FIXTURE_HPP
