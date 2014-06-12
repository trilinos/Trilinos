/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef unit_test_UnitTestRingFixture_hpp
#define unit_test_UnitTestRingFixture_hpp

#include <stk_mesh/fixtures/RingFixture.hpp>

// must be in same namespace as BulkData to make friending w/out exposure possible

namespace stk_classic {
namespace unit_test {

/**
 * TODO - Document what this does
 */
void test_shift_ring( stk_classic::mesh::fixtures::RingFixture& ring,
                      bool generate_aura=true );

} // namespace unit_test
} // namespace stk_classic

#endif
