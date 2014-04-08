/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <ostream>                      // for basic_ostream::operator<<
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/util/Foreach.hpp>    // for stk_foreach
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector
#include "boost/foreach.hpp"            // for auto_any_base
#include "boost/range/iterator_range_core.hpp"  // for operator==, etc
#include "boost/range/iterator_range_io.hpp"  // for operator<<
#include "boost/range/sub_range.hpp"    // for sub_range, operator==
#include "gtest/gtest.h"                // for AssertHelper



STKUNIT_UNIT_TEST( UnitTestRange, range )
{
  typedef std::vector<int> IntVector;
  typedef boost::sub_range<IntVector> IntRange;

  IntVector array(10,0);

  IntRange a_range = std::make_pair(array.begin(),array.end());
  IntRange b_range(std::make_pair(array.begin(),array.end()));
  IntRange c_range(array.begin(),array.end() );
  IntRange d_range(array);

  STKUNIT_EXPECT_EQ( a_range, b_range );
  STKUNIT_EXPECT_EQ( a_range, c_range );
  STKUNIT_EXPECT_EQ( a_range, d_range );

  stk_foreach( int & i, a_range)
  {
    STKUNIT_EXPECT_EQ( i, 0 );
    i = 1;
  }

  stk_foreach( int i, c_range)
  {
    STKUNIT_EXPECT_EQ( i, 1 );
  }

  STKUNIT_EXPECT_TRUE( a_range == array);

}


