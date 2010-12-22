/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <unistd.h>

#include <iostream>

#include <stk_util/environment/WallTime.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

STKUNIT_UNIT_TEST(UnitTestWallTime, UnitTest)
{
  double wall_now = stk::wall_time();
  
  ::sleep(1);
  
  double wall_delta = stk::wall_time() - wall_now;
  
  STKUNIT_ASSERT(wall_delta >= 1.0 && wall_delta <= 2.0);

  double wall_delta2 = stk::wall_dtime(wall_now);

  STKUNIT_ASSERT(wall_delta2 >= 1.0 && wall_delta2 <= 2.0);
}
