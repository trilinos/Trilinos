/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_search/IdentProc.hpp>

#include <iostream>
#include <sstream>

namespace  {

STKUNIT_UNIT_TEST(stk_search, ident_proc)
{

  stk::search::IdentProc<int,int> a(1,0), b;
  b = a;
  stk::search::IdentProc<int,int> c(a), d(1,1), e(0,0);

  STKUNIT_EXPECT_EQUAL(c.proc(), 0);
  STKUNIT_EXPECT_EQUAL(c.id(), 1);

  STKUNIT_EXPECT_EQUAL((a == b),true);
  STKUNIT_EXPECT_EQUAL((a != d),true);
  STKUNIT_EXPECT_EQUAL((a <  d),true);
  STKUNIT_EXPECT_EQUAL((a >  e),true);
  STKUNIT_EXPECT_EQUAL((a <= b),true);
  STKUNIT_EXPECT_EQUAL((a <= d),true);
  STKUNIT_EXPECT_EQUAL((a >= b),true);
  STKUNIT_EXPECT_EQUAL((a >= e),true);

  STKUNIT_EXPECT_EQUAL((a == d),false);
  STKUNIT_EXPECT_EQUAL((a != b),false);
  STKUNIT_EXPECT_EQUAL((a <  b),false);
  STKUNIT_EXPECT_EQUAL((a >  b),false);
  STKUNIT_EXPECT_EQUAL((a <= e),false);
  STKUNIT_EXPECT_EQUAL((a >= d),false);

  {
    std::ostringstream out;
    out << a;
    STKUNIT_EXPECT_EQUAL( out.str(), std::string("{id:1,proc:0}"));
  }

}

} // namespace

