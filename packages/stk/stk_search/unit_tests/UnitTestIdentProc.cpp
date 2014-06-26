/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <gtest/gtest.h>

#include <stk_search/IdentProc.hpp>

#include <iostream>
#include <sstream>

namespace  {

TEST(stk_search, ident_proc)
{

  stk::search::IdentProc<int,int> a(1,0), b;
  b = a;
  stk::search::IdentProc<int,int> c(a), d(1,1), e(0,0);

  EXPECT_EQ(c.proc(), 0);
  EXPECT_EQ(c.id(), 1);

  EXPECT_EQ((a == b),true);
  EXPECT_EQ((a != d),true);
  EXPECT_EQ((a <  d),true);
  EXPECT_EQ((a >  e),true);
  EXPECT_EQ((a <= b),true);
  EXPECT_EQ((a <= d),true);
  EXPECT_EQ((a >= b),true);
  EXPECT_EQ((a >= e),true);

  EXPECT_EQ((a == d),false);
  EXPECT_EQ((a != b),false);
  EXPECT_EQ((a <  b),false);
  EXPECT_EQ((a >  b),false);
  EXPECT_EQ((a <= e),false);
  EXPECT_EQ((a >= d),false);

  {
    std::ostringstream out;
    out << a;
    EXPECT_EQ( out.str(), std::string("{id:1,proc:0}"));
  }

}

} // namespace

