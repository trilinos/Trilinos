/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "stk_util/registry/VersionNumber.hpp"

namespace stk
{
namespace util
{

TEST(VersionNumber, LessThan)
{
  const stk::util::VersionNumber a(4, 46);
  const stk::util::VersionNumber b(4, 48);
  const stk::util::VersionNumber c(3, 48);

  EXPECT_TRUE(a < b);
  EXPECT_FALSE(b < a);

  EXPECT_TRUE(c < b);
  EXPECT_FALSE(b < c);

  EXPECT_TRUE(c < a);
  EXPECT_FALSE(a < c);

  EXPECT_FALSE(a < a);
}

TEST(VersionNumber, SetCurrentFromString)
{
  stk::util::VersionNumber::set_current_version("4.46");
  EXPECT_EQ(stk::util::VersionNumber(4, 46), stk::util::VersionNumber::current_version());
  stk::util::VersionNumber::set_current_version("4.45.1");
  EXPECT_EQ(stk::util::VersionNumber(4, 45), stk::util::VersionNumber::current_version());
  stk::util::VersionNumber::set_current_version("4.47.5");
  EXPECT_EQ(stk::util::VersionNumber(4, 47), stk::util::VersionNumber::current_version());
  stk::util::VersionNumber::set_current_version("3.49.5-216-g400c8f6b-modified");
  EXPECT_EQ(stk::util::VersionNumber(3, 49), stk::util::VersionNumber::current_version());
  stk::util::VersionNumber::set_current_version("badstring");
  EXPECT_EQ(stk::util::VersionNumber(-1, -1), stk::util::VersionNumber::current_version());
}
}
}
