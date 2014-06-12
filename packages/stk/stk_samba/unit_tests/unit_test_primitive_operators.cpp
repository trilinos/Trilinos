#include <gtest/gtest.h>

#include <samba/utility.hpp>

#include <iostream>
#include <sstream>

//put inside of samba namespace to test
//argument dependant lookup (ADL)
//for operators
namespace samba {

SAMBA_MAKE_PRIMITIVE(A,A,int);

} //namespace samba

//the following macros must be called outside any
//namespace without a trailing semi-colon
SAMBA_IS_PRIMITIVE(samba::A)


TEST(samba, compare_primitives)
{
  using namespace samba;

  A a = {0}, b = {1}, c={0};

  EXPECT_EQ(a,c);
  EXPECT_NE(a,b);

  EXPECT_LT(a,b);
  EXPECT_LE(a,b);

  EXPECT_GT(b,c);
  EXPECT_GE(b,c);
}


TEST(samba, compare_primitive_against_value_type)
{
  using namespace samba;

  A a = {0};
  int b = 1;
  int c = 0;


  EXPECT_EQ(a,c);
  EXPECT_NE(a,b);

  EXPECT_LT(a,b);
  EXPECT_LE(a,b);

  EXPECT_GT(b,a);
  EXPECT_GE(b,a);
}

TEST(samba, assign_primitive_to_value)
{
  using namespace samba;

  {
    A a;
    a = 0;
  }

  {
    int i=0;
    for (A a={0},e={10}; a<e; ++a,++i)
    {
      EXPECT_EQ(a(),i);
    }
  }

}



