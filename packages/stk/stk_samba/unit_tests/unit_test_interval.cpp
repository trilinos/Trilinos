#include <gtest/gtest.h>

#include <samba/utility.hpp>

#include <iostream>
#include <sstream>

//put inside of samba namespace to test
//argument dependant lookup (ADL)
//for operators
namespace samba {
SAMBA_MAKE_PRIMITIVE(A,A,int);
}

//the following macros must be declared outside of any namespace
//and called without the trailing semicolon
SAMBA_IS_PRIMITIVE(samba::A)


TEST(samba, interval)
{
  using namespace samba;

  typedef interval<A> interval;

  interval i;

  EXPECT_TRUE( i.empty() );
  EXPECT_EQ( i.size(), 0u);

  i = interval(0,10);

  EXPECT_FALSE( i.empty() );
  EXPECT_EQ( i.size(), 10u);

  {
    std::ostringstream out;
    out << i;
    EXPECT_EQ(out.str(),"{A:[0,10)}");
  }

  {
    int count = 0;
    for ( interval::const_iterator b=i.begin(), e=i.end();
        b != e;
        ++b, ++count )
    {
      EXPECT_EQ( *b, A::create(count));
    }
  }

}

