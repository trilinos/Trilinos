#include <gtest/gtest.h>

#include <samba/utility.hpp>
#include <samba/utility/debug_message.hpp>

#include <boost/type_traits.hpp>

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

TEST(samba, interval_set)
{
  using namespace samba;

  typedef interval<A> interval;
  typedef interval_set<A> interval_set;

  interval_set s;

  EXPECT_TRUE( s.empty() );
  EXPECT_EQ( s.size(), 0u);

  s += interval(0,10);
  s -= A::create(5);

  EXPECT_FALSE( s.empty() );
  EXPECT_EQ( s.size(), 9u);

  {
    std::ostringstream out;
    out << s;
    EXPECT_EQ(out.str(),"{A:[0,5),[6,10),\b}");
  }


  {
    int count=0;
    for ( interval_set::const_iterator b=s.begin(), e=s.end();
        b != e;
        ++b,++count)
    {
      if(count==5) ++count;
      A a = {count};
      EXPECT_EQ(*b,a);
    }
  }

  {
    int count=0;
    for ( interval_set::const_interval_iterator b=s.interval_begin(), e=s.interval_end();
        b != e;
        ++b,++count)
    {
      if (count==0) EXPECT_EQ( interval(0,5), *b );
      if (count==1) EXPECT_EQ( interval(6,10), *b );
    }
  }

}
