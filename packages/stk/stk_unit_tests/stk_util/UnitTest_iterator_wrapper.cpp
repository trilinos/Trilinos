
#include <gtest/gtest.h>
#include <stk_util/util/random_access_iterator_wrapper.hpp>
#include "boost/iterator/iterator_facade.hpp"  // for iterator_facade, etc

TEST( random_access_iterator_wrapper, basic )
{
  typedef stk_util::random_access_iterator_wrapper<int> iterator;
  typedef stk_util::random_access_iterator_wrapper<const int> const_iterator;

  int a[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  {
    iterator i(a);
    const_iterator end(a+10);
    for (int tmp=1; i!=end; ++i, ++tmp) {
      EXPECT_EQ(*i,tmp);
    }
  }

  {
    iterator i(a);
    for (int n=0; n<10; ++n) {
      EXPECT_EQ(i[n],a[n]);
    }
  }



}
