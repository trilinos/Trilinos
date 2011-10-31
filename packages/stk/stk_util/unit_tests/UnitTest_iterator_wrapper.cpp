#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/util/random_access_iterator_wrapper.hpp>

#include <iostream>

STKUNIT_UNIT_TEST( random_access_iterator_wrapper, basic )
{
  typedef stk_util::random_access_iterator_wrapper<int> iterator;
  typedef stk_util::random_access_iterator_wrapper<const int> const_iterator;

  int a[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  {
    iterator i(a);
    const_iterator end(a+10);
    for (int tmp=1; i!=end; ++i, ++tmp) {
      STKUNIT_EXPECT_EQ(*i,tmp);
    }
  }

  {
    iterator i(a);
    for (int n=0; n<10; ++n) {
      STKUNIT_EXPECT_EQ(i[n],a[n]);
    }
  }



}
