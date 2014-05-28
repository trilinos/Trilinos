#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/util/nested_iterator.hpp>
#include <stk_util/util/nested_range.hpp>

#include <vector>
#include <list>
#include <iostream>

namespace vector_vector_int {

  typedef std::vector< std::vector<int> > nested_type;

  typedef stk::util::nested_iterator< std::vector<std::vector<int> >,
                                      std::vector<int>,
                                      stk::util::details::identity<std::vector<int> >
                                    > nested_iterator;

  typedef stk::util::nested_iterator< const std::vector<std::vector<int> >,
                                      std::vector<int>,
                                      stk::util::details::identity<std::vector<int> >
                                    > const_nested_iterator;
}

STKUNIT_UNIT_TEST ( nested_iterator, vector_vector_int)
{
  using namespace vector_vector_int;

  const int OUTER = 10;
  const int INNER = 10;

  nested_type a(OUTER);

  {
    int count = 0;
    for (int i=0; i<OUTER; ++i) {
      a[i].resize(INNER);
      for (int j=0; j<INNER; ++j) {
        a[i][j] = ++count;
      }
    }
  }

  nested_iterator itr(a);
  const nested_iterator end;

  STKUNIT_EXPECT_EQUAL( OUTER*INNER, std::distance(itr,end));

  {
    int count = 0;
    for(; itr != end; ++itr) {
      *itr *= 3;
      ++count;
      STKUNIT_EXPECT_EQUAL(count*3,*itr);
    }
  }

}

STKUNIT_UNIT_TEST ( nested_iterator, vector_vector_int_nonconst_to_const)
{
  using namespace vector_vector_int;

  const int OUTER = 10;
  const int INNER = 10;

  nested_type a(OUTER);

  {
    int count = 0;
    for (int i=0; i<OUTER; ++i) {
      a[i].resize(INNER);
      for (int j=0; j<INNER; ++j) {
        a[i][j] = ++count;
      }
    }
  }

  nested_iterator itr(a);

  const_nested_iterator const_itr = itr;

  const_nested_iterator end;

  STKUNIT_EXPECT_EQUAL( OUTER*INNER, std::distance(const_itr,end));

  {
    int count = 0;
    for(; const_itr != end; ++const_itr) {
      ++count;
      STKUNIT_EXPECT_EQUAL(count,*const_itr);
    }
  }

}

namespace list_vector_int {

  typedef std::list< std::vector<int> > nested_type;

  typedef stk::util::nested_iterator< std::list<std::vector<int> >,
                                      std::vector<int>,
                                      stk::util::details::identity<std::vector<int> >
                                    > nested_iterator;
}

STKUNIT_UNIT_TEST ( nested_iterator, list_vector_int)
{
  using namespace list_vector_int;

  const int OUTER = 10;
  const int INNER = 10;

  nested_type a(OUTER);

  {
    int count = 0;
    for (int i=0; i<OUTER; ++i) {
      std::vector<int> tmp;
      a.push_back(tmp);
      for (int j=0; j<INNER; ++j) {
        a.back().push_back(++count);
      }
    }
  }

  nested_iterator itr(a);
  const nested_iterator end;

  STKUNIT_EXPECT_EQUAL( OUTER*INNER, std::distance(itr,end));

  {
    int count = 0;
    for(; itr != end; ++itr) {
      *itr *= 3;
      ++count;
      STKUNIT_EXPECT_EQUAL(count*3,*itr);
    }
  }

}

