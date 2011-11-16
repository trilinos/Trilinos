#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/util/nested_iterator.hpp>

#include <vector>
#include <list>
#include <iostream>

namespace vector_vector_int {

  typedef std::vector< std::vector<int> > nested_type;

  struct to_inner_range {
    typedef std::vector<int> & value_type;
    typedef boost::iterator_range< std::vector<int>::iterator> result_type;

    result_type operator()(value_type value) const {
      return result_type(boost::begin(value),boost::end(value));
    }
  };

  typedef stk::util::nested_iterator< std::vector<std::vector<int> >::iterator,
                                      std::vector<int>::iterator,
                                      to_inner_range
                                    > nested_iterator;
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

  nested_iterator itr(a,to_inner_range());
  const nested_iterator end;

  STKUNIT_EXPECT_EQUAL( OUTER*INNER, std::distance(itr,end));

  {
    int count = 0;
    for(; itr != end; ++itr) {
      *itr *= 3;
      STKUNIT_EXPECT_EQUAL(++count*3,*itr);
    }
  }

}

namespace list_vector_int {

  typedef std::list< std::vector<int> > nested_type;

  struct to_inner_range {
    typedef std::vector<int> & value_type;
    typedef boost::iterator_range< std::vector<int>::iterator> result_type;

    result_type operator()(value_type value) const {
      return result_type(boost::begin(value),boost::end(value));
    }
  };

  typedef stk::util::nested_iterator< std::list<std::vector<int> >::iterator,
                                      std::vector<int>::iterator,
                                      to_inner_range
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

  nested_iterator itr(a,to_inner_range());
  const nested_iterator end;

  STKUNIT_EXPECT_EQUAL( OUTER*INNER, std::distance(itr,end));

  {
    int count = 0;
    for(; itr != end; ++itr) {
      *itr *= 3;
      STKUNIT_EXPECT_EQUAL(++count*3,*itr);
    }
  }

}
