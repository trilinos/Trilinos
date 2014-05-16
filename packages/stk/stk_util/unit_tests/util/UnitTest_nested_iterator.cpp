
#include <iterator>                     // for distance
#include <list>                         // for list
#include <gtest/gtest.h>
#include <stk_util/util/nested_iterator.hpp>  // for nested_iterator
#include <stk_util/util/nested_range.hpp>  // for identity
#include <vector>                       // for vector
#include "boost/iterator/iterator_facade.hpp"  // for iterator_facade, etc
#include "boost/optional/optional.hpp"  // for operator==, operator!=


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
/// srk 12/20/12 - these tests seem to hang on boost 1.50 / Trilinos build
#if defined(STK_BUILT_IN_SIERRA)
TEST ( nested_iterator, vector_vector_int)
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

  EXPECT_EQ( OUTER*INNER, std::distance(itr,end));

  {
    int count = 0;
    for(; itr != end; ++itr) {
      *itr *= 3;
      EXPECT_EQ(++count*3,*itr);
    }
  }

}

TEST ( nested_iterator, vector_vector_int_nonconst_to_const)
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

  EXPECT_EQ( OUTER*INNER, std::distance(const_itr,end));

  {
    int count = 0;
    for(; const_itr != end; ++const_itr) {
      EXPECT_EQ(++count,*const_itr);
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

TEST ( nested_iterator, list_vector_int)
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

  EXPECT_EQ( OUTER*INNER, std::distance(itr,end));

  {
    int count = 0;
    for(; itr != end; ++itr) {
      *itr *= 3;
      EXPECT_EQ(++count*3,*itr);
    }
  }

}

#endif
