#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/util/nested_range.hpp>
#include <boost/foreach.hpp>

#include <vector>
#include <list>
#include <iostream>

namespace vector_vector_int {

typedef std::vector< std::vector<int> > nested_type;
typedef std::vector< std::pair<int*,int*> > nested_pair_type;

typedef stk_classic::util::nested_range< nested_type > range;
typedef stk_classic::util::nested_range< nested_pair_type > pair_range;


}

STKUNIT_UNIT_TEST ( nested_range, basic)
{
  using namespace vector_vector_int;

  const int OUTER = 10;
  const int INNER = 10;

  nested_type a(OUTER);
  nested_pair_type ap(OUTER);

  {
    int count = 0;
    for (int i=0; i<OUTER; ++i) {
      a[i].resize(INNER);
      int* int_ptr = NULL;
      ap[i] = std::make_pair(int_ptr,int_ptr);
      for (int j=0; j<INNER; ++j) {
        a[i][j] = ++count;
      }
    }
  }

  range rng(a);
rng.begin();
//  BOOST_FOREACH(int i, rng) {
//   std::cout<<i<<std::endl;
//  }

  pair_range prng(ap);
prng.begin();
//  BOOST_FOREACH(int i, prng) {
//   std::cout<<i<<std::endl;
//  }

}

//namespace vector_vector_int {
//
//struct to_inner_range {
//  typedef std::vector<int> result_type;
//  typedef std::vector<int>* value_type;
//
//  result_type& operator()(value_type& r) const { return *r; }
//  const result_type& operator()(const value_type& r) const { return *r; }
//};
//
//typedef std::vector< std::vector<int>* > nested_ptr_type;
//
//typedef stk_classic::util::nested_range< nested_ptr_type, std::vector<int>, to_inner_range > ptr_range;
//
//}
//
//STKUNIT_UNIT_TEST ( nested_range, nested_ptr)
//{
//  using namespace vector_vector_int;
//
//  const int OUTER = 10;
//  const int INNER = 10;
//
//  nested_ptr_type a(OUTER);
//
//  {
//    int count = 0;
//    for (int i=0; i<OUTER; ++i) {
//      a[i] = new std::vector<int>(INNER);
//      for (int j=0; j<INNER; ++j) {
//        (*a[i])[j] = ++count;
//      }
//    }
//  }
//
//  BOOST_FOREACH(std::vector<int>* vecptr,a) {
//    std::vector<int>& vec = *vecptr;
//    BOOST_FOREACH(int i,vec) {
//      std::cout<<i<<std::endl;
//    }
//  }
//
//  ptr_range rng(a,to_inner_range());
//  BOOST_FOREACH(int i, rng) {
//   std::cout<<i<<std::endl;
//  }
//
//}

