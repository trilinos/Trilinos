
#include <stddef.h>                     // for size_t
#include <boost/iterator/filter_iterator.hpp>  // for filter_iterator
#include <iostream>                     // for basic_ostream::operator<<, etc
#include <gtest/gtest.h>
#include <vector>                       // for vector, vector<>::iterator
#include "boost/iterator/iterator_facade.hpp"  // for iterator_facade, etc


namespace vector_vector_int {

struct vec_filter {
 vec_filter(size_t sz=0) : m_sz(sz) {}

 bool operator()(const std::vector<int>& vec) const { return vec.size() > m_sz; }

 private:
  size_t m_sz;
};

typedef std::vector< std::vector<int> > nested_type;

}

TEST ( filter_iterator, vector_vector_int)
{
  using namespace vector_vector_int;

  const int OUTER = 10;

  nested_type a(OUTER);

  {
    int count = 0;
    for (int i=0; i<OUTER; ++i) {
      a[i].resize(i+1);
      for (int j=0; j<i+1; ++j) {
        a[i][j] = ++count;
      }
    }
  }

  boost::filter_iterator<vec_filter,std::vector<std::vector<int> >::iterator> itr(vec_filter(4),a.begin(),a.end());
  boost::filter_iterator<vec_filter,std::vector<std::vector<int> >::iterator> itr_end(a.end(), a.end());

  for(; itr!=itr_end; ++itr)
  {
    std::vector<int>& v = *itr;
    for(std::vector<int>::iterator vit=v.begin(),vend=v.end(); vit!=vend; ++vit) {
      std::cout<<*vit<<" ";
    }
    std::cout<<std::endl;
  }

}
