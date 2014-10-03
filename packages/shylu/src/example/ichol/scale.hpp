#pragma once
#ifndef __SCALE_HPP__
#define __SCALE_HPP__

namespace Example { 

  using namespace std;
  
  template<typename CrsMatrixView>
  inline int
  scale(const typename CrsMatrixView::value_type alpha,
        CrsMatrixView &A) {
    typedef typename CrsMatrixView::ordinal_type ordinal_type;

    for (ordinal_type i=0;i<A.NumRows();++i) {
      auto row = A.extractRow(i);
      for (ordinal_type j=0;j<row.NumNonZeros();++j)
        row.Value(j) *= alpha;
    }

    return 0;
  }

}

#endif
