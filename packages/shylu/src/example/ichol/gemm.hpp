#pragma once
#ifndef __GEMM_HPP__
#define __GEMM_HPP__

namespace Example { 

  using namespace std;
  
  template<typename CrsMatViewType>
  inline int
  gemm_nt_t(const typename CrsMatViewType::value_type alpha,
            const CrsMatViewType A,
            const CrsMatViewType X,
            const typename CrsMatViewType::value_type beta,
            const CrsMatViewType Y) {
    typedef typename CrsMatViewType::ordinal_type ordinal_type;
    typedef typename CrsMatViewType::value_type   value_type;
    
    // case that X.transpose, A.no_transpose, Y.no_transpose

    for (ordinal_type j=0;j<X.NumRows();++j) {
      auto x = X.extractRow(j);
      if (x.NumNonZeros()) {
        for (ordinal_type i=0;i<Y.NumRows();++i) {
          auto y = Y.extractRow(i);
          auto a = A.extractRow(i);

          if (y.NumNonZeros() && a.NumNonZeros()) {
            ordinal_type id = y.Index(j);
            if (id >= 0) {
              value_type &upsilon = y.Value(id);
              upsilon = beta*upsilon + alpha*dot(a, x);
            }
          }
        }
      }
    } 
    
    return 0;
  }

}

#endif
