#pragma once
#ifndef __GEMV_HPP__
#define __GEMV_HPP__

namespace Example { 

  using namespace std;
  
  template<typename CrsMatViewType>
  inline int
  gemv_nt_t(const typename CrsMatViewType::value_type alpha,
            const CrsMatViewType A,
            const CrsMatViewType x,
            const typename CrsMatViewType::value_type beta,
            const CrsMatViewType y) {
    typedef typename CrsMatViewType::ordinal_type ordinal_type;
    typedef typename CrsMatViewType::value_type   value_type;

    // case that x is x.transpose, A.no_transpose, y.no_transpose

    auto xx = x.extractRow(0);
    if (xx.NumNonZeros()) {
      for (ordinal_type i=0;i<y.NumRows();++i) {
        auto yy = y.extractRow(i);
        auto aa = A.extractRow(i);

        // grep the scalar located at index 0 in the row
        if (yy.NumNonZeros() && aa.NumNonZeros()) {
          ordinal_type id = yy.Index(0);
          if (id >= 0) {
            value_type &upsilon = yy.Value(id);
            upsilon = beta*upsilon + alpha*dot(aa, xx);
          }
        }
      }
    } 

    return 0;
  }

}

#endif
