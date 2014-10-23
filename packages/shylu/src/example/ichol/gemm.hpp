#pragma once
#ifndef __GEMM_HPP__
#define __GEMM_HPP__

/// \file gemm.hpp
/// \brief Sparse matrix-matrix multiplication on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;
  
  template<typename ScalarType, 
           typename CrsMatViewType>
  KOKKOS_INLINE_FUNCTION 
  int
  gemm_nt_t(const ScalarType alpha,
            const CrsMatViewType A,
            const CrsMatViewType X,
            const ScalarType beta,
            const CrsMatViewType Y) {
    typedef typename CrsMatViewType::ordinal_type  ordinal_type;
    typedef typename CrsMatViewType::value_type    value_type;
    typedef typename CrsMatViewType::row_view_type row_view_type;

    // case that X.transpose, A.no_transpose, Y.no_transpose

    for (ordinal_type j=0;j<X.NumRows();++j) {
      row_view_type x = X.extractRow(j);
      if (x.NumNonZeros()) {
        for (ordinal_type i=0;i<Y.NumRows();++i) {
          row_view_type y = Y.extractRow(i);
          row_view_type a = A.extractRow(i);

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
