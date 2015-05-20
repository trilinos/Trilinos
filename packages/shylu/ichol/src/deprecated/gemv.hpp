#pragma once
#ifndef __GEMV_HPP__
#define __GEMV_HPP__

/// \file gemv.hpp
/// \brief Sparse matrix-vector multiplication on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;
  
  template<typename ScalarType, 
           typename CrsMatViewType>
  KOKKOS_INLINE_FUNCTION 
  int
  gemv_nt_t(const ScalarType alpha,
            const CrsMatViewType A,
            const CrsMatViewType x,
            const ScalarType beta,
            const CrsMatViewType y) {
    typedef typename CrsMatViewType::ordinal_type  ordinal_type;
    typedef typename CrsMatViewType::value_type    value_type;
    typedef typename CrsMatViewType::row_view_type row_view_type;    

    // case that x is x.transpose, A.no_transpose, y.no_transpose

    row_view_type xx = x.extractRow(0);
    if (xx.NumNonZeros()) {
      for (ordinal_type i=0;i<y.NumRows();++i) {
        row_view_type yy = y.extractRow(i);
        row_view_type aa = A.extractRow(i);

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
