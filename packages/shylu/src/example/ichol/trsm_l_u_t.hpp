#pragma once
#ifndef __TRSM_L_U_T_HPP__
#define __TRSM_L_U_T_HPP__

/// \file trsm_l_u_t.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///

namespace Example { 

  using namespace std;
  
  template<>
  template<typename ScalarType,
           typename CrsMatViewType>
  KOKKOS_INLINE_FUNCTION 
  int
  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
       AlgoTrsm::ForRightBlocked>
  ::invoke(const int diag,
           const ScalarType alpha,
           const CrsMatViewType A,
           const CrsMatViewType B) {
    typedef typename CrsMatViewType::ordinal_type  ordinal_type;
    typedef typename CrsMatViewType::value_type    value_type;
    typedef typename CrsMatViewType::row_view_type row_view_type;

    row_view_type a, b1, b2;

    for (ordinal_type k=0;k<A.NumRows();++k) {
      // pick a diag
      a.setView(A, k);
      const value_type diag = a.Value(0);

      // invert
      b1.setView(B, k);

      const ordinal_type nnz_b1 = b1.NumNonZeros();
      for (ordinal_type j=0;j<nnz_b1;++j) 
        b1.Value(j) /= diag;
      
      // update 
      const ordinal_type nnz_a = a.NumNonZeros();
      for (ordinal_type i=1;i<nnz_a;++i) {
        const ordinal_type row_at_i = a.Col(i);
        const value_type   val_at_i = conj(a.Value(i));
          
        b2.setView(B, row_at_i);
        ordinal_type prev = 0;

        for (ordinal_type j=0;j<nnz_b1;++j) {
          ordinal_type col_at_j = b1.Col(j);
          value_type   val_at_j = b1.Value(j);

          ordinal_type idx = b2.Index(col_at_j, prev);
          if (idx >= 0) {
            b2.Value(idx) += alpha*val_at_i*val_at_j;
            prev = idx;
          }
        }
      }
    }

    return 0;
  }

}

#endif
