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
           typename CrsMatViewType,
           typename ParallelForType>
  KOKKOS_INLINE_FUNCTION 
  int
  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
       AlgoTrsm::ForRightBlocked>
  ::invoke(const ParallelForType::member_type member,
           const int diag,
           const ScalarType alpha,
           const CrsMatViewType A,
           const CrsMatViewType B) {
    typedef typename CrsMatViewType::ordinal_type  ordinal_type;
    typedef typename CrsMatViewType::value_type    value_type;
    typedef typename CrsMatViewType::row_view_type row_view_type;

    //row_view_type a, b1, b2;

    for (ordinal_type k=0;k<A.NumRows();++k) {
      // pick a diag
      //a.setView(A, k);
      row_view_type &a = A.RowView(k);
      const value_type diag = a.Value(0);

      // invert
      //b1.setView(B, k);
      row_view_type &b1 = B.RowView(k);

      const ordinal_type nnz_b1 = b1.NumNonZeros();

      ParallelFor(ParallelForType::TeamThreadLoop(member, 0, nnz_b1),
                  [&](const ordinal_type j) {
                    b1.Value(j) /= diag;
                  });
      
      
      // update 
      const ordinal_type nnz_a = a.NumNonZeros();
      ParallelFor(ParallelForType::TeamThreadLoop(member, 1, nnz_a),
                  [&](const ordinal_type i) {
                    const ordinal_type row_at_i = a.Col(i);
                    const value_type   val_at_i = conj(a.Value(i));
                    
                    //b2.setView(B, row_at_i);
                    row_view_type &b2 = B.RowView(row_at_i);
                    
                    ordinal_type idx = 0;
                    for (ordinal_type j=0;j<nnz_b1 && (idx > -2);++j) {
                      ordinal_type col_at_j = b1.Col(j);
                      value_type   val_at_j = b1.Value(j);
                      
                      idx = b2.Index(col_at_j, idx);
                      if (idx >= 0) 
                        b2.Value(idx) += alpha*val_at_i*val_at_j;
                    }
                  });
    }

    return 0;
  }

}

#endif
