#pragma once
#ifndef __HERK_U_T_HPP__
#define __HERK_U_T_HPP__

/// \file herk_u_t.hpp
/// \brief Sparse hermitian rank one update on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;
  
  template<>
  template<typename ScalarType, 
           typename CrsMatViewType>
  KOKKOS_INLINE_FUNCTION 
  int
  Herk<Uplo::Upper,Trans::ConjTranspose,
       AlgoHerk::ForRightBlocked>
  ::invoke(const ScalarType alpha,
           const CrsMatViewType A,
           const ScalarType beta,
           const CrsMatViewType C) {
    typedef typename CrsMatViewType::ordinal_type  ordinal_type;
    typedef typename CrsMatViewType::value_type    value_type;
    typedef typename CrsMatViewType::row_view_type row_view_type;
    typedef DotTraits<value_type> dot_traits;

    scale(beta, C);

    row_view_type a, c;
    for (ordinal_type k=0;k<A.NumRows();++k) {
      a.setView(A, k);
      const ordinal_type nnz = a.NumNonZeros();

      for (ordinal_type i=0;i<nnz;++i) {
        const ordinal_type row_at_i = a.Col(i);
        const value_type   val_at_i = conj(a.Value(i));
        
        c.setView(C, row_at_i);
        ordinal_type prev = 0;
        
        for (ordinal_type j=0;j<nnz;++j) {
          ordinal_type col_at_j = a.Col(j);
          value_type   val_at_j = a.Value(j);
          
          ordinal_type idx = c.Index(col_at_j, prev);
          if (idx >= 0) {
            c.Value(idx) += alpha*val_at_i*val_at_j;
            prev = idx;
          }
        }
      }
    }

    return 0;
  }

}

#endif
