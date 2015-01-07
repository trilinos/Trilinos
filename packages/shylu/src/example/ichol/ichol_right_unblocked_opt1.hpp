#pragma once
#ifndef __ICHOL_RIGHT_UNBLOCKED_OPT1_HPP__
#define __ICHOL_RIGHT_UNBLOCKED_OPT1_HPP__

/// \file ichol_right_unblocked_opt1.hpp
/// \brief Unblocked incomplete Chloesky factorization.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "util.hpp"
#include "partition.hpp"

namespace Example { 

  using namespace std;

  template<>
  template<typename CrsMatViewType>
  KOKKOS_INLINE_FUNCTION
  int 
  IChol<Uplo::Upper,AlgoIChol::RightUnblockedOpt1>::invoke(const CrsMatViewType A) {
    
    typedef typename CrsMatViewType::value_type    value_type;
    typedef typename CrsMatViewType::ordinal_type  ordinal_type;
    typedef typename CrsMatViewType::row_view_type row_view_type;

    row_view_type r1t, r2t;

    for (ordinal_type i=0;i<A.NumRows();++i) {
      r1t.setView(A, i);

      // extract diagonal from alpha11
      value_type &alpha = r1t.Value(0);

      // if encounter null diag or wrong index, return -(row + 1)
      if (abs(alpha) == 0.0 || r1t.Col(0) != i)
        return -(i + 1);
      
      // sqrt on diag
      alpha = sqrt(real(alpha));

      // inverse scale 
      ordinal_type nnz_r1t = r1t.NumNonZeros();
      for (ordinal_type j=1;j<nnz_r1t;++j) 
        r1t.Value(j) /= alpha;
      
      // hermitian rank update
      for (ordinal_type i=1;i<r1t.NumNonZeros();++i) {
        const ordinal_type row_at_i = r1t.Col(i);
        const value_type   val_at_i = r1t.Value(i);

        r2t.setView(A, row_at_i);
        ordinal_type prev = 0;
        
        for (ordinal_type j=1;j<nnz_r1t;++j) {
          const ordinal_type col_at_j = r1t.Col(j);
          const value_type   val_at_j = r1t.Value(j);

          const ordinal_type idx = r2t.Index(col_at_j, prev);
          if (idx >= 0) {
            r2t.Value(idx) -= val_at_i*conj(val_at_j);
            prev = idx;
          }
        }
      }
    }
    
    return 0;
  }

}

#endif
