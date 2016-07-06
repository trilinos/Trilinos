#pragma once
#ifndef __GEMM_NT_T_HPP__
#define __GEMM_NT_T_HPP__

/// \file gemm_nt_t.hpp
/// \brief Sparse matrix-matrix multiplication on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;
  
  template<>
  template<typename ScalarType, 
           typename CrsMatViewType>
  KOKKOS_INLINE_FUNCTION 
  int
  Gemm<Trans::NoTranspose,Trans::ConjTranspose,
       AlgoGemm::ForLeftBlocked>
  ::invoke(const ScalarType alpha,
           const CrsMatViewType A,
           const CrsMatViewType B,
           const ScalarType beta,
           const CrsMatViewType C) {
    typedef typename CrsMatViewType::ordinal_type  ordinal_type;
    typedef typename CrsMatViewType::value_type    value_type;
    typedef typename CrsMatViewType::row_view_type row_view_type;

    scale(beta, C);

    row_view_type a, b, c;
    for (ordinal_type j=0;j<B.NumRows();++j) {
      b.setView(B, j);
      if (b.NumNonZeros()) {
        for (ordinal_type i=0;i<C.NumRows();++i) {
          c.setView(C, i);
          a.setView(A, i);

          if (c.NumNonZeros() && a.NumNonZeros()) {
            ordinal_type id = c.Index(j);
            if (id >= 0) 
              c.Value(id) += alpha*dot(a, b);
          }
        }
      }
    } 
    
    return 0;
  }
  
}

#endif
