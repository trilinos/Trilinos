#pragma once
#ifndef __GEMM_NT_NT_HPP__
#define __GEMM_NT_NT_HPP__

/// \file gemm_nt_nt.hpp
/// \brief Sparse matrix-matrix multiplication on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example {

  using namespace std;

  // B and C are dense matrices
  template<>
  template<typename ParallelForType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename DenseExecViewTypeB,
           typename DenseExecViewTypeC>
  KOKKOS_INLINE_FUNCTION
  int
  Gemm<Trans::NoTranspose,Trans::NoTranspose,
       AlgoGemm::ForTriSolveBlocked>
  ::invoke(const typename CrsExecViewTypeA::policy_type::member_type &member,
           const ScalarType alpha,
           const CrsExecViewTypeA &A,
           const DenseExecViewTypeB &B,
           const ScalarType beta,
           const DenseExecViewTypeC &C) {
    typedef typename CrsExecViewTypeA::ordinal_type      ordinal_type;
    typedef typename CrsExecViewTypeA::value_type        value_type;
    typedef typename CrsExecViewTypeA::row_view_type     row_view_type;
    typedef typename CrsExecViewTypeA::team_factory_type team_factory_type;

    scaleDenseMatrix<ParallelForType,ScalarType,DenseExecViewTypeC>(member, beta, C);
    
    for (ordinal_type i=0;i<A.NumRows();++i) {
      row_view_type &a = A.RowView(i);
      const ordinal_type nnz_a = a.NumNonZeros();

      ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, nnz_a),
                      [&](const ordinal_type k) {
                        const ordinal_type col_at_k = a.Col(k);
                        const value_type   val_at_k = a.Value(k);
                        
                        for (ordinal_type j=0;j<B.NumCols();++j) {
                          const value_type val_at_j = B.Value(k, j);
                          
                          C.Value(i, j) += alpha*val_at_k*val_at_j;
                        }
                      });
    }

    return 0;
  }

}

#endif
