#pragma once
#ifndef __GEMM_CT_NT_HPP__
#define __GEMM_CT_NT_HPP__

/// \file gemm_t_nt.hpp
/// \brief Sparse matrix-matrix multiplication on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example {

  using namespace std;

  // 
  template<>
  template<typename ParallelForType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename CrsExecViewTypeB,
           typename CrsExecViewTypeC>
  KOKKOS_INLINE_FUNCTION
  int
  Gemm<Trans::ConjTranspose,Trans::NoTranspose,
       AlgoGemm::ForFactorRightBlocked>
  ::invoke(const typename CrsExecViewTypeA::policy_type::member_type &member,
           const ScalarType alpha,
           const CrsExecViewTypeA &A,
           const CrsExecViewTypeB &B,
           const ScalarType beta,
           const CrsExecViewTypeC &C) {
    typedef typename CrsExecViewTypeA::ordinal_type      ordinal_type;
    typedef typename CrsExecViewTypeA::value_type        value_type;
    typedef typename CrsExecViewTypeA::row_view_type     row_view_type;
    typedef typename CrsExecViewTypeA::team_factory_type team_factory_type;

    scaleCrsMatrix<ParallelForType,ScalarType,CrsExecViewTypeC>(member, beta, C);

    //row_view_type a, b, c;
    for (ordinal_type k=0;k<A.NumRows();++k) {
      //a.setView(A, k);
      row_view_type &a = A.RowView(k);
      const ordinal_type nnz_a = a.NumNonZeros();

      //b.setView(B, k);
      row_view_type &b = B.RowView(k);
      const ordinal_type nnz_b = b.NumNonZeros();

      ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, nnz_a),
                      [&](const ordinal_type i) {
                        const ordinal_type row_at_i = a.Col(i);
                        const value_type   val_at_i = conj(a.Value(i));

                        //c.setView(C, row_at_i);
                        row_view_type &c = C.RowView(row_at_i);

                        ordinal_type idx = 0;
                        for (ordinal_type j=0;j<nnz_b && (idx > -2);++j) {
                          ordinal_type col_at_j = b.Col(j);
                          value_type   val_at_j = b.Value(j);

                          idx = c.Index(col_at_j, idx);
                          if (idx >= 0)
                            c.Value(idx) += alpha*val_at_i*val_at_j;
                        }
                      });
    }

    return 0;
  }


  // B and C are dense matrices
  template<>
  template<typename ParallelForType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename DenseExecViewTypeB,
           typename DenseExecViewTypeC>
  KOKKOS_INLINE_FUNCTION
  int
  Gemm<Trans::ConjTranspose,Trans::NoTranspose,
       AlgoGemm::ForSolveDenseMatrix>
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

    for (ordinal_type k=0;k<A.NumRows();++k) {
      row_view_type &a = A.RowView(k);
      const ordinal_type nnz_a = a.NumNonZeros();

      ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, nnz_a),
                      [&](const ordinal_type i) {
                        const ordinal_type row_at_i = a.Col(i);
                        const value_type   val_at_i = conj(a.Value(i));
                        
                        for (ordinal_type j=0;j<B.NumCols();++j) {
                          const value_type val_at_j = B.Value(k, j);

                          C.Value(row_at_i, j) += alpha*val_at_i*val_at_j;
                        }
                      });
    }

    return 0;
  }

}

#endif
