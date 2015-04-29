#pragma once
#ifndef __TRSM_L_U_NT_HPP__
#define __TRSM_L_U_NT_HPP__

/// \file trsm_l_u_nt.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///

namespace Example {

  using namespace std;

  template<>
  template<typename ParallelForType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename DenseExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,
       AlgoTrsm::ForTriSolveBlocked>
  ::invoke(const typename CrsExecViewTypeA::policy_type::member_type &member,
           const int diag,
           const ScalarType alpha,
           const CrsExecViewTypeA &A,
           const DenseExecViewTypeB &B) {
    typedef typename CrsExecViewTypeA::ordinal_type      ordinal_type;
    typedef typename CrsExecViewTypeA::value_type        value_type;
    typedef typename CrsExecViewTypeA::row_view_type     row_view_type;
    typedef typename CrsExecViewTypeA::team_factory_type team_factory_type;

    scaleDenseMatrix<ParallelForType,ScalarType,DenseExecViewTypeB>(member, alpha, B);

    for (ordinal_type k=A.NumRows()-1;k>=0;--k) {
      row_view_type &a = A.RowView(k);
      const value_type diag = a.Value(0);
      
      const ordinal_type nnz_a = a.NumNonZeros();
      const ordinal_type nB = B.NumCols();

      // update
      ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, nB),
                      [&](const ordinal_type j) {
                        for (ordinal_type ja=1;ja<nnz_a;++ja) {
                          ordinal_type col_at_ja = a.Col(ja);
                          value_type   val_at_ja = a.Value(ja);

                          B.Value(k, col_at_ja) -= val_at_ja*B.Value(col_at_ja, j);
                        }
                      });

      // invert
      ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, nB),
                      [&](const ordinal_type j) {
                        B.Value(k, j) /= diag;
                      });
    }
    
    return 0;
  }

}

#endif
