#pragma once
#ifndef __TRSM_L_U_CT_HPP__
#define __TRSM_L_U_CT_HPP__

/// \file trsm_l_u_t.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///

namespace Example {

  using namespace std;

  template<>
  template<typename ParallelForType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename CrsExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
       AlgoTrsm::ForFactorBlocked>
  ::invoke(const typename CrsExecViewTypeA::policy_type::member_type &member,
           const int diag,
           const ScalarType alpha,
           const CrsExecViewTypeA &A,
           const CrsExecViewTypeB &B) {
    typedef typename CrsExecViewTypeA::ordinal_type      ordinal_type;
    typedef typename CrsExecViewTypeA::value_type        value_type;
    typedef typename CrsExecViewTypeA::row_view_type     row_view_type;
    typedef typename CrsExecViewTypeA::team_factory_type team_factory_type;

    scaleCrsMatrix<ParallelForType,ScalarType,CrsExecViewTypeB>(member, alpha, B);

    for (ordinal_type k=0;k<A.NumRows();++k) {
      // pick a diag
      row_view_type &a = A.RowView(k);
      const value_type diag = a.Value(0);

      // invert
      row_view_type &b1 = B.RowView(k);

      const ordinal_type nnz_b1 = b1.NumNonZeros();

      ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, nnz_b1),
                      [&](const ordinal_type j) {
                        b1.Value(j) /= diag;
                      });
      

      // update
      const ordinal_type nnz_a = a.NumNonZeros();
      ParallelForType(team_factory_type::createThreadLoopRegion(member, 1, nnz_a),
                      [&](const ordinal_type i) {
                        const ordinal_type row_at_i = a.Col(i);
                        const value_type   val_at_i = conj(a.Value(i));
                        
                        row_view_type &b2 = B.RowView(row_at_i);
                        
                        ordinal_type idx = 0;
                        for (ordinal_type j=0;j<nnz_b1 && (idx > -2);++j) {
                          ordinal_type col_at_j = b1.Col(j);
                          value_type   val_at_j = b1.Value(j);
                          
                          idx = b2.Index(col_at_j, idx);
                          if (idx >= 0)
                            b2.Value(idx) -= val_at_i*val_at_j;
                        }
                      });
    }
    
    return 0;
  }

}

#endif
