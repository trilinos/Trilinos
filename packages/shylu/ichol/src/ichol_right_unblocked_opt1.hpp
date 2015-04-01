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
  template<typename ParallelForType,
           typename CrsExecViewType>
  KOKKOS_INLINE_FUNCTION
  int
  IChol<Uplo::Upper,AlgoIChol::RightUnblockedOpt1>
  ::invoke(const typename CrsExecViewType::policy_type::member_type &member,
           const CrsExecViewType &A) {

    typedef typename CrsExecViewType::value_type        value_type;
    typedef typename CrsExecViewType::ordinal_type      ordinal_type;
    typedef typename CrsExecViewType::row_view_type     row_view_type;
    typedef typename CrsExecViewType::team_factory_type team_factory_type;

    // row_view_type r1t, r2t;

    for (ordinal_type k=0;k<A.NumRows();++k) {
      //r1t.setView(A, k);
      row_view_type &r1t = A.RowView(k);

      // extract diagonal from alpha11
      value_type &alpha = r1t.Value(0);

      if (member.team_rank() == 0) {
        // if encounter null diag or wrong index, return -(row + 1)
        if (abs(alpha) == 0.0 || r1t.Col(0) != k)
          return -(k + 1);

        // error handling should be more carefully designed

        // sqrt on diag
        alpha = sqrt(real(alpha));
      }
      member.team_barrier();

      const ordinal_type nnz_r1t = r1t.NumNonZeros();

      if (nnz_r1t) {
        // inverse scale
        ParallelForType(team_factory_type::createThreadLoopRegion(member, 1, nnz_r1t),
                        [&](const ordinal_type j) {
                          r1t.Value(j) /= alpha;
                        });

        member.team_barrier();

        // hermitian rank update
        ParallelForType(team_factory_type::createThreadLoopRegion(member, 1, nnz_r1t),
                        [&](const ordinal_type i) {
                          const ordinal_type row_at_i = r1t.Col(i);
                          const value_type   val_at_i = conj(r1t.Value(i));

                          //r2t.setView(A, row_at_i);
                          row_view_type &r2t = A.RowView(row_at_i);
                          ordinal_type idx = 0;

                          for (ordinal_type j=i;j<nnz_r1t && (idx > -2);++j) {
                            const ordinal_type col_at_j = r1t.Col(j);
                            const value_type   val_at_j = r1t.Value(j);

                            idx = r2t.Index(col_at_j, idx);
                            if (idx >= 0)
                              r2t.Value(idx) -= val_at_i*val_at_j;
                          }
                        });
      }
    }
    return 0;
  }

}

#endif
