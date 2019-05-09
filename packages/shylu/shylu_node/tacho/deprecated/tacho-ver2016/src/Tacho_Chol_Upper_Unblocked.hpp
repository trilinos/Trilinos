#ifndef __TACHO_CHOL_UPPER_UNBLOCKED_HPP__
#define __TACHO_CHOL_UPPER_UNBLOCKED_HPP__

/// \file Tacho_Chol_Upper_Unblocked.hpp
/// \brief Unblocked incomplete Chloesky factorization.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  /// Native unblocked sparse Chol
  /// ============================
  template<>
  template<typename PolicyType,
           typename MemberType,
           typename CrsExecViewTypeA>
  KOKKOS_INLINE_FUNCTION
  int
  Chol<Uplo::Upper,
       AlgoChol::Unblocked,Variant::One>
  ::invoke(PolicyType &policy,
           MemberType &member,
           CrsExecViewTypeA &A) {

    typedef typename CrsExecViewTypeA::value_type    value_type;
    typedef typename CrsExecViewTypeA::ordinal_type  ordinal_type;
    typedef typename CrsExecViewTypeA::row_view_type row_view_type;

    // row_view_type r1t, r2t;

    for (ordinal_type k=0;k<A.NumRows();++k) {
      //r1t.setView(A, k);
      row_view_type &r1t = A.RowView(k);

      // extract diagonal from alpha11
      value_type &alpha = r1t.Value(0);

      if (member.team_rank() == 0) {
        // if encounter null diag or wrong index, return -(row + 1)
        TACHO_TEST_FOR_ABORT( r1t.Col(0) != k, "Chol::Unblocked:: Diagonal does not exist");        
        if (Util::real(alpha) <= 0.0) {
          // warning message
          fprintf(stderr, "   diagonal = %f, local col = %d, global col = %d\n", 
                  Util::real(alpha), k, r1t.OffsetCols() + k);
          // proceed with epsilon; for incomplete factorization, Cholesky factor may not exit
          alpha = 1.0e-8;

          //TACHO_TEST_FOR_ABORT( true, "Chol::Unblocked:: Diagonal is negative");
          //return -(k + 1);
        }

        // error handling should be more carefully designed

        // sqrt on diag
        alpha = sqrt(Util::real(alpha));
      }
      member.team_barrier();

      const ordinal_type nnz_r1t = r1t.NumNonZeros();

      if (nnz_r1t) {
        // inverse scale
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 1, nnz_r1t),
                             [&](const ordinal_type j) {
                               r1t.Value(j) /= alpha;
                             });

        member.team_barrier();

        // hermitian rank update
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 1, nnz_r1t),
                             [&](const ordinal_type i) {
                               const ordinal_type row_at_i = r1t.Col(i);
                               const value_type   val_at_i = Util::conj(r1t.Value(i));

                               //r2t.setView(A, row_at_i);
                               row_view_type &r2t = A.RowView(row_at_i);
                               ordinal_type idx = 0;

                               for (ordinal_type j=i;j<nnz_r1t && (idx > -2);++j) {
                                 const ordinal_type col_at_j = r1t.Col(j);
                                 idx = r2t.Index(col_at_j, idx);

                                 if (idx >= 0) {
                                   const value_type val_at_j = r1t.Value(j);
                                   r2t.Value(idx) -= val_at_i*val_at_j;
                                 }
                               }
                             });
      }
    }
    return 0;
  }

}

#endif
