#pragma once
#ifndef __TRSM_L_U_CT_HPP__
#define __TRSM_L_U_CT_HPP__

/// \file trsm_l_u_ct.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///

namespace Example {

  using namespace std;

  // Trsm used in the factorization phase: data parallel on b1t
  // ==========================================================
  template<>
  template<typename ParallelForType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename CrsExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
       AlgoTrsm::ForFactorBlockedVar1>
  ::invoke(const typename CrsExecViewTypeA::policy_type::member_type &member,
           const int diagA,
           const ScalarType alpha,
           const CrsExecViewTypeA &A,
           const CrsExecViewTypeB &B) {
    typedef typename CrsExecViewTypeA::ordinal_type      ordinal_type;
    typedef typename CrsExecViewTypeA::value_type        value_type;
    typedef typename CrsExecViewTypeA::row_view_type     row_view_type;
    typedef typename CrsExecViewTypeA::team_factory_type team_factory_type;

    // scale the matrix B with alpha
    scaleCrsMatrix<ParallelForType,ScalarType,CrsExecViewTypeB>(member, alpha, B);

    // Solve a system: AX = B -> B := inv(A) B
    const ordinal_type mA = A.NumRows();
    const ordinal_type nB = B.NumCols();

    if (nB > 0) {
      for (ordinal_type k=0;k<mA;++k) {
        row_view_type &a = A.RowView(k);
        const value_type cdiag = conj(a.Value(0));

        // invert
        row_view_type &b1 = B.RowView(k);
        const ordinal_type nnz_b1 = b1.NumNonZeros();

        if (diagA != Diag::Unit && nnz_b1 > 0) {
          // b1t = b1t / conj(diag)
          ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, nnz_b1),
                          [&](const ordinal_type j) {
                            b1.Value(j) /= cdiag;
                          });
        }

        // update
        const ordinal_type nnz_a = a.NumNonZeros();
        if (nnz_a > 0) {
          // B2 = B2 - trans(conj(a12t)) b1t
          ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, nnz_b1),
                          [&](const ordinal_type j) {
                            // grab b1
                            const ordinal_type col_at_j = b1.Col(j);
                            const value_type   val_at_j = b1.Value(j);

                            for (ordinal_type i=1;i<nnz_a;++i) {
                              // grab a12t
                              const ordinal_type row_at_i = a.Col(i);
                              const value_type   val_at_i = conj(a.Value(i));

                              // grab b2
                              row_view_type &b2 = B.RowView(row_at_i);

                              // check and update
                              ordinal_type idx = 0;
                              idx = b2.Index(col_at_j, idx);
                              if (idx >= 0)
                                b2.Value(idx) -= val_at_i*val_at_j;
                            }
                          });
        }
        member.team_barrier();
      }
    }

    return 0;
  }

  // Trsm used in the factorization phase: data parallel on a1t
  // ==========================================================
  template<>
  template<typename ParallelForType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename CrsExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
       AlgoTrsm::ForFactorBlockedVar2>
  ::invoke(const typename CrsExecViewTypeA::policy_type::member_type &member,
           const int diagA,
           const ScalarType alpha,
           const CrsExecViewTypeA &A,
           const CrsExecViewTypeB &B) {
    typedef typename CrsExecViewTypeA::ordinal_type      ordinal_type;
    typedef typename CrsExecViewTypeA::value_type        value_type;
    typedef typename CrsExecViewTypeA::row_view_type     row_view_type;
    typedef typename CrsExecViewTypeA::team_factory_type team_factory_type;

    // scale the matrix B with alpha
    scaleCrsMatrix<ParallelForType,ScalarType,CrsExecViewTypeB>(member, alpha, B);

    // Solve a system: AX = B -> B := inv(A) B
    const ordinal_type mA = A.NumRows();
    const ordinal_type nB = B.NumCols();

    if (nB > 0) {
      for (ordinal_type k=0;k<mA;++k) {
        row_view_type &a = A.RowView(k);
        const value_type cdiag = conj(a.Value(0));

        // invert
        row_view_type &b1 = B.RowView(k);
        const ordinal_type nnz_b1 = b1.NumNonZeros();

        if (diagA != Diag::Unit && nnz_b1 > 0) {
          // b1t = b1t / conj(diag)
          ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, nnz_b1),
                          [&](const ordinal_type j) {
                            b1.Value(j) /= cdiag;
                          });
          member.team_barrier();
        }

        // update
        const ordinal_type nnz_a = a.NumNonZeros();
        if (nnz_a > 0) {
          // B2 = B2 - trans(conj(a12t)) b1t
          ParallelForType(team_factory_type::createThreadLoopRegion(member, 1, nnz_a),
                          [&](const ordinal_type i) {
                            // grab a12t
                            const ordinal_type row_at_i = a.Col(i);
                            const value_type   val_at_i = conj(a.Value(i));

                            // grab b2
                            row_view_type &b2 = B.RowView(row_at_i);

                            ordinal_type idx = 0;
                            for (ordinal_type j=0;j<nnz_b1 && (idx > -2);++j) {
                              // grab b1
                              const ordinal_type col_at_j = b1.Col(j);
                              const value_type   val_at_j = b1.Value(j);

                              // check and update
                              idx = b2.Index(col_at_j, idx);
                              if (idx >= 0)
                                b2.Value(idx) -= val_at_i*val_at_j;
                            }
                          });
          member.team_barrier();
        }
      }
    }

    return 0;
  }

  // Trsm used in the tri-solve phase: Multiple RHS
  // ==============================================
  template<>
  template<typename ParallelForType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename DenseExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
       AlgoTrsm::ForTriSolveBlockedVar1>
  ::invoke(const typename CrsExecViewTypeA::policy_type::member_type &member,
           const int diagA,
           const ScalarType alpha,
           const CrsExecViewTypeA &A,
           const DenseExecViewTypeB &B) {
    typedef typename CrsExecViewTypeA::ordinal_type      ordinal_type;
    typedef typename CrsExecViewTypeA::value_type        value_type;
    typedef typename CrsExecViewTypeA::row_view_type     row_view_type;
    typedef typename CrsExecViewTypeA::team_factory_type team_factory_type;

    // scale the matrix B with alpha
    scaleDenseMatrix<ParallelForType,ScalarType,DenseExecViewTypeB>(member, alpha, B);

    // Solve a system: AX = B -> B := inv(A) B
    const ordinal_type mA = A.NumRows();
    const ordinal_type nB = B.NumCols();

    if (nB > 0) {
      for (ordinal_type k=0;k<mA;++k) {
        row_view_type &a = A.RowView(k);
        const value_type cdiag = conj(a.Value(0));

        // invert
        if (diagA != Diag::Unit) {
          // b1t = b1t / conj(diag);
          ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, nB),
                          [&](const ordinal_type j) {
                            B.Value(k, j) /= cdiag;
                          });
        }

        // update
        const ordinal_type nnz_a = a.NumNonZeros();
        if (nnz_a > 0) {
          // B2 = B2 - trans(conj(a12t)) b1t
          ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, nB),
                          [&](const ordinal_type j) {
                            // grab b1t
                            const value_type val_at_j = B.Value(k, j);
                            
                            for (ordinal_type i=1;i<nnz_a;++i) {
                              // grab a12t
                              const ordinal_type row_at_i = a.Col(i);
                              const value_type   val_at_i = conj(a.Value(i));
                              
                              // update B2
                              B.Value(row_at_i, j) -= val_at_i*val_at_j;
                            }
                          });
        }
        member.team_barrier();
      }
    }

    return 0;
  }


  // Trsm used in the tri-solve phase: Single RHS
  // ============================================
  template<>
  template<typename ParallelForType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename DenseExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
       AlgoTrsm::ForTriSolveBlockedVar2>
  ::invoke(const typename CrsExecViewTypeA::policy_type::member_type &member,
           const int diagA,
           const ScalarType alpha,
           const CrsExecViewTypeA &A,
           const DenseExecViewTypeB &B) {
    typedef typename CrsExecViewTypeA::ordinal_type      ordinal_type;
    typedef typename CrsExecViewTypeA::value_type        value_type;
    typedef typename CrsExecViewTypeA::row_view_type     row_view_type;
    typedef typename CrsExecViewTypeA::team_factory_type team_factory_type;

    // scale the matrix B with alpha
    scaleDenseMatrix<ParallelForType,ScalarType,DenseExecViewTypeB>(member, alpha, B);

    // Solve a system: AX = B -> B := inv(A) B
    const ordinal_type mA = A.NumRows();
    const ordinal_type nB = B.NumCols();

    if (nB > 0) {
      for (ordinal_type k=0;k<mA;++k) {
        row_view_type &a = A.RowView(k);
        const value_type cdiag = conj(a.Value(0));

        // invert
        if (diagA != Diag::Unit) {
          // b1t = b1t / conj(diag);
          ParallelForType(team_factory_type::createThreadLoopRegion(member, 0, nB),
                          [&](const ordinal_type j) {
                            B.Value(k, j) /= cdiag;
                          });
          member.team_barrier();
        }

        // update
        const ordinal_type nnz_a = a.NumNonZeros();
        if (nnz_a > 0) {
          // B2 = B2 - trans(conj(a12t)) b1t
          ParallelForType(team_factory_type::createThreadLoopRegion(member, 1, nnz_a),
                          [&](const ordinal_type i) {
                            // grab a12t
                            const ordinal_type row_at_i = a.Col(i);
                            const value_type   val_at_i = conj(a.Value(i));

                            for (ordinal_type j=0;j<nB;++j) {
                              // grab b1t
                              const ordinal_type col_at_j = j;
                              const value_type   val_at_j = B.Value(k, j);

                              // update B2
                              B.Value(row_at_i, col_at_j) -= val_at_i*val_at_j;
                            }
                          });
          member.team_barrier();
        }
      }
    }
    return 0;
  }
  
}

#endif
