#ifndef __HERK_UPPER_CONJTRANS_FOR_FACTORIZATION_HPP__
#define __HERK_UPPER_CONJTRANS_FOR_FACTORIZATION_HPP__

/// \file Tacho_Herk_Upper_ConjTrans_ForFactorization.hpp
/// \brief Sparse hermitian rank one update on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  // Herk used in the factorization phase
  // ====================================
  template<>
  template<typename PolicyType,
           typename MemberType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename CrsExecViewTypeC>
  KOKKOS_INLINE_FUNCTION
  int
  Herk<Uplo::Upper,Trans::ConjTranspose,
       AlgoHerk::ForFactorization,Variant::One>
  ::invoke(PolicyType &policy,
           const MemberType &member,
           const ScalarType alpha,
           CrsExecViewTypeA &A,
           const ScalarType beta,
           CrsExecViewTypeC &C) {

    typedef typename CrsExecViewTypeA::ordinal_type  ordinal_type;
    typedef typename CrsExecViewTypeA::value_type    value_type;
    typedef typename CrsExecViewTypeA::row_view_type row_view_type;

    // scale the matrix C with beta
    ScaleCrsMatrix::invoke(policy, member,
                           beta, C);

    // KJ :: fix later
    //
    // the following implementation is not stable 
    // c.Value possible very large and the rank update
    // may be so small. k-loop should goes inside 
    // but it will be expensive as A is row-major.
    // a potential fix is to use small workspace and run
    // the loop with batched computations.
    // 
    // FYI: this is in particular a problem for Cholesky.
    //      at least diagonal should be reconsidered.

    // C(i,j) += alpha*A'(i,k)*A(k,j)
    const ordinal_type mA = A.NumRows();
    for (ordinal_type k=0;k<A.NumRows();++k) {
      row_view_type &a = A.RowView(k);
      const ordinal_type nnz = a.NumNonZeros();
      
      if (nnz > 0) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, nnz),
                             [&](const ordinal_type i) {
                               const ordinal_type row_at_i  = a.Col(i);
                               const value_type   val_at_ik = Util::conj(a.Value(i));

                               row_view_type &c = C.RowView(row_at_i);

                               ordinal_type idx = 0;
                               for (ordinal_type j=i;j<nnz && (idx > -2);++j) {
                                 const ordinal_type col_at_j  = a.Col(j);
                                 const value_type   val_at_kj = a.Value(j);

                                 idx = c.Index(col_at_j, idx);
                                 if (idx >= 0)
                                   c.Value(idx) += alpha*val_at_ik*val_at_kj;
                               }
                             });
        member.team_barrier();
      }
    }

    return 0;
  }

}

#endif
