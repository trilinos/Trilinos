#ifndef __HERK_UPPER_CONJTRANS_SPARSE_SPARSE_SUPERNODES_HPP__
#define __HERK_UPPER_CONJTRANS_SPARSE_SPARSE_SUPERNODES_HPP__

/// \file Tacho_Herk_Upper_ConjTrans_SparseSparseSuperNodes.hpp
/// \brief Hermitian rank-k update for supernodal factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  template<typename MT>
  class DenseMatrixView;

  // Herk used in the supernodal factorization
  // =========================================
  template<>
  template<typename PolicyType,
           typename MemberType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename CrsExecViewTypeC>
  KOKKOS_INLINE_FUNCTION
  int
  Herk<Uplo::Upper,Trans::ConjTranspose,
       AlgoHerk::SparseSparseSuperNodes,Variant::One>
  ::invoke(PolicyType &policy,
           const MemberType &member,
           const ScalarType alpha,
           CrsExecViewTypeA &A,
           const ScalarType beta,
           CrsExecViewTypeC &C) {

    typedef typename CrsExecViewTypeA::ordinal_type ordinal_type;

    if (member.team_rank() == 0) {
      DenseMatrixView<typename CrsExecViewTypeA::flat_mat_base_type> AA; //(A.Flat());
      DenseMatrixView<typename CrsExecViewTypeA::flat_mat_base_type> CC; //(C.Flat());
      
      {
        ordinal_type tr, br, lc, rc;

        A.getDataRegion(tr, br, lc, rc);
        const ordinal_type offm = tr, offn = lc, m = br - tr + 1, n = rc - lc + 1;

        AA.setView(A.Flat(), 
                   offm, m,
                   offn, n);

        CC.setView(C.Flat(),
                   offn, n,
                   offn, n);
      }

      Herk<Uplo::Upper,Trans::ConjTranspose,
        AlgoHerk::ExternalBlas,Variant::One>
        ::invoke(policy, member,
                 alpha, AA, beta, CC);
    }

    return 0;
  }

}

#endif
