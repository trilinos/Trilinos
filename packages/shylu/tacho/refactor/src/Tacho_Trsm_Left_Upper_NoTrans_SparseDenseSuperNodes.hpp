#ifndef __TRSM_LEFT_UPPER_NOTRANS_SPARSE_DENSE_SUPERNODES_HPP__
#define __TRSM_LEFT_UPPER_NOTRANS_SPARSE_DENSE_SUPERNODES_HPP__

/// \file Tacho_Trsm_Left_Upper_NoTrans_SparseDenseSuperNodes.hpp
/// \brief triangular solve for supernodal factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  template<typename MT>
  class DenseMatrixView;

  // Trsm for supernodal factorization
  // =================================
  template<>
  template<typename PolicyType,
           typename MemberType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename DenseExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,
       AlgoTrsm::SparseDenseSuperNodes,Variant::One>
  ::invoke(PolicyType &policy,
           MemberType &member,
           const int diagA,
           const ScalarType alpha,
           CrsExecViewTypeA &A,
           DenseExecViewTypeB &B) {

    if (member.team_rank() == 0) {
      DenseMatrixView<typename CrsExecViewTypeA::flat_mat_base_type> AA(A.Flat());

      Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,
        AlgoTrsm::ExternalBlas,Variant::One>
        ::invoke(policy, member, diagA, alpha, AA, B);
    }

    return 0;
  }

}

#endif
