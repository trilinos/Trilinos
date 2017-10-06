#ifndef __GEMM_NOTRANS_NOTRANS_SPARSE_DENSE_SUPERNODES_HPP__
#define __GEMM_NOTRANS_NOTRANS_SPARSE_DENSE_SUPERNODES_HPP__

/// \file Tacho_Gemm_NoTrans_NoTrans_SparseDenseSuperNodes.hpp
/// \brief matrix-matrix multiplication for supernodal factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  template<typename MT>
  class DenseMatrixView;

  /// Gemm for supernodal factorization
  /// =================================
  template<>
  template<typename PolicyType,
           typename MemberType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename DenseExecViewTypeB,
           typename DenseExecViewTypeC>
  KOKKOS_INLINE_FUNCTION
  int
  Gemm<Trans::NoTranspose,Trans::NoTranspose,
       AlgoGemm::SparseDenseSuperNodes,Variant::One>
  ::invoke(PolicyType &policy,
           MemberType &member,
           const ScalarType alpha,
           CrsExecViewTypeA &A,
           DenseExecViewTypeB &B,
           const ScalarType beta,
           DenseExecViewTypeC &C) {

    if (member.team_rank() == 0) {
      DenseMatrixView<typename CrsExecViewTypeA::flat_mat_base_type> AA(A.Flat());
      
      Gemm<Trans::NoTranspose,Trans::NoTranspose,
        AlgoGemm::ExternalBlas,Variant::One>
        ::invoke(policy, member,
                 alpha, AA, B, beta, C);
    }

    return 0;
  }

}

#endif
