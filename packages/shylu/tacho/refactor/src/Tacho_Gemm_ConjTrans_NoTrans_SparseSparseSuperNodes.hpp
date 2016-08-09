#ifndef __GEMM_CONJTRANS_NOTRANS_SPARSE_SPARSE_SUPERNODES_HPP__
#define __GEMM_CONJTRANS_NOTRANS_SPARSE_SPARSE_SUPERNODES_HPP__

/// \file Tacho_Gemm_ConjTrans_NoTrans_SparseSparseSuperNodes.hpp
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
           typename CrsExecViewTypeB,
           typename CrsExecViewTypeC>
  KOKKOS_INLINE_FUNCTION
  int
  Gemm<Trans::ConjTranspose,Trans::NoTranspose,
       AlgoGemm::SparseSparseSuperNodes,Variant::One>
  ::invoke(PolicyType &policy,
           MemberType &member,
           const ScalarType alpha,
           CrsExecViewTypeA &A,
           CrsExecViewTypeB &B,
           const ScalarType beta,
           CrsExecViewTypeC &C) {

    if (member.team_rank() == 0) {
      DenseMatrixView<typename CrsExecViewTypeA::flat_mat_base_type> AA(A.Flat());
      DenseMatrixView<typename CrsExecViewTypeA::flat_mat_base_type> BB(B.Flat());
      DenseMatrixView<typename CrsExecViewTypeA::flat_mat_base_type> CC(C.Flat());
      
      Gemm<Trans::ConjTranspose,Trans::NoTranspose,
        AlgoGemm::ExternalBlas,Variant::One>
        ::invoke(policy, member,
                 alpha, AA, BB, beta, CC);
    }

    return 0;
  }

}

#endif
