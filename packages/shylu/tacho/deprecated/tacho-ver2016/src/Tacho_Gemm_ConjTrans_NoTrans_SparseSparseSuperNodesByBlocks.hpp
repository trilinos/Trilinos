#ifndef __GEMM_CONJTRANS_NOTRANS_SPARSE_SPARSE_SUPERNODES_BY_BLOCKS_HPP__
#define __GEMM_CONJTRANS_NOTRANS_SPARSE_SPARSE_SUPERNODES_BY_BLOCKS_HPP__

/// \file Tacho_Gemm_ConjTrans_NoTrans_SparseSparseSuperNodesByBlocks.hpp
/// \brief matrix-matrix multiplication for supernodal factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  class Util;

  template<typename MT>
  class DenseMatrixView;

  /// Gemm for supernodal factorization
  /// =================================
  template<>
  template<typename ScalarType,
           typename CrsExecViewTypeA,
           typename CrsExecViewTypeB,
           typename CrsExecViewTypeC>
  inline
  Stat
  Gemm<Trans::ConjTranspose,Trans::NoTranspose,
       AlgoGemm::SparseSparseSuperNodesByBlocks,Variant::One>
  ::stat(const ScalarType alpha,
         CrsExecViewTypeA &A,
         CrsExecViewTypeB &B,
         const ScalarType beta,
         CrsExecViewTypeC &C) {
    DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> AA(A.Hier());
    DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> BB(B.Hier());
    DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> CC(C.Hier());
    
    return Gemm<Trans::ConjTranspose,Trans::NoTranspose,
      AlgoGemm::DenseByBlocks,Variant::One>
      ::stat(alpha, AA, BB, beta, CC);
  }

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
       AlgoGemm::SparseSparseSuperNodesByBlocks,Variant::One>
  ::invoke(PolicyType &policy,
           MemberType &member,
           const ScalarType alpha,
           CrsExecViewTypeA &A,
           CrsExecViewTypeB &B,
           const ScalarType beta,
           CrsExecViewTypeC &C) {

    if (member.team_rank() == 0) {
      DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> AA(A.Hier());
      DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> BB(B.Hier());
      DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> CC(C.Hier());
      
      Gemm<Trans::ConjTranspose,Trans::NoTranspose,
        AlgoGemm::DenseByBlocks,Variant::One>
        ::invoke(policy, member,
                 alpha, AA, BB, beta, CC);
    }

    return 0;
  }

}

#endif
