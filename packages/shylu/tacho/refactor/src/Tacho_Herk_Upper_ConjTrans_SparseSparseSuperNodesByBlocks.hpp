#ifndef __HERK_UPPER_CONJTRANS_SPARSE_SPARSE_SUPERNODES_BY_BLOCKS_HPP__
#define __HERK_UPPER_CONJTRANS_SPARSE_SPARSE_SUPERNODES_BY_BLOCKS_HPP__

/// \file Tacho_Herk_Upper_ConjTrans_SparseSparseSuperNodesByBlocks.hpp
/// \brief Hermitian rank-k update for supernodal factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  class Util;

  template<typename MT>
  class DenseMatrixView;

  template<>
  template<typename ScalarType,
           typename CrsExecViewTypeA,
           typename CrsExecViewTypeC>
  inline
  Stat
  Herk<Uplo::Upper,Trans::ConjTranspose,
       AlgoHerk::SparseSparseSuperNodesByBlocks,Variant::One>
  ::stat(const ScalarType alpha,
         CrsExecViewTypeA &A,
         const ScalarType beta,
         CrsExecViewTypeC &C) {
    DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> AA(A.Hier());
    DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> CC(C.Hier());
    
    return Herk<Uplo::Upper,Trans::ConjTranspose,
      AlgoHerk::DenseByBlocks,Variant::One>
      ::stat(alpha, AA, beta, CC);
  }

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
       AlgoHerk::SparseSparseSuperNodesByBlocks,Variant::One>
  ::invoke(PolicyType &policy,
           MemberType &member,
           const ScalarType alpha,
           CrsExecViewTypeA &A,
           const ScalarType beta,
           CrsExecViewTypeC &C) {



    if (member.team_rank() == 0) {
      DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> AA(A.Hier());
      DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> CC(C.Hier());
      
      Herk<Uplo::Upper,Trans::ConjTranspose,
        AlgoHerk::DenseByBlocks,Variant::One>
        ::invoke(policy, member,
                 alpha, AA, beta, CC);
    }

    return 0;
  }

}

#endif
