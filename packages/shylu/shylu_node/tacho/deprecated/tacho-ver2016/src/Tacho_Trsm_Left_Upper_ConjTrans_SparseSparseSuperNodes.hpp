#ifndef __TRSM_LEFT_UPPER_CONJTRANS_SPARSE_SPARSE_SUPERNODES_HPP__
#define __TRSM_LEFT_UPPER_CONJTRANS_SPARSE_SPARSE_SUPERNODES_HPP__

/// \file Tacho_Trsm_Left_Upper_ConjTrans_SparseSparseSuperNodes.hpp
/// \brief triangular solve for supernodal factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  template<typename MT>
  class DenseMatrixView;

  // Trsm for supernodal factorization
  // =================================
  template<>
  template<typename ScalarType,
           typename CrsExecViewTypeA,
           typename CrsExecViewTypeB>
  inline
  Stat
  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
       AlgoTrsm::SparseSparseSuperNodes,Variant::One>
  ::stat(const int diagA,
         const ScalarType alpha,
         CrsExecViewTypeA &A,
         CrsExecViewTypeB &B) {
    DenseMatrixView<typename CrsExecViewTypeA::flat_mat_base_type> AA(A.Flat());
    DenseMatrixView<typename CrsExecViewTypeA::flat_mat_base_type> BB(B.Flat());
    
    return Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
      AlgoTrsm::ExternalBlas,Variant::One>
      ::stat(diagA, alpha, AA, BB);
  }

  // Trsm for supernodal factorization
  // =================================
  template<>
  template<typename PolicyType,
           typename MemberType,
           typename ScalarType,
           typename CrsExecViewTypeA,
           typename CrsExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
       AlgoTrsm::SparseSparseSuperNodes,Variant::One>
  ::invoke(PolicyType &policy,
           MemberType &member,
           const int diagA,
           const ScalarType alpha,
           CrsExecViewTypeA &A,
           CrsExecViewTypeB &B) {

    if (member.team_rank() == 0) {
      DenseMatrixView<typename CrsExecViewTypeA::flat_mat_base_type> AA(A.Flat());
      DenseMatrixView<typename CrsExecViewTypeA::flat_mat_base_type> BB(B.Flat());

      // all diagonal blocks are supposed and assumed to be full matrix
      // B matrix dimensions should match to A
      Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
        AlgoTrsm::ExternalBlas,Variant::One>
        ::invoke(policy, member, diagA, alpha, AA, BB);
    }

    return 0;
  }

}

#endif
