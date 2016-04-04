#ifndef __TRSM_LEFT_UPPER_CONJTRANS_SPARSE_SPARSE_SUPERNODES_BY_BLOCKS_HPP__
#define __TRSM_LEFT_UPPER_CONJTRANS_SPARSE_SPARSE_SUPERNODES_BY_BLOCKS_HPP__

/// \file Tacho_Trsm_Left_Upper_ConjTrans_SparseSparseSuperNodesByBlocks.hpp
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
           typename CrsExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
       AlgoTrsm::SparseSparseSuperNodesByBlocks,Variant::One>
  ::invoke(PolicyType &policy,
           const MemberType &member,
           const int diagA,
           const ScalarType alpha,
           CrsExecViewTypeA &A,
           CrsExecViewTypeB &B) {

    if (member.team_rank() == 0) {
      DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> AA(A.Hier());
      DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> BB(B.Hier());

      // {
      //   ordinal_type tr, br, lc, rc;

      //   B.getDataRegion(tr, br, lc, rc);
      //   const ordinal_type offm = tr, offn = lc, m = br - tr + 1, n = rc - lc + 1;

      //   BB.setView(B.Flat(),
      //              offm, m,
      //              offn, n);

      //   AA.setView(A.Flat(),
      //              offm, m,
      //              offm, m);
      // }

      // all diagonal blocks are supposed and assumed to be full matrix
      // B matrix dimensions should match to A
      Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
        AlgoTrsm::DenseByBlocks,Variant::One>
        ::invoke(policy, member, diagA, alpha, AA, BB);
    }

    return 0;
  }

}

#endif
