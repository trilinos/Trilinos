#ifndef __GEMM_CONJTRANS_NOTRANS_SPARSE_SPARSE_SUPERNODES_BY_BLOCKS_HPP__
#define __GEMM_CONJTRANS_NOTRANS_SPARSE_SPARSE_SUPERNODES_BY_BLOCKS_HPP__

/// \file Tacho_Gemm_ConjTrans_NoTrans_SparseSparseSuperNodesByBlocks.hpp
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
       AlgoGemm::SparseSparseSuperNodesByBlocks,Variant::One>
  ::invoke(PolicyType &policy,
           const MemberType &member,
           const ScalarType alpha,
           CrsExecViewTypeA &A,
           CrsExecViewTypeB &B,
           const ScalarType beta,
           CrsExecViewTypeC &C) {

    //typedef typename CrsExecViewTypeA::ordinal_type ordinal_type;

    if (member.team_rank() == 0) {
      DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> AA(A.Hier());
      DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> BB(B.Hier());
      DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> CC(C.Hier());
      
      // be careful for scaling on C as we narrow computation region for factorization
      // {
      //   ordinal_type tr, br, lc, rc;

      //   A.getDataRegion(tr, br, lc, rc);
      //   const ordinal_type offmA = tr, offnA = lc, mA = br - tr + 1, nA = rc - lc + 1;

      //   B.getDataRegion(tr, br, lc, rc);
      //   const ordinal_type offmB = tr, offnB = lc, mB = br - tr + 1, nB = rc - lc + 1;

      //   C.getDataRegion(tr, br, lc, rc);
      //   const ordinal_type offmC = tr, offnC = lc, mC = br - tr + 1, nC = rc - lc + 1;
        
      //   AA.setView(A.Flat(),
      //              Util::max(offmA, offmB), Util::min(mA, mB),
      //              Util::max(offnA, offmC), Util::min(nA, mC));

      //   BB.setView(B.Flat(),
      //              Util::max(offmA, offmB), Util::min(mA, mB),
      //              Util::max(offnB, offnC), Util::min(nB, nC));

      //   CC.setView(C.Flat(),
      //              AA.OffsetCols(), AA.NumCols(),
      //              BB.OffsetCols(), BB.NumCols());
      // }

      Gemm<Trans::ConjTranspose,Trans::NoTranspose,
        AlgoGemm::DenseByBlocks,Variant::One>
        ::invoke(policy, member,
                 alpha, AA, BB, beta, CC);
    }

    return 0;
  }

}

#endif
