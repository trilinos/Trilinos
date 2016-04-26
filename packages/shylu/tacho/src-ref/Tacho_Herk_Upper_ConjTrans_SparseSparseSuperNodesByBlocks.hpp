#ifndef __HERK_UPPER_CONJTRANS_SPARSE_SPARSE_SUPERNODES_BY_BLOCKS_HPP__
#define __HERK_UPPER_CONJTRANS_SPARSE_SPARSE_SUPERNODES_BY_BLOCKS_HPP__

/// \file Tacho_Herk_Upper_ConjTrans_SparseSparseSuperNodesByBlocks.hpp
/// \brief Hermitian rank-k update for supernodal factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  class Util;

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
       AlgoHerk::SparseSparseSuperNodesByBlocks,Variant::One>
  ::invoke(PolicyType &policy,
           const MemberType &member,
           const ScalarType alpha,
           CrsExecViewTypeA &A,
           const ScalarType beta,
           CrsExecViewTypeC &C) {



    if (member.team_rank() == 0) {
      DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> AA; //(A.Hier());
      DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> CC; //(C.Hier());
      
      {
        typedef typename CrsExecViewTypeA::ordinal_type ordinal_type;
        const ordinal_type blksize = Util::max(A.Hier().Value(0,0).NumRows(),
                                               A.Hier().Value(0,0).NumCols());
                                               
        ordinal_type tr, br, lc, rc;
        
        A.getDataRegion(tr, br, lc, rc);
        const ordinal_type 
          offm = tr/blksize, m = br/blksize - offm + 1, 
          offn = lc/blksize, n = rc/blksize - offn + 1;

        AA.setView(A.Hier(), 
                   offm, m,
                   offn, n);

        CC.setView(C.Hier(),
                   offn, n,
                   offn, n);
      }

      Herk<Uplo::Upper,Trans::ConjTranspose,
        AlgoHerk::DenseByBlocks,Variant::One>
        ::invoke(policy, member,
                 alpha, AA, beta, CC);
    }

    return 0;
  }

}

#endif
