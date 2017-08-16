#ifndef __TACHO_CHOL_UPPER_SUPERNODES_BY_BLOCKS_HPP__
#define __TACHO_CHOL_UPPER_SUPERNODES_BY_BLOCKS_HPP__

/// \file Tacho_Chol_Upper_SuperNodesByBlocks.hpp
/// \brief Supernodal Chloesky factorization.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  template<typename MT>
  class DenseMatrixView;

  template<>
  template<typename CrsExecViewTypeA>
  inline
  Stat
  Chol<Uplo::Upper,
       AlgoChol::SuperNodesByBlocks,Variant::One>
  ::stat(CrsExecViewTypeA &A) {

    DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> AA(A.Hier());
      
    return Chol<Uplo::Upper,
      AlgoChol::DenseByBlocks,Variant::One>
      ::stat(AA);
  }

  /// Supernodal Cholesky
  /// ===================
  template<>
  template<typename PolicyType,
           typename MemberType,
           typename CrsExecViewTypeA>
  KOKKOS_INLINE_FUNCTION
  int
  Chol<Uplo::Upper,
       AlgoChol::SuperNodesByBlocks,Variant::One>
  ::invoke(PolicyType &policy,
           MemberType &member,
           CrsExecViewTypeA &A) {

    int r_val = 0;
    if (member.team_rank() == 0) {
      DenseMatrixView<typename CrsExecViewTypeA::hier_mat_base_type> AA(A.Hier());
      
      // all diagonal blocks are supposed and assumed to be full matrix
      r_val = Chol<Uplo::Upper,
        AlgoChol::DenseByBlocks,Variant::One>
        ::invoke(policy, member,
                 AA);
    }
    
    return r_val;
  }

}

#endif
