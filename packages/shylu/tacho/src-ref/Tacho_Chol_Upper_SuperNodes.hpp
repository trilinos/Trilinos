#ifndef __TACHO_CHOL_UPPER_SUPERNODES_HPP__
#define __TACHO_CHOL_UPPER_SUPERNODES_HPP__

/// \file Tacho_Chol_Upper_SuperNodes.hpp
/// \brief Supernodal Chloesky factorization.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  template<typename MT>
  class DenseMatrixView;

  /// Supernodal Cholesky
  /// ===================
  template<>
  template<typename PolicyType,
           typename MemberType,
           typename CrsExecViewTypeA>
  KOKKOS_INLINE_FUNCTION
  int
  Chol<Uplo::Upper,
       AlgoChol::SuperNodes,Variant::One>
  ::invoke(PolicyType &policy,
           const MemberType &member,
           CrsExecViewTypeA &A) {

    int r_val = 0;
    if (member.team_rank() == 0) {
      DenseMatrixView<typename CrsExecViewTypeA::flat_mat_base_type> AA(A.Flat());
      
      // all diagonal blocks are supposed and assumed to be full matrix
      r_val = Chol<Uplo::Upper,
        AlgoChol::ExternalLapack,Variant::One>
        ::invoke(policy, member,
                 AA);
    }
    
    return r_val;
  }

}

#endif
