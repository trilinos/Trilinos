#ifndef __TACHO_TRI_SOLVE_UPPER_NOTRANS_SUPERNODES_HPP__
#define __TACHO_TRI_SOLVE_UPPER_NOTRANS_SUPERNODES_HPP__


/// \file Tacho_TriSolve_Upper_NoTrans_SuperNodes.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///

namespace Tacho {

  template<>
  template<typename PolicyType,
           typename MemberType,
           typename CrsExecViewTypeA,
           typename DenseExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  TriSolve<Uplo::Upper,Trans::NoTranspose,
           AlgoTriSolve::SuperNodes,Variant::One>
    ::invoke(PolicyType &policy,
             MemberType &member,
             const int diagA,
             CrsExecViewTypeA &A,
             DenseExecViewTypeB &B) {

    if (member.team_rank() == 0) {
      typedef typename CrsExecViewTypeA::flat_mat_base_type flat_mat_base_type;
      DenseMatrixView<flat_mat_base_type> AA(A.Flat());
      Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,
      AlgoTrsm::SparseDenseSuperNodes,Variant::One>
      ::invoke(policy, member, diagA, 1.0, AA, B);
    }
    return 0;
  }

}

#endif
