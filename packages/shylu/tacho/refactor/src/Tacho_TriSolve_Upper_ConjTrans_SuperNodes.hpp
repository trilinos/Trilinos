#ifndef __TACHO_TRI_SOLVE_UPPER_CONJTRANS_SUPERNODES_HPP__
#define __TACHO_TRI_SOLVE_UPPER_CONJTRANS_SUPERNODES_HPP__

/// \file Tacho_TriSolve_Upper_ConjTrans_SuperNodes.hpp
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
  TriSolve<Uplo::Upper,Trans::ConjTranspose,
           AlgoTriSolve::SuperNodes,Variant::One>
  ::invoke(PolicyType &policy,
           MemberType &member,
           const int diagA,
           CrsExecViewTypeA &A,
           DenseExecViewTypeB &B) {
    return 
      Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
      AlgoTrsm::SparseDenseSuperNodes,Variant::One>
      ::invoke(policy, member, diagA, 1.0, A, B);
  }

}

#endif
