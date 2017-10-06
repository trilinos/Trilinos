#ifndef __TACHO_TRI_SOLVE_UPPER_NOTRANS_UNBLOCKED_HPP__
#define __TACHO_TRI_SOLVE_UPPER_NOTRANS_UNBLOCKED_HPP__


/// \file Tacho_TriSolve_Upper_NoTrans_Unblocked.hpp
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
           AlgoTriSolve::Unblocked,Variant::One>
    ::invoke(PolicyType &policy,
             const MemberType &member,
             const int diagA,
             CrsExecViewTypeA &A,
             DenseExecViewTypeB &B) {
    return 
      Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,
      AlgoTrsm::SparseDenseUnblocked,Variant::One>
      ::invoke(policy, member, diagA, 1.0, A, B);
  }

}

#endif
