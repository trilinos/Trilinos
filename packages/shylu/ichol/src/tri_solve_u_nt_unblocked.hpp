#pragma once
#ifndef __TRI_SOLVE_U_NT_UNBLOCKED_HPP__
#define __TRI_SOLVE_U_NT_UNBLOCKED_HPP__

/// \file trsm_l_u_nt.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///

namespace Example {

  using namespace std;

  template<>
  template<typename ParallelForType,
           typename CrsExecViewTypeA,
           typename DenseExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  TriSolve<Uplo::Upper,Trans::NoTranspose,
           AlgoTriSolve::Unblocked>
  ::invoke(const typename CrsExecViewTypeA::policy_type::member_type &member,
           const int diagA,
           CrsExecViewTypeA &A,
           DenseExecViewTypeB &B) {

    return Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,AlgoTrsm::ForTriSolveBlocked>
      ::invoke<ParallelForType>(member, diagA, 1.0, A, B);
  }

}

#endif
