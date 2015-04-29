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
           typename CrsExecViewType,
           typename DenseExecViewType>
  KOKKOS_INLINE_FUNCTION
  int
  TriSolve<Uplo::Upper,Trans::NoTranspose,
           AlgoTriSolve::Unblocked>
  ::invoke(const typename CrsExecViewType::policy_type::member_type &member,
           const int diagA,
           const CrsExecViewType &A,
           const DenseExecViewType &B) {

    Trsm::<Side::Left,Uplo::Upper,Trans::NoTranspose,AlgoTrsm::ForTriSolveBlocked>
      ::invoke<ParallelForType,CrsExecViewType,DenseExecViewType>(member, diagA, A, B);

    return 0;
  }

}

#endif
