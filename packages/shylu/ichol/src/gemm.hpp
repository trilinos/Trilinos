#pragma once
#ifndef __GEMM_HPP__
#define __GEMM_HPP__

/// \file gemm.hpp
/// \brief Sparse matrix-matrix multiplication on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example {

  using namespace std;

  template<int ArgTransA, int ArgTransB, int ArgAlgo>
  struct Gemm {

    // data-parallel interface
    // =======================
    template<typename ParallelForType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB,
             typename ExecViewTypeC>
    KOKKOS_INLINE_FUNCTION
    static int invoke(const typename ExecViewTypeA::policy_type::member_type &member,
                      const ScalarType alpha,
                      ExecViewTypeA &A,
                      ExecViewTypeB &B,
                      const ScalarType beta,
                      ExecViewTypeC &C);

    // task-data parallel interface
    // ============================
    template<typename ParallelForType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB,
             typename ExecViewTypeC>
    class TaskFunctor {
    private:
      ScalarType _alpha, _beta;
      ExecViewTypeA _A;
      ExecViewTypeB _B;
      ExecViewTypeC _C;

    public:
      typedef typename ExecViewTypeA::policy_type policy_type;
      typedef typename policy_type::member_type member_type;
      typedef int value_type;

      TaskFunctor(const ScalarType alpha,
                  const ExecViewTypeA A,
                  const ExecViewTypeB B,
                  const ScalarType beta,
                  const ExecViewTypeC C)
        : _alpha(alpha),
          _beta(beta),
          _A(A),
          _B(B),
          _C(C)
      { }

      string Label() const { return "Gemm"; }

      // task execution
      void apply(value_type &r_val) {
        r_val = Gemm::invoke<ParallelForType,ScalarType,
          ExecViewTypeA, ExecViewTypeB, ExecViewTypeC>(policy_type::member_null(), _alpha, _A, _B, _beta, _C);
      }

      // task-data execution
      void apply(const member_type &member, value_type &r_val) {
        r_val = Gemm::invoke<ParallelForType,ScalarType,
          ExecViewTypeA, ExecViewTypeB, ExecViewTypeC>(member, _alpha, _A, _B, _beta, _C);
      }

    };

  };

}


#include "gemm_nt_nt.hpp"
#include "gemm_ct_nt.hpp"

#endif
