#pragma once
#ifndef __TRI_SOLVE_HPP__
#define __TRI_SOLVE_HPP__

/// \file tri_solve.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example {

  using namespace std;

  template<int ArgUplo, int ArgTrans, int ArgAlgo>
  struct TriSolve {
    static int blocksize;

    // data-parallel interface
    // =======================
    template<typename ParallelForType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    KOKKOS_INLINE_FUNCTION
    static int invoke(const typename ExecViewTypeA::policy_type::member_type &member,
                      const int diagA,
                      ExecViewTypeA &A,
                      ExecViewTypeB &B);

    // task-data parallel interface
    // ============================
    template<typename ParallelForType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    class TaskFunctor {
    private:
      int _diagA;
      ExecViewTypeA _A;
      ExecViewTypeB _B;

    public:
      typedef typename ExecViewTypeA::policy_type policy_type;
      typedef typename policy_type::member_type member_type;
      typedef int value_type;

      TaskFunctor(const int diagA,
                  const ExecViewTypeA A,
                  const ExecViewTypeB B)
        : _diagA(diagA),
          _A(A),
          _B(B)
      { }

      string Label() const { return "TriSolve"; }

      // task execution
      void apply(value_type &r_val) {
        r_val = TriSolve::invoke<ParallelForType>(policy_type::member_null(), _diagA, _A, _B);
      }

      // task-data execution
      void apply(const member_type &member, value_type &r_val) {
        r_val = TriSolve::invoke<ParallelForType>(member, _diagA, _A, _B);
      }

    };
  };

  template<int ArgUplo, int ArgTrans, int ArgAlgo> int TriSolve<ArgUplo,ArgTrans,ArgAlgo>::blocksize = 32;
}

// basic utils
#include "util.hpp"
#include "partition.hpp"

// unblocked version blas operations
#include "scale.hpp"

// blocked version blas operations
#include "gemm.hpp"
#include "trsm.hpp"
#include "herk.hpp"

// triangular solve
#include "tri_solve_u_ct_unblocked.hpp"
#include "tri_solve_u_ct_blocked.hpp"
#include "tri_solve_u_ct_by_blocks.hpp"

#include "tri_solve_u_nt_unblocked.hpp"
#include "tri_solve_u_nt_blocked.hpp"
#include "tri_solve_u_nt_by_blocks.hpp"

#endif
