#pragma once
#ifndef __TRI_SOLVE_HPP__
#define __TRI_SOLVE_HPP__

/// \file trsm.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example {

  using namespace std;

  template<int ArgUplo, int ArgTrans, int ArgAlgo>
  struct Trsm {

    // data-parallel interface
    // =======================
    template<typename ParallelForType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    KOKKOS_INLINE_FUNCTION
    static int invoke(const typename ExecViewTypeA::policy_type::member_type &member,
                      const int diag,
                      const ExecViewTypeA &A,
                      const ExecViewTypeB &B);

    // task-data parallel interface
    // ============================
    template<typename ParallelForType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    class TaskFunctor {
    private:
      int _diag;
      ExecViewTypeA _A;
      ExecViewTypeB _B;

    public:
      typedef typename ExecViewTypeA::policy_type policy_type;
      typedef typename policy_type::member_type member_type;
      typedef int value_type;

      TaskFunctor(const int diag,
                  const ExecViewTypeA A,
                  const ExecViewTypeB B)
        : _diag(diag),
          _A(A),
          _B(B)
      { }

      string Label() const { return "TriSolve"; }

      // task execution
      void apply(value_type &r_val) {
        r_val = TriSolve::invoke<ParallelForType,ScalarType,
          ExecViewTypeA,ExecViewTypeB>(policy_type::member_null(), _diag, _A, _B);
      }

      // task-data execution
      void apply(const member_type &member, value_type &r_val) const {
        r_val = TriSolve::invoke<ParallelForType,ScalarType,
          ExecViewTypeA,ExecViewTypeB>(member, _diag,_A, _B);
      }

    };
  };

}


//#include "tri_solve_u_ct_unblocked.hpp"
//#include "tri_solve_u_ct_blocked.hpp"
#include "tri_solve_u_ct_by_blocks.hpp"

//#include "tri_solve_u_nt_unblocked.hpp"
//#include "tri_solve_u_ct_blocked.hpp"
#include "tri_solve_u_ct_by_blocks.hpp"

#endif
