#pragma once
#ifndef __CHOL_HPP__
#define __CHOL_HPP__

/// \file chol.hpp
/// \brief Incomplete Cholesky factorization front interface.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

// basic utils
#include "util.hpp"
#include "control.hpp"
#include "partition.hpp"

namespace Tacho { 

  using namespace std;

  // polymophic function body accordingto different algorithms and control trees
  // ===========================================================================
  template<int ArgUplo, int ArgAlgo, 
           template<int> class ControlType,
           typename ParallelForType,
           typename ExecViewType>
  class FunctionChol {
  public:
    KOKKOS_INLINE_FUNCTION
    int getReturnValue() const { return -1; } 

    KOKKOS_INLINE_FUNCTION
    FunctionChol(typename ExecViewType::policy_type &policy, 
                 const typename ExecViewType::policy_type::member_type &member, 
                 ExecViewType &A) { 
      // specialized versions over-ride this
      ERROR(MSG_INVALID_TEMPLATE_ARGS);
    }
  };    

  template<int ArgUplo, int ArgAlgo, template<int> class ControlType = Control>
  class Chol {
  public:

    // function interface
    // ==================
    template<typename ParallelForType,
             typename ExecViewType>
    KOKKOS_INLINE_FUNCTION
    static int invoke(typename ExecViewType::policy_type &policy, 
                      const typename ExecViewType::policy_type::member_type &member, 
                      ExecViewType &A) {
      return FunctionChol
        <ArgUplo,ArgAlgo,ControlType,ParallelForType,ExecViewType>
        (policy, member, A).getReturnValue();
    }

    // task-data parallel interface
    // ============================
    template<typename ParallelForType,
             typename ExecViewType>
    class TaskFunctor {
    public:
      typedef typename ExecViewType::policy_type policy_type;
      typedef typename policy_type::member_type member_type;
      typedef int value_type;

    private:
      ExecViewType _A;
      
      policy_type &_policy;

    public:
      TaskFunctor(const ExecViewType A)
        : _A(A),
          _policy(ExecViewType::task_factory_type::Policy())
      { } 

      string Label() const { return "Chol"; }
      
      // task execution
      void apply(value_type &r_val) {
        r_val = Chol::invoke<ParallelForType>(_policy, _policy.member_single(), 
                                               _A);
      }

      // task-data execution
      void apply(const member_type &member, value_type &r_val) {
        r_val = Chol::invoke<ParallelForType>(_policy, member, 
                                               _A);
      }

    };

  };
}


// unblocked version blas operations
#include "scale.hpp"

// blocked version blas operations
#include "gemm.hpp"
#include "trsm.hpp"
#include "herk.hpp"

// cholesky
#include "chol_u.hpp"

#endif
