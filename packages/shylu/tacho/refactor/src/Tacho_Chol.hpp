#ifndef __TACHO_CHOL_HPP__
#define __TACHO_CHOL_HPP__

/// \file Tacho_Trsm.hpp
/// \brief Front interface for Trsm operators
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"
#include "Tacho_Control.hpp"
#include "Tacho_Partition.hpp"

namespace Tacho {
  
  template<int ArgUplo,
           int ArgAlgo, int ArgVariant,
           template<int,int> class ControlType = Control>
  class Chol {
  public:
    // data-parallel interface with nested task generation
    // ===================================================
    template<typename PolicyType,
             typename MemberType,
             typename ExecViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      MemberType &member,
                      ExecViewTypeA &A) {
      fprintf(stderr, ">> Template Args - Uplo %d, Algo %d, Variant %d\n", 
              ArgUplo, ArgAlgo, ArgVariant);           
      TACHO_TEST_FOR_ABORT( true, MSG_INVALID_TEMPLATE_ARGS );
      return -1;
    }

    // task-data parallel interface
    // ===================\=========
    template<typename PolicyType,
             typename ExecViewTypeA>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;

    private:
      ExecViewTypeA _A;

      PolicyType _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor() {}

      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const ExecViewTypeA &A)
        : _A(A),
          _policy(policy)
      { }
      
      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "Chol"; }

      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        const int ierr = Chol::invoke(_policy, member,
                                      _A);
        
        // return for only team leader
        if (!member.team_rank()) { 
          _A.setFuture(typename ExecViewTypeA::future_type());
          r_val = ierr; 
        }
      }

    };

    template<typename PolicyType,
             typename ExecViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static 
    TaskFunctor<PolicyType,ExecViewTypeA>
    createTaskFunctor(const PolicyType &policy,
                      const ExecViewTypeA &A) {
      return TaskFunctor<PolicyType,ExecViewTypeA>
        (policy, A);
    }
    
  };
  
}

#include "Tacho_Chol_Upper.hpp"

#endif
