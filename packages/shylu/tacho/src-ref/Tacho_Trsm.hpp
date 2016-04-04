#ifndef __TACHO_TRSM_HPP__
#define __TACHO_TRSM_HPP__

/// \file Tacho_Trsm.hpp
/// \brief Front interface for Trsm operators
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"
#include "Tacho_Control.hpp"
#include "Tacho_Partition.hpp"

namespace Tacho {
  
  template<int ArgSide,int ArgUplo, int ArgTrans,
           int ArgAlgo, int ArgVariant,
           template<int,int> class ControlType = Control>
  class Trsm {
  public:
    // data-parallel interface with nested task generation
    // ===================================================
    template<typename PolicyType,
             typename MemberType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      const MemberType &member,
                      const int diagA,
                      const ScalarType alpha,
                      ExecViewTypeA &A,
                      ExecViewTypeB &B) {
      fprintf(stderr, ">> Template Args - Side %d, Uplo %d, Trans %d, Algo %d, Variant %d\n", 
              ArgSide, ArgUplo, ArgTrans, ArgAlgo, ArgVariant);           
      TACHO_TEST_FOR_ABORT( true, MSG_INVALID_TEMPLATE_ARGS );
      return 0;
    }

    // task-data parallel interface
    // ===================\=========
    template<typename PolicyType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;

    private:
      int _diagA;
      ScalarType _alpha;

      ExecViewTypeA _A;
      ExecViewTypeB _B;

      PolicyType _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const int diagA,
                  const ScalarType alpha,
                  const ExecViewTypeA &A,
                  const ExecViewTypeB &B)
        : _diagA(diagA),
          _alpha(alpha),
          _A(A),
          _B(B),
          _policy(policy)
      { }

      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "Trsm"; }

      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = Trsm::invoke(_policy, _policy.member_single(),
                             _diagA, _alpha, _A, _B);
        _B.setFuture(typename ExecViewTypeB::future_type());
      }

      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {
        const int ierr = Trsm::invoke(_policy, member,
                                      _diagA, _alpha, _A, _B);
        
        // return for only team leader
        if (!member.team_rank()) { 
          _B.setFuture(typename ExecViewTypeB::future_type());
          r_val = ierr; 
        }
      }

    };

    template<typename PolicyType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    KOKKOS_INLINE_FUNCTION
    static 
    TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeB>
    createTaskFunctor(const PolicyType &policy,
                      const int diagA,
                      const ScalarType alpha,
                      const ExecViewTypeA &A,
                      const ExecViewTypeB &B) {
      return TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeB>
        (policy, diagA, alpha, A, B);
    }
    
  };
  
}

#include "Tacho_Trsm_Left_Upper_ConjTrans.hpp"

#endif
