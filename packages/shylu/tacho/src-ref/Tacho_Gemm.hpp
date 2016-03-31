#ifndef __TACHO_GEMM_HPP__
#define __TACHO_GEMM_HPP__

/// \file Tacho_Gemm.hpp
/// \brief Front interface for Gemm operators
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"
#include "Tacho_Control.hpp"
#include "Tacho_Partition.hpp"

#include "Tacho_ScaleCrsMatrix.hpp"
#include "Tacho_ScaleDenseMatrix.hpp"

namespace Tacho {
  
  template<int ArgTransA, int ArgTransB, 
           int ArgAlgo, int ArgVariant,
           template<int,int> class ControlType = Control>
  struct Gemm {

    // data-parallel interface with nested task generation
    // ===================================================
    template<typename PolicyType,
             typename MemberType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB,
             typename ExecViewTypeC>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      const MemberType &member,
                      const ScalarType alpha,
                      ExecViewTypeA &A,
                      ExecViewTypeB &B,
                      const ScalarType beta,
                      ExecViewTypeC &C) { 
      fprintf(stderr, ">> Template Args - TransA %d, TransB %d, Algo %d, Variant %d\n", 
              ArgTransA, ArgTransB, ArgAlgo, ArgVariant);  
      TACHO_TEST_FOR_ABORT( true, MSG_INVALID_TEMPLATE_ARGS );
      return 0;
    }

    // task-data parallel interface
    // ===================\=========
    template<typename PolicyType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB,
             typename ExecViewTypeC>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;

    private:
      ScalarType _alpha, _beta;
      ExecViewTypeA _A;
      ExecViewTypeB _B;
      ExecViewTypeC _C;

      PolicyType _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const ScalarType alpha,
                  const ExecViewTypeA &A,
                  const ExecViewTypeB &B,
                  const ScalarType beta,
                  const ExecViewTypeC &C)
        : _alpha(alpha),
          _beta(beta),
          _A(A),
          _B(B),
          _C(C),
          _policy(policy)
      { }

      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "Gemm"; }

      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = Gemm::invoke(_policy, _policy.member_single(),
                             _alpha, _A, _B, _beta, _C);
        _C.setFuture(typename ExecViewTypeC::future_type());
      }

      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {
        const int ierr = Gemm::invoke(_policy, member,
                                      _alpha, _A, _B, _beta, _C);
        
        // return for only team leader
        if (!member.team_rank()) { 
          _C.setFuture(typename ExecViewTypeC::future_type());
          r_val = ierr; 
        }
      }

    };

    template<typename PolicyType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB,
             typename ExecViewTypeC>
    KOKKOS_INLINE_FUNCTION
    static 
    TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeB,ExecViewTypeC>
    createTaskFunctor(const PolicyType &policy,
                      const ScalarType alpha,
                      const ExecViewTypeA &A,
                      const ExecViewTypeB &B,
                      const ScalarType beta,
                      const ExecViewTypeC &C) {
      return TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeB,ExecViewTypeC>
        (policy, alpha, A, B, beta, C);
    }
    
  };
  
}

#include "Tacho_Gemm_NoTrans_NoTrans.hpp"
#include "Tacho_Gemm_ConjTrans_NoTrans.hpp"

#endif
