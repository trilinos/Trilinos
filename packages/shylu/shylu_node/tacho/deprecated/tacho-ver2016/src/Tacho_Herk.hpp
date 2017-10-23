#ifndef __TACHO_HERK_HPP__
#define __TACHO_HERK_HPP__

/// \file Tacho_Herk.hpp
/// \brief Front interface for Herk operators
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"
#include "Tacho_Control.hpp"
#include "Tacho_Partition.hpp"

namespace Tacho {
  
  template<int ArgUplo, int ArgTrans,
           int ArgAlgo, int ArgVariant,
           template<int,int> class ControlType = Control>
  class Herk {
  public:
    // statistics
    // ==========
    template<typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeC>
    inline
    static Stat stat(const ScalarType alpha,
                     ExecViewTypeA &A,
                     const ScalarType beta,
                     ExecViewTypeC &C) { 
      printf(">> Template Args - Uplo %d, Trans %d, Algo %d, Variant %d\n", 
             ArgUplo, ArgTrans, ArgAlgo, ArgVariant);  
      TACHO_TEST_FOR_ABORT( true, MSG_INVALID_TEMPLATE_ARGS );
      return Stat();
    }

    // data-parallel interface with nested task generation
    // ===================================================
    template<typename PolicyType,
             typename MemberType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeC>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      MemberType &member,
                      const ScalarType alpha,
                      ExecViewTypeA &A,
                      const ScalarType beta,
                      ExecViewTypeC &C) { 
      printf(">> Template Args - Uplo %d, Trans %d, Algo %d, Variant %d\n", 
             ArgUplo, ArgTrans, ArgAlgo, ArgVariant);  
      TACHO_TEST_FOR_ABORT( true, MSG_INVALID_TEMPLATE_ARGS );
      return -1;
    }

    // task-data parallel interface
    // ===================\=========
    template<typename PolicyType,
             typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeC>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;

    private:
      ScalarType _alpha, _beta;
      ExecViewTypeA _A;
      ExecViewTypeC _C;

      PolicyType _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const ScalarType alpha,
                  const ExecViewTypeA &A,
                  const ScalarType beta,
                  const ExecViewTypeC &C)
        : _alpha(alpha),
          _beta(beta),
          _A(A),
          _C(C),
          _policy(policy)
      { }

      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "Herk"; }

      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        const int ierr = Herk::invoke(_policy, member,
                                      _alpha, _A, _beta, _C);
        
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
             typename ExecViewTypeC>
    KOKKOS_INLINE_FUNCTION
    static 
    TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeC>
    createTaskFunctor(const PolicyType &policy,
                      const ScalarType alpha,
                      const ExecViewTypeA &A,
                      const ScalarType beta,
                      const ExecViewTypeC &C) {
      return TaskFunctor<PolicyType,ScalarType,ExecViewTypeA,ExecViewTypeC>
        (policy, alpha, A, beta, C);
    }
    
  };
  
}

#include "Tacho_Herk_Upper_ConjTrans.hpp"

#endif
