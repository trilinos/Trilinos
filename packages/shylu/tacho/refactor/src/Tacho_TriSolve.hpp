#ifndef __TACHO_TRI_SOLVE_HPP__
#define __TACHO_TRI_SOLVE_HPP__

/// \file Tacho_TriSolve.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  template<int ArgUplo, int ArgTrans, 
           int ArgAlgo, int ArgVariant,
           template<int,int> class ControlType = Control>
  class TriSolve {
  public:

    // data-parallel interface
    // =======================
    template<typename PolicyType,
             typename MemberType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      const MemberType &member,
                      const int diagA,
                      ExecViewTypeA &A,
                      ExecViewTypeB &B) {
      fprintf(stderr, ">> Template Args - Uplo %d, Trans %d, Algo %d, Variant %d\n",
              ArgUplo, ArgTrans, ArgAlgo, ArgVariant);
      TACHO_TEST_FOR_ABORT( true, MSG_INVALID_TEMPLATE_ARGS );
      return -1;
    }

    // task-data parallel interface
    // ============================
    template<typename PolicyType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;

    private:
      int _diagA;
      ExecViewTypeA _A;
      ExecViewTypeB _B;

      PolicyType _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const int diagA,
                  const ExecViewTypeA &A,
                  const ExecViewTypeB &B)
        : _diagA(diagA),
          _A(A),
          _B(B),
          _policy(policy)
      { }

      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "TriSolve"; }

      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = TriSolve::invoke(_policy, _policy.member_single(), 
                                 _diagA, _A, _B);
        _B.setFuture(typename ExecViewTypeB::future_type());
      }

      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {
        const int ierr = TriSolve::invoke(_policy, member, 
                                          _diagA, _A, _B);
        
        // return for only team leader
        if (!member.team_rank()) {
          _B.setFuture(typename ExecViewTypeB::future_type());
          r_val = ierr;
        }
      }

    };

    template<typename PolicyType,
             typename ExecViewTypeA,
             typename ExecViewTypeB>
    KOKKOS_INLINE_FUNCTION
    static
    TaskFunctor<PolicyType,ExecViewTypeA,ExecViewTypeB>
    createTaskFunctor(const PolicyType &policy,
                      const int diagA,
                      const ExecViewTypeA &A,
                      const ExecViewTypeA &B) {
      return TaskFunctor<PolicyType,ExecViewTypeA,ExecViewTypeB>
        (policy, diagA, A, B);
    }

  };

}

//#include "Tacho_TriSolve_Upper_ConjTrans_Unblocked.hpp"
//#include "Tacho_TriSolve_Upper_NoTrans_Unblocked.hpp"

#include "Tacho_TriSolve_Upper_ConjTrans_SuperNodes.hpp"
#include "Tacho_TriSolve_Upper_NoTrans_SuperNodes.hpp"

//#include "Tacho_TriSolve_Upper_ConjTrans_SuperNodesByBlocks.hpp"
//#include "Tacho_TriSolve_Upper_NoTrans_SuperNodesByBlocks.hpp"

#include "Tacho_TriSolve_Upper_ConjTrans_ByBlocks.hpp"
#include "Tacho_TriSolve_Upper_NoTrans_ByBlocks.hpp"

#endif
