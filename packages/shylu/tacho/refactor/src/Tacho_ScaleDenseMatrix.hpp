#ifndef __TACHO_SCALE_DENSEMATRIX_HPP__
#define __TACHO_SCALE_DENSEMATRIX_HPP__

/// \file Tacho_ScaleDenseMatrix.hpp
/// \brief Front interface for Trsm operators
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"
#include "Tacho_Control.hpp"
#include "Tacho_Partition.hpp"

namespace Tacho {
  
  class ScaleDenseMatrix {
  public:

    // data-parallel interface with nested task generation
    // ===================================================
    template<typename PolicyType,
             typename MemberType,
             typename ScalarType,
             typename DenseExecViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      const MemberType &member,
                      const ScalarType alpha,
                      DenseExecViewTypeA &A) {
      typedef typename DenseExecViewTypeA::ordinal_type  ordinal_type;

      if (alpha == 1) {
        // do nothing
      } else {
        if (A.BaseObject().ColStride() > A.BaseObject().RowStride()) {
          const ordinal_type nA = A.NumCols();
          if (nA > 0) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, nA),
                                 [&](const ordinal_type j) {
                                   const ordinal_type mA = A.NumRows();
                                   for (ordinal_type i=0;i<mA;++i)
                                     A.Value(i, j) *= alpha;
                                 });
          }
        } else {
          const ordinal_type mA = A.NumRows();
          if (mA > 0) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, mA),
                                 [&](const ordinal_type i) {
                                   const ordinal_type nA = A.NumCols();
                                   for (ordinal_type j=0;j<nA;++j)
                                     A.Value(i, j) *= alpha;
                                 });
          }
        }
        member.team_barrier();
      }

      return 0;
    }

    // task-data parallel interface
    // ===================\=========
    template<typename PolicyType,
             typename ScalarType,
             typename DenseExecViewTypeA>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;

    private:
      ScalarType _alpha;
      DenseExecViewTypeA _A;

      PolicyType _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const ScalarType alpha,
                  const DenseExecViewTypeA &A)
        : _alpha(alpha),
          _A(A),
          _policy(policy)
      { }
      
      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "Scal"; }

      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = ScaleDenseMatrix::invoke(_policy, _policy.member_single(),
                                         _alpha, _A);
        _A.setFuture(typename DenseExecViewTypeA::future_type());
      }

      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {
        const int ierr = ScaleDenseMatrix::invoke(_policy, member,
                                                  _alpha, _A);
        
        // return for only team leader
        if (!member.team_rank()) { 
          _A.setFuture(typename DenseExecViewTypeA::future_type());
          r_val = ierr; 
        }
      }

    };

    template<typename PolicyType,
             typename ScalarType,
             typename DenseExecViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static 
    TaskFunctor<PolicyType,ScalarType,DenseExecViewTypeA>
    createTaskFunctor(const PolicyType &policy,
                      const ScalarType alpha,
                      const DenseExecViewTypeA &A) {
      return TaskFunctor<PolicyType,ScalarType,DenseExecViewTypeA>
        (policy, alpha, A);
    }
    
  };
  
}

#endif
