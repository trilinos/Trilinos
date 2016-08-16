#ifndef __TACHO_SCALE_CRSMATRIX_HPP__
#define __TACHO_SCALE_CRSMATRIX_HPP__

/// \file Tacho_ScaleCrsMatrix.hpp
/// \brief Front interface for Trsm operators
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"
#include "Tacho_Control.hpp"
#include "Tacho_Partition.hpp"

namespace Tacho {
  
  class ScaleCrsMatrix {
  public:

    // data-parallel interface with nested task generation
    // ===================================================
    template<typename PolicyType,
             typename MemberType,
             typename ScalarType,
             typename CrsExecViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static int invoke(PolicyType &policy,
                      const MemberType &member,
                      const ScalarType alpha,
                      CrsExecViewTypeA &A) {
      typedef typename CrsExecViewTypeA::ordinal_type  ordinal_type;
      typedef typename CrsExecViewTypeA::row_view_type row_view_type;

      if (alpha == 1) {
        // do nothing
      } else {
        const ordinal_type mA = A.NumRows();
        if (mA > 0) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, mA),
                               [&](const ordinal_type i) {
                                 row_view_type &row = A.RowView(i);
                                 const ordinal_type nnz = row.NumNonZeros();
                                 for (ordinal_type j=0;j<nnz;++j)
                                   row.Value(j) *= alpha;
                               });
          member.team_barrier();
        }
      }
      return 0;
    }

    // task-data parallel interface
    // ===================\=========
    template<typename PolicyType,
             typename ScalarType,
             typename CrsExecViewTypeA>
    class TaskFunctor {
    public:
      typedef typename PolicyType::member_type member_type;
      typedef int value_type;

    private:
      ScalarType _alpha;
      CrsExecViewTypeA _A;

      PolicyType _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor() {}
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const PolicyType &policy,
                  const ScalarType alpha,
                  const CrsExecViewTypeA &A)
        : _alpha(alpha),
          _A(A),
          _policy(policy)
      { }
      
      KOKKOS_INLINE_FUNCTION
      const char* Label() const { return "Scal"; }

      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = ScaleCrsMatrix::invoke(_policy, _policy.member_single(),
                                       _alpha, _A);
        _A.setFuture(typename CrsExecViewTypeA::future_type());
      }

      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {
        const int ierr = ScaleCrsMatrix::invoke(_policy, member,
                                                _alpha, _A);
        
        // return for only team leader
        if (!member.team_rank()) { 
          _A.setFuture(typename CrsExecViewTypeA::future_type());
          r_val = ierr; 
        }
      }

    };

    template<typename PolicyType,
             typename ScalarType,
             typename CrsExecViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static 
    TaskFunctor<PolicyType,ScalarType,CrsExecViewTypeA>
    createTaskFunctor(const PolicyType &policy,
                      const ScalarType alpha,
                      const CrsExecViewTypeA &A) {
      return TaskFunctor<PolicyType,ScalarType,CrsExecViewTypeA>
        (policy, alpha, A);
    }
    
  };
  
}

#endif
