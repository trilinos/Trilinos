#pragma once
#ifndef __HERK_HPP__
#define __HERK_HPP__

/// \file herk.hpp
/// \brief Sparse hermitian rank one update on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example {

  using namespace std;

  template<int ArgUplo, int ArgTrans, int ArgAlgo>
  struct Herk {

    // data-parallel interface
    // =======================
    template<typename ScalarType,
             typename CrsExecViewType,
             typename ParallelForType>
    KOKKOS_INLINE_FUNCTION
    static int invoke(const typename CrsExecViewType::policy_type::member_type &member,
                      const ScalarType alpha,
                      const CrsExecViewType &A,
                      const ScalarType beta,
                      const CrsExecViewType &C);

    // task-data parallel interface
    // ============================
    template<typename ScalarType,
             typename CrsExecViewType,
             typename ParallelForType>
    class TaskFunctor {
    private:
      ScalarType _alpha, _beta;
      CrsExecViewType _A, _C;

    public:
      typedef typename CrsExecViewType::policy_type::member_type member_type;
      typedef int value_type;

      TaskFunctor(const ScalarType alpha,
                  const CrsExecViewType A,
                  const ScalarType beta,
                  const CrsExecViewType C)
        : _alpha(alpha),
          _beta(beta),
          _A(A),
          _C(C)
      { }

      string Label() const { return "Herk"; }

      // task execution
      void apply(value_type &r_val) {
        r_val = Herk::invoke<ScalarType,CrsExecViewType,ParallelForType>(member_type(), _alpha, _A, _beta, _C);
      }

      // task-data execution
      void apply(const member_type &member, value_type &r_val) const {
        r_val = Herk::invoke<ScalarType,CrsExecViewType,ParallelForType>(member, _alpha, _A, _beta, _C);
      }

    };

  };

}

#include "herk_u_t.hpp"

#endif
