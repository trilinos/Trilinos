#pragma once
#ifndef __TRSM_HPP__
#define __TRSM_HPP__

/// \file trsm.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example {

  using namespace std;

  template<int ArgSide,int ArgUplo, int ArgTrans, int ArgAlgo>
  struct Trsm {

    // data-parallel interface
    // =======================
    template<typename ScalarType,
             typename CrsExecViewType,
             typename ParallelForType>
    KOKKOS_INLINE_FUNCTION
    static int invoke(const typename CrsExecViewType::policy_type::member_type &member,
                      const int diag,
                      const ScalarType alpha,
                      const CrsExecViewType &A,
                      const CrsExecViewType &B);

    // task-data parallel interface
    // ============================
    template<typename ScalarType,
             typename CrsExecViewType,
             typename ParallelForType>
    class TaskFunctor {
    private:
      int _diag;
      ScalarType _alpha;
      CrsExecViewType _A, _B;

    public:
      typedef typename CrsExecViewType::policy_type policy_type;
      typedef typename policy_type::member_type member_type;
      typedef int value_type;

      TaskFunctor(const int diag,
                  const ScalarType alpha,
                  const CrsExecViewType A,
                  const CrsExecViewType B)
        : _diag(diag),
          _alpha(alpha),
          _A(A),
          _B(B)
      { }

      string Label() const { return "Trsm"; }

      // task execution
      void apply(value_type &r_val) {
        r_val = Trsm::invoke<ScalarType,CrsExecViewType,ParallelForType>(policy_type::member_null(),
                                                                         _diag, _alpha, _A, _B);
      }

      // task-data execution
      void apply(const member_type &member, value_type &r_val) const {
        r_val = Trsm::invoke<ScalarType,CrsExecViewType,ParallelForType>(member,
                                                                        _diag, _alpha, _A, _B);
      }

    };
  };

}

//#include "trsm_r_l_t.hpp"
#include "trsm_l_u_t.hpp"

#endif
