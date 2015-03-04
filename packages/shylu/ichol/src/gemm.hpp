#pragma once
#ifndef __GEMM_HPP__
#define __GEMM_HPP__

/// \file gemm.hpp
/// \brief Sparse matrix-matrix multiplication on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  template<int ArgTransA, int ArgTransB, int ArgAlgo>
  struct Gemm {
    template<typename ScalarType, 
             typename CrsMatViewType,
             typename ParallelForType>
    KOKKOS_INLINE_FUNCTION
    static int invoke(const ScalarType alpha,
                      const CrsMatViewType A,
                      const CrsMatViewType B,
                      const ScalarType beta,
                      const CrsMatViewType C);

    template<typename ScalarType, 
             typename CrsMatViewType,
             typename ParallelForType>
    class TaskFunctor {
    private:
      ScalarType _alpha, _beta;
      CrsMatViewType _A, _B, _C;
      
    public:
      TaskFunctor(const ScalarType alpha,
                  const CrsMatViewType A,
                  const CrsMatViewType B,
                  const ScalarType beta,
                  const CrsMatViewType C) 
        : _alpha(alpha),
          _beta(beta),
          _A(A),
          _B(B),
          _C(C) 
      { }

      string Label() const { return "Gemm"; }

      // task execution
      typedef int value_type;      
      void apply(value_type &r_val) {
        r_val = Gemm::invoke<ScalarType,CrsMatViewType,ParallelForType>(_alpha, _A, _B, _beta, _C);
      }

      // task-data execution
      void operator()(const member_type &member) const {
        Gemm::invoke<ScalarType,CrsMatViewType,ParallelForType>(_alpha, _A, _B, _beta, _C);
      }
        
    };

  };

}


// #include "gemm_nt_t.hpp"
#include "gemm_t_nt.hpp"

#endif
