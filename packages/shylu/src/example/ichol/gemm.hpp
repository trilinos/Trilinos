#pragma once
#ifndef __GEMM_HPP__
#define __GEMM_HPP__

/// \file gemm.hpp
/// \brief Sparse matrix-matrix multiplication on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  template<int ArgTransA, int ArgTransX>
  struct Gemm {
    template<typename ScalarType, 
             typename CrsMatViewType>
    KOKKOS_INLINE_FUNCTION
    static int invoke(const ScalarType alpha,
                      const CrsMatViewType A,
                      const CrsMatViewType X,
                      const ScalarType beta,
                      const CrsMatViewType Y);

    template<typename ScalarType, 
             typename CrsMatViewType>
    class TaskFunctor {
    private:
      ScalarType _alpha, _beta;
      CrsMatViewType _A, _X, _Y;
      
    public:
      TaskFunctor(const ScalarType alpha,
                  const CrsMatViewType A,
                  const CrsMatViewType X,
                  const ScalarType beta,
                  const CrsMatViewType Y) 
        : _alpha(alpha),
          _beta(beta),
          _A(A),
          _X(X),
          _Y(Y) 
      { }

      typedef int value_type;      
      void apply(value_type &r_val) {
        r_val = Gemm::invoke(_alpha, _A, _X, _beta, _Y);
      }
    };

  };

}

#include "gemm_nt_t.hpp"

#endif
