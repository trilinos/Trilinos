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
    template<typename ScalarType,
             typename CrsMatViewType>
    KOKKOS_INLINE_FUNCTION
    static int invoke(const int diag,
                      const ScalarType alpha,
                      const CrsMatViewType A,
                      const CrsMatViewType B);

    template<typename ScalarType,
             typename CrsMatViewType>
    class TaskFunctor {
    private:
      int _diag;
      ScalarType _alpha;
      CrsMatViewType _A, _B;

    public:
      TaskFunctor(const int diag,
                  const ScalarType alpha,
                  const CrsMatViewType A,
                  const CrsMatViewType B)
        : _diag(diag),
          _alpha(alpha),
          _A(A),
          _B(B) 
      { } 

      string Label() const { return "Trsm"; }
      
      typedef int value_type;
      void apply(value_type &r_val) {
        r_val = Trsm::invoke(_diag, _alpha, _A, _B);
      }
    };
  };
  
}

#include "trsm_r_l_t.hpp"
#include "trsm_l_u_t.hpp"

#endif
