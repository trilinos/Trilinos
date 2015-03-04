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
             typename CrsMatViewType,
             typename ParallelForType>
    KOKKOS_INLINE_FUNCTION
    static int invoke(const ParallelForType::member_type member,
                      const int diag,
                      const ScalarType alpha,
                      const CrsMatViewType A,
                      const CrsMatViewType B);

    template<typename ScalarType,
             typename CrsMatViewType,
             typename ParallelForType>
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
      
      // task execution
      typedef int value_type;
      void apply(value_type &r_val) {
        r_val = Trsm::invoke<ScalarType,CrsMatViewType,ParallelForType>(ParallelForType::Root,
                                                                        _diag, _alpha, _A, _B);
      }

      // task-data execution
      void operator()(const member_type &member) const {
        Trsm::invoke<ScalarType,CrsMatViewType,ParallelForType>(member, 
                                                                _diag, _alpha, _A, _B);
      }

    };
  };
  
}

#include "trsm_r_l_t.hpp"
#include "trsm_l_u_t.hpp"

#endif
