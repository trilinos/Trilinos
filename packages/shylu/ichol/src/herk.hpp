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
    template<typename ScalarType, 
             typename CrsMatViewType,
             typename ParallelForType>
    KOKKOS_INLINE_FUNCTION
    static int invoke(const ParallelForType::member_type member,
                      const ScalarType alpha,
                      const CrsMatViewType A,
                      const ScalarType beta,
                      const CrsMatViewType C);

    template<typename ScalarType, 
             typename CrsMatViewType,
             typename ParallelForType>
    class TaskFunctor {
    private:
      ScalarType _alpha, _beta;
      CrsMatViewType _A, _C;
      
    public:
      TaskFunctor(const ScalarType alpha,
                  const CrsMatViewType A,
                  const ScalarType beta,
                  const CrsMatViewType C) 
        : _alpha(alpha),
          _beta(beta),
          _A(A),
          _C(C)
      { }

      string Label() const { return "Herk"; }

      // task execution
      typedef int value_type;      
      void apply(value_type &r_val) {
        r_val = Herk::invoke<ScalarType,CrsMatViewType,ParallelForType>(ParallelForType::Root, _alpha, _A, _beta, _C);
      }
      // task-data execution
      void operator()(const member_type &member) const {
        Herk::invoke<ScalarType,CrsMatViewType,ParallelForType>(member, _alpha, _A, _beta, _C);
      }

    };

  };

}

#include "herk_u_t.hpp"

#endif
