#pragma once
#ifndef __ICHOL_HPP__
#define __ICHOL_HPP__

/// \file ichol.hpp
/// \brief Incomplete Cholesky factorization front interface.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  template<int ArgUplo, 
           int ArgAlgo = Algo::LeftUnblocked>
  class IChol {
  public:
    static int blocksize;

    // flat matrix interface
    // -------------------------------------------------------------------
    template<typename CrsMatViewType>
    KOKKOS_INLINE_FUNCTION
    static int invoke(const CrsMatViewType A);

    template<typename CrsMatViewType>
    class TaskFunctor {
    private:
      CrsMatViewType _A;
      
    public:
      TaskFunctor(const CrsMatViewType A)
        : _A(A)
        { } 

      string Label() const { return "IChol"; }
      
      typedef int value_type;
      void apply(value_type &r_val) {
        r_val = IChol::invoke(_A);
      }
    };

  };

  template<int ArgUplo,int ArgAlgo> int IChol<ArgUplo,ArgAlgo>::blocksize = 32;
}

#include "ichol_left_unblocked.hpp"
#include "ichol_left_blocked.hpp"

#endif
