#pragma once
#ifndef __ICHOL_HPP__
#define __ICHOL_HPP__

/// \file ichol.hpp
/// \brief Incomplete Cholesky factorization front interface.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  template<int ArgUplo, int ArgAlgo>
  class IChol {
  public:
    static int blocksize;

    // data-parallel interface
    // =======================
    template<typename ParallelForType,
             typename CrsExecViewType>
    KOKKOS_INLINE_FUNCTION
    static int invoke(const typename CrsExecViewType::policy_type::member_type &member, 
                      const CrsExecViewType &A);

    // task-data parallel interface
    // ============================
    template<typename ParallelForType,
             typename CrsExecViewType>
    class TaskFunctor {
    private:
      CrsExecViewType _A;
      
    public:
      typedef typename CrsExecViewType::policy_type policy_type;
      typedef typename policy_type::member_type member_type;
      typedef int value_type;

      TaskFunctor(const CrsExecViewType A)
        : _A(A)
      { } 

      string Label() const { return "IChol"; }
      
      // task execution
      void apply(value_type &r_val) {
        r_val = IChol::invoke<ParallelForType,CrsExecViewType>(policy_type::member_null(), _A);
      }

      // task-data execution
      void apply(const member_type &member, value_type &r_val) const {
        r_val = IChol::invoke<ParallelForType,CrsExecViewType>(member, _A);
      }

    };

  };

  template<int ArgUplo, int ArgAlgo> int IChol<ArgUplo,ArgAlgo>::blocksize = 32;
}

// basic utils
#include "util.hpp"
#include "partition.hpp"

// unblocked version blas operations
//#include "dot.hpp"
#include "scale.hpp"

// blocked version blas operations
#include "gemm.hpp"
#include "trsm.hpp"
#include "herk.hpp"

// left looking: only for testing 
//#include "ichol_left_unblocked.hpp"
//#include "ichol_left_blocked.hpp"
//#include "ichol_left_by_blocks.hpp"

// right looking: performance with CRS
//#include "ichol_right_unblocked.hpp"
#include "ichol_right_unblocked_opt1.hpp"
#include "ichol_right_blocked.hpp"

// task / task-data parallel
#include "ichol_right_by_blocks.hpp"

#endif
