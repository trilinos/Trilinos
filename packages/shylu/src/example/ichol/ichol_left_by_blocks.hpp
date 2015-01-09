#pragma once
#ifndef __ICHOL_LEFT_BY_BLOCKS_HPP__
#define __ICHOL_LEFT_BY_BLOCKS_HPP__

/// \file ichol_left_by_blocks.hpp
/// \brief Incomplete Cholesky factorization by blocks.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "partition.hpp"

#include "scale.hpp"
#include "dot.hpp"

#include "gemv.hpp"
#include "gemm.hpp"

#include "trsv.hpp"
#include "trsm.hpp"

#include "ichol.hpp"

namespace Example { 

  using namespace std;

  template<int ArgUplo>
  class ICholLeftByBlocks {
  public:

    // hierarchically partitioned matrix interface
    // -------------------------------------------------------------------
    template<typename CrsTaskViewType>
    KOKKOS_INLINE_FUNCTION
    static int invoke(typename CrsTaskViewType::policy_type &policy, 
                      const CrsTaskViewType A);

    template<typename CrsTaskViewType>
    class TaskFunctor {
    private:
      CrsTaskViewType _A;
      
    public:
      TaskFunctor(const CrsTaskViewType A)
        : _A(A)
        { } 
      
      typedef int value_type;
      void apply(value_type &r_val) {
        typename CrsTaskViewType::policy_type policy;
        r_val = ICholLeftByBlocks::invoke(policy, _A);
      }
    };
    
    template<typename CrsTaskViewType>
    KOKKOS_INLINE_FUNCTION
    static int genScalarTask(typename CrsTaskViewType::policy_type &policy, 
                             const CrsTaskViewType A);

    template<typename CrsTaskViewType>
    KOKKOS_INLINE_FUNCTION
    static int genGemmTasks(typename CrsTaskViewType::policy_type &policy, 
                            const CrsTaskViewType A,
                            const CrsTaskViewType X,
                            const CrsTaskViewType Y);
    
    template<typename CrsTaskViewType>
    KOKKOS_INLINE_FUNCTION
    static int genTrsmTasks(typename CrsTaskViewType::policy_type &policy,
                            const CrsTaskViewType A,
                            const CrsTaskViewType B);
  };

}

#include "ichol_left_by_blocks_lower.hpp"
#include "ichol_left_by_blocks_lower_var1.hpp"

//#include "ichol_upper_left_by_blocks.hpp"

#endif
