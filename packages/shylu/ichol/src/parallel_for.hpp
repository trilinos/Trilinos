#pragma once
#ifndef __PARALLEL_FOR_HPP__
#define __PARALLEL_FOR_HPP__

/// \file parallel_for.hpp
/// \brief A wrapper of "for" loop.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  // Parallel for loop
  // =================
  class ParallelFor {
  public:
    template<typename TeamThreadLoopRegionType, typename LambdaType>
    ParallelFor(const TeamThreadLoopRegionType &loop, 
                const LambdaType &lambda) {
      Kokkos::parallel_for(loop, lambda);
    }
  };
  
}

#endif
