#ifndef KOKKOS_THRUSTGPUNODE_CUH_
#define KOKKOS_THRUSTGPUNODE_CUH_

#include <thrust/for_each.h>
#include <thrust/transform_reduce.h>

// must define this before including any kernels
#define KERNEL_PREFIX __device__ __host__

// MUST define this to prevent bringing in implementation of CUDANodeMemoryModel (and therefore, half of Teuchos)
#define KOKKOS_NO_INCLUDE_INSTANTIATIONS
#include <Kokkos_ThrustGPUNode.hpp>

namespace Kokkos {

  template <class WDP>
  void ThrustGPUNode::parallel_for(int begin, int end, WDP wd) {
    // wrap in Thrust and hand to thrust::for_each
  }

  template <class WDP>
  typename WDP::ReductionType
  ThrustGPUNode::parallel_reduce(int begin, int end, WDP wd) 
  {
    // wrap in Thrust and hand to thrust::transform_reduce
    return 0;
  }

}

#endif
