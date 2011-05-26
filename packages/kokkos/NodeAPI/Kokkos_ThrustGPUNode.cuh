#ifndef KOKKOS_THRUSTGPUNODE_CUH_
#define KOKKOS_THRUSTGPUNODE_CUH_

#include <thrust/for_each.h>
#include <thrust/transform_reduce.h>
#include <thrust/iterator/counting_iterator.h>

// must define this before including any kernels
#define KERNEL_PREFIX __device__ __host__

// MUST define this to prevent bringing in implementation of CUDANodeMemoryModel (and therefore, half of Teuchos)
#define KOKKOS_NO_INCLUDE_INSTANTIATIONS
#include <Kokkos_ThrustGPUNode.hpp>

namespace Kokkos {

  template <class WDPin> 
  struct ThrustExecuteWrapper {
    mutable WDPin wd;

    inline ThrustExecuteWrapper(WDPin in) : wd(in) {}

    __device__ __host__ inline void operator()(int i) const {
      wd.execute(i);
    }
  };

  template <class WDPin> 
  struct ThrustReduceWrapper {
    mutable WDPin wd;
    inline ThrustReduceWrapper (WDPin in) : wd(in) {}

    __device__ __host__ inline 
    typename WDPin::ReductionType 
    operator()(typename WDPin::ReductionType x, typename WDPin::ReductionType y) {
      return wd.reduce(x,y);
    }
  };

  template <class WDPin>
  struct ThrustGenerateWrapper {
    mutable WDPin wd;
    inline ThrustGenerateWrapper (WDPin in) : wd(in) {}
 
    __device__ __host__ inline 
    typename WDPin::ReductionType
    operator()(int i) {
      return wd.generate(i);
    }
  };

  template <class WDP>
  void ThrustGPUNode::parallel_for(int begin, int end, WDP wd) {
#ifdef HAVE_KOKKOS_DEBUG
    cudaError_t err = cudaGetLastError();
    TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error, 
        "Kokkos::ThrustGPUNode::" << __FUNCTION__ << ": " 
        << "cudaGetLastError() returned error before function call:\n"
        << cudaGetErrorString(err) );
#endif
    // wrap in Thrust and hand to thrust::for_each
    ThrustExecuteWrapper<WDP> body(wd);  
    thrust::counting_iterator<int,thrust::device_space_tag> bit(begin),
                                                            eit(end);
    thrust::for_each( bit, eit, body );
#ifdef HAVE_KOKKOS_DEBUG
    err = cudaThreadSynchronize();
    TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error, 
        "Kokkos::ThrustGPUNode::" << __FUNCTION__ << ": " 
        << "cudaThreadSynchronize() returned error after function call:\n"
        << cudaGetErrorString(err) );
#endif
  };

  template <class WDP>
  typename WDP::ReductionType
  ThrustGPUNode::parallel_reduce(int begin, int end, WDP wd) 
  {
#ifdef HAVE_KOKKOS_DEBUG
    cudaError_t err = cudaGetLastError();
    TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error, 
        "Kokkos::ThrustGPUNode::" << __FUNCTION__ << ": " 
        << "cudaGetLastError() returned error before function call:\n"
        << cudaGetErrorString(err) );
#endif
    // wrap in Thrust and hand to thrust::transform_reduce
    thrust::counting_iterator<int,thrust::device_space_tag> bit(begin),
                                                            eit(end);
    ThrustReduceWrapper<WDP> ROp(wd);
    ThrustGenerateWrapper<WDP> TOp(wd);
    typename WDP::ReductionType init = wd.identity(), ret;
    ret = thrust::transform_reduce( bit, eit, TOp, init, ROp );
#ifdef HAVE_KOKKOS_DEBUG
    err = cudaThreadSynchronize();
    TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error, 
        "Kokkos::ThrustGPUNode::" << __FUNCTION__ << ": " 
        << "cudaThreadSynchronize() returned error after function call:\n"
        << cudaGetErrorString(err) );
#endif
    return ret;
  }

}

#endif
