// include for ThrustGPUNode method implementations
#include "Kokkos_ThrustGPUNode.cuh"

// includes for all operators
#include "TestOps.hpp"
#include <thrust/generate.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>

template void Kokkos::ThrustGPUNode::parallel_for<InitOp<int> >(int, int, InitOp<int>);
template SumOp<float>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce<SumOp<float> >(int, int, SumOp<float>);
template void Kokkos::ThrustGPUNode::parallel_for<InitOp<float> >(int, int, InitOp<float>);
template NullOp::ReductionType Kokkos::ThrustGPUNode::parallel_reduce<NullOp>(int, int, NullOp);
template SumOp<int>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce<SumOp<int> >(int, int, SumOp<int>);

struct thrust_test_constant_float
{
  __host__ __device__
  float operator()() {return 1.0f;}
};
struct thrust_test_constant_int
{
  __host__ __device__
  int operator()() {return 1;}
};

void thrust_float_init(thrust::device_vector<float> &buff) {
  thrust::generate( buff.begin(), buff.end(), thrust_test_constant_float() );
}
void thrust_int_init(thrust::device_vector<int> &buff) {
  thrust::generate( buff.begin(), buff.end(), thrust_test_constant_int() );
}

float thrust_float_sum(const thrust::device_vector<float> &buff) {
  return thrust::reduce( buff.begin(), buff.end(), 0.0f, thrust::plus<float>() );
}
int thrust_int_sum(const thrust::device_vector<int> &buff) {
  return thrust::reduce( buff.begin(), buff.end(), 0,    thrust::plus<int>() );
}
