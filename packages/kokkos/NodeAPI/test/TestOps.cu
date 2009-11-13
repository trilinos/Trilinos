// include for ThrustGPUNode method implementations
#include <Kokkos_ThrustGPUNode.cuh>

// includes for all operators
#include "TestOps.hpp"

template void Kokkos::ThrustGPUNode::parallel_for<InitOp<int> >(int, int, InitOp<int>);
template SumOp<float>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce<SumOp<float> >(int, int, SumOp<float>);
template void Kokkos::ThrustGPUNode::parallel_for<InitOp<float> >(int, int, InitOp<float>);
template NullOp::ReductionType Kokkos::ThrustGPUNode::parallel_reduce<NullOp>(int, int, NullOp);
template SumOp<int>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce<SumOp<int> >(int, int, SumOp<int>);
