// include for CudaNode method implementations
#include <Kokkos_CudaNode.cuh>

// includes for all operators
#include "TestOps.hpp"

template void Kokkos::CUDANode::parallel_for<InitOp<int> >(int, int, InitOp<int>);
template SumOp<float>::ReductionType Kokkos::CUDANode::parallel_reduce<SumOp<float> >(int, int, SumOp<float>);
template void Kokkos::CUDANode::parallel_for<InitOp<float> >(int, int, InitOp<float>);
template NullOp::ReductionType Kokkos::CUDANode::parallel_reduce<NullOp>(int, int, NullOp);
template SumOp<int>::ReductionType Kokkos::CUDANode::parallel_reduce<SumOp<int> >(int, int, SumOp<int>);
