// include for ThrustGPUNode method implementations
#include "Kokkos_ThrustGPUNode.cuh"

// includes for all ops
#include "Kokkos_DefaultSparseMultiplyKernelOps.hpp"

#define INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ORDINAL, SCALAR) \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseMultiplyOp1         <SCALAR, ORDINAL, SCALAR, SCALAR> >(int, int, Kokkos::DefaultSparseMultiplyOp1         <SCALAR, ORDINAL, SCALAR, SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeMultiplyOp1<SCALAR, ORDINAL, SCALAR, SCALAR> >(int, int, Kokkos::DefaultSparseTransposeMultiplyOp1<SCALAR, ORDINAL, SCALAR, SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseMultiplyOp2         <SCALAR, ORDINAL, SCALAR, SCALAR> >(int, int, Kokkos::DefaultSparseMultiplyOp2         <SCALAR, ORDINAL, SCALAR, SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeMultiplyOp2<SCALAR, ORDINAL, SCALAR, SCALAR> >(int, int, Kokkos::DefaultSparseTransposeMultiplyOp2<SCALAR, ORDINAL, SCALAR, SCALAR> );

INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(int,float)
typedef short int ShortInt; INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,float)
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(int,int)
typedef short int ShortInt; INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,int)
