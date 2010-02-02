// include for ThrustGPUNode method implementations
#include "Kokkos_ThrustGPUNode.cuh"

// includes for all ops
#include "Kokkos_DefaultSparseMultiplyKernelOps.hpp"

#define INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ORDINAL, SCALAR) \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseMultiplyOp1         <SCALAR, ORDINAL, SCALAR, SCALAR, 0> >(int, int, Kokkos::DefaultSparseMultiplyOp1         <SCALAR, ORDINAL, SCALAR, SCALAR, 0> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseMultiplyOp2         <SCALAR, ORDINAL, SCALAR, SCALAR, 0> >(int, int, Kokkos::DefaultSparseMultiplyOp2         <SCALAR, ORDINAL, SCALAR, SCALAR, 0> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeMultiplyOp1<SCALAR, ORDINAL, SCALAR, SCALAR, 0> >(int, int, Kokkos::DefaultSparseTransposeMultiplyOp1<SCALAR, ORDINAL, SCALAR, SCALAR, 0> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeMultiplyOp2<SCALAR, ORDINAL, SCALAR, SCALAR, 0> >(int, int, Kokkos::DefaultSparseTransposeMultiplyOp2<SCALAR, ORDINAL, SCALAR, SCALAR, 0> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseMultiplyOp1         <SCALAR, ORDINAL, SCALAR, SCALAR, 1> >(int, int, Kokkos::DefaultSparseMultiplyOp1         <SCALAR, ORDINAL, SCALAR, SCALAR, 1> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseMultiplyOp2         <SCALAR, ORDINAL, SCALAR, SCALAR, 1> >(int, int, Kokkos::DefaultSparseMultiplyOp2         <SCALAR, ORDINAL, SCALAR, SCALAR, 1> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeMultiplyOp1<SCALAR, ORDINAL, SCALAR, SCALAR, 1> >(int, int, Kokkos::DefaultSparseTransposeMultiplyOp1<SCALAR, ORDINAL, SCALAR, SCALAR, 1> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeMultiplyOp2<SCALAR, ORDINAL, SCALAR, SCALAR, 1> >(int, int, Kokkos::DefaultSparseTransposeMultiplyOp2<SCALAR, ORDINAL, SCALAR, SCALAR, 1> );

#define INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR_SCALAR(ORDINAL, SCALAR, DRSCALAR) \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseMultiplyOp1         <SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 0> >(int, int, Kokkos::DefaultSparseMultiplyOp1         <SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 0> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseMultiplyOp2         <SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 0> >(int, int, Kokkos::DefaultSparseMultiplyOp2         <SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 0> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeMultiplyOp1<SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 0> >(int, int, Kokkos::DefaultSparseTransposeMultiplyOp1<SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 0> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeMultiplyOp2<SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 0> >(int, int, Kokkos::DefaultSparseTransposeMultiplyOp2<SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 0> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseMultiplyOp1         <SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 1> >(int, int, Kokkos::DefaultSparseMultiplyOp1         <SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 1> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseMultiplyOp2         <SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 1> >(int, int, Kokkos::DefaultSparseMultiplyOp2         <SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 1> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeMultiplyOp1<SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 1> >(int, int, Kokkos::DefaultSparseTransposeMultiplyOp1<SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 1> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeMultiplyOp2<SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 1> >(int, int, Kokkos::DefaultSparseTransposeMultiplyOp2<SCALAR, ORDINAL, DRSCALAR, DRSCALAR, 1> );

#ifdef HAVE_KOKKOS_CUDA_FLOAT
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(int,float)
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR_SCALAR(int,int,float)
typedef short int ShortInt; INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,float)
#endif

#ifdef HAVE_KOKKOS_CUDA_DOUBLE
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(int,double)
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR_SCALAR(int,int,double)
typedef short int ShortInt; INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,double)
#endif

INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(int,int)
typedef short int ShortInt; INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,int)
