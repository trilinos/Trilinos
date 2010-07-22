// include for ThrustGPUNode method implementations
#include "Kokkos_ThrustGPUNode.cuh"

// includes for all ops
#include "Kokkos_DefaultBlockSparseMultiplyKernelOps.hpp"

#define INSTANTIATE_BLOCKSPARSEMULTIPLY_ORDINAL_SCALAR(ORDINAL, SCALAR) \
  template void Kokkos::ThrustGPUNode::parallel_for<Kokkos::DefaultBlockSparseMultiplyOp1         <SCALAR, ORDINAL, SCALAR, SCALAR> >(int, int, Kokkos::DefaultBlockSparseMultiplyOp1         <SCALAR, ORDINAL, SCALAR, SCALAR>); \
  template void Kokkos::ThrustGPUNode::parallel_for<Kokkos::DefaultBlockSparseMultiplyOp1Transpose<SCALAR, ORDINAL, SCALAR, SCALAR> >(int, int, Kokkos::DefaultBlockSparseMultiplyOp1Transpose<SCALAR, ORDINAL, SCALAR, SCALAR>);

#define INSTANTIATE_BLOCKSPARSEMULTIPLY_ORDINAL_SCALAR_SCALAR(ORDINAL, SCALAR, DRSCALAR) \
  template void Kokkos::ThrustGPUNode::parallel_for<Kokkos::DefaultBlockSparseMultiplyOp1         <SCALAR, ORDINAL, DRSCALAR, DRSCALAR> >(int, int, Kokkos::DefaultBlockSparseMultiplyOp1         <SCALAR, ORDINAL, DRSCALAR, DRSCALAR>); \
  template void Kokkos::ThrustGPUNode::parallel_for<Kokkos::DefaultBlockSparseMultiplyOp1Transpose<SCALAR, ORDINAL, DRSCALAR, DRSCALAR> >(int, int, Kokkos::DefaultBlockSparseMultiplyOp1Transpose<SCALAR, ORDINAL, DRSCALAR, DRSCALAR>);

typedef short int ShortInt; 
#ifdef HAVE_KOKKOS_CUDA_FLOAT
  INSTANTIATE_BLOCKSPARSEMULTIPLY_ORDINAL_SCALAR(int,float)
  INSTANTIATE_BLOCKSPARSEMULTIPLY_ORDINAL_SCALAR_SCALAR(int,int,float)
  INSTANTIATE_BLOCKSPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,float)
#endif

#ifdef HAVE_KOKKOS_CUDA_DOUBLE
  INSTANTIATE_BLOCKSPARSEMULTIPLY_ORDINAL_SCALAR(int,double)
  INSTANTIATE_BLOCKSPARSEMULTIPLY_ORDINAL_SCALAR_SCALAR(int,int,double)
  INSTANTIATE_BLOCKSPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,double)
#endif

#if defined(HAVE_KOKKOS_CUDA_DOUBLE) && defined(HAVE_KOKKOS_CUDA_FLOAT)
  INSTANTIATE_BLOCKSPARSEMULTIPLY_ORDINAL_SCALAR_SCALAR(int,float,double)
  INSTANTIATE_BLOCKSPARSEMULTIPLY_ORDINAL_SCALAR_SCALAR(int,double,float)
#endif

INSTANTIATE_BLOCKSPARSEMULTIPLY_ORDINAL_SCALAR(int,int)
INSTANTIATE_BLOCKSPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,int)
