// include for ThrustGPUNode method implementations
#include "Kokkos_ThrustGPUNode.cuh"

// includes for all ops
#include "Kokkos_DefaultRelaxationKernelOps.hpp"

#define INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ORDINAL, SCALAR) \
  template void Kokkos::ThrustGPUNode::parallel_for<Kokkos::DefaultChebyshevOp1<SCALAR, ORDINAL>                    >(int, int, Kokkos::DefaultChebyshevOp1<SCALAR, ORDINAL>                   ); \
  template void Kokkos::ThrustGPUNode::parallel_for<Kokkos::DefaultJacobiOp1<SCALAR, ORDINAL>                       >(int, int, Kokkos::DefaultJacobiOp1<SCALAR, ORDINAL>                      ); \
  template void Kokkos::ThrustGPUNode::parallel_for<Kokkos::DefaultJacobiOp2<SCALAR, ORDINAL>                       >(int, int, Kokkos::DefaultJacobiOp2<SCALAR, ORDINAL>                      ); \
  template void Kokkos::ThrustGPUNode::parallel_for<Kokkos::DefaultCoarseGrainHybridGaussSeidelOp1<SCALAR, ORDINAL> >(int, int, Kokkos::DefaultCoarseGrainHybridGaussSeidelOp1<SCALAR, ORDINAL>); \
  template void Kokkos::ThrustGPUNode::parallel_for<Kokkos::DefaultCoarseGrainHybridGaussSeidelOp2<SCALAR, ORDINAL> >(int, int, Kokkos::DefaultCoarseGrainHybridGaussSeidelOp2<SCALAR, ORDINAL>); \
  template void Kokkos::ThrustGPUNode::parallel_for<Kokkos::DefaultFineGrainHybridGaussSeidelOp1<SCALAR, ORDINAL>   >(int, int, Kokkos::DefaultFineGrainHybridGaussSeidelOp1<SCALAR, ORDINAL>  ); \
  template void Kokkos::ThrustGPUNode::parallel_for<Kokkos::DefaultFineGrainHybridGaussSeidelOp2<SCALAR, ORDINAL>   >(int, int, Kokkos::DefaultFineGrainHybridGaussSeidelOp2<SCALAR, ORDINAL>  ); \
  template void Kokkos::ThrustGPUNode::parallel_for<Kokkos::ExtractDiagonalOp1<SCALAR, ORDINAL>                     >(int, int, Kokkos::ExtractDiagonalOp1<SCALAR, ORDINAL>                    ); \
  template void Kokkos::ThrustGPUNode::parallel_for<Kokkos::ExtractDiagonalOp2<SCALAR, ORDINAL>                     >(int, int, Kokkos::ExtractDiagonalOp2<SCALAR, ORDINAL>                    );

#ifdef HAVE_KOKKOS_CUDA_FLOAT
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(int,float)
typedef short int ShortInt; INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,float)
#endif

#ifdef HAVE_KOKKOS_CUDA_DOUBLE
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(int,double)
typedef short int ShortInt; INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,double)
#endif

INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(int,int)
typedef short int ShortInt; INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,int)

