// include for ThrustGPUNode method implementations
#include "Kokkos_ThrustGPUNode.cuh"

// includes for all ops
#include "Kokkos_MultiVectorKernelOps.hpp"


#define INSTANTIATE_MULTIVECTOR_SCALAR(SCALAR) \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::InitOp       <SCALAR> >(int, int, Kokkos::InitOp       <SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::AssignOp     <SCALAR> >(int, int, Kokkos::AssignOp     <SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::SingleScaleOp<SCALAR> >(int, int, Kokkos::SingleScaleOp<SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::MVScaleOp    <SCALAR> >(int, int, Kokkos::MVScaleOp    <SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::AbsOp        <SCALAR> >(int, int, Kokkos::AbsOp        <SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::RecipOp      <SCALAR> >(int, int, Kokkos::RecipOp      <SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::GESUMOp      <SCALAR> >(int, int, Kokkos::GESUMOp      <SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::GESUMOp3     <SCALAR> >(int, int, Kokkos::GESUMOp3     <SCALAR> ); \
  template Kokkos::SumAbsOp    <SCALAR>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::SumAbsOp    <SCALAR> >(int, int, Kokkos::SumAbsOp    <SCALAR> ); \
  template Kokkos::WeightNormOp<SCALAR>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::WeightNormOp<SCALAR> >(int, int, Kokkos::WeightNormOp<SCALAR> ); \
  template Kokkos::SumOp       <SCALAR>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::SumOp       <SCALAR> >(int, int, Kokkos::SumOp       <SCALAR> ); \
  template Kokkos::MaxAbsOp    <SCALAR>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::MaxAbsOp    <SCALAR> >(int, int, Kokkos::MaxAbsOp    <SCALAR> ); \
  template Kokkos::DotOp1      <SCALAR>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::DotOp1      <SCALAR> >(int, int, Kokkos::DotOp1      <SCALAR> ); \
  template Kokkos::DotOp2      <SCALAR>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::DotOp2      <SCALAR> >(int, int, Kokkos::DotOp2      <SCALAR> );

INSTANTIATE_MULTIVECTOR_SCALAR(float)
INSTANTIATE_MULTIVECTOR_SCALAR(int)
