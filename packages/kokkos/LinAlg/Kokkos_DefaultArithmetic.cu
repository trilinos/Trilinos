// include for ThrustGPUNode method implementations
#include "Kokkos_ThrustGPUNode.cuh"

// includes for all ops
#include "Kokkos_MultiVectorKernelOps.hpp"

template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::InitOp       <float> >(int, int, Kokkos::InitOp       <float> );
template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::AssignOp     <float> >(int, int, Kokkos::AssignOp     <float> );
template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::SingleScaleOp<float> >(int, int, Kokkos::SingleScaleOp<float> );
template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::MVScaleOp    <float> >(int, int, Kokkos::MVScaleOp    <float> );
template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::AbsOp        <float> >(int, int, Kokkos::AbsOp        <float> );
template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::RecipOp      <float> >(int, int, Kokkos::RecipOp      <float> );
template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::GESUMOp      <float> >(int, int, Kokkos::GESUMOp      <float> );
template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::GESUMOp3     <float> >(int, int, Kokkos::GESUMOp3     <float> );

template Kokkos::SumAbsOp    <float>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::SumAbsOp    <float> >(int, int, Kokkos::SumAbsOp    <float> );
template Kokkos::WeightNormOp<float>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::WeightNormOp<float> >(int, int, Kokkos::WeightNormOp<float> );
template Kokkos::SumOp       <float>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::SumOp       <float> >(int, int, Kokkos::SumOp       <float> );
template Kokkos::MaxAbsOp    <float>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::MaxAbsOp    <float> >(int, int, Kokkos::MaxAbsOp    <float> );
template Kokkos::DotOp1      <float>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::DotOp1      <float> >(int, int, Kokkos::DotOp1      <float> );
template Kokkos::DotOp2      <float>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::DotOp2      <float> >(int, int, Kokkos::DotOp2      <float> );
