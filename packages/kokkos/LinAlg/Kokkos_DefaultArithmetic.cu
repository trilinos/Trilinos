/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <cuComplex.h>

// include for ThrustGPUNode method implementations
#include "Kokkos_ThrustGPUNode.cuh"

// includes for all ops
#include "Kokkos_MultiVectorKernelOps.hpp"

#ifdef HAVE_KOKKOS_CUSP
#include <cusp/complex.h>
#endif


#define INSTANTIATE_MULTIVECTOR_SCALAR(SCALAR) \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::InitOp       <SCALAR> >(int, int, Kokkos::InitOp       <SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::AssignOp     <SCALAR> >(int, int, Kokkos::AssignOp     <SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::SingleScaleOp<SCALAR> >(int, int, Kokkos::SingleScaleOp<SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::MVScaleOp    <SCALAR> >(int, int, Kokkos::MVScaleOp    <SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::AbsOp        <SCALAR> >(int, int, Kokkos::AbsOp        <SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::RecipOp      <SCALAR> >(int, int, Kokkos::RecipOp      <SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::GESUMOp      <SCALAR> >(int, int, Kokkos::GESUMOp      <SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::GESUMOp3     <SCALAR> >(int, int, Kokkos::GESUMOp3     <SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::MVElemMultOp <SCALAR> >(int, int, Kokkos::MVElemMultOp <SCALAR> ); \
  template Kokkos::SumAbsOp    <SCALAR>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::SumAbsOp    <SCALAR> >(int, int, Kokkos::SumAbsOp    <SCALAR> ); \
  template Kokkos::WeightNormOp<SCALAR>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::WeightNormOp<SCALAR> >(int, int, Kokkos::WeightNormOp<SCALAR> ); \
  template Kokkos::SumOp       <SCALAR>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::SumOp       <SCALAR> >(int, int, Kokkos::SumOp       <SCALAR> ); \
  template Kokkos::MaxAbsOp    <SCALAR>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::MaxAbsOp    <SCALAR> >(int, int, Kokkos::MaxAbsOp    <SCALAR> ); \
  template Kokkos::DotOp1      <SCALAR>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::DotOp1      <SCALAR> >(int, int, Kokkos::DotOp1      <SCALAR> ); \
  template Kokkos::DotOp2      <SCALAR>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce< Kokkos::DotOp2      <SCALAR> >(int, int, Kokkos::DotOp2      <SCALAR> );

INSTANTIATE_MULTIVECTOR_SCALAR(int)
INSTANTIATE_MULTIVECTOR_SCALAR(long)

#ifdef HAVE_KOKKOS_CUDA_FLOAT
INSTANTIATE_MULTIVECTOR_SCALAR(float)
#endif
#ifdef HAVE_KOKKOS_CUDA_DOUBLE
INSTANTIATE_MULTIVECTOR_SCALAR(double)
#endif
//typedef cusp::complex<float>  ComplexFloat;
//typedef cusp::complex<double> ComplexDouble;
//#ifdef HAVE_KOKKOS_CUDA_COMPLEX_FLOAT
//INSTANTIATE_MULTIVECTOR_SCALAR(ComplexFloat)
//#endif
//#ifdef HAVE_KOKKOS_CUDA_COMPLEX_DOUBLE
//INSTANTIATE_MULTIVECTOR_SCALAR(ComplexDouble)
//#endif
