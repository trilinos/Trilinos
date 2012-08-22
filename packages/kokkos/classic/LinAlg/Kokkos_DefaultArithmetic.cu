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

#ifdef HAVE_KOKKOSCLASSIC_CUSP
#include <cusp/complex.h>
#endif


#define INSTANTIATE_MULTIVECTOR_SCALAR(SCALAR) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_FOR( Kokkos::InitOp       <SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_FOR( Kokkos::AssignOp     <SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_FOR( Kokkos::SingleScaleOp<SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_FOR( Kokkos::MVScaleOp    <SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_FOR( Kokkos::AbsOp        <SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_FOR( Kokkos::RecipOp      <SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_FOR( Kokkos::GESUMOp      <SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_FOR( Kokkos::GESUMOp3     <SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_FOR( Kokkos::MVElemMultOp <SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_RED( Kokkos::SumAbsOp     <SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_RED( Kokkos::WeightNormOp <SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_RED( Kokkos::SumOp        <SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_RED( Kokkos::MaxAbsOp     <SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_RED( Kokkos::DotOp1       <SCALAR> ) \
  KOKKOS_INSTANT_THRUSTGPUNODE_PARALLEL_RED( Kokkos::DotOp2       <SCALAR> )

INSTANTIATE_MULTIVECTOR_SCALAR(int)
INSTANTIATE_MULTIVECTOR_SCALAR(long)

#ifdef HAVE_KOKKOSCLASSIC_CUDA_FLOAT
INSTANTIATE_MULTIVECTOR_SCALAR(float)
#endif
#ifdef HAVE_KOKKOSCLASSIC_CUDA_DOUBLE
INSTANTIATE_MULTIVECTOR_SCALAR(double)
#endif
//typedef cusp::complex<float>  ComplexFloat;
//typedef cusp::complex<double> ComplexDouble;
//#ifdef HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_FLOAT
//INSTANTIATE_MULTIVECTOR_SCALAR(ComplexFloat)
//#endif
//#ifdef HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_DOUBLE
//INSTANTIATE_MULTIVECTOR_SCALAR(ComplexDouble)
//#endif
