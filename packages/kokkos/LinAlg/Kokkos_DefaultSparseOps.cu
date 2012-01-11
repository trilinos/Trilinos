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


// include for ThrustGPUNode method implementations
#include "Kokkos_ThrustGPUNode.cuh"

// includes for all ops
#include "Kokkos_DefaultSparseMultiplyKernelOps.hpp"
#include "Kokkos_DefaultSparseSolveKernelOps.hpp"

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

typedef short int ShortInt;

#ifdef HAVE_KOKKOS_CUDA_FLOAT
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(int,float)
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR_SCALAR(int,int,float)
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,float)
#endif

#ifdef HAVE_KOKKOS_CUDA_DOUBLE
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(int,double)
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR_SCALAR(int,int,double)
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,double)
#endif

#if defined(HAVE_KOKKOS_CUDA_DOUBLE) && defined(HAVE_KOKKOS_CUDA_FLOAT)
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR_SCALAR(int,float,double)
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR_SCALAR(int,double,float)
#endif

INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(int,int)
INSTANTIATE_SPARSEMULTIPLY_ORDINAL_SCALAR(ShortInt,int)

#define INSTANTIATE_SPARSESOLVE_ORDINAL_SCALAR(ORDINAL, SCALAR) \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseSolveOp1            <SCALAR, ORDINAL, SCALAR, SCALAR> >(int, int, Kokkos::DefaultSparseSolveOp1            <SCALAR, ORDINAL, SCALAR, SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseSolveOp2            <SCALAR, ORDINAL, SCALAR, SCALAR> >(int, int, Kokkos::DefaultSparseSolveOp2            <SCALAR, ORDINAL, SCALAR, SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeSolveOp1   <SCALAR, ORDINAL, SCALAR, SCALAR> >(int, int, Kokkos::DefaultSparseTransposeSolveOp1   <SCALAR, ORDINAL, SCALAR, SCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeSolveOp2   <SCALAR, ORDINAL, SCALAR, SCALAR> >(int, int, Kokkos::DefaultSparseTransposeSolveOp2   <SCALAR, ORDINAL, SCALAR, SCALAR> );

#define INSTANTIATE_SPARSESOLVE_ORDINAL_SCALAR_SCALAR(ORDINAL, SCALAR, DRSCALAR) \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseSolveOp1            <SCALAR, ORDINAL, DRSCALAR, DRSCALAR> >(int, int, Kokkos::DefaultSparseSolveOp1            <SCALAR, ORDINAL, DRSCALAR, DRSCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseSolveOp2            <SCALAR, ORDINAL, DRSCALAR, DRSCALAR> >(int, int, Kokkos::DefaultSparseSolveOp2            <SCALAR, ORDINAL, DRSCALAR, DRSCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeSolveOp1   <SCALAR, ORDINAL, DRSCALAR, DRSCALAR> >(int, int, Kokkos::DefaultSparseTransposeSolveOp1   <SCALAR, ORDINAL, DRSCALAR, DRSCALAR> ); \
  template void Kokkos::ThrustGPUNode::parallel_for< Kokkos::DefaultSparseTransposeSolveOp2   <SCALAR, ORDINAL, DRSCALAR, DRSCALAR> >(int, int, Kokkos::DefaultSparseTransposeSolveOp2   <SCALAR, ORDINAL, DRSCALAR, DRSCALAR> );

#ifdef HAVE_KOKKOS_CUDA_FLOAT
INSTANTIATE_SPARSESOLVE_ORDINAL_SCALAR(int,float)
INSTANTIATE_SPARSESOLVE_ORDINAL_SCALAR_SCALAR(int,int,float)
typedef short int ShortInt; INSTANTIATE_SPARSESOLVE_ORDINAL_SCALAR(ShortInt,float)
#endif

#ifdef HAVE_KOKKOS_CUDA_DOUBLE
INSTANTIATE_SPARSESOLVE_ORDINAL_SCALAR(int,double)
INSTANTIATE_SPARSESOLVE_ORDINAL_SCALAR_SCALAR(int,int,double)
typedef short int ShortInt; INSTANTIATE_SPARSESOLVE_ORDINAL_SCALAR(ShortInt,double)
#endif

#if defined(HAVE_KOKKOS_CUDA_DOUBLE) && defined(HAVE_KOKKOS_CUDA_FLOAT)
INSTANTIATE_SPARSESOLVE_ORDINAL_SCALAR_SCALAR(int,float,double)
INSTANTIATE_SPARSESOLVE_ORDINAL_SCALAR_SCALAR(int,double,float)
#endif

INSTANTIATE_SPARSESOLVE_ORDINAL_SCALAR(int,int)
typedef short int ShortInt; INSTANTIATE_SPARSESOLVE_ORDINAL_SCALAR(ShortInt,int)
