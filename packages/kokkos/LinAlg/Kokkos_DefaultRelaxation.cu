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

