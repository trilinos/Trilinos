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

#include <cuComplex.h>

#include "KokkosClassic_config.h"
#include "Kokkos_CuspOps.cuh"

// includes for all ops
#include "Kokkos_MultiVectorKernelOps.hpp"

// cusp doesn't currently support mixed precision, but maybe it will one day...
#define INSTANTIATE_CUSP_ORDINAL_SCALAR(OFFSET,ORDINAL,SCALARA,SCALARX,SCALARY)                       \
  template void Kokkos::Cuspdetails::cuspCrsMultiply<OFFSET,ORDINAL,SCALARA,SCALARX,SCALARY>          \
                       ( ORDINAL, ORDINAL, ORDINAL, const OFFSET *, const ORDINAL *, const SCALARA *, \
                         ORDINAL, const SCALARX *, ORDINAL, SCALARY *, ORDINAL );                     \
  template void Kokkos::Cuspdetails::cuspCrsTranspose<OFFSET,ORDINAL,SCALARA>(ORDINAL, ORDINAL, ORDINAL, \
                        const OFFSET *, const ORDINAL *, const SCALARA *,                                \
                        OFFSET *, ORDINAL *, SCALARA *);

#ifdef HAVE_KOKKOSCLASSIC_CUDA_FLOAT
INSTANTIATE_CUSP_ORDINAL_SCALAR(short,short,float,float,float)
INSTANTIATE_CUSP_ORDINAL_SCALAR(int,int,float,float,float)
#endif
#ifdef HAVE_KOKKOSCLASSIC_CUDA_DOUBLE
INSTANTIATE_CUSP_ORDINAL_SCALAR(short,short,double,double,double)
INSTANTIATE_CUSP_ORDINAL_SCALAR(int,int,double,double,double)
#endif
//typedef cusp::complex<float>  ComplexFloat;
//typedef cusp::complex<double> ComplexDouble;
//#ifdef HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_FLOAT
//INSTANTIATE_CUSP_ORDINAL_SCALAR(short,ComplexFloat)
//INSTANTIATE_CUSP_ORDINAL_SCALAR(int,ComplexFloat)
//#endif
//#ifdef HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_DOUBLE
//INSTANTIATE_CUSP_ORDINAL_SCALAR(short,ComplexDouble)
//INSTANTIATE_CUSP_ORDINAL_SCALAR(int,ComplexDouble)
//#endif
