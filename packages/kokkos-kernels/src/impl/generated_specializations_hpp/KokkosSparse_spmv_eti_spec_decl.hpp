#ifndef KOKKOSSPARSE_SPMV_ETI_SPEC_DECL_HPP_
#define KOKKOSSPARSE_SPMV_ETI_SPEC_DECL_HPP_
/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

namespace KokkosSparse {
namespace Impl {

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(double, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(float, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_CUDA) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HOSTSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_OPENMP) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_THREADS) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Threads, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, int, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
 && defined (KOKKOSKERNELS_INST_EXECSPACE_SERIAL) \
 && defined (KOKKOSKERNELS_INST_MEMSPACE_HBWSPACE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) 
 KOKKOSSPARSE_SPMV_ETI_SPEC_DECL(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::Experimental::HBWSpace)
#endif
} // Impl
} // KokkosSparse
#endif // KOKKOSSPARSE_SPMV_ETI_SPEC_DECL_HPP_
