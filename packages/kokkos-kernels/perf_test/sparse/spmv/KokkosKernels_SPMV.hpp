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

#ifndef KOKKOSKERNELS_SPMV_HPP_
#define KOKKOSKERNELS_SPMV_HPP_

#ifdef MAKE_BUILD
#ifdef KOKKOS_ENABLE_CUDA
  #define KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS()  \
        typedef Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space> Kokkos_Device0Kokkos_Cuda_Kokkos_CudaSpace0; \
        typedef Kokkos::complex<double> Kokkos_complex0double0; \
        typedef long long longlong;
#else
  #ifdef KOKKOS_ENABLE_OPENMP
    #define KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS()  \
        typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::OpenMP::memory_space> Kokkos_Device0Kokkos_OpenMP_Kokkos_HostSpace0; \
        typedef Kokkos::complex<double> Kokkos_complex0double0; \
        typedef long long longlong;
  #else
    #ifdef KOKKOS_ENABLE_THREADS
      #define KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS()  \
        typedef Kokkos::Device<Kokkos::Threads, Kokkos::Threads::memory_space> Kokkos_Device0Kokkos_Threads_Kokkos_HostSpace0; \
        typedef Kokkos::complex<double> Kokkos_complex0double0; \
        typedef long long longlong;
    #else
      #define KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS()  \
        typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> Kokkos_Device0Kokkos_OpenMP_Kokkos_HostSpace0; \
        typedef Kokkos::complex<double> Kokkos_complex0double0; \
        typedef long long longlong;
    #endif
  #endif
#endif

#endif

#include <KokkosBlas.hpp>
#include <KokkosSparse_spmv.hpp>

#ifdef HAVE_KK_KERNELS


template<typename AType, typename XType, typename YType>
void kokkoskernels_matvec(AType A, XType x, YType y, int rows_per_thread, int team_size, int vector_length) {
  KokkosSparse::spmv (KokkosSparse::NoTranspose,1.0,A,x,0.0,y);
}
#endif

#endif /* KOKKOSKERNELS_SPMV_HPP_ */
