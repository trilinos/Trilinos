//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSKERNELS_SPMV_HPP_
#define KOKKOSKERNELS_SPMV_HPP_

#ifdef MAKE_BUILD
#ifdef KOKKOS_ENABLE_CUDA
#define KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS()                                                                   \
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space> Kokkos_Device0Kokkos_Cuda_Kokkos_CudaSpace0; \
  typedef Kokkos::complex<double> Kokkos_complex0double0;                                                       \
  typedef long long longlong;
#else
#ifdef KOKKOS_ENABLE_OPENMP
#define KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS()                                                                         \
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::OpenMP::memory_space> Kokkos_Device0Kokkos_OpenMP_Kokkos_HostSpace0; \
  typedef Kokkos::complex<double> Kokkos_complex0double0;                                                             \
  typedef long long longlong;
#else
#ifdef KOKKOS_ENABLE_THREADS
#define KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS()                            \
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::Threads::memory_space> \
      Kokkos_Device0Kokkos_Threads_Kokkos_HostSpace0;                    \
  typedef Kokkos::complex<double> Kokkos_complex0double0;                \
  typedef long long longlong;
#else
#define KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS()                                                              \
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> Kokkos_Device0Kokkos_OpenMP_Kokkos_HostSpace0; \
  typedef Kokkos::complex<double> Kokkos_complex0double0;                                                  \
  typedef long long longlong;
#endif
#endif
#endif

#endif

#include <KokkosBlas.hpp>
#include <KokkosSparse_spmv.hpp>

#ifdef HAVE_KK_KERNELS

template <typename AType, typename XType, typename YType>
void kokkoskernels_matvec(AType A, XType x, YType y, int rows_per_thread, int team_size, int vector_length) {
  KokkosSparse::spmv(KokkosSparse::NoTranspose, 1.0, A, x, 0.0, y);
}
#endif

#endif /* KOKKOSKERNELS_SPMV_HPP_ */
