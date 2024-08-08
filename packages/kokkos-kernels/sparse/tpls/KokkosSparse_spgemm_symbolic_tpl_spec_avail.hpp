/*
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
*/

#ifndef KOKKOSPARSE_SPGEMM_SYMBOLIC_TPL_SPEC_AVAIL_HPP_
#define KOKKOSPARSE_SPGEMM_SYMBOLIC_TPL_SPEC_AVAIL_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include <mkl.h>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class b_size_view_t_, class b_lno_view_t,
          class c_size_view_t_>
struct spgemm_symbolic_tpl_spec_avail {
  enum : bool { value = false };
};

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
// NOTE: all versions of cuSPARSE 10.x and 11.x support exactly the same matrix
// types, so there is no ifdef'ing on versions needed in avail. Offset and
// Ordinal must both be 32-bit. Even though the "generic" API lets you specify
// offsets and ordinals independently as either 16, 32 or 64-bit, SpGEMM will
// just fail at runtime if you don't use 32 for both.

#define SPGEMM_SYMBOLIC_AVAIL_CUSPARSE(SCALAR, MEMSPACE)                                                           \
  template <>                                                                                                      \
  struct spgemm_symbolic_tpl_spec_avail<                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const int, const int, const SCALAR, Kokkos::Cuda, MEMSPACE, \
                                                       MEMSPACE>,                                                  \
      Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                      \
      Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                      \
      Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                      \
      Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                      \
      Kokkos::View<int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                                   \
    enum : bool { value = true };                                                                                  \
  };

#define SPGEMM_SYMBOLIC_AVAIL_CUSPARSE_S(SCALAR)            \
  SPGEMM_SYMBOLIC_AVAIL_CUSPARSE(SCALAR, Kokkos::CudaSpace) \
  SPGEMM_SYMBOLIC_AVAIL_CUSPARSE(SCALAR, Kokkos::CudaUVMSpace)

#if (CUDA_VERSION < 11000) || (CUDA_VERSION >= 11040)
SPGEMM_SYMBOLIC_AVAIL_CUSPARSE_S(float)
SPGEMM_SYMBOLIC_AVAIL_CUSPARSE_S(double)
SPGEMM_SYMBOLIC_AVAIL_CUSPARSE_S(Kokkos::complex<float>)
SPGEMM_SYMBOLIC_AVAIL_CUSPARSE_S(Kokkos::complex<double>)
#endif
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
#define SPGEMM_SYMBOLIC_AVAIL_ROCSPARSE(SCALAR)                                                         \
  template <>                                                                                           \
  struct spgemm_symbolic_tpl_spec_avail<                                                                \
      KokkosKernels::Experimental::KokkosKernelsHandle<const int, const int, const SCALAR, Kokkos::HIP, \
                                                       Kokkos::HIPSpace, Kokkos::HIPSpace>,             \
      Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                           \
      Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                           \
      Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                           \
      Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                           \
      Kokkos::View<int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                        \
    enum : bool { value = true };                                                                       \
  };

SPGEMM_SYMBOLIC_AVAIL_ROCSPARSE(float)
SPGEMM_SYMBOLIC_AVAIL_ROCSPARSE(double)
SPGEMM_SYMBOLIC_AVAIL_ROCSPARSE(Kokkos::complex<float>)
SPGEMM_SYMBOLIC_AVAIL_ROCSPARSE(Kokkos::complex<double>)
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#define SPGEMM_SYMBOLIC_AVAIL_MKL(SCALAR, EXEC)                                                          \
  template <>                                                                                            \
  struct spgemm_symbolic_tpl_spec_avail<                                                                 \
      KokkosKernels::Experimental::KokkosKernelsHandle<const MKL_INT, const MKL_INT, const SCALAR, EXEC, \
                                                       Kokkos::HostSpace, Kokkos::HostSpace>,            \
      Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                            \
      Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                            \
      Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                            \
      Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                            \
      Kokkos::View<MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                         \
    enum : bool { value = true };                                                                        \
  };

#define SPGEMM_SYMBOLIC_AVAIL_MKL_E(EXEC)                 \
  SPGEMM_SYMBOLIC_AVAIL_MKL(float, EXEC)                  \
  SPGEMM_SYMBOLIC_AVAIL_MKL(double, EXEC)                 \
  SPGEMM_SYMBOLIC_AVAIL_MKL(Kokkos::complex<float>, EXEC) \
  SPGEMM_SYMBOLIC_AVAIL_MKL(Kokkos::complex<double>, EXEC)

#ifdef KOKKOS_ENABLE_SERIAL
SPGEMM_SYMBOLIC_AVAIL_MKL_E(Kokkos::Serial)
#endif
#ifdef KOKKOS_ENABLE_OPENMP
SPGEMM_SYMBOLIC_AVAIL_MKL_E(Kokkos::OpenMP)
#endif
#endif  // KOKKOSKERNELS_ENABLE_TPL_MKL

}  // namespace Impl
}  // namespace KokkosSparse

#endif
