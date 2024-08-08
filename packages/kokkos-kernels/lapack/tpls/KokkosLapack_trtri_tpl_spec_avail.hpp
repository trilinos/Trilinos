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

#ifndef KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_HPP_
#define KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_HPP_

namespace KokkosLapack {
namespace Impl {

// Specialization struct which defines whether a specialization exists
template <class RVT, class AVT>
struct trtri_tpl_spec_avail {
  enum : bool { value = false };
};

// Generic Host side LAPACK (could be MKL or whatever)
#define KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL(SCALAR, LAYOUTA, MEMSPACE)                                       \
  template <class ExecSpace>                                                                               \
  struct trtri_tpl_spec_avail<                                                                             \
      Kokkos::View<int, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEMSPACE>,                                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                           \
    enum : bool { value = true };                                                                          \
  };

#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK
#define KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_LAPACK(SCALAR, LAYOUTA, MEMSPACE) \
  KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL(SCALAR, LAYOUTA, MEMSPACE)
#else
#define KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_LAPACK(SCALAR, LAYOUTA, MEMSPACE)
#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK

#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#define KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(SCALAR, LAYOUTA, MEMSPACE) \
  KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL(SCALAR, LAYOUTA, MEMSPACE)
#else
#define KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(SCALAR, LAYOUTA, MEMSPACE)
#endif  // KOKKOSKERNELS_ENABLE_TPL_MAGMA

KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_LAPACK(double, Kokkos::LayoutLeft, Kokkos::HostSpace)
#ifdef KOKKOS_ENABLE_CUDA
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(double, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_LAPACK(float, Kokkos::LayoutLeft, Kokkos::HostSpace)
#ifdef KOKKOS_ENABLE_CUDA
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(float, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HostSpace)
#ifdef KOKKOS_ENABLE_CUDA
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HostSpace)
#ifdef KOKKOS_ENABLE_CUDA
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_LAPACK(double, Kokkos::LayoutRight, Kokkos::HostSpace)
#ifdef KOKKOS_ENABLE_CUDA
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(double, Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(double, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
#endif

KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_LAPACK(float, Kokkos::LayoutRight, Kokkos::HostSpace)
#ifdef KOKKOS_ENABLE_CUDA
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(float, Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(float, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
#endif

KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::HostSpace)
#ifdef KOKKOS_ENABLE_CUDA
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
#endif

KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_LAPACK(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::HostSpace)
#ifdef KOKKOS_ENABLE_CUDA
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSLAPACK_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
#endif

}  // namespace Impl
}  // namespace KokkosLapack

#endif  // KOKKOSLAPACKy_TRTRI_TPL_SPEC_AVAIL_HPP_
