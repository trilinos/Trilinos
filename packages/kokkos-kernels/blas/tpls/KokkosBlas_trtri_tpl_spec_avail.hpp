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

#ifndef KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {

// Specialization struct which defines whether a specialization exists
template <class RVT, class AVT>
struct trtri_tpl_spec_avail {
  enum : bool { value = false };
};

// Generic Host side LAPACK (could be MKL or whatever)
#define KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL(SCALAR, LAYOUTA, MEMSPACE)         \
  template <class ExecSpace>                                               \
  struct trtri_tpl_spec_avail<                                             \
      Kokkos::View<int, Kokkos::LayoutRight, Kokkos::HostSpace,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
      Kokkos::View<SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {           \
    enum : bool { value = true };                                          \
  };

#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#define KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUTA, MEMSPACE) \
  KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL(SCALAR, LAYOUTA, MEMSPACE)
#else
#define KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUTA, MEMSPACE)
#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#define KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(SCALAR, LAYOUTA, MEMSPACE) \
  KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL(SCALAR, LAYOUTA, MEMSPACE)
#else
#define KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(SCALAR, LAYOUTA, MEMSPACE)
#endif  // KOKKOSKERNELS_ENABLE_TPL_MAGMA

KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft,
                                     Kokkos::HostSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(double, Kokkos::LayoutLeft,
                                      Kokkos::CudaSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(double, Kokkos::LayoutLeft,
                                      Kokkos::CudaUVMSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft,
                                     Kokkos::HostSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(float, Kokkos::LayoutLeft,
                                      Kokkos::CudaSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(float, Kokkos::LayoutLeft,
                                      Kokkos::CudaUVMSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>,
                                     Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<double>,
                                      Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<double>,
                                      Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft,
                                     Kokkos::HostSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<float>,
                                      Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<float>,
                                      Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)

KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutRight,
                                     Kokkos::HostSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(double, Kokkos::LayoutRight,
                                      Kokkos::CudaSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(double, Kokkos::LayoutRight,
                                      Kokkos::CudaUVMSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutRight,
                                     Kokkos::HostSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(float, Kokkos::LayoutRight,
                                      Kokkos::CudaSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(float, Kokkos::LayoutRight,
                                      Kokkos::CudaUVMSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>,
                                     Kokkos::LayoutRight, Kokkos::HostSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<double>,
                                      Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<double>,
                                      Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>,
                                     Kokkos::LayoutRight, Kokkos::HostSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<float>,
                                      Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_MAGMA(Kokkos::complex<float>,
                                      Kokkos::LayoutRight, Kokkos::CudaUVMSpace)

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS_TRTRI_TPL_SPEC_AVAIL_HPP_
