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

#ifndef KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class AV, class XMV, class YMV, int Xrank = XMV::rank, int Yrank = YMV::rank>
struct dot_tpl_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

namespace KokkosBlas {
namespace Impl {

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
// double
#define KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUT, MEMSPACE)                                                  \
  template <class ExecSpace>                                                                                           \
  struct dot_tpl_spec_avail<ExecSpace,                                                                                 \
                            Kokkos::View<SCALAR, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                            Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                   \
                                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                    \
                            Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                   \
                                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                    \
                            1, 1> {                                                                                    \
    enum : bool { value = true };                                                                                      \
  };

KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft, Kokkos::HostSpace)

// TODO: we met difficuties in FindTPLMKL.cmake to set the BLAS library properly
// such that the test in CheckHostBlasReturnComplex.cmake could not be
// compiled and run to give a correct answer on KK_BLAS_RESULT_AS_POINTER_ARG.
// This resulted in segfault in dot() with MKL and complex.
// So we just temporarily disable it until FindTPLMKL.cmake is fixed.
#if !defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HostSpace)
#endif

#endif

#define KOKKOSBLAS1_DOT_TPL_SPEC(SCALAR, LAYOUT, EXECSPACE, MEMSPACE)                                                  \
  template <>                                                                                                          \
  struct dot_tpl_spec_avail<EXECSPACE,                                                                                 \
                            Kokkos::View<SCALAR, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                            Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                   \
                                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                    \
                            Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                   \
                                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                    \
                            1, 1> {                                                                                    \
    enum : bool { value = true };                                                                                      \
  };

#define KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL(LAYOUT, EXECSPACE, MEMSPACE)             \
  KOKKOSBLAS1_DOT_TPL_SPEC(float, LAYOUT, EXECSPACE, MEMSPACE)                  \
  KOKKOSBLAS1_DOT_TPL_SPEC(double, LAYOUT, EXECSPACE, MEMSPACE)                 \
  KOKKOSBLAS1_DOT_TPL_SPEC(Kokkos::complex<float>, LAYOUT, EXECSPACE, MEMSPACE) \
  KOKKOSBLAS1_DOT_TPL_SPEC(Kokkos::complex<double>, LAYOUT, EXECSPACE, MEMSPACE)

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
// Note BMK: CUBLAS dot is consistently slower than our native dot
// (measured 11.2, 11.8, 12.0 using perf test, and all are similar)
// If a future version improves performance, re-enable it here and
// in the tpl_spec_decl file.
#if 0
KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL(Kokkos::LayoutLeft, Kokkos::Cuda,
                               Kokkos::CudaSpace)
#endif
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace)
#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_MKL) && defined(KOKKOS_ENABLE_SYCL)
KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL(Kokkos::LayoutLeft, Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLDeviceUSMSpace)
#endif
}  // namespace Impl
}  // namespace KokkosBlas
#endif
