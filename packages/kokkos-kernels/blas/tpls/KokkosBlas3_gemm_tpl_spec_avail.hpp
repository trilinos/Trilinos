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

#ifndef KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class AT, class XT, class YT>
struct gemm_tpl_spec_avail {
  enum : bool { value = false };
};

// Generic Host side BLAS (could be MKL or whatever)
#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)

#define KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUTA, LAYOUTB, LAYOUTC, MEMSPACE)               \
  template <class ExecSpace>                                                                            \
  struct gemm_tpl_spec_avail<ExecSpace,                                                                 \
                             Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEMSPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
                             Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEMSPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
                             Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<ExecSpace, MEMSPACE>,       \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                 \
    enum : bool { value = true };                                                                       \
  };

KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                     Kokkos::HostSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                     Kokkos::HostSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                     Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                     Kokkos::HostSpace)

KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                     Kokkos::HostSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                     Kokkos::HostSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                     Kokkos::LayoutRight, Kokkos::HostSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                     Kokkos::LayoutRight, Kokkos::HostSpace)

#endif

// cuBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)

#define KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(SCALAR, LAYOUTA, LAYOUTB, LAYOUTC, MEMSPACE)             \
  template <class ExecSpace>                                                                            \
  struct gemm_tpl_spec_avail<ExecSpace,                                                                 \
                             Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEMSPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
                             Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEMSPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
                             Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<ExecSpace, MEMSPACE>,       \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                 \
    enum : bool { value = true };                                                                       \
  };

KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       Kokkos::CudaSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       Kokkos::CudaUVMSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       Kokkos::CudaSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       Kokkos::CudaUVMSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)

KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       Kokkos::CudaSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       Kokkos::CudaUVMSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       Kokkos::CudaSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       Kokkos::CudaUVMSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       Kokkos::LayoutRight, Kokkos::CudaUVMSpace)

#endif

// rocBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)

#define KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(SCALAR, LAYOUT, MEMSPACE)                              \
  template <class ExecSpace>                                                                           \
  struct gemm_tpl_spec_avail<ExecSpace,                                                                \
                             Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                   \
                             Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                   \
                             Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,       \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                \
    enum : bool { value = true };                                                                      \
  };

KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(double, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(float, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPSpace)

KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(double, Kokkos::LayoutRight, Kokkos::HIPSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(float, Kokkos::LayoutRight, Kokkos::HIPSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::HIPSpace)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::HIPSpace)

#endif
}  // namespace Impl
}  // namespace KokkosBlas

#endif
