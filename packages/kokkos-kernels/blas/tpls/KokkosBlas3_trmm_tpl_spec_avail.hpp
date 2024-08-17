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

#ifndef KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {

// Specialization struct which defines whether a specialization exists
template <class execution_space, class AVT, class BVT>
struct trmm_tpl_spec_avail {
  enum : bool { value = false };
};

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS

#define KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUTA, LAYOUTB, MEMSPACE)                        \
  template <class ExecSpace>                                                                            \
  struct trmm_tpl_spec_avail<ExecSpace,                                                                 \
                             Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEMSPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
                             Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEMSPACE>,       \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                 \
    enum : bool { value = true };                                                                       \
  };

KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace)

KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                     Kokkos::HostSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                     Kokkos::HostSpace)

#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#define KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(SCALAR, LAYOUTA, LAYOUTB, MEMSPACE)                      \
  template <class ExecSpace>                                                                            \
  struct trmm_tpl_spec_avail<Kokkos::Cuda,                                                              \
                             Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEMSPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
                             Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEMSPACE>,       \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                 \
    enum : bool { value = true };                                                                       \
  };

KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       Kokkos::CudaSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       Kokkos::CudaUVMSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       Kokkos::CudaSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       Kokkos::CudaUVMSpace)

KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       Kokkos::CudaSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       Kokkos::CudaUVMSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       Kokkos::CudaSpace)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       Kokkos::CudaUVMSpace)

#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS
}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_HPP_
