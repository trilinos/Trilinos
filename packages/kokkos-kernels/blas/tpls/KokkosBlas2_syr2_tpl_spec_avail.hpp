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

#ifndef KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class EXEC_SPACE, class XT, class YT, class AT>
struct syr2_tpl_spec_avail {
  enum : bool { value = false };
};

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS

#define KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                     \
  template <>                                                                                           \
  struct syr2_tpl_spec_avail<EXEC_SPACE,                                                                \
                             Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
                             Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
                             Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,      \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                 \
    enum : bool { value = true };                                                                       \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)

KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)

KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#define KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                   \
  template <>                                                                                           \
  struct syr2_tpl_spec_avail<EXEC_SPACE,                                                                \
                             Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
                             Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
                             Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,      \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                 \
    enum : bool { value = true };                                                                       \
  };

KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)

KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)

KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace)

KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace)

#endif

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS

#define KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_ROCBLAS(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE)                  \
  template <>                                                                                           \
  struct syr2_tpl_spec_avail<EXEC_SPACE,                                                                \
                             Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
                             Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                    \
                             Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,      \
                                          Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                 \
    enum : bool { value = true };                                                                       \
  };

KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_ROCBLAS(double, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_ROCBLAS(float, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace)

KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_ROCBLAS(double, Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_ROCBLAS(float, Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace)
KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace)

#endif
}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS2_SYR2_TPL_SPEC_AVAIL_HPP_
