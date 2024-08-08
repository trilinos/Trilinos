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

#ifndef KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class RV, class XMV, int Xrank = XMV::rank>
struct iamax_tpl_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

namespace KokkosBlas {
namespace Impl {

// Generic Host side BLAS (could be MKL or whatever)
#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)
// double
#define KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS(INDEX_TYPE, SCALAR, LAYOUT, MEMSPACE)                             \
  template <class ExecSpace>                                                                                    \
  struct iamax_tpl_spec_avail<                                                                                  \
      ExecSpace, Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                   \
      1> {                                                                                                      \
    enum : bool { value = true };                                                                               \
  };

KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS(unsigned long, double, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS(unsigned long, float, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS(unsigned long, Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS(unsigned long, Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HostSpace)

#endif

// cuBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
// double
#define KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(INDEX_TYPE, SCALAR, LAYOUT, MEMSPACE)                              \
  template <>                                                                                                      \
  struct iamax_tpl_spec_avail<                                                                                     \
      Kokkos::Cuda, Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                      \
      1> {                                                                                                         \
    enum : bool { value = true };                                                                                  \
  };                                                                                                               \
  template <>                                                                                                      \
  struct iamax_tpl_spec_avail<Kokkos::Cuda,                                                                        \
                              Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,             \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                              \
                              Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,          \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                              \
                              1> {                                                                                 \
    enum : bool { value = true };                                                                                  \
  };

KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned long, double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned int, double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned long, float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned int, float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned long, Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned int, Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned long, Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned int, Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)

KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned long, double, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned int, double, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned long, float, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned int, float, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned long, Kokkos::complex<double>, Kokkos::LayoutLeft,
                                        Kokkos::CudaUVMSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned int, Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned long, Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned int, Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)

#endif

// rocBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)

#define KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(INDEX_TYPE, SCALAR, LAYOUT, MEMSPACE)                            \
  template <>                                                                                                     \
  struct iamax_tpl_spec_avail<                                                                                    \
      Kokkos::HIP, Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>,                                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                     \
      1> {                                                                                                        \
    enum : bool { value = true };                                                                                 \
  };                                                                                                              \
  template <>                                                                                                     \
  struct iamax_tpl_spec_avail<Kokkos::HIP,                                                                        \
                              Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>,             \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                             \
                              Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>,          \
                                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                             \
                              1> {                                                                                \
    enum : bool { value = true };                                                                                 \
  };

KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned long, double, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned int, double, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned long, float, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned int, float, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned long, Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned int, Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned long, Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned int, Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIPSpace)

#endif

}  // namespace Impl
}  // namespace KokkosBlas
#endif
