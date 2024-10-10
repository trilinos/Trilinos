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

#ifndef KOKKOSPARSE_SPMV_TPL_SPEC_AVAIL_HPP_
#define KOKKOSPARSE_SPMV_TPL_SPEC_AVAIL_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include <mkl.h>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector>
struct spmv_tpl_spec_avail {
  enum : bool { value = false };
};

// cuSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

// These versions of cuSPARSE require the ordinal and offset types to be the
// same. For KokkosKernels, this means int/int only.

#define KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(SCALAR, ORDINAL, OFFSET, XL, YL, MEMSPACE)                        \
  template <>                                                                                                       \
  struct spmv_tpl_spec_avail<                                                                                       \
      Kokkos::Cuda, KokkosSparse::Impl::SPMVHandleImpl<Kokkos::Cuda, MEMSPACE, SCALAR, OFFSET, ORDINAL>,            \
      KokkosSparse::CrsMatrix<const SCALAR, const ORDINAL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                  \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET>,                               \
      Kokkos::View<const SCALAR*, XL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                 \
      Kokkos::View<SCALAR*, YL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> { \
    enum : bool { value = true };                                                                                   \
  };

#if (9000 <= CUDA_VERSION)

KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                          Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                          Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                          Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                          Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                          Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                          Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                          Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                          Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                          Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                          Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                          Kokkos::CudaUVMSpace)

// CUDA_VERSION by itself cannot determine whether the generic cuSPARSE API is
// available: cuSPARSE version 10.1.105 does not have the generic API, but it
// comes with the same CUDA_VERSION (10010) as 10.1.243 which does.
#if defined(CUSPARSE_VERSION) && (10300 <= CUSPARSE_VERSION)

// Can enable int64/size_t.
// TODO: if Nvidia ever supports int/size_t, add that too.

KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                          Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                          Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int64_t, size_t, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                          Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int64_t, size_t, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                          Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                          Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                          Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int64_t, size_t, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                          Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int64_t, size_t, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                          Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutLeft,
                                          Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutLeft,
                                          Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutRight,
                                          Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutRight,
                                          Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutLeft,
                                          Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutLeft,
                                          Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutRight,
                                          Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutRight,
                                          Kokkos::LayoutRight, Kokkos::CudaUVMSpace)

#endif  // CUSPARSE >= 10.3 (nested, implies >= 9.0)
#endif  // CUDA/CUSPARSE >= 9.0?
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

#undef KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE

#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)

#define KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ROCSPARSE(SCALAR, LAYOUT)                                              \
  template <>                                                                                                   \
  struct spmv_tpl_spec_avail<                                                                                   \
      Kokkos::HIP,                                                                                              \
      KokkosSparse::Impl::SPMVHandleImpl<Kokkos::HIP, Kokkos::HIPSpace, SCALAR, rocsparse_int, rocsparse_int>,  \
      KokkosSparse::CrsMatrix<const SCALAR, const rocsparse_int, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const rocsparse_int>,                    \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                             \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,                              \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                  \
    enum : bool { value = true };                                                                               \
  };

KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ROCSPARSE(double, Kokkos::LayoutLeft)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ROCSPARSE(float, Kokkos::LayoutLeft)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<double>, Kokkos::LayoutLeft)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<float>, Kokkos::LayoutLeft)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ROCSPARSE(double, Kokkos::LayoutRight)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ROCSPARSE(float, Kokkos::LayoutRight)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<double>, Kokkos::LayoutRight)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<float>, Kokkos::LayoutRight)

#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#define KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_MKL(SCALAR, EXECSPACE)                                              \
  template <>                                                                                                \
  struct spmv_tpl_spec_avail<                                                                                \
      EXECSPACE, KokkosSparse::Impl::SPMVHandleImpl<EXECSPACE, Kokkos::HostSpace, SCALAR, MKL_INT, MKL_INT>, \
      KokkosSparse::CrsMatrix<const SCALAR, const MKL_INT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,     \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const MKL_INT>,                       \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                          \
      Kokkos::View<SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                               \
    enum : bool { value = true };                                                                            \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_MKL(float, Kokkos::Serial)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_MKL(double, Kokkos::Serial)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_MKL(Kokkos::complex<float>, Kokkos::Serial)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_MKL(Kokkos::complex<double>, Kokkos::Serial)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_MKL(float, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_MKL(double, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_MKL(Kokkos::complex<float>, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_MKL(Kokkos::complex<double>, Kokkos::OpenMP)
#endif

#if defined(KOKKOS_ENABLE_SYCL)
#define KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ONEMKL(SCALAR, ORDINAL, MEMSPACE)                                       \
  template <>                                                                                                    \
  struct spmv_tpl_spec_avail<                                                                                    \
      Kokkos::Experimental::SYCL,                                                                                \
      KokkosSparse::Impl::SPMVHandleImpl<Kokkos::Experimental::SYCL, MEMSPACE, SCALAR, ORDINAL, ORDINAL>,        \
      KokkosSparse::CrsMatrix<const SCALAR, const ORDINAL, Kokkos::Device<Kokkos::Experimental::SYCL, MEMSPACE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const ORDINAL>,                           \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Experimental::SYCL, MEMSPACE>,      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                              \
      Kokkos::View<SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Experimental::SYCL, MEMSPACE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                                   \
    enum : bool { value = true };                                                                                \
  };

// intel-oneapi-mkl/2023.2.0: spmv with complex data types produce:
// oneapi::mkl::sparse::gemv: unimplemented functionality: currently does not
// support complex data types.
// TODO: Revisit with later versions and selectively enable this if it's
// working.

KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ONEMKL(float, std::int32_t, Kokkos::Experimental::SYCLDeviceUSMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ONEMKL(double, std::int32_t, Kokkos::Experimental::SYCLDeviceUSMSpace)
/*
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ONEMKL(
    Kokkos::complex<float>, std::int32_t,
    Kokkos::Experimental::SYCLDeviceUSMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ONEMKL(
    Kokkos::complex<double>, std::int32_t,
    Kokkos::Experimental::SYCLDeviceUSMSpace)
*/

KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ONEMKL(float, std::int64_t, Kokkos::Experimental::SYCLDeviceUSMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ONEMKL(double, std::int64_t, Kokkos::Experimental::SYCLDeviceUSMSpace)
/*
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ONEMKL(
    Kokkos::complex<float>, std::int64_t,
    Kokkos::Experimental::SYCLDeviceUSMSpace)
KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_ONEMKL(
    Kokkos::complex<double>, std::int64_t,
    Kokkos::Experimental::SYCLDeviceUSMSpace)
*/
#endif

#endif  // KOKKOSKERNELS_ENABLE_TPL_MKL

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSPARSE_SPMV_TPL_SPEC_AVAIL_HPP_
