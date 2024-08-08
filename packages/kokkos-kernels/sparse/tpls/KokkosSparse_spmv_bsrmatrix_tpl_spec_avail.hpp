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

#ifndef KOKKOSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_HPP_
#define KOKKOSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include <mkl.h>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector>
struct spmv_bsrmatrix_tpl_spec_avail {
  enum : bool { value = false };
};

// cuSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

// These versions of cuSPARSE require the ordinal and offset types to be the
// same. For KokkosKernels, this means int/int only.

#define KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(SCALAR, ORDINAL, OFFSET, XL, YL, MEMSPACE)              \
  template <>                                                                                                       \
  struct spmv_bsrmatrix_tpl_spec_avail<                                                                             \
      Kokkos::Cuda, KokkosSparse::Impl::SPMVHandleImpl<Kokkos::Cuda, MEMSPACE, SCALAR, OFFSET, ORDINAL>,            \
      ::KokkosSparse::Experimental::BsrMatrix<const SCALAR, const ORDINAL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,  \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET>,               \
      Kokkos::View<const SCALAR*, XL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                 \
      Kokkos::View<SCALAR*, YL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> { \
    enum : bool { value = true };                                                                                   \
  };

#if (9000 <= CUDA_VERSION)

KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                                    Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                                    Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                                    Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                                    Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                                    Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                                    Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                                    Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                                    Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                                                    Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft,
                                                    Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight,
                                                    Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight,
                                                    Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                                                    Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft,
                                                    Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight,
                                                    Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight,
                                                    Kokkos::LayoutRight, Kokkos::CudaUVMSpace)

#endif  // CUDA/CUSPARSE >= 9.0?
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#define KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(SCALAR, EXECSPACE)                                    \
  template <>                                                                                                \
  struct spmv_bsrmatrix_tpl_spec_avail<                                                                      \
      EXECSPACE, KokkosSparse::Impl::SPMVHandleImpl<EXECSPACE, Kokkos::HostSpace, SCALAR, MKL_INT, MKL_INT>, \
      ::KokkosSparse::Experimental::BsrMatrix<const SCALAR, const MKL_INT,                                   \
                                              Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                  \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const MKL_INT>,       \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                          \
      Kokkos::View<SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                               \
    enum : bool { value = true };                                                                            \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(float, Kokkos::Serial)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(double, Kokkos::Serial)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(Kokkos::complex<float>, Kokkos::Serial)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(Kokkos::complex<double>, Kokkos::Serial)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(float, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(double, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(Kokkos::complex<float>, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(Kokkos::complex<double>, Kokkos::OpenMP)
#endif

#endif

// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector,
          const bool integerScalarType = std::is_integral_v<typename AMatrix::non_const_value_type>>
struct spmv_mv_bsrmatrix_tpl_spec_avail {
  enum : bool { value = false };
};

// cuSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

// These versions of cuSPARSE require the ordinal and offset types to be the
// same. For KokkosKernels, this means int/int only.
// cuSparse level 3 does not currently support LayoutRight
#define KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(SCALAR, ORDINAL, OFFSET, LAYOUT, MEMSPACE)              \
  template <>                                                                                                          \
  struct spmv_mv_bsrmatrix_tpl_spec_avail<                                                                             \
      Kokkos::Cuda, KokkosSparse::Impl::SPMVHandleImpl<Kokkos::Cuda, MEMSPACE, SCALAR, OFFSET, ORDINAL>,               \
      ::KokkosSparse::Experimental::BsrMatrix<const SCALAR, const ORDINAL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,     \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET>,                  \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                                     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                    \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      false> {                                                                                                         \
    enum : bool { value = true };                                                                                      \
  };

#if (9000 <= CUDA_VERSION)

KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                                                       Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft,
                                                       Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                                                       Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft,
                                                       Kokkos::CudaUVMSpace)

#endif  // CUDA/CUSPARSE >= 9.0?
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#define KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(SCALAR, EXECSPACE)                                         \
  template <>                                                                                                        \
  struct spmv_mv_bsrmatrix_tpl_spec_avail<                                                                           \
      EXECSPACE, KokkosSparse::Impl::SPMVHandleImpl<EXECSPACE, Kokkos::HostSpace, SCALAR, MKL_INT, MKL_INT>,         \
      ::KokkosSparse::Experimental::BsrMatrix<const SCALAR, const int, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>, \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const int>,                   \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                  \
      Kokkos::View<SCALAR*, Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                         \
      true> {                                                                                                        \
    enum : bool { value = true };                                                                                    \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(float, Kokkos::Serial)
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(double, Kokkos::Serial)
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(Kokkos::complex<float>, Kokkos::Serial)
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(Kokkos::complex<double>, Kokkos::Serial)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(float, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(double, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(Kokkos::complex<float>, Kokkos::OpenMP)
KOKKOSSPARSE_SPMV_MV_BSRMATRIX_TPL_SPEC_AVAIL_MKL(Kokkos::complex<double>, Kokkos::OpenMP)
#endif

#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)

#include "KokkosSparse_Utils_rocsparse.hpp"

#define KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(SCALAR, ORDINAL, OFFSET, LAYOUT, MEMSPACE)                \
  template <>                                                                                                          \
  struct spmv_bsrmatrix_tpl_spec_avail<                                                                                \
      Kokkos::HIP, KokkosSparse::Impl::SPMVHandleImpl<Kokkos::HIP, MEMSPACE, SCALAR, OFFSET, ORDINAL>,                 \
      ::KokkosSparse::Experimental::BsrMatrix<const SCALAR, const ORDINAL, Kokkos::Device<Kokkos::HIP, MEMSPACE>,      \
                                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET>,                  \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                    \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::HIP, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> { \
    enum : bool { value = true };                                                                                      \
  };

#if KOKKOSSPARSE_IMPL_ROCM_VERSION >= 50200

KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(float, rocsparse_int, rocsparse_int, Kokkos::LayoutLeft,
                                                     Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(double, rocsparse_int, rocsparse_int, Kokkos::LayoutLeft,
                                                     Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(float, rocsparse_int, rocsparse_int, Kokkos::LayoutRight,
                                                     Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(double, rocsparse_int, rocsparse_int, Kokkos::LayoutRight,
                                                     Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<float>, rocsparse_int, rocsparse_int,
                                                     Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<double>, rocsparse_int, rocsparse_int,
                                                     Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<float>, rocsparse_int, rocsparse_int,
                                                     Kokkos::LayoutRight, Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<double>, rocsparse_int, rocsparse_int,
                                                     Kokkos::LayoutRight, Kokkos::HIPSpace)

#endif  // KOKKOSSPARSE_IMPL_ROCM_VERSION >= 50200

#undef KOKKOSSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_ROCSPARSE

#endif  // defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSPARSE_SPMV_BSRMATRIX_TPL_SPEC_AVAIL_HPP_
