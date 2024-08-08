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

#ifndef KOKKOSPARSE_SPMV_MV_TPL_SPEC_AVAIL_HPP_
#define KOKKOSPARSE_SPMV_MV_TPL_SPEC_AVAIL_HPP_

namespace KokkosSparse {
namespace Impl {

// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class Handle, class AMatrix, class XVector, class YVector,
          const bool integerScalarType = std::is_integral_v<typename AMatrix::non_const_value_type>>
struct spmv_mv_tpl_spec_avail {
  enum : bool { value = false };
};

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#define KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(SCALAR, ORDINAL, OFFSET, XL, YL, MEMSPACE)                      \
  template <>                                                                                                        \
  struct spmv_mv_tpl_spec_avail<                                                                                     \
      Kokkos::Cuda, KokkosSparse::Impl::SPMVHandleImpl<Kokkos::Cuda, MEMSPACE, SCALAR, OFFSET, ORDINAL>,             \
      KokkosSparse::CrsMatrix<const SCALAR, const ORDINAL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                   \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET>,                                \
      Kokkos::View<const SCALAR**, XL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                  \
      Kokkos::View<SCALAR**, YL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> { \
    enum : bool { value = true };                                                                                    \
  };

/* CUSPARSE_VERSION 10300 and lower seem to have a bug in cusparseSpMM
non-transpose that produces incorrect result. This is cusparse distributed with
CUDA 10.1.243. The bug seems to be resolved by CUSPARSE 10301 (present by
CUDA 10.2.89) */

/* cusparseSpMM also produces incorrect results for some inputs in CUDA 11.6.1.
 * (CUSPARSE_VERSION 11702).
 * ALG1 and ALG3 produce completely incorrect results for one set of inputs.
 * ALG2 works for that case, but has low numerical accuracy in another case.
 */
#if defined(CUSPARSE_VERSION) && (10301 <= CUSPARSE_VERSION) && (CUSPARSE_VERSION != 11702)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::Experimental::half_t, int, int, Kokkos::LayoutLeft,
                                             Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::Experimental::half_t, int, int, Kokkos::LayoutRight,
                                             Kokkos::LayoutLeft, Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::Experimental::half_t, int, int, Kokkos::LayoutLeft,
                                             Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::Experimental::half_t, int, int, Kokkos::LayoutRight,
                                             Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)

#endif
#endif  // defined(CUSPARSE_VERSION) && (10300 <= CUSPARSE_VERSION)
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
#define KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_ROCSPARSE(SCALAR, XL, YL, MEMSPACE)                                     \
  template <>                                                                                                       \
  struct spmv_mv_tpl_spec_avail<                                                                                    \
      Kokkos::HIP, KokkosSparse::Impl::SPMVHandleImpl<Kokkos::HIP, MEMSPACE, SCALAR, rocsparse_int, rocsparse_int>, \
      KokkosSparse::CrsMatrix<const SCALAR, const rocsparse_int, Kokkos::Device<Kokkos::HIP, MEMSPACE>,             \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const rocsparse_int>,                        \
      Kokkos::View<const SCALAR**, XL, Kokkos::Device<Kokkos::HIP, MEMSPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                 \
      Kokkos::View<SCALAR**, YL, Kokkos::Device<Kokkos::HIP, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> { \
    enum : bool { value = true };                                                                                   \
  };

#define AVAIL_ROCSPARSE_SCALAR_MEMSPACE(SCALAR, MEMSPACE)                                                  \
  KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_ROCSPARSE(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutLeft, MEMSPACE)  \
  KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_ROCSPARSE(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutRight, MEMSPACE) \
  KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_ROCSPARSE(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutLeft, MEMSPACE) \
  KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_ROCSPARSE(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutRight, MEMSPACE)

#define AVAIL_ROCSPARSE_SCALAR(SCALAR)                      \
  AVAIL_ROCSPARSE_SCALAR_MEMSPACE(SCALAR, Kokkos::HIPSpace) \
  AVAIL_ROCSPARSE_SCALAR_MEMSPACE(SCALAR, Kokkos::HIPManagedSpace)

AVAIL_ROCSPARSE_SCALAR(float)
AVAIL_ROCSPARSE_SCALAR(double)
AVAIL_ROCSPARSE_SCALAR(Kokkos::complex<float>)
AVAIL_ROCSPARSE_SCALAR(Kokkos::complex<double>)

#undef AVAIL_ROCSPARSE_SCALAR_MEMSPACE
#undef AVAIL_ROCSPARSE_SCALAR
#undef KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_ROCSPARSE

#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSPARSE_SPMV_MV_TPL_SPEC_AVAIL_HPP_
