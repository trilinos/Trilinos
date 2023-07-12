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
template <class AT, class AO, class AD, class AM, class AS, class XT, class XL,
          class XD, class XM, class YT, class YL, class YD, class YM,
          const bool integerScalarType =
              std::is_integral<typename std::decay<AT>::type>::value>
struct spmv_mv_tpl_spec_avail {
  enum : bool { value = false };
};

#define KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(SCALAR, ORDINAL, OFFSET, \
                                                     XL, YL, MEMSPACE)        \
  template <>                                                                 \
  struct spmv_mv_tpl_spec_avail<                                              \
      const SCALAR, const ORDINAL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,    \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET, const SCALAR**,  \
      XL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                             \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>,         \
      SCALAR**, YL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                   \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > {                             \
    enum : bool { value = true };                                             \
  };

/* CUSPARSE_VERSION 10300 and lower seem to have a bug in cusparseSpMM
non-transpose that produces incorrect result. This is cusparse distributed with
CUDA 10.1.243. The bug seems to be resolved by CUSPARSE 10301 (present by
CUDA 10.2.89) */
#if defined(CUSPARSE_VERSION) && (10301 <= CUSPARSE_VERSION)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int,
                                             Kokkos::LayoutLeft,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int,
                                             Kokkos::LayoutRight,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int,
                                             Kokkos::LayoutLeft,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int,
                                             Kokkos::LayoutRight,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int,
                                             Kokkos::LayoutLeft,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int,
                                             Kokkos::LayoutRight,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int,
                                             Kokkos::LayoutLeft,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int,
                                             Kokkos::LayoutRight,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int,
                                             Kokkos::LayoutLeft,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int,
                                             Kokkos::LayoutRight,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int,
                                             Kokkos::LayoutLeft,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int,
                                             Kokkos::LayoutRight,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int,
                                             Kokkos::LayoutLeft,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int,
                                             Kokkos::LayoutRight,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int,
                                             Kokkos::LayoutLeft,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int,
                                             Kokkos::LayoutRight,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::Experimental::half_t, int,
                                             int, Kokkos::LayoutLeft,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::Experimental::half_t, int,
                                             int, Kokkos::LayoutRight,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::Experimental::half_t, int,
                                             int, Kokkos::LayoutLeft,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::Experimental::half_t, int,
                                             int, Kokkos::LayoutRight,
                                             Kokkos::LayoutLeft,
                                             Kokkos::CudaUVMSpace)

#endif
#endif  // defined(CUSPARSE_VERSION) && (10300 <= CUSPARSE_VERSION)

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSPARSE_SPMV_MV_TPL_SPEC_AVAIL_HPP_
