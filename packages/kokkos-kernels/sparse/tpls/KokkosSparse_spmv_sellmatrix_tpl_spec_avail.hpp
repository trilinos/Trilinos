// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_HPP
#define KOKKOSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_HPP

namespace KokkosSparse {
namespace Impl {

// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class AMatrix, class XVector, class YVector>
struct spmv_sellmatrix_tpl_spec_avail {
  enum : bool { value = false };
};

// cuSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

#define KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(SCALAR, ORDINAL, OFFSET, XL, YL, MEMSPACE)             \
  template <>                                                                                                       \
  struct spmv_sellmatrix_tpl_spec_avail<                                                                            \
      Kokkos::Cuda,                                                                                                 \
      KokkosSparse::Experimental::SellMatrix<const SCALAR, const ORDINAL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,   \
                                             Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET>,                \
      Kokkos::View<const SCALAR*, XL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                 \
      Kokkos::View<SCALAR*, YL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> { \
    enum : bool { value = true };                                                                                   \
  };

KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                                     Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                                     Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutRight,
                                                     Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutRight,
                                                     Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                                     Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                                     Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                                     Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                                     Kokkos::CudaUVMSpace)

KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                                     Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                                     Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutRight,
                                                     Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutRight,
                                                     Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                                     Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                                     Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                                     Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                                     Kokkos::CudaUVMSpace)

KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                                                     Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                                                     Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                                                     Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                                                     Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight,
                                                     Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight,
                                                     Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight,
                                                     Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight,
                                                     Kokkos::LayoutRight, Kokkos::CudaUVMSpace)

KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft,
                                                     Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft,
                                                     Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft,
                                                     Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft,
                                                     Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight,
                                                     Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight,
                                                     Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight,
                                                     Kokkos::LayoutRight, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight,
                                                     Kokkos::LayoutRight, Kokkos::CudaUVMSpace)

#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

// cuSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE

#define KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(SCALAR, ORDINAL, OFFSET, XL, YL, MEMSPACE)           \
  template <>                                                                                                      \
  struct spmv_sellmatrix_tpl_spec_avail<                                                                           \
      Kokkos::HIP,                                                                                                 \
      KokkosSparse::Experimental::SellMatrix<const SCALAR, const ORDINAL, Kokkos::Device<Kokkos::HIP, MEMSPACE>,   \
                                             Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET>,               \
      Kokkos::View<const SCALAR*, XL, Kokkos::Device<Kokkos::HIP, MEMSPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                \
      Kokkos::View<SCALAR*, YL, Kokkos::Device<Kokkos::HIP, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> { \
    enum : bool { value = true };                                                                                  \
  };

KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                                      Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutRight,
                                                      Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                                      Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                                      Kokkos::HIPSpace)

KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                                      Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutRight,
                                                      Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                                      Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight,
                                                      Kokkos::HIPSpace)

KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                                                      Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                                                      Kokkos::LayoutRight, Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight,
                                                      Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight,
                                                      Kokkos::LayoutRight, Kokkos::HIPSpace)

KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft,
                                                      Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft,
                                                      Kokkos::LayoutRight, Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight,
                                                      Kokkos::LayoutLeft, Kokkos::HIPSpace)
KOKKOSSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_ROCSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight,
                                                      Kokkos::LayoutRight, Kokkos::HIPSpace)
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSPARSE_SPMV_SELLMATRIX_TPL_SPEC_AVAIL_HPP
