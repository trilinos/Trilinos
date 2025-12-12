// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class RV, class AV, class XV, int Xrank = XV::rank>
struct scal_tpl_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

namespace KokkosBlas {
namespace Impl {

// Generic Host side BLAS (could be MKL or whatever)
#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)
// double
#define KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUT, MEMSPACE)                                              \
  template <class ExecSpace>                                                                                        \
  struct scal_tpl_spec_avail<                                                                                       \
      ExecSpace,                                                                                                    \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      SCALAR,                                                                                                       \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                       \
      1> {                                                                                                          \
    enum : bool { value = true };                                                                                   \
  };

KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HostSpace)

#endif

// cuBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
// double
#define KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(SCALAR, LAYOUT, EXECSPACE, MEMSPACE)                                 \
  template <>                                                                                                       \
  struct scal_tpl_spec_avail<                                                                                       \
      EXECSPACE,                                                                                                    \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      SCALAR,                                                                                                       \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                       \
      1> {                                                                                                          \
    enum : bool { value = true };                                                                                   \
  };

KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)

KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace)

#endif

// rocBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)

#define KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(SCALAR, LAYOUT, EXECSPACE, MEMSPACE)                                \
  template <>                                                                                                       \
  struct scal_tpl_spec_avail<                                                                                       \
      EXECSPACE,                                                                                                    \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      SCALAR,                                                                                                       \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                       \
      1> {                                                                                                          \
    enum : bool { value = true };                                                                                   \
  };

KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(double, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(float, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace)

#endif

}  // namespace Impl
}  // namespace KokkosBlas
#endif
