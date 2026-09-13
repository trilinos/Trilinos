// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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

#define KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUT)                                                \
  template <typename ExecSpace>                                                                             \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                           \
  struct gemm_tpl_spec_avail<                                                                               \
      ExecSpace, Kokkos::View<const SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,            \
      Kokkos::View<SCALAR**, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {               \
    enum : bool { value = true };                                                                           \
  };

KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutLeft)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft)

KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutRight)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutRight)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutRight)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutRight)

#endif

// cuBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)

#define KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(SCALAR, LAYOUT)                                                    \
  template <>                                                                                                     \
  struct gemm_tpl_spec_avail<                                                                                     \
      Kokkos::Cuda, Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,               \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                  \
    enum : bool { value = true };                                                                                 \
  };

KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft)

KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutRight)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutRight)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutRight)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutRight)

#endif

// rocBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)

#define KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(SCALAR, LAYOUT)                                                 \
  template <>                                                                                                   \
  struct gemm_tpl_spec_avail<                                                                                   \
      Kokkos::HIP, Kokkos::View<const SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,              \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                 \
    enum : bool { value = true };                                                                               \
  };

KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(double, Kokkos::LayoutLeft)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(float, Kokkos::LayoutLeft)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft)

KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(double, Kokkos::LayoutRight)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(float, Kokkos::LayoutRight)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<double>, Kokkos::LayoutRight)
KOKKOSBLAS3_GEMM_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<float>, Kokkos::LayoutRight)

#endif
}  // namespace Impl
}  // namespace KokkosBlas

#endif
