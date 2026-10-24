// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
#define KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS(INDEX_TYPE, SCALAR)                                                \
  template <typename ExecSpace>                                                                                  \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                \
  struct iamax_tpl_spec_avail<                                                                                   \
      ExecSpace,                                                                                                 \
      Kokkos::View<INDEX_TYPE, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1> { \
    enum : bool { value = true };                                                                                \
  };

KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS(unsigned long, double)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS(unsigned long, float)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS(unsigned long, Kokkos::complex<double>)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS(unsigned long, Kokkos::complex<float>)

#endif

// cuBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
// double
#define KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(INDEX_TYPE, SCALAR)                                                 \
  template <>                                                                                                       \
  struct iamax_tpl_spec_avail<                                                                                      \
      Kokkos::Cuda,                                                                                                 \
      Kokkos::View<INDEX_TYPE, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1> { \
    enum : bool { value = true };                                                                                   \
  };                                                                                                                \
  template <>                                                                                                       \
  struct iamax_tpl_spec_avail<                                                                                      \
      Kokkos::Cuda,                                                                                                 \
      Kokkos::View<INDEX_TYPE, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1> { \
    enum : bool { value = true };                                                                                   \
  };

KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned long, double)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned int, double)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned long, float)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned int, float)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned long, Kokkos::complex<double>)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned int, Kokkos::complex<double>)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned long, Kokkos::complex<float>)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS(unsigned int, Kokkos::complex<float>)

#endif

// rocBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)

#define KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(INDEX_TYPE, SCALAR)                                               \
  template <>                                                                                                      \
  struct iamax_tpl_spec_avail<                                                                                     \
      Kokkos::HIP,                                                                                                 \
      Kokkos::View<INDEX_TYPE, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,   \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1> { \
    enum : bool { value = true };                                                                                  \
  };                                                                                                               \
  template <>                                                                                                      \
  struct iamax_tpl_spec_avail<                                                                                     \
      Kokkos::HIP,                                                                                                 \
      Kokkos::View<INDEX_TYPE, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1> { \
    enum : bool { value = true };                                                                                  \
  };

KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned long, double)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned int, double)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned long, float)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned int, float)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned long, Kokkos::complex<double>)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned int, Kokkos::complex<double>)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned long, Kokkos::complex<float>)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_ROCBLAS(unsigned int, Kokkos::complex<float>)

#endif

}  // namespace Impl
}  // namespace KokkosBlas
#endif
