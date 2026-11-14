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

#define KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(SCALAR)                                                                \
  template <typename ExecSpace>                                                                                     \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                   \
  struct scal_tpl_spec_avail<                                                                                       \
      ExecSpace, Kokkos::View<SCALAR*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      SCALAR, Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      1> {                                                                                                          \
    enum : bool { value = true };                                                                                   \
  };

KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(double)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(float)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>)

#endif

// cuBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
// double
#define KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(SCALAR)                                                                 \
  template <>                                                                                                          \
  struct scal_tpl_spec_avail<                                                                                          \
      Kokkos::Cuda, Kokkos::View<SCALAR*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      SCALAR, Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      1> {                                                                                                             \
    enum : bool { value = true };                                                                                      \
  };

KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(double)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(float)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>)

#endif

// rocBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)

#define KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(SCALAR)                                                               \
  template <>                                                                                                         \
  struct scal_tpl_spec_avail<                                                                                         \
      Kokkos::HIP, Kokkos::View<SCALAR*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,  \
      SCALAR, Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      1> {                                                                                                            \
    enum : bool { value = true };                                                                                     \
  };

KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(double)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(float)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<double>)
KOKKOSBLAS1_SCAL_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<float>)

#endif

}  // namespace Impl
}  // namespace KokkosBlas
#endif
