// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class ExecutionSpace, class VectorView, class MagnitudeView, class ScalarView>
struct rot_tpl_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

namespace KokkosBlas {
namespace Impl {

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#define KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_BLAS(SCALAR)                                                              \
  template <typename ExecSpace>                                                                                  \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                \
  struct rot_tpl_spec_avail<                                                                                     \
      ExecSpace,                                                                                                 \
      Kokkos::View<SCALAR*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,     \
      Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, Kokkos::LayoutLeft, Kokkos::HostSpace, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                     \
      Kokkos::View<SCALAR, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {    \
    enum : bool { value = true };                                                                                \
  };

KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_BLAS(double)
KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_BLAS(float)
KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>)
KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>)
#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#define KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_CUBLAS(SCALAR)                                                                 \
  template <>                                                                                                         \
  struct rot_tpl_spec_avail<                                                                                          \
      Kokkos::Cuda, Kokkos::View<SCALAR*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<typename KokkosKernels::ArithTraits<SCALAR>::mag_type, Kokkos::LayoutLeft, Kokkos::Cuda,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<SCALAR, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {              \
    enum : bool { value = true };                                                                                     \
  };

KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_CUBLAS(double)
KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_CUBLAS(float)
KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>)
KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>)
#endif

// rocBLAS
/*
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#define KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_ROCBLAS(SCALAR)                  \
  template <>                                                             \
  struct rot_tpl_spec_avail<                                              \
      Kokkos::HIP,                                                          \
      Kokkos::View<SCALAR*, Kokkos::LayoutLeft, Kokkos::HIP,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,              \
      Kokkos::View<SCALAR, Kokkos::LayoutLeft, Kokkos::HIP,   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {            \
    enum : bool { value = true };                                         \
  };

KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_ROCBLAS(double)
KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_ROCBLAS(float)
KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<double>)
KOKKOSBLAS1_ROT_TPL_SPEC_AVAIL_ROCBLAS(Kokkos::complex<float>)
#endif
*/

}  // namespace Impl
}  // namespace KokkosBlas
#endif
