// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class DXView, class YView, class PView>
struct rotmg_tpl_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

namespace KokkosBlas {
namespace Impl {

// Generic Host side BLAS (could be MKL or whatever)
// ARMPL is disabled as it does not detect some corner
// cases correctly which leads to failing unit-tests
#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS) && !defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL)
#define KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUT)                                                         \
  template <typename ExecSpace>                                                                                       \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                     \
  struct rotmg_tpl_spec_avail<ExecSpace,                                                                              \
                              Kokkos::View<SCALAR, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,       \
                              Kokkos::View<const SCALAR, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
                              Kokkos::View<SCALAR[5], LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {  \
    enum : bool { value = true };                                                                                     \
  };

KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft)
KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft)
KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutRight)
KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutRight)
#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#define KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_CUBLAS(SCALAR, LAYOUT)                                          \
  template <>                                                                                            \
  struct rotmg_tpl_spec_avail<                                                                           \
      Kokkos::Cuda, Kokkos::View<SCALAR, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<const SCALAR, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,         \
      Kokkos::View<SCALAR[5], LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {          \
    enum : bool { value = true };                                                                        \
  };

KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft)
KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft)
KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutRight)
KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutRight)
#endif

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#define KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_ROCBLAS(SCALAR, LAYOUT)                                       \
  template <>                                                                                          \
  struct rotmg_tpl_spec_avail<                                                                         \
      Kokkos::HIP, Kokkos::View<SCALAR, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<const SCALAR, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,        \
      Kokkos::View<SCALAR[5], LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {         \
    enum : bool { value = true };                                                                      \
  };

// Turning off use of rocBLAS as it returns false results in some of the
// param  values as well as x1 and d1. These values can e Inf which really
// should be checked by rocBLAS : (
// KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_ROCBLAS(double, Kokkos::LayoutLeft)
// KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_ROCBLAS(float, Kokkos::LayoutLeft)
// KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_ROCBLAS(double, Kokkos::LayoutRight)
// KOKKOSBLAS1_ROTMG_TPL_SPEC_AVAIL_ROCBLAS(float, Kokkos::LayoutRight)
#endif

}  // namespace Impl
}  // namespace KokkosBlas
#endif
