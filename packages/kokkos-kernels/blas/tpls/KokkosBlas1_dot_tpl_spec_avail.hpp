// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class AV, class XMV, class YMV, int Xrank = XMV::rank, int Yrank = YMV::rank>
struct dot_tpl_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

namespace KokkosBlas {
namespace Impl {

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS

#define KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_BLAS(SCALAR)                                                                 \
  template <typename ExecSpace>                                                                                     \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                   \
  struct dot_tpl_spec_avail<                                                                                        \
      ExecSpace,                                                                                                    \
      Kokkos::View<SCALAR, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,        \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, 1> { \
    enum : bool { value = true };                                                                                   \
  };

KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_BLAS(double)
KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_BLAS(float)

// TODO: we met difficuties in FindTPLMKL.cmake to set the BLAS library properly
// such that the test in CheckHostBlasReturnComplex.cmake could not be
// compiled and run to give a correct answer on KK_BLAS_RESULT_AS_POINTER_ARG.
// This resulted in segfault in dot() with MKL and complex.
// So we just temporarily disable it until FindTPLMKL.cmake is fixed.
#if !defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>)
KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>)
#endif

#endif

#define KOKKOSBLAS1_DOT_TPL_SPEC(SCALAR, EXECSPACE)                                                                 \
  template <>                                                                                                       \
  struct dot_tpl_spec_avail<                                                                                        \
      EXECSPACE,                                                                                                    \
      Kokkos::View<SCALAR, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,        \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, EXECSPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
      Kokkos::View<const SCALAR*, Kokkos::LayoutLeft, EXECSPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, 1> { \
    enum : bool { value = true };                                                                                   \
  };

#define KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL(EXECSPACE)             \
  KOKKOSBLAS1_DOT_TPL_SPEC(float, EXECSPACE)                  \
  KOKKOSBLAS1_DOT_TPL_SPEC(double, EXECSPACE)                 \
  KOKKOSBLAS1_DOT_TPL_SPEC(Kokkos::complex<float>, EXECSPACE) \
  KOKKOSBLAS1_DOT_TPL_SPEC(Kokkos::complex<double>, EXECSPACE)

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
// Note BMK: CUBLAS dot is consistently slower than our native dot
// (measured 11.2, 11.8, 12.0 using perf test, and all are similar)
// If a future version improves performance, re-enable it here and
// in the tpl_spec_decl file.
#if 0
KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL(Kokkos::Cuda)
#endif
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL(Kokkos::HIP)
#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_MKL) && defined(KOKKOS_ENABLE_SYCL)
KOKKOSBLAS1_DOT_TPL_SPEC_AVAIL(Kokkos::SYCL)
#endif
}  // namespace Impl
}  // namespace KokkosBlas
#endif
