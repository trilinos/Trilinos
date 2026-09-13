// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {

// Specialization struct which defines whether a specialization exists
template <class execution_space, class AVT, class BVT>
struct trmm_tpl_spec_avail {
  enum : bool { value = false };
};

// Generic Host side BLAS (could be MKL or whatever)
#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)

#define KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUTA, LAYOUTB)                                       \
  template <typename ExecSpace>                                                                              \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                            \
  struct trmm_tpl_spec_avail<                                                                                \
      ExecSpace, Kokkos::View<const SCALAR**, LAYOUTA, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<SCALAR**, LAYOUTB, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {               \
    enum : bool { value = true };                                                                            \
  };

KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft, Kokkos::LayoutLeft)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft, Kokkos::LayoutLeft)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::LayoutLeft)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::LayoutLeft)

KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutRight, Kokkos::LayoutRight)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutRight, Kokkos::LayoutRight)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::LayoutRight)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::LayoutRight)

#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#define KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(SCALAR, LAYOUTA, LAYOUTB)                                           \
  template <>                                                                                                      \
  struct trmm_tpl_spec_avail<                                                                                      \
      Kokkos::Cuda, Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      Kokkos::View<SCALAR**, LAYOUTB, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                  \
    enum : bool { value = true };                                                                                  \
  };

KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft, Kokkos::LayoutLeft)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft, Kokkos::LayoutLeft)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::LayoutLeft)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::LayoutLeft)

KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutRight, Kokkos::LayoutRight)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutRight, Kokkos::LayoutRight)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::LayoutRight)
KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::LayoutRight)

#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS
}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS3_TRMM_TPL_SPEC_AVAIL_HPP_
