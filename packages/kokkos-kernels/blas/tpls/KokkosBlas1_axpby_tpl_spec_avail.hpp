// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_AXPBY_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS1_AXPBY_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class execution_space, class AV, class XMV, class BV, class YMV, int rank = YMV::rank>
struct axpby_tpl_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosBlas

namespace KokkosBlas {
namespace Impl {

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS

#define KOKKOSBLAS1_AXPBY_TPL_SPEC_AVAIL_BLAS(SCALAR, LAYOUT, MEMSPACE)                                             \
  template <class ExecSpace>                                                                                        \
  struct axpby_tpl_spec_avail<                                                                                      \
      ExecSpace, SCALAR,                                                                                            \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                       \
      SCALAR,                                                                                                       \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      1> {                                                                                                          \
    enum : bool { value = true };                                                                                   \
  };

KOKKOSBLAS1_AXPBY_TPL_SPEC_AVAIL_BLAS(double, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_AXPBY_TPL_SPEC_AVAIL_BLAS(float, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_AXPBY_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HostSpace)
KOKKOSBLAS1_AXPBY_TPL_SPEC_AVAIL_BLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HostSpace)

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#define KOKKOSBLAS1_AXPBY_TPL_SPEC_AVAIL_CUBLAS(SCALAR, LAYOUT, MEMSPACE)                                           \
  template <class ExecSpace>                                                                                        \
  struct axpby_tpl_spec_avail<                                                                                      \
      ExecSpace, SCALAR,                                                                                            \
      Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                       \
      SCALAR,                                                                                                       \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      1> {                                                                                                          \
    enum : bool { value = true };                                                                                   \
  };

KOKKOSBLAS1_AXPBY_TPL_SPEC_AVAIL_CUBLAS(double, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_AXPBY_TPL_SPEC_AVAIL_CUBLAS(float, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_AXPBY_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_AXPBY_TPL_SPEC_AVAIL_CUBLAS(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)

#endif
}  // namespace Impl
}  // namespace KokkosBlas
#endif
