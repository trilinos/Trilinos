// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HPP_
#define KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HPP_

namespace KokkosLapack {
namespace Impl {
// Backend specializations will opt in only for supported view types.
template <class ExecutionSpace, class AMatrix, class IpivView, class BMatrix, class InfoView>
struct getrs_tpl_spec_avail {
  enum : bool { value = false };
};

#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) || defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
#define KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST(SCALAR, EXEC)                                    \
  template <>                                                                                   \
  struct getrs_tpl_spec_avail<                                                                  \
      EXEC,                                                                                     \
      Kokkos::View<const SCALAR**, Kokkos::LayoutLeft, Kokkos::Device<EXEC, Kokkos::HostSpace>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                    \
      Kokkos::View<const int*, Kokkos::LayoutLeft, Kokkos::Device<EXEC, Kokkos::HostSpace>,     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                    \
      Kokkos::View<SCALAR**, Kokkos::LayoutLeft, Kokkos::Device<EXEC, Kokkos::HostSpace>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                    \
      Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::Device<EXEC, Kokkos::HostSpace>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {                                  \
    enum : bool { value = true };                                                               \
  };
#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST(float, Kokkos::Serial)
KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST(double, Kokkos::Serial)
KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST(Kokkos::complex<float>, Kokkos::Serial)
KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST(Kokkos::complex<double>, Kokkos::Serial)
#endif
#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST(float, Kokkos::OpenMP)
KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST(double, Kokkos::OpenMP)
KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST(Kokkos::complex<float>, Kokkos::OpenMP)
KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST(Kokkos::complex<double>, Kokkos::OpenMP)
#endif
#ifdef KOKKOS_ENABLE_THREADS
KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST(float, Kokkos::Threads)
KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST(double, Kokkos::Threads)
KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST(Kokkos::complex<float>, Kokkos::Threads)
KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST(Kokkos::complex<double>, Kokkos::Threads)
#endif
#undef KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HOST
#endif
}  // namespace Impl
}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_GETRS_TPL_SPEC_AVAIL_HPP_
