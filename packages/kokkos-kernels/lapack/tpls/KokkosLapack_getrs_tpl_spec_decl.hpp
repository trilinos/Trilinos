// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSLAPACK_GETRS_TPL_SPEC_DECL_HPP_
#define KOKKOSLAPACK_GETRS_TPL_SPEC_DECL_HPP_

#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) || defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
#include <KokkosLapack_Host_tpl.hpp>
#include <KokkosKernels_ArithTraits.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <limits>
#include <complex>

namespace KokkosLapack {
namespace Impl {

template <class AView, class PivotView, class BView, class InfoView>
void lapackGetrsWrapper(const char trans[], const AView& A, const PivotView& Ipiv, const BView& B,
                        const InfoView& Info) {
  const auto max_int = static_cast<size_t>(std::numeric_limits<int>::max());
  if (A.extent(0) > max_int || B.extent(1) > max_int || A.stride(1) > max_int || B.stride(1) > max_int) {
    KokkosKernels::Impl::throw_runtime_exception("KokkosLapack::getrs: dimensions and strides must fit LAPACK int.");
  }
  const int n    = static_cast<int>(A.extent(0));
  const int nrhs = static_cast<int>(B.extent(1));
  const int lda  = static_cast<int>(A.stride(1));
  const int ldb  = static_cast<int>(B.stride(1));
  using Scalar   = typename AView::non_const_value_type;
  if constexpr (KokkosKernels::ArithTraits<Scalar>::is_complex) {
    using HostScalar = std::complex<typename KokkosKernels::ArithTraits<Scalar>::mag_type>;
    HostLapack<HostScalar>::getrs(trans[0], n, nrhs, reinterpret_cast<const HostScalar*>(A.data()), lda, Ipiv.data(),
                                  reinterpret_cast<HostScalar*>(B.data()), ldb, Info.data());
  } else {
    HostLapack<Scalar>::getrs(trans[0], n, nrhs, A.data(), lda, Ipiv.data(), B.data(), ldb, Info.data());
  }
}

#define KOKKOSLAPACK_GETRS_LAPACK(SCALAR, LAYOUT, EXECSPACE, MEM_SPACE)                                                \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct GETRS<                                                                                                        \
      EXECSPACE,                                                                                                       \
      Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const int*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<int*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, \
      ETI_SPEC_AVAIL> {                                                                                                \
    using AViewType = Kokkos::View<const SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>,                       \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                           \
    using BViewType =                                                                                                  \
        Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using IpivView = Kokkos::View<const int*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>,                            \
                                  Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                            \
    using InfoViewType =                                                                                               \
        Kokkos::View<int*, LAYOUT, Kokkos::Device<EXECSPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
                                                                                                                       \
    static void getrs(const EXECSPACE& space, const char trans[], const AViewType& A, const IpivView& Ipiv,            \
                      const BViewType& B, const InfoViewType& Info) {                                                  \
      Kokkos::Profiling::ScopedRegion region("KokkosLapack::getrs[TPL_LAPACK," #SCALAR "]");                           \
      space.fence("KokkosLapack::getrs: synchronize inputs before host LAPACK");                                       \
      lapackGetrsWrapper(trans, A, Ipiv, B, Info);                                                                     \
    }                                                                                                                  \
  };

#if defined(KOKKOS_ENABLE_SERIAL)
KOKKOSLAPACK_GETRS_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GETRS_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GETRS_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
KOKKOSLAPACK_GETRS_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
KOKKOSLAPACK_GETRS_LAPACK(float, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GETRS_LAPACK(double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GETRS_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
KOKKOSLAPACK_GETRS_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace)
#endif

#if defined(KOKKOS_ENABLE_THREADS)
KOKKOSLAPACK_GETRS_LAPACK(float, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GETRS_LAPACK(double, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GETRS_LAPACK(Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
KOKKOSLAPACK_GETRS_LAPACK(Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace)
#endif

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK

#endif  // KOKKOSLAPACK_GETRS_TPL_SPEC_DECL_HPP_
