//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSBLAS1_NRMINF_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_NRMINF_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
template <class RV, class XV>
inline void nrminf_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
  printf("KokkosBlas1::nrminf<> TPL Blas specialization for < %s , %s >\n", typeid(RV).name(), typeid(XV).name());
#endif
}
}  // namespace
}  // namespace Impl
}  // namespace KokkosBlas

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DNRMINF_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                     \
  template <class ExecSpace>                                                                                         \
  struct NrmInf<ExecSpace, Kokkos::View<double, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
                Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                             \
                             Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                               \
                1, true, ETI_SPEC_AVAIL> {                                                                           \
    typedef Kokkos::View<double, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> RV;             \
    typedef Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                    \
        XV;                                                                                                          \
    typedef typename XV::size_type size_type;                                                                        \
    typedef Kokkos::Details::InnerProductSpaceTraits<double> IPT;                                                    \
                                                                                                                     \
    static void nrminf(const ExecSpace& space, RV& R, const XV& X) {                                                 \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrminf[TPL_BLAS,double]");                                          \
      const size_type numElems = X.extent(0);                                                                        \
      if (numElems == 0) {                                                                                           \
        R() = 0.0;                                                                                                   \
        return;                                                                                                      \
      }                                                                                                              \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                              \
        nrminf_print_specialization<RV, XV>();                                                                       \
        int N   = numElems;                                                                                          \
        int one = 1;                                                                                                 \
        int idx = HostBlas<double>::iamax(N, X.data(), one) - 1;                                                     \
        R()     = IPT::norm(X(idx));                                                                                 \
      } else {                                                                                                       \
        NrmInf<ExecSpace, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrminf(space, R, X);                                    \
      }                                                                                                              \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS1_SNRMINF_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                    \
  template <class ExecSpace>                                                                                        \
  struct NrmInf<ExecSpace, Kokkos::View<float, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
                Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                             \
                             Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
                1, true, ETI_SPEC_AVAIL> {                                                                          \
    typedef Kokkos::View<float, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> RV;             \
    typedef Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                   \
        XV;                                                                                                         \
    typedef typename XV::size_type size_type;                                                                       \
    typedef Kokkos::Details::InnerProductSpaceTraits<float> IPT;                                                    \
                                                                                                                    \
    static void nrminf(const ExecSpace& space, RV& R, const XV& X) {                                                \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrminf[TPL_BLAS,float]");                                          \
      const size_type numElems = X.extent(0);                                                                       \
      if (numElems == 0) {                                                                                          \
        R() = 0.0f;                                                                                                 \
        return;                                                                                                     \
      }                                                                                                             \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                             \
        nrminf_print_specialization<RV, XV>();                                                                      \
        int N   = numElems;                                                                                         \
        int one = 1;                                                                                                \
        int idx = HostBlas<float>::iamax(N, X.data(), one) - 1;                                                     \
        R()     = IPT::norm(X(idx));                                                                                \
      } else {                                                                                                      \
        NrmInf<ExecSpace, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrminf(space, R, X);                                   \
      }                                                                                                             \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_ZNRMINF_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                     \
  template <class ExecSpace>                                                                                         \
  struct NrmInf<ExecSpace, Kokkos::View<double, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
                Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,            \
                             Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                               \
                1, true, ETI_SPEC_AVAIL> {                                                                           \
    typedef Kokkos::View<double, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> RV;             \
    typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                    \
        XV;                                                                                                          \
    typedef typename XV::size_type size_type;                                                                        \
    typedef Kokkos::Details::InnerProductSpaceTraits<Kokkos::complex<double>> IPT;                                   \
                                                                                                                     \
    static void nrminf(const ExecSpace& space, RV& R, const XV& X) {                                                 \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrminf[TPL_BLAS,complex<double>]");                                 \
      const size_type numElems = X.extent(0);                                                                        \
      if (numElems == 0) {                                                                                           \
        R() = 0.0;                                                                                                   \
        return;                                                                                                      \
      }                                                                                                              \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                              \
        nrminf_print_specialization<RV, XV>();                                                                       \
        int N   = numElems;                                                                                          \
        int one = 1;                                                                                                 \
        int idx =                                                                                                    \
            HostBlas<std::complex<double>>::iamax(N, reinterpret_cast<const std::complex<double>*>(X.data()), one) - \
            1;                                                                                                       \
        R() = IPT::norm(X(idx));                                                                                     \
      } else {                                                                                                       \
        NrmInf<ExecSpace, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrminf(space, R, X);                                    \
      }                                                                                                              \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS1_CNRMINF_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                      \
  template <class ExecSpace>                                                                                          \
  struct NrmInf<ExecSpace, Kokkos::View<float, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
                Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,              \
                             Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                \
                1, true, ETI_SPEC_AVAIL> {                                                                            \
    typedef Kokkos::View<float, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> RV;               \
    typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                  \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>                                                     \
        XV;                                                                                                           \
    typedef typename XV::size_type size_type;                                                                         \
    typedef Kokkos::Details::InnerProductSpaceTraits<Kokkos::complex<float>> IPT;                                     \
                                                                                                                      \
    static void nrminf(const ExecSpace& space, RV& R, const XV& X) {                                                  \
      Kokkos::Profiling::pushRegion("KokkosBlas::nrminf[TPL_BLAS,complex<float>]");                                   \
      const size_type numElems = X.extent(0);                                                                         \
      if (numElems == 0) {                                                                                            \
        R() = 0.0f;                                                                                                   \
        return;                                                                                                       \
      }                                                                                                               \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                               \
        nrminf_print_specialization<RV, XV>();                                                                        \
        int N   = numElems;                                                                                           \
        int one = 1;                                                                                                  \
        int idx =                                                                                                     \
            HostBlas<std::complex<float>>::iamax(N, reinterpret_cast<const std::complex<float>*>(X.data()), one) - 1; \
        R() = IPT::norm(X(idx));                                                                                      \
      } else {                                                                                                        \
        NrmInf<ExecSpace, RV, XV, 1, false, ETI_SPEC_AVAIL>::nrminf(space, R, X);                                     \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

KOKKOSBLAS1_DNRMINF_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_DNRMINF_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_SNRMINF_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_SNRMINF_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_ZNRMINF_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_ZNRMINF_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_CNRMINF_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_CNRMINF_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#endif
