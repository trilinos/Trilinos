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

#ifndef KOKKOSBLAS1_DOT_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
template <class RV, class XV, class YV>
inline void dot_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
  printf("KokkosBlas1::dot<> TPL Blas specialization for < %s , %s , %s >\n",
         typeid(RV).name(), typeid(XV).name(), typeid(YV).name());
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

#define KOKKOSBLAS1_DDOT_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)  \
  template <class ExecSpace>                                                   \
  struct Dot<                                                                  \
      ExecSpace,                                                               \
      Kokkos::View<double, LAYOUT, Kokkos::HostSpace,                          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      1, 1, true, ETI_SPEC_AVAIL> {                                            \
    typedef Kokkos::View<double, LAYOUT, Kokkos::HostSpace,                    \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        RV;                                                                    \
    typedef Kokkos::View<const double*, LAYOUT,                                \
                         Kokkos::Device<ExecSpace, MEMSPACE>,                  \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XV;                                                                    \
    typedef typename XV::size_type size_type;                                  \
                                                                               \
    static void dot(const ExecSpace& space, RV& R, const XV& X, const XV& Y) { \
      Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_BLAS,double]");       \
      const size_type numElems = X.extent(0);                                  \
      if (numElems < static_cast<size_type>(INT_MAX)) {                        \
        dot_print_specialization<RV, XV, XV>();                                \
        int N   = numElems;                                                    \
        int one = 1;                                                           \
        R()     = HostBlas<double>::dot(N, X.data(), one, Y.data(), one);      \
      } else {                                                                 \
        Dot<ExecSpace, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R, \
                                                                     X, Y);    \
      }                                                                        \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS1_SDOT_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)  \
  template <class ExecSpace>                                                   \
  struct Dot<                                                                  \
      ExecSpace,                                                               \
      Kokkos::View<float, LAYOUT, Kokkos::HostSpace,                           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      1, 1, true, ETI_SPEC_AVAIL> {                                            \
    typedef Kokkos::View<float, LAYOUT, Kokkos::HostSpace,                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        RV;                                                                    \
    typedef Kokkos::View<const float*, LAYOUT,                                 \
                         Kokkos::Device<ExecSpace, MEMSPACE>,                  \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XV;                                                                    \
    typedef typename XV::size_type size_type;                                  \
                                                                               \
    static void dot(const ExecSpace& space, RV& R, const XV& X, const XV& Y) { \
      Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_BLAS,float]");        \
      const size_type numElems = X.extent(0);                                  \
      if (numElems < static_cast<size_type>(INT_MAX)) {                        \
        dot_print_specialization<RV, XV, XV>();                                \
        int N   = numElems;                                                    \
        int one = 1;                                                           \
        R()     = HostBlas<float>::dot(N, X.data(), one, Y.data(), one);       \
      } else {                                                                 \
        Dot<ExecSpace, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R, \
                                                                     X, Y);    \
      }                                                                        \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS1_ZDOT_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)    \
  template <class ExecSpace>                                                     \
  struct Dot<ExecSpace,                                                          \
             Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::HostSpace,    \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT,                \
                          Kokkos::Device<ExecSpace, MEMSPACE>,                   \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT,                \
                          Kokkos::Device<ExecSpace, MEMSPACE>,                   \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
             1, 1, true, ETI_SPEC_AVAIL> {                                       \
    typedef Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::HostSpace,     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >               \
        RV;                                                                      \
    typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT,                 \
                         Kokkos::Device<ExecSpace, MEMSPACE>,                    \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >               \
        XV;                                                                      \
    typedef typename XV::size_type size_type;                                    \
                                                                                 \
    static void dot(const ExecSpace& space, RV& R, const XV& X, const XV& Y) {   \
      Kokkos::Profiling::pushRegion(                                             \
          "KokkosBlas::dot[TPL_BLAS,complex<double>]");                          \
      const size_type numElems = X.extent(0);                                    \
      if (numElems < static_cast<size_type>(INT_MAX)) {                          \
        dot_print_specialization<RV, XV, XV>();                                  \
        int N   = numElems;                                                      \
        int one = 1;                                                             \
        R()     = HostBlas<std::complex<double> >::dot(                          \
            N, reinterpret_cast<const std::complex<double>*>(X.data()), one, \
            reinterpret_cast<const std::complex<double>*>(Y.data()), one);   \
      } else {                                                                   \
        Dot<ExecSpace, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R,   \
                                                                     X, Y);      \
      }                                                                          \
      Kokkos::Profiling::popRegion();                                            \
    }                                                                            \
  };

#define KOKKOSBLAS1_CDOT_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)   \
  template <class ExecSpace>                                                    \
  struct Dot<ExecSpace,                                                         \
             Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::HostSpace,    \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,            \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT,                \
                          Kokkos::Device<ExecSpace, MEMSPACE>,                  \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,            \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT,                \
                          Kokkos::Device<ExecSpace, MEMSPACE>,                  \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,            \
             1, 1, true, ETI_SPEC_AVAIL> {                                      \
    typedef Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::HostSpace,     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >              \
        RV;                                                                     \
    typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT,                 \
                         Kokkos::Device<ExecSpace, MEMSPACE>,                   \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >              \
        XV;                                                                     \
    typedef typename XV::size_type size_type;                                   \
                                                                                \
    static void dot(const ExecSpace& space, RV& R, const XV& X, const XV& Y) {  \
      Kokkos::Profiling::pushRegion(                                            \
          "KokkosBlas::dot[TPL_BLAS,complex<float>]");                          \
      const size_type numElems = X.extent(0);                                   \
      if (numElems < static_cast<size_type>(INT_MAX)) {                         \
        dot_print_specialization<RV, XV, XV>();                                 \
        int N   = numElems;                                                     \
        int one = 1;                                                            \
        R()     = HostBlas<std::complex<float> >::dot(                          \
            N, reinterpret_cast<const std::complex<float>*>(X.data()), one, \
            reinterpret_cast<const std::complex<float>*>(Y.data()), one);   \
      } else {                                                                  \
        Dot<ExecSpace, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R,  \
                                                                     X, Y);     \
      }                                                                         \
      Kokkos::Profiling::popRegion();                                           \
    }                                                                           \
  };

KOKKOSBLAS1_DDOT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_DDOT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace,
                                    false)

KOKKOSBLAS1_SDOT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_SDOT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace,
                                    false)

KOKKOSBLAS1_ZDOT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_ZDOT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace,
                                    false)

KOKKOSBLAS1_CDOT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_CDOT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace,
                                    false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DDOT_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE,     \
                                              ETI_SPEC_AVAIL)                  \
  template <>                                                                  \
  struct Dot<                                                                  \
      EXECSPACE,                                                               \
      Kokkos::View<double, LAYOUT, Kokkos::HostSpace,                          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      1, 1, true, ETI_SPEC_AVAIL> {                                            \
    typedef Kokkos::View<double, LAYOUT, Kokkos::HostSpace,                    \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        RV;                                                                    \
    typedef Kokkos::View<const double*, LAYOUT,                                \
                         Kokkos::Device<EXECSPACE, MEMSPACE>,                  \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XV;                                                                    \
    typedef typename XV::size_type size_type;                                  \
                                                                               \
    static void dot(const EXECSPACE& space, RV& R, const XV& X, const XV& Y) { \
      Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_CUBLAS,double]");     \
      const size_type numElems = X.extent(0);                                  \
      if (numElems < static_cast<size_type>(INT_MAX)) {                        \
        dot_print_specialization<RV, XV, XV>();                                \
        const int N       = static_cast<int>(numElems);                        \
        constexpr int one = 1;                                                 \
        KokkosBlas::Impl::CudaBlasSingleton& s =                               \
            KokkosBlas::Impl::CudaBlasSingleton::singleton();                  \
        cublasDdot(s.handle, N, X.data(), one, Y.data(), one, &R());           \
      } else {                                                                 \
        Dot<EXECSPACE, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R, \
                                                                     X, Y);    \
      }                                                                        \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS1_SDOT_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE,     \
                                              ETI_SPEC_AVAIL)                  \
  template <>                                                                  \
  struct Dot<                                                                  \
      EXECSPACE,                                                               \
      Kokkos::View<float, LAYOUT, Kokkos::HostSpace,                           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      1, 1, true, ETI_SPEC_AVAIL> {                                            \
    typedef Kokkos::View<float, LAYOUT, Kokkos::HostSpace,                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        RV;                                                                    \
    typedef Kokkos::View<const float*, LAYOUT,                                 \
                         Kokkos::Device<EXECSPACE, MEMSPACE>,                  \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XV;                                                                    \
    typedef typename XV::size_type size_type;                                  \
                                                                               \
    static void dot(const EXECSPACE& space, RV& R, const XV& X, const XV& Y) { \
      Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_CUBLAS,float]");      \
      const size_type numElems = X.extent(0);                                  \
      if (numElems < static_cast<size_type>(INT_MAX)) {                        \
        dot_print_specialization<RV, XV, XV>();                                \
        const int N       = static_cast<int>(numElems);                        \
        constexpr int one = 1;                                                 \
        KokkosBlas::Impl::CudaBlasSingleton& s =                               \
            KokkosBlas::Impl::CudaBlasSingleton::singleton();                  \
        cublasSdot(s.handle, N, X.data(), one, Y.data(), one, &R());           \
      } else {                                                                 \
        Dot<EXECSPACE, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R, \
                                                                     X, Y);    \
      }                                                                        \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS1_ZDOT_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE,     \
                                              ETI_SPEC_AVAIL)                  \
  template <>                                                                  \
  struct Dot<EXECSPACE,                                                        \
             Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::HostSpace,  \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,           \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT,              \
                          Kokkos::Device<EXECSPACE, MEMSPACE>,                 \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,           \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT,              \
                          Kokkos::Device<EXECSPACE, MEMSPACE>,                 \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,           \
             1, 1, true, ETI_SPEC_AVAIL> {                                     \
    typedef Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::HostSpace,   \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        RV;                                                                    \
    typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT,               \
                         Kokkos::Device<EXECSPACE, MEMSPACE>,                  \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XV;                                                                    \
    typedef typename XV::size_type size_type;                                  \
                                                                               \
    static void dot(const EXECSPACE& space, RV& R, const XV& X, const XV& Y) { \
      Kokkos::Profiling::pushRegion(                                           \
          "KokkosBlas::dot[TPL_CUBLAS,complex<double>]");                      \
      const size_type numElems = X.extent(0);                                  \
      if (numElems < static_cast<size_type>(INT_MAX)) {                        \
        dot_print_specialization<RV, XV, XV>();                                \
        const int N       = static_cast<int>(numElems);                        \
        constexpr int one = 1;                                                 \
        KokkosBlas::Impl::CudaBlasSingleton& s =                               \
            KokkosBlas::Impl::CudaBlasSingleton::singleton();                  \
        cublasZdotc(s.handle, N,                                               \
                    reinterpret_cast<const cuDoubleComplex*>(X.data()), one,   \
                    reinterpret_cast<const cuDoubleComplex*>(Y.data()), one,   \
                    reinterpret_cast<cuDoubleComplex*>(&R()));                 \
      } else {                                                                 \
        Dot<EXECSPACE, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R, \
                                                                     X, Y);    \
      }                                                                        \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS1_CDOT_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE,     \
                                              ETI_SPEC_AVAIL)                  \
  template <>                                                                  \
  struct Dot<EXECSPACE,                                                        \
             Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::HostSpace,   \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,           \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT,               \
                          Kokkos::Device<EXECSPACE, MEMSPACE>,                 \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,           \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT,               \
                          Kokkos::Device<EXECSPACE, MEMSPACE>,                 \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,           \
             1, 1, true, ETI_SPEC_AVAIL> {                                     \
    typedef Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::HostSpace,    \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        RV;                                                                    \
    typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT,                \
                         Kokkos::Device<EXECSPACE, MEMSPACE>,                  \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XV;                                                                    \
    typedef typename XV::size_type size_type;                                  \
                                                                               \
    static void dot(const EXECSPACE& space, RV& R, const XV& X, const XV& Y) { \
      Kokkos::Profiling::pushRegion(                                           \
          "KokkosBlas::dot[TPL_CUBLAS,complex<float>]");                       \
      const size_type numElems = X.extent(0);                                  \
      if (numElems < static_cast<size_type>(INT_MAX)) {                        \
        dot_print_specialization<RV, XV, XV>();                                \
        const int N       = static_cast<int>(numElems);                        \
        constexpr int one = 1;                                                 \
        KokkosBlas::Impl::CudaBlasSingleton& s =                               \
            KokkosBlas::Impl::CudaBlasSingleton::singleton();                  \
        cublasCdotc(s.handle, N, reinterpret_cast<const cuComplex*>(X.data()), \
                    one, reinterpret_cast<const cuComplex*>(Y.data()), one,    \
                    reinterpret_cast<cuComplex*>(&R()));                       \
      } else {                                                                 \
        Dot<EXECSPACE, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R, \
                                                                     X, Y);    \
      }                                                                        \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

KOKKOSBLAS1_DDOT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda,
                                      Kokkos::CudaSpace, true)
KOKKOSBLAS1_DDOT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda,
                                      Kokkos::CudaSpace, false)

KOKKOSBLAS1_SDOT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda,
                                      Kokkos::CudaSpace, true)
KOKKOSBLAS1_SDOT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda,
                                      Kokkos::CudaSpace, false)

KOKKOSBLAS1_ZDOT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda,
                                      Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZDOT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda,
                                      Kokkos::CudaSpace, false)

KOKKOSBLAS1_CDOT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda,
                                      Kokkos::CudaSpace, true)
KOKKOSBLAS1_CDOT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda,
                                      Kokkos::CudaSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#endif
