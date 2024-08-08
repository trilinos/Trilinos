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

#ifndef KOKKOSBLAS1_AXPBY_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_AXPBY_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
template <class AV, class XMV, class BV, class YMV>
inline void axpby_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
  printf(
      "KokkosBlas1::axpby<> TPL Blas specialization for < %s , %s , %s , %s "
      ">\n",
      typeid(AV).name(), typeid(XMV).name(), typeid(BV).name(), typeid(YMV).name());
#endif
}
}  // namespace
}  // namespace Impl
}  // namespace KokkosBlas

#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#include "KokkosBlas_Host_tpl.hpp"
namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DAXPBY_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                                      \
  template <class ExecSpace>                                                                                           \
  struct Axpby<                                                                                                        \
      ExecSpace, double,                                                                                               \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      double,                                                                                                          \
      Kokkos::View<double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    typedef double AV;                                                                                                 \
    typedef double BV;                                                                                                 \
    typedef Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                   \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        XV;                                                                                                            \
    typedef Kokkos::View<double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                         \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        YV;                                                                                                            \
                                                                                                                       \
    static void axpby(const ExecSpace& space, const AV& alpha, const XV& X, const BV& beta, const YV& Y) {             \
      Kokkos::Profiling::pushRegion("KokkosBlas::axpby[TPL_BLAS,double]");                                             \
      if ((X.extent(0) < INT_MAX) && (beta == 1.0)) {                                                                  \
        axpby_print_specialization<AV, XV, BV, YV>();                                                                  \
        int N   = X.extent(0);                                                                                         \
        int one = 1;                                                                                                   \
        HostBlas<double>::axpy(N, alpha, X.data(), one, Y.data(), one);                                                \
      } else                                                                                                           \
        Axpby<ExecSpace, AV, XV, BV, YV, YV::rank, false, ETI_SPEC_AVAIL>::axpby(space, alpha, X, beta, Y);            \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS1_SAXPBY_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                                     \
  template <class ExecSpace>                                                                                          \
  struct Axpby<                                                                                                       \
      ExecSpace, float,                                                                                               \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      float,                                                                                                          \
      Kokkos::View<float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    typedef float AV;                                                                                                 \
    typedef float BV;                                                                                                 \
    typedef Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                   \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        XV;                                                                                                           \
    typedef Kokkos::View<float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                         \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        YV;                                                                                                           \
                                                                                                                      \
    static void axpby(const ExecSpace& space, const AV& alpha, const XV& X, const BV& beta, const YV& Y) {            \
      Kokkos::Profiling::pushRegion("KokkosBlas::axpby[TPL_BLAS,float]");                                             \
      if ((X.extent(0) < INT_MAX) && (beta == 1.0f)) {                                                                \
        axpby_print_specialization<AV, XV, BV, YV>();                                                                 \
        int N   = X.extent(0);                                                                                        \
        int one = 1;                                                                                                  \
        HostBlas<float>::axpy(N, alpha, X.data(), one, Y.data(), one);                                                \
      } else                                                                                                          \
        Axpby<ExecSpace, AV, XV, BV, YV, YV::rank, false, ETI_SPEC_AVAIL>::axpby(space, alpha, X, beta, Y);           \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS1_ZAXPBY_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                                    \
  template <class ExecSpace>                                                                                         \
  struct Axpby<ExecSpace, Kokkos::complex<double>,                                                                   \
               Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,             \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                               \
               Kokkos::complex<double>,                                                                              \
               Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                   \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                               \
               1, true, ETI_SPEC_AVAIL> {                                                                            \
    typedef Kokkos::complex<double> AV;                                                                              \
    typedef Kokkos::complex<double> BV;                                                                              \
    typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                   \
        XV;                                                                                                          \
    typedef Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                      \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                   \
        YV;                                                                                                          \
                                                                                                                     \
    static void axpby(const ExecSpace& space, const AV& alpha, const XV& X, const BV& beta, const YV& Y) {           \
      Kokkos::Profiling::pushRegion("KokkosBlas::axpby[TPL_BLAS,complex<double>]");                                  \
      if ((X.extent(0) < INT_MAX) && (beta == 1.0f)) {                                                               \
        axpby_print_specialization<AV, XV, BV, YV>();                                                                \
        int N                                = X.extent(0);                                                          \
        int one                              = 1;                                                                    \
        const std::complex<double> alpha_val = alpha;                                                                \
        HostBlas<std::complex<double> >::axpy(N, alpha_val, reinterpret_cast<const std::complex<double>*>(X.data()), \
                                              one, reinterpret_cast<std::complex<double>*>(Y.data()), one);          \
      } else                                                                                                         \
        Axpby<ExecSpace, AV, XV, BV, YV, YV::rank, false, ETI_SPEC_AVAIL>::axpby(space, alpha, X, beta, Y);          \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS1_CAXPBY_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                                  \
  template <class ExecSpace>                                                                                       \
  struct Axpby<ExecSpace, Kokkos::complex<float>,                                                                  \
               Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,            \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                             \
               Kokkos::complex<float>,                                                                             \
               Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                  \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                             \
               1, true, ETI_SPEC_AVAIL> {                                                                          \
    typedef Kokkos::complex<float> AV;                                                                             \
    typedef Kokkos::complex<float> BV;                                                                             \
    typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,               \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                 \
        XV;                                                                                                        \
    typedef Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                 \
        YV;                                                                                                        \
                                                                                                                   \
    static void axpby(const ExecSpace& space, const AV& alpha, const XV& X, const BV& beta, const YV& Y) {         \
      Kokkos::Profiling::pushRegion("KokkosBlas::axpby[TPL_BLAS,complex<float>]");                                 \
      if ((X.extent(0) < INT_MAX) && (beta == 1.0f)) {                                                             \
        axpby_print_specialization<AV, XV, BV, YV>();                                                              \
        int N                               = X.extent(0);                                                         \
        int one                             = 1;                                                                   \
        const std::complex<float> alpha_val = alpha;                                                               \
        HostBlas<std::complex<float> >::axpy(N, alpha_val, reinterpret_cast<const std::complex<float>*>(X.data()), \
                                             one, reinterpret_cast<std::complex<float>*>(Y.data()), one);          \
      } else                                                                                                       \
        Axpby<ExecSpace, AV, XV, BV, YV, YV::rank, false, ETI_SPEC_AVAIL>::axpby(space, alpha, X, beta, Y);        \
      Kokkos::Profiling::popRegion();                                                                              \
    }                                                                                                              \
  };

KOKKOSBLAS1_DAXPBY_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_DAXPBY_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_SAXPBY_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_SAXPBY_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_ZAXPBY_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_ZAXPBY_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_CAXPBY_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_CAXPBY_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

#undef KOKKOSBLAS1_DAXPBY_BLAS
#undef KOKKOSBLAS1_SAXPBY_BLAS
#undef KOKKOSBLAS1_ZAXPBY_BLAS
#undef KOKKOSBLAS1_CAXPBY_BLAS
}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DAXPBY_CUBLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                                    \
  template <class ExecSpace>                                                                                           \
  struct Axpby<                                                                                                        \
      ExecSpace, double,                                                                                               \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      double,                                                                                                          \
      Kokkos::View<double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    typedef double AV;                                                                                                 \
    typedef double BV;                                                                                                 \
    typedef Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                   \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        XV;                                                                                                            \
    typedef Kokkos::View<double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                         \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        YV;                                                                                                            \
    typedef typename XV::size_type size_type;                                                                          \
                                                                                                                       \
    static void axpby(const ExecSpace& space, const AV& alpha, const XV& X, const BV& beta, const YV& Y) {             \
      Kokkos::Profiling::pushRegion("KokkosBlas::axpby[TPL_CUBLAS,double]");                                           \
      const size_type numElems = X.extent(0);                                                                          \
      if ((numElems < static_cast<size_type>(INT_MAX)) && (beta == 1.0)) {                                             \
        axpby_print_specialization<AV, XV, BV, YV>();                                                                  \
        const int N                            = static_cast<int>(numElems);                                           \
        constexpr int one                      = 1;                                                                    \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                     \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                  \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasDaxpy(s.handle, N, &alpha, X.data(), one, Y.data(), one));                  \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                                 \
      } else                                                                                                           \
        Axpby<ExecSpace, AV, XV, BV, YV, YV::rank, false, ETI_SPEC_AVAIL>::axpby(space, alpha, X, beta, Y);            \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS1_SAXPBY_CUBLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                                   \
  template <class ExecSpace>                                                                                          \
  struct Axpby<                                                                                                       \
      ExecSpace, float,                                                                                               \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      float,                                                                                                          \
      Kokkos::View<float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    typedef float AV;                                                                                                 \
    typedef float BV;                                                                                                 \
    typedef Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                   \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        XV;                                                                                                           \
    typedef Kokkos::View<float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                                         \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        YV;                                                                                                           \
    typedef typename XV::size_type size_type;                                                                         \
                                                                                                                      \
    static void axpby(const ExecSpace& space, const AV& alpha, const XV& X, const BV& beta, const YV& Y) {            \
      Kokkos::Profiling::pushRegion("KokkosBlas::axpby[TPL_CUBLAS,float]");                                           \
      const size_type numElems = X.extent(0);                                                                         \
      if ((numElems < static_cast<size_type>(INT_MAX)) && (beta == 1.0f)) {                                           \
        axpby_print_specialization<AV, XV, BV, YV>();                                                                 \
        const int N                            = static_cast<int>(numElems);                                          \
        constexpr int one                      = 1;                                                                   \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                    \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                 \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSaxpy(s.handle, N, &alpha, X.data(), one, Y.data(), one));                 \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                                \
      } else                                                                                                          \
        Axpby<ExecSpace, AV, XV, BV, YV, YV::rank, false, ETI_SPEC_AVAIL>::axpby(space, alpha, X, beta, Y);           \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS1_ZAXPBY_CUBLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                             \
  template <class ExecSpace>                                                                                    \
  struct Axpby<ExecSpace, Kokkos::complex<double>,                                                              \
               Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,        \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                          \
               Kokkos::complex<double>,                                                                         \
               Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,              \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                          \
               1, true, ETI_SPEC_AVAIL> {                                                                       \
    typedef Kokkos::complex<double> AV;                                                                         \
    typedef Kokkos::complex<double> BV;                                                                         \
    typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,           \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                              \
        XV;                                                                                                     \
    typedef Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                              \
        YV;                                                                                                     \
    typedef typename XV::size_type size_type;                                                                   \
                                                                                                                \
    static void axpby(const ExecSpace& space, const AV& alpha, const XV& X, const BV& beta, const YV& Y) {      \
      Kokkos::Profiling::pushRegion("KokkosBlas::axpby[TPL_CUBLAS,complex<double>]");                           \
      const size_type numElems = X.extent(0);                                                                   \
      if ((numElems < static_cast<size_type>(INT_MAX)) && (beta == 1.0f)) {                                     \
        axpby_print_specialization<AV, XV, BV, YV>();                                                           \
        const int N                            = static_cast<int>(numElems);                                    \
        constexpr int one                      = 1;                                                             \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();              \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                           \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasZaxpy(s.handle, N, reinterpret_cast<const cuDoubleComplex*>(&alpha), \
                                                 reinterpret_cast<const cuDoubleComplex*>(X.data()), one,       \
                                                 reinterpret_cast<cuDoubleComplex*>(Y.data()), one));           \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                          \
      } else                                                                                                    \
        Axpby<ExecSpace, AV, XV, BV, YV, YV::rank, false, ETI_SPEC_AVAIL>::axpby(space, alpha, X, beta, Y);     \
      Kokkos::Profiling::popRegion();                                                                           \
    }                                                                                                           \
  };

#define KOKKOSBLAS1_CAXPBY_CUBLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                         \
  template <class ExecSpace>                                                                                \
  struct Axpby<ExecSpace, Kokkos::complex<float>,                                                           \
               Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,     \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                      \
               Kokkos::complex<float>,                                                                      \
               Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,           \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                      \
               1, true, ETI_SPEC_AVAIL> {                                                                   \
    typedef Kokkos::complex<float> AV;                                                                      \
    typedef Kokkos::complex<float> BV;                                                                      \
    typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,        \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                          \
        XV;                                                                                                 \
    typedef Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,              \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                          \
        YV;                                                                                                 \
    typedef typename XV::size_type size_type;                                                               \
                                                                                                            \
    static void axpby(const ExecSpace& space, const AV& alpha, const XV& X, const BV& beta, const YV& Y) {  \
      Kokkos::Profiling::pushRegion("KokkosBlas::axpby[TPL_CUBLAS,complex<float>]");                        \
      const size_type numElems = X.extent(0);                                                               \
      if ((numElems < static_cast<size_type>(INT_MAX)) && (beta == 1.0f)) {                                 \
        axpby_print_specialization<AV, XV, BV, YV>();                                                       \
        const int N                            = static_cast<int>(numElems);                                \
        constexpr int one                      = 1;                                                         \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();          \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                       \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasCaxpy(s.handle, N, reinterpret_cast<const cuComplex*>(&alpha),   \
                                                 reinterpret_cast<const cuComplex*>(X.data()), one,         \
                                                 reinterpret_cast<cuComplex*>(Y.data()), one));             \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                      \
      } else                                                                                                \
        Axpby<ExecSpace, AV, XV, BV, YV, YV::rank, false, ETI_SPEC_AVAIL>::axpby(space, alpha, X, beta, Y); \
      Kokkos::Profiling::popRegion();                                                                       \
    }                                                                                                       \
  };

KOKKOSBLAS1_DAXPBY_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DAXPBY_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_SAXPBY_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SAXPBY_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_ZAXPBY_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZAXPBY_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_CAXPBY_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_CAXPBY_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

#undef KOKKOSBLAS1_DAXPBY_CUBLAS
#undef KOKKOSBLAS1_SAXPBY_CUBLAS
#undef KOKKOSBLAS1_ZAXPBY_CUBLAS
#undef KOKKOSBLAS1_CAXPBY_CUBLAS
}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#endif
