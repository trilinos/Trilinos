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

#ifndef KOKKOSBLAS1_ROTG_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_ROTG_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
template <class Scalar, class ExecutionSpace>
inline void rotg_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
  printf("KokkosBlas1::rotg<> TPL Blas specialization for < %s, %s >\n", typeid(Scalar).name(),
         typeid(ExecutionSpace).name);
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

#define KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_BLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                           \
  template <>                                                                                                       \
  struct Rotg<                                                                                                      \
      EXECSPACE,                                                                                                    \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      true, ETI_SPEC_AVAIL> {                                                                                       \
    using SViewType =                                                                                               \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using MViewType =                                                                                               \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    static void rotg(EXECSPACE const, SViewType const& a, SViewType const& b, MViewType const& c,                   \
                     SViewType const& s) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_BLAS,double]");                                           \
      HostBlas<double>::rotg(a.data(), b.data(), c.data(), s.data());                                               \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_BLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                              \
  template <>                                                                                                          \
  struct Rotg<                                                                                                         \
      EXECSPACE,                                                                                                       \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,       \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, \
      ETI_SPEC_AVAIL> {                                                                                                \
    using SViewType =                                                                                                  \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    using MViewType =                                                                                                  \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    static void rotg(EXECSPACE const, SViewType const& a, SViewType const& b, MViewType const& c,                      \
                     SViewType const& s) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_BLAS,float]");                                               \
      HostBlas<float>::rotg(a.data(), b.data(), c.data(), s.data());                                                   \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_BLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                           \
  template <>                                                                                                       \
  struct Rotg<                                                                                                      \
      EXECSPACE,                                                                                                    \
      Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                        \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      true, ETI_SPEC_AVAIL> {                                                                                       \
    using SViewType = Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,            \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                        \
    using MViewType =                                                                                               \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    static void rotg(EXECSPACE const, SViewType const& a, SViewType const& b, MViewType const& c,                   \
                     SViewType const& s) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_BLAS,complex<double>]");                                  \
      HostBlas<std::complex<double>>::rotg(reinterpret_cast<std::complex<double>*>(a.data()),                       \
                                           reinterpret_cast<std::complex<double>*>(b.data()), c.data(),             \
                                           reinterpret_cast<std::complex<double>*>(s.data()));                      \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_BLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                              \
  template <>                                                                                                          \
  struct Rotg<                                                                                                         \
      EXECSPACE,                                                                                                       \
      Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, \
      ETI_SPEC_AVAIL> {                                                                                                \
    using SViewType = Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                           \
    using MViewType =                                                                                                  \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    static void rotg(EXECSPACE const, SViewType const& a, SViewType const& b, MViewType const& c,                      \
                     SViewType const& s) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_BLAS,complex<float>]");                                      \
      HostBlas<std::complex<float>>::rotg(reinterpret_cast<std::complex<float>*>(a.data()),                            \
                                          reinterpret_cast<std::complex<float>*>(b.data()), c.data(),                  \
                                          reinterpret_cast<std::complex<float>*>(s.data()));                           \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)

KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)

KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)

KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)

KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)

KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)

KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)
#endif

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                         \
  template <>                                                                                                       \
  struct Rotg<                                                                                                      \
      EXECSPACE,                                                                                                    \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      true, ETI_SPEC_AVAIL> {                                                                                       \
    using SViewType =                                                                                               \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using MViewType =                                                                                               \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    static void rotg(EXECSPACE const& space, SViewType const& a, SViewType const& b, MViewType const& c,            \
                     SViewType const& s) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_CUBLAS,double]");                                         \
      rotg_print_specialization<double, EXECSPACE>();                                                               \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();            \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(singleton.handle, space.cuda_stream()));                         \
      cublasPointerMode_t pointer_mode;                                                                             \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasGetPointerMode(singleton.handle, &pointer_mode));                          \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));             \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasDrotg(singleton.handle, a.data(), b.data(), c.data(), s.data()));          \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetPointerMode(singleton.handle, pointer_mode));                           \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                            \
  template <>                                                                                                          \
  struct Rotg<                                                                                                         \
      EXECSPACE,                                                                                                       \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,       \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, \
      ETI_SPEC_AVAIL> {                                                                                                \
    using SViewType =                                                                                                  \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    using MViewType =                                                                                                  \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    static void rotg(EXECSPACE const& space, SViewType const& a, SViewType const& b, MViewType const& c,               \
                     SViewType const& s) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_CUBLAS,float]");                                             \
      rotg_print_specialization<float, EXECSPACE>();                                                                   \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();               \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(singleton.handle, space.cuda_stream()));                            \
      cublasPointerMode_t pointer_mode;                                                                                \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasGetPointerMode(singleton.handle, &pointer_mode));                             \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));                \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSrotg(singleton.handle, a.data(), b.data(), c.data(), s.data()));             \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetPointerMode(singleton.handle, pointer_mode));                              \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                         \
  template <>                                                                                                       \
  struct Rotg<                                                                                                      \
      EXECSPACE,                                                                                                    \
      Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                        \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      true, ETI_SPEC_AVAIL> {                                                                                       \
    using SViewType = Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,            \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                        \
    using MViewType =                                                                                               \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    static void rotg(EXECSPACE const& space, SViewType const& a, SViewType const& b, MViewType const& c,            \
                     SViewType const& s) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_CUBLAS,complex<double>]");                                \
      rotg_print_specialization<Kokkos::complex<double>, EXECSPACE>();                                              \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();            \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(singleton.handle, space.cuda_stream()));                         \
      cublasPointerMode_t pointer_mode;                                                                             \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasGetPointerMode(singleton.handle, &pointer_mode));                          \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));             \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasZrotg(singleton.handle, reinterpret_cast<cuDoubleComplex*>(a.data()),      \
                                               reinterpret_cast<cuDoubleComplex*>(b.data()), c.data(),              \
                                               reinterpret_cast<cuDoubleComplex*>(s.data())));                      \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetPointerMode(singleton.handle, pointer_mode));                           \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                            \
  template <>                                                                                                          \
  struct Rotg<                                                                                                         \
      EXECSPACE,                                                                                                       \
      Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, \
      ETI_SPEC_AVAIL> {                                                                                                \
    using SViewType = Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                           \
    using MViewType =                                                                                                  \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    static void rotg(EXECSPACE const& space, SViewType const& a, SViewType const& b, MViewType const& c,               \
                     SViewType const& s) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_CUBLAS,complex<float>]");                                    \
      rotg_print_specialization<Kokkos::complex<float>, EXECSPACE>();                                                  \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();               \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(singleton.handle, space.cuda_stream()));                            \
      cublasPointerMode_t pointer_mode;                                                                                \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasGetPointerMode(singleton.handle, &pointer_mode));                             \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));                \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasCrotg(singleton.handle, reinterpret_cast<cuComplex*>(a.data()),               \
                                               reinterpret_cast<cuComplex*>(b.data()), c.data(),                       \
                                               reinterpret_cast<cuComplex*>(s.data())));                               \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetPointerMode(singleton.handle, pointer_mode));                              \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                        \
  template <>                                                                                                       \
  struct Rotg<                                                                                                      \
      EXECSPACE,                                                                                                    \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      true, ETI_SPEC_AVAIL> {                                                                                       \
    using SViewType =                                                                                               \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using MViewType =                                                                                               \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    static void rotg(EXECSPACE const& space, SViewType const& a, SViewType const& b, MViewType const& c,            \
                     SViewType const& s) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_ROCBLAS,double]");                                        \
      rotg_print_specialization<double, EXECSPACE>();                                                               \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();              \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(singleton.handle, space.hip_stream()));                      \
      rocblas_pointer_mode pointer_mode;                                                                            \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_get_pointer_mode(singleton.handle, &pointer_mode));                     \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(singleton.handle, rocblas_pointer_mode_device));       \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_drotg(singleton.handle, a.data(), b.data(), c.data(), s.data()));       \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(singleton.handle, pointer_mode));                      \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                           \
  template <>                                                                                                          \
  struct Rotg<                                                                                                         \
      EXECSPACE,                                                                                                       \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,       \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, \
      ETI_SPEC_AVAIL> {                                                                                                \
    using SViewType =                                                                                                  \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    using MViewType =                                                                                                  \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    static void rotg(EXECSPACE const& space, SViewType const& a, SViewType const& b, MViewType const& c,               \
                     SViewType const& s) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_ROCBLAS,float]");                                            \
      rotg_print_specialization<float, EXECSPACE>();                                                                   \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();                 \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(singleton.handle, space.hip_stream()));                         \
      rocblas_pointer_mode pointer_mode;                                                                               \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_get_pointer_mode(singleton.handle, &pointer_mode));                        \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(singleton.handle, rocblas_pointer_mode_device));          \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_srotg(singleton.handle, a.data(), b.data(), c.data(), s.data()));          \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(singleton.handle, pointer_mode));                         \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                        \
  template <>                                                                                                       \
  struct Rotg<                                                                                                      \
      EXECSPACE,                                                                                                    \
      Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                        \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      true, ETI_SPEC_AVAIL> {                                                                                       \
    using SViewType = Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,            \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                        \
    using MViewType =                                                                                               \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    static void rotg(EXECSPACE const& space, SViewType const& a, SViewType const& b, MViewType const& c,            \
                     SViewType const& s) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_ROCBLAS,complex<double>]");                               \
      rotg_print_specialization<Kokkos::complex<double>, EXECSPACE>();                                              \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();              \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(singleton.handle, space.hip_stream()));                      \
      rocblas_pointer_mode pointer_mode;                                                                            \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_get_pointer_mode(singleton.handle, &pointer_mode));                     \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(singleton.handle, rocblas_pointer_mode_device));       \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_zrotg(singleton.handle,                                                 \
                                                  reinterpret_cast<rocblas_double_complex*>(a.data()),              \
                                                  reinterpret_cast<rocblas_double_complex*>(b.data()), c.data(),    \
                                                  reinterpret_cast<rocblas_double_complex*>(s.data())));            \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(singleton.handle, pointer_mode));                      \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                           \
  template <>                                                                                                          \
  struct Rotg<                                                                                                         \
      EXECSPACE,                                                                                                       \
      Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                                \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, \
      ETI_SPEC_AVAIL> {                                                                                                \
    using SViewType = Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                \
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                           \
    using MViewType =                                                                                                  \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    static void rotg(EXECSPACE const& space, SViewType const& a, SViewType const& b, MViewType const& c,               \
                     SViewType const& s) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_ROCBLAS,complex<float>]");                                   \
      rotg_print_specialization<Kokkos::complex<float>, EXECSPACE>();                                                  \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();                 \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(singleton.handle, space.hip_stream()));                         \
      rocblas_pointer_mode pointer_mode;                                                                               \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_get_pointer_mode(singleton.handle, &pointer_mode));                        \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(singleton.handle, rocblas_pointer_mode_device));          \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_crotg(singleton.handle,                                                    \
                                                  reinterpret_cast<rocblas_float_complex*>(a.data()),                  \
                                                  reinterpret_cast<rocblas_float_complex*>(b.data()), c.data(),        \
                                                  reinterpret_cast<rocblas_float_complex*>(s.data())));                \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(singleton.handle, pointer_mode));                         \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#endif
