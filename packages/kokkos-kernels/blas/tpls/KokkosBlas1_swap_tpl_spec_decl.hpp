/*
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
*/

#ifndef KOKKOSBLAS1_SWAP_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_SWAP_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
template <class ExecutionSpace, class XVector, class YVector>
inline void swap_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
  printf("KokkosBlas::swap<> TPL Blas specialization for < %s, %s, %s >\n", typeid(XVector).name(),
         typeid(YVector).name(), typeid(ExecutionSpace).name);
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

#define KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_BLAS(LAYOUT, EXECSPACE, ETI_SPEC_AVAIL)                 \
  template <>                                                                                   \
  struct Swap<EXECSPACE,                                                                        \
              Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,       \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                            \
              Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,       \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                            \
              true, ETI_SPEC_AVAIL> {                                                           \
    using XVector = Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>, \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                      \
    using YVector = Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>, \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                      \
    static void swap(EXECSPACE const& /*space*/, XVector const& X, YVector const& Y) {          \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_BLAS,double]");                       \
      HostBlas<double>::swap(X.extent_int(0), X.data(), 1, Y.data(), 1);                        \
      Kokkos::Profiling::popRegion();                                                           \
    }                                                                                           \
  };

#define KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_BLAS(LAYOUT, EXECSPACE, ETI_SPEC_AVAIL)                \
  template <>                                                                                  \
  struct Swap<EXECSPACE,                                                                       \
              Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,       \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                           \
              Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,       \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                           \
              true, ETI_SPEC_AVAIL> {                                                          \
    using XVector = Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>, \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                     \
    using YVector = Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>, \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                     \
    static void swap(EXECSPACE const& /*space*/, XVector const& X, YVector const& Y) {         \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_BLAS,float]");                       \
      HostBlas<float>::swap(X.extent_int(0), X.data(), 1, Y.data(), 1);                        \
      Kokkos::Profiling::popRegion();                                                          \
    }                                                                                          \
  };

#define KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_BLAS(LAYOUT, EXECSPACE, ETI_SPEC_AVAIL)                                   \
  template <>                                                                                                     \
  struct Swap<EXECSPACE,                                                                                          \
              Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,        \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
              Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,        \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                              \
              true, ETI_SPEC_AVAIL> {                                                                             \
    using XVector = Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,  \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                        \
    using YVector = Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,  \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                        \
    static void swap(EXECSPACE const& /*space*/, XVector const& X, YVector const& Y) {                            \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_BLAS,complex<double>]");                                \
      HostBlas<std::complex<double>>::swap(X.extent_int(0), reinterpret_cast<std::complex<double>*>(X.data()), 1, \
                                           reinterpret_cast<std::complex<double>*>(Y.data()), 1);                 \
      Kokkos::Profiling::popRegion();                                                                             \
    }                                                                                                             \
  };

#define KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_BLAS(LAYOUT, EXECSPACE, ETI_SPEC_AVAIL)                                 \
  template <>                                                                                                   \
  struct Swap<EXECSPACE,                                                                                        \
              Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,       \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                            \
              Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,       \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                            \
              true, ETI_SPEC_AVAIL> {                                                                           \
    using XVector = Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>, \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                      \
    using YVector = Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>, \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                      \
    static void swap(EXECSPACE const& /*space*/, XVector const& X, YVector const& Y) {                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_BLAS,complex<float>]");                               \
      HostBlas<std::complex<float>>::swap(X.extent_int(0), reinterpret_cast<std::complex<float>*>(X.data()), 1, \
                                          reinterpret_cast<std::complex<float>*>(Y.data()), 1);                 \
      Kokkos::Profiling::popRegion();                                                                           \
    }                                                                                                           \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, true)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, false)

KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, true)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, false)

KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, true)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, false)

KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, true)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, false)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, true)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, false)

KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, true)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, false)

KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, true)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, false)

KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, true)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, false)
#endif

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                          \
  template <>                                                                                                        \
  struct Swap<                                                                                                       \
      EXECSPACE,                                                                                                     \
      Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      true, ETI_SPEC_AVAIL> {                                                                                        \
    using XVector =                                                                                                  \
        Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using YVector =                                                                                                  \
        Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    static void swap(EXECSPACE const& space, XVector const& X, YVector const& Y) {                                   \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_CUBLAS,double]");                                          \
      swap_print_specialization<EXECSPACE, XVector, YVector>();                                                      \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();             \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(singleton.handle, space.cuda_stream()));                          \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasDswap(singleton.handle, X.extent_int(0), X.data(), 1, Y.data(), 1));        \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                         \
  template <>                                                                                                       \
  struct Swap<                                                                                                      \
      EXECSPACE,                                                                                                    \
      Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      true, ETI_SPEC_AVAIL> {                                                                                       \
    using XVector =                                                                                                 \
        Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using YVector =                                                                                                 \
        Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    static void swap(EXECSPACE const& space, XVector const& X, YVector const& Y) {                                  \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_CUBLAS,float]");                                          \
      swap_print_specialization<EXECSPACE, XVector, YVector>();                                                     \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();            \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(singleton.handle, space.cuda_stream()));                         \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSswap(singleton.handle, X.extent_int(0), X.data(), 1, Y.data(), 1));       \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)              \
  template <>                                                                                            \
  struct Swap<EXECSPACE,                                                                                 \
              Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,        \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                     \
              Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,        \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                     \
              true, ETI_SPEC_AVAIL> {                                                                    \
    using XVector = Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,  \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                               \
    using YVector = Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,  \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                               \
    static void swap(EXECSPACE const& space, XVector const& X, YVector const& Y) {                       \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_CUBLAS,complex<double>]");                     \
      swap_print_specialization<EXECSPACE, XVector, YVector>();                                          \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(singleton.handle, space.cuda_stream()));              \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasZswap(singleton.handle, X.extent_int(0),                        \
                                               reinterpret_cast<cuDoubleComplex*>(X.data()), 1,          \
                                               reinterpret_cast<cuDoubleComplex*>(Y.data()), 1));        \
      Kokkos::Profiling::popRegion();                                                                    \
    }                                                                                                    \
  };

#define KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)              \
  template <>                                                                                            \
  struct Swap<EXECSPACE,                                                                                 \
              Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,         \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                     \
              Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,         \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                     \
              true, ETI_SPEC_AVAIL> {                                                                    \
    using XVector = Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,   \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                               \
    using YVector = Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,   \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                               \
    static void swap(EXECSPACE const& space, XVector const& X, YVector const& Y) {                       \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_CUBLAS,complex<float>]");                      \
      swap_print_specialization<EXECSPACE, XVector, YVector>();                                          \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(singleton.handle, space.cuda_stream()));              \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasCswap(singleton.handle, X.extent_int(0),                        \
                                               reinterpret_cast<cuComplex*>(X.data()), 1,                \
                                               reinterpret_cast<cuComplex*>(Y.data()), 1));              \
      Kokkos::Profiling::popRegion();                                                                    \
    }                                                                                                    \
  };

KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                         \
  template <>                                                                                                        \
  struct Swap<                                                                                                       \
      EXECSPACE,                                                                                                     \
      Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      true, ETI_SPEC_AVAIL> {                                                                                        \
    using XVector =                                                                                                  \
        Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using YVector =                                                                                                  \
        Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    static void swap(EXECSPACE const& space, XVector const& X, YVector const& Y) {                                   \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_ROCBLAS,double]");                                         \
      swap_print_specialization<EXECSPACE, XVector, YVector>();                                                      \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();               \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(singleton.handle, space.hip_stream()));                       \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_dswap(singleton.handle, X.extent_int(0), X.data(), 1, Y.data(), 1));     \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                        \
  template <>                                                                                                       \
  struct Swap<                                                                                                      \
      EXECSPACE,                                                                                                    \
      Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      true, ETI_SPEC_AVAIL> {                                                                                       \
    using XVector =                                                                                                 \
        Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using YVector =                                                                                                 \
        Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    static void swap(EXECSPACE const& space, XVector const& X, YVector const& Y) {                                  \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_ROCBLAS,float]");                                         \
      swap_print_specialization<EXECSPACE, XVector, YVector>();                                                     \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();              \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(singleton.handle, space.hip_stream()));                      \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_sswap(singleton.handle, X.extent_int(0), X.data(), 1, Y.data(), 1));    \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                \
  template <>                                                                                               \
  struct Swap<EXECSPACE,                                                                                    \
              Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,           \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
              Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,           \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
              true, ETI_SPEC_AVAIL> {                                                                       \
    using XVector = Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,     \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                  \
    using YVector = Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,     \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                  \
    static void swap(EXECSPACE const& space, XVector const& X, YVector const& Y) {                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_ROCBLAS,complex_double]");                        \
      swap_print_specialization<EXECSPACE, XVector, YVector>();                                             \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();      \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(singleton.handle, space.hip_stream()));              \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_zswap(singleton.handle, X.extent_int(0),                        \
                                                  reinterpret_cast<rocblas_double_complex*>(X.data()), 1,   \
                                                  reinterpret_cast<rocblas_double_complex*>(Y.data()), 1)); \
      Kokkos::Profiling::popRegion();                                                                       \
    }                                                                                                       \
  };

#define KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)               \
  template <>                                                                                              \
  struct Swap<EXECSPACE,                                                                                   \
              Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,           \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
              Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,           \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
              true, ETI_SPEC_AVAIL> {                                                                      \
    using XVector = Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,     \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                 \
    using YVector = Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,     \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                 \
    static void swap(EXECSPACE const& space, XVector const& X, YVector const& Y) {                         \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_ROCBLAS,complex_float]");                        \
      swap_print_specialization<EXECSPACE, XVector, YVector>();                                            \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();     \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(singleton.handle, space.hip_stream()));             \
      KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_cswap(singleton.handle, X.extent_int(0),                       \
                                                  reinterpret_cast<rocblas_float_complex*>(X.data()), 1,   \
                                                  reinterpret_cast<rocblas_float_complex*>(Y.data()), 1)); \
      Kokkos::Profiling::popRegion();                                                                      \
    }                                                                                                      \
  };

KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)
}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCBLAS

#endif
