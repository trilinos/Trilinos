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
  printf("KokkosBlas1::dot<> TPL Blas specialization for < %s , %s , %s >\n", typeid(RV).name(), typeid(XV).name(),
         typeid(YV).name());
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
#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_BLAS(LAYOUT, KOKKOS_TYPE, TPL_TYPE, MEMSPACE, ETI_SPEC_AVAIL)                \
  template <class ExecSpace>                                                                                       \
  struct Dot<ExecSpace,                                                                                            \
             Kokkos::View<KOKKOS_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,       \
             Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                         \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                               \
             Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                         \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                               \
             1, 1, true, ETI_SPEC_AVAIL> {                                                                         \
    typedef Kokkos::View<KOKKOS_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV;     \
    typedef Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                          \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                 \
        XV;                                                                                                        \
    typedef typename XV::size_type size_type;                                                                      \
                                                                                                                   \
    static void dot(const ExecSpace& space, RV& R, const XV& X, const XV& Y) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_BLAS," + Kokkos::ArithTraits<KOKKOS_TYPE>::name() + "]"); \
      const size_type numElems = X.extent(0);                                                                      \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                            \
        dot_print_specialization<RV, XV, XV>();                                                                    \
        int N   = numElems;                                                                                        \
        int one = 1;                                                                                               \
        R()     = HostBlas<TPL_TYPE>::dot(N, reinterpret_cast<const TPL_TYPE*>(X.data()), one,                     \
                                          reinterpret_cast<const TPL_TYPE*>(Y.data()), one);                       \
      } else {                                                                                                     \
        Dot<ExecSpace, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R, X, Y);                              \
      }                                                                                                            \
      Kokkos::Profiling::popRegion();                                                                              \
    }                                                                                                              \
  };

#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_BLAS_EXT(ETI_SPEC_AVAIL)                                              \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, float, float, Kokkos::HostSpace, ETI_SPEC_AVAIL)   \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, double, double, Kokkos::HostSpace, ETI_SPEC_AVAIL) \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::complex<float>, std::complex<float>,       \
                                     Kokkos::HostSpace, ETI_SPEC_AVAIL)                                     \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::complex<double>, std::complex<double>,     \
                                     Kokkos::HostSpace, ETI_SPEC_AVAIL)

KOKKOSBLAS1_DOT_TPL_SPEC_DECL_BLAS_EXT(true)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_BLAS_EXT(false)
}  // namespace Impl
}  // namespace KokkosBlas
#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
// Disabled because native has better performance.
// See tpl_spec_avail file for more details
#if 0
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {
#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_CUBLAS(LAYOUT, KOKKOS_TYPE, TPL_TYPE, EXECSPACE, MEMSPACE, TPL_DOT,            \
                                             ETI_SPEC_AVAIL)                                                         \
  template <>                                                                                                        \
  struct Dot<EXECSPACE,                                                                                              \
             Kokkos::View<KOKKOS_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
             Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                           \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
             Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                           \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
             1, 1, true, ETI_SPEC_AVAIL> {                                                                           \
    typedef Kokkos::View<KOKKOS_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV;       \
    typedef Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                            \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                   \
        XV;                                                                                                          \
    typedef typename XV::size_type size_type;                                                                        \
                                                                                                                     \
    static void dot(const EXECSPACE& space, RV& R, const XV& X, const XV& Y) {                                       \
      Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_CUBLAS," + Kokkos::ArithTraits<KOKKOS_TYPE>::name() + "]"); \
      const size_type numElems = X.extent(0);                                                                        \
      /* TODO: CUDA-12's 64-bit indices allow larger numElems */                                                     \
      if (numElems <= static_cast<size_type>(std::numeric_limits<int>::max())) {                                     \
        dot_print_specialization<RV, XV, XV>();                                                                      \
        const int N                            = static_cast<int>(numElems);                                         \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                   \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(TPL_DOT(s.handle, N, reinterpret_cast<const TPL_TYPE*>(X.data()), 1,            \
                                             reinterpret_cast<const TPL_TYPE*>(Y.data()), 1,                         \
                                             reinterpret_cast<TPL_TYPE*>(&R())));                                    \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                               \
      } else {                                                                                                       \
        Dot<EXECSPACE, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R, X, Y);                                \
      }                                                                                                              \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_CUBLAS_EXT(ETI_SPEC_AVAIL)                                                      \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, float, float, Kokkos::Cuda, Kokkos::CudaSpace, cublasSdot, \
                                       ETI_SPEC_AVAIL)                                                                \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, double, double, Kokkos::Cuda, Kokkos::CudaSpace,           \
                                       cublasDdot, ETI_SPEC_AVAIL)                                                    \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::complex<float>, cuComplex, Kokkos::Cuda,           \
                                       Kokkos::CudaSpace, cublasCdotc, ETI_SPEC_AVAIL)                                \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::complex<double>, cuDoubleComplex, Kokkos::Cuda,    \
                                       Kokkos::CudaSpace, cublasZdotc, ETI_SPEC_AVAIL)

KOKKOSBLAS1_DOT_TPL_SPEC_DECL_CUBLAS_EXT(true)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_CUBLAS_EXT(false)
}  // namespace Impl
}  // namespace KokkosBlas
#endif
#endif

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {
#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ROCBLAS(LAYOUT, KOKKOS_TYPE, TPL_TYPE, EXECSPACE, MEMSPACE, TPL_DOT,            \
                                              ETI_SPEC_AVAIL)                                                         \
  template <>                                                                                                         \
  struct Dot<EXECSPACE,                                                                                               \
             Kokkos::View<KOKKOS_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,          \
             Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                            \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                            \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             1, 1, true, ETI_SPEC_AVAIL> {                                                                            \
    typedef Kokkos::View<KOKKOS_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV;        \
    typedef Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                             \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        XV;                                                                                                           \
    typedef typename XV::size_type size_type;                                                                         \
                                                                                                                      \
    static void dot(const EXECSPACE& space, RV& R, const XV& X, const XV& Y) {                                        \
      Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_ROCBLAS," + Kokkos::ArithTraits<KOKKOS_TYPE>::name() + "]"); \
      const size_type numElems = X.extent(0);                                                                         \
      if (numElems <= static_cast<size_type>(std::numeric_limits<rocblas_int>::max())) {                              \
        dot_print_specialization<RV, XV, XV>();                                                                       \
        const rocblas_int N                   = static_cast<rocblas_int>(numElems);                                   \
        KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                      \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                              \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(TPL_DOT(s.handle, N, reinterpret_cast<const TPL_TYPE*>(X.data()), 1,            \
                                              reinterpret_cast<const TPL_TYPE*>(Y.data()), 1,                         \
                                              reinterpret_cast<TPL_TYPE*>(&R())));                                    \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, NULL));                                            \
      } else {                                                                                                        \
        Dot<EXECSPACE, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R, X, Y);                                 \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ROCBLAS_EXT(ETI_SPEC_AVAIL)                                                      \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, float, float, Kokkos::HIP, Kokkos::HIPSpace, rocblas_sdot, \
                                        ETI_SPEC_AVAIL)                                                                \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, double, double, Kokkos::HIP, Kokkos::HIPSpace,             \
                                        rocblas_ddot, ETI_SPEC_AVAIL)                                                  \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::complex<float>, rocblas_float_complex,             \
                                        Kokkos::HIP, Kokkos::HIPSpace, rocblas_cdotc, ETI_SPEC_AVAIL)                  \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::complex<double>, rocblas_double_complex,           \
                                        Kokkos::HIP, Kokkos::HIPSpace, rocblas_zdotc, ETI_SPEC_AVAIL)

KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ROCBLAS_EXT(true)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ROCBLAS_EXT(false)
}  // namespace Impl
}  // namespace KokkosBlas
#endif

// ONEMKL
#if defined(KOKKOSKERNELS_ENABLE_TPL_MKL) && defined(KOKKOS_ENABLE_SYCL)
#include <mkl.h>
#include <oneapi/mkl/blas.hpp>
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {
#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ONEMKL(LAYOUT, KOKKOS_TYPE, TPL_TYPE, EXECSPACE, MEMSPACE, TPL_DOT,            \
                                             ETI_SPEC_AVAIL)                                                         \
  template <>                                                                                                        \
  struct Dot<EXECSPACE,                                                                                              \
             Kokkos::View<KOKKOS_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
             Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                           \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
             Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                           \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                 \
             1, 1, true, ETI_SPEC_AVAIL> {                                                                           \
    typedef Kokkos::View<KOKKOS_TYPE, LAYOUT, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV;       \
    typedef Kokkos::View<const KOKKOS_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                            \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                   \
        XV;                                                                                                          \
    typedef typename XV::size_type size_type;                                                                        \
                                                                                                                     \
    static void dot(const EXECSPACE& exec, RV& R, const XV& X, const XV& Y) {                                        \
      Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_ONEMKL," + Kokkos::ArithTraits<KOKKOS_TYPE>::name() + "]"); \
      const size_type numElems = X.extent(0);                                                                        \
      if (numElems <= static_cast<size_type>(std::numeric_limits<std::int64_t>::max())) {                            \
        dot_print_specialization<RV, XV, XV>();                                                                      \
        const std::int64_t N = static_cast<std::int64_t>(numElems);                                                  \
        TPL_DOT(exec.sycl_queue(), N, reinterpret_cast<const TPL_TYPE*>(X.data()), 1,                                \
                reinterpret_cast<const TPL_TYPE*>(Y.data()), 1, reinterpret_cast<TPL_TYPE*>(&R()));                  \
      } else {                                                                                                       \
        Dot<EXECSPACE, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(exec, R, X, Y);                                 \
      }                                                                                                              \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ONEMKL_EXT(ETI_SPEC_AVAIL)                                                    \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ONEMKL(Kokkos::LayoutLeft, float, float, Kokkos::Experimental::SYCL,                \
                                       Kokkos::Experimental::SYCLDeviceUSMSpace, oneapi::mkl::blas::row_major::dot, \
                                       ETI_SPEC_AVAIL)                                                              \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ONEMKL(Kokkos::LayoutLeft, double, double, Kokkos::Experimental::SYCL,              \
                                       Kokkos::Experimental::SYCLDeviceUSMSpace, oneapi::mkl::blas::row_major::dot, \
                                       ETI_SPEC_AVAIL)                                                              \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ONEMKL(Kokkos::LayoutLeft, Kokkos::complex<float>, std::complex<float>,             \
                                       Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLDeviceUSMSpace,        \
                                       oneapi::mkl::blas::row_major::dotc, ETI_SPEC_AVAIL)                          \
  KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ONEMKL(Kokkos::LayoutLeft, Kokkos::complex<double>, std::complex<double>,           \
                                       Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLDeviceUSMSpace,        \
                                       oneapi::mkl::blas::row_major::dotc, ETI_SPEC_AVAIL)

KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ONEMKL_EXT(true)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ONEMKL_EXT(false)
}  // namespace Impl
}  // namespace KokkosBlas
#endif

#endif
