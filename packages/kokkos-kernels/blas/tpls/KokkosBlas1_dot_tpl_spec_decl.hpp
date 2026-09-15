// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_BLAS(KOKKOS_TYPE, TPL_TYPE)                                                      \
  template <typename ExecSpace, bool ETI_SPEC_AVAIL>                                                                   \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                      \
  struct Dot<                                                                                                          \
      ExecSpace,                                                                                                       \
      Kokkos::View<KOKKOS_TYPE, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,      \
      Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,       \
      Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, 1, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    typedef Kokkos::View<KOKKOS_TYPE, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
        RV;                                                                                                            \
    typedef Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >  \
        XV;                                                                                                            \
    typedef typename XV::size_type size_type;                                                                          \
                                                                                                                       \
    static void dot(const ExecSpace& space, RV& R, const XV& X, const XV& Y) {                                         \
      Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_BLAS," + KokkosKernels::ArithTraits<KOKKOS_TYPE>::name() +    \
                                    "]");                                                                              \
      const size_type numElems = X.extent(0);                                                                          \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                                \
        dot_print_specialization<RV, XV, XV>();                                                                        \
        int N   = numElems;                                                                                            \
        int one = 1;                                                                                                   \
        R()     = HostBlas<TPL_TYPE>::dot(N, reinterpret_cast<const TPL_TYPE*>(X.data()), one,                         \
                                          reinterpret_cast<const TPL_TYPE*>(Y.data()), one);                           \
      } else {                                                                                                         \
        Dot<ExecSpace, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R, X, Y);                                  \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_DOT_TPL_SPEC_DECL_BLAS(float, float)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_BLAS(double, double)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_BLAS(Kokkos::complex<float>, std::complex<float>)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_BLAS(Kokkos::complex<double>, std::complex<double>)

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
#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_CUBLAS(KOKKOS_TYPE, TPL_TYPE, TPL_DOT)                                           \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct Dot<                                                                                                          \
      Kokkos::Cuda,                                                                                                    \
      Kokkos::View<KOKKOS_TYPE, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,      \
      Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, \
      1, true, ETI_SPEC_AVAIL> {                                                                                       \
    typedef Kokkos::View<KOKKOS_TYPE, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
        RV;                                                                                                            \
    typedef Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::Cuda,                                         \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        XV;                                                                                                            \
    typedef typename XV::size_type size_type;                                                                          \
                                                                                                                       \
    static void dot(const Kokkos::Cuda& space, RV& R, const XV& X, const XV& Y) {                                      \
      Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_CUBLAS," + KokkosKernels::ArithTraits<KOKKOS_TYPE>::name() +  \
                                    "]");                                                                              \
      const size_type numElems = X.extent(0);                                                                          \
      /* TODO: CUDA-12's 64-bit indices allow larger numElems */                                                       \
      if (numElems <= static_cast<size_type>(std::numeric_limits<int>::max())) {                                       \
        dot_print_specialization<RV, XV, XV>();                                                                        \
        const int N                            = static_cast<int>(numElems);                                           \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                     \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                              \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(TPL_DOT(s.handle, N, reinterpret_cast<const TPL_TYPE*>(X.data()), 1,          \
                                                 reinterpret_cast<const TPL_TYPE*>(Y.data()), 1,                       \
                                                 reinterpret_cast<TPL_TYPE*>(&R())));                                  \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));                                             \
      } else {                                                                                                         \
        Dot<Kokkos::Cuda, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R, X, Y);                               \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_DOT_TPL_SPEC_DECL_CUBLAS(float, float, cublasSdot)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_CUBLAS(double, double, cublasDdot)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<float>, cuComplex, cublasCdotc)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<double>, cuDoubleComplex, cublasZdotc)

}  // namespace Impl
}  // namespace KokkosBlas
#endif
#endif

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {
#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ROCBLAS(KOKKOS_TYPE, TPL_TYPE, TPL_DOT)                                          \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct Dot<                                                                                                          \
      Kokkos::HIP,                                                                                                     \
      Kokkos::View<KOKKOS_TYPE, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,      \
      Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,     \
      Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1,  \
      1, true, ETI_SPEC_AVAIL> {                                                                                       \
    typedef Kokkos::View<KOKKOS_TYPE, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
        RV;                                                                                                            \
    typedef Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::HIP,                                          \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        XV;                                                                                                            \
    typedef typename XV::size_type size_type;                                                                          \
                                                                                                                       \
    static void dot(const Kokkos::HIP& space, RV& R, const XV& X, const XV& Y) {                                       \
      Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_ROCBLAS," + KokkosKernels::ArithTraits<KOKKOS_TYPE>::name() + \
                                    "]");                                                                              \
      const size_type numElems = X.extent(0);                                                                          \
      if (numElems <= static_cast<size_type>(std::numeric_limits<rocblas_int>::max())) {                               \
        dot_print_specialization<RV, XV, XV>();                                                                        \
        const rocblas_int N                   = static_cast<rocblas_int>(numElems);                                    \
        KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                       \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                           \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(TPL_DOT(s.handle, N, reinterpret_cast<const TPL_TYPE*>(X.data()), 1,         \
                                                  reinterpret_cast<const TPL_TYPE*>(Y.data()), 1,                      \
                                                  reinterpret_cast<TPL_TYPE*>(&R())));                                 \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                         \
      } else {                                                                                                         \
        Dot<Kokkos::HIP, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(space, R, X, Y);                                \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ROCBLAS(float, float, rocblas_sdot)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ROCBLAS(double, double, rocblas_ddot)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<float>, rocblas_float_complex, rocblas_cdotc)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<double>, rocblas_double_complex, rocblas_zdotc)

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
#define KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ONEMKL(KOKKOS_TYPE, TPL_TYPE, TPL_DOT)                                           \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct Dot<                                                                                                          \
      Kokkos::SYCL,                                                                                                    \
      Kokkos::View<KOKKOS_TYPE, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,      \
      Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::SYCL, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,    \
      Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::SYCL, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, \
      1, true, ETI_SPEC_AVAIL> {                                                                                       \
    typedef Kokkos::View<KOKKOS_TYPE, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
        RV;                                                                                                            \
    typedef Kokkos::View<const KOKKOS_TYPE*, Kokkos::LayoutLeft, Kokkos::SYCL,                                         \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        XV;                                                                                                            \
    typedef typename XV::size_type size_type;                                                                          \
                                                                                                                       \
    static void dot(const Kokkos::SYCL& exec, RV& R, const XV& X, const XV& Y) {                                       \
      Kokkos::Profiling::pushRegion("KokkosBlas::dot[TPL_ONEMKL," + KokkosKernels::ArithTraits<KOKKOS_TYPE>::name() +  \
                                    "]");                                                                              \
      const size_type numElems = X.extent(0);                                                                          \
      if (numElems <= static_cast<size_type>(std::numeric_limits<std::int64_t>::max())) {                              \
        dot_print_specialization<RV, XV, XV>();                                                                        \
        const std::int64_t N = static_cast<std::int64_t>(numElems);                                                    \
        Kokkos::View<typename RV::value_type, Kokkos::SYCL> device_result("device_result");                            \
        TPL_DOT(exec.sycl_queue(), N, reinterpret_cast<const TPL_TYPE*>(X.data()), 1,                                  \
                reinterpret_cast<const TPL_TYPE*>(Y.data()), 1, reinterpret_cast<TPL_TYPE*>(device_result.data()));    \
        Kokkos::deep_copy(R, device_result);                                                                           \
      } else {                                                                                                         \
        Dot<Kokkos::SYCL, RV, XV, XV, 1, 1, false, ETI_SPEC_AVAIL>::dot(exec, R, X, Y);                                \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ONEMKL(float, float, oneapi::mkl::blas::row_major::dot)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ONEMKL(double, double, oneapi::mkl::blas::row_major::dot)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ONEMKL(Kokkos::complex<float>, std::complex<float>, oneapi::mkl::blas::row_major::dotc)
KOKKOSBLAS1_DOT_TPL_SPEC_DECL_ONEMKL(Kokkos::complex<double>, std::complex<double>, oneapi::mkl::blas::row_major::dotc)

}  // namespace Impl
}  // namespace KokkosBlas
#endif

#endif
