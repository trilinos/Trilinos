// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_SWAP_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_SWAP_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
template <class ExecutionSpace, class XVector, class YVector>
inline void swap_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
  printf("KokkosBlas::swap<> TPL Blas specialization for < %s, %s, %s >\n", typeid(XVector).name(),
         typeid(YVector).name(), typeid(ExecutionSpace).name());
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

template <typename ExecSpace, bool ETI_SPEC_AVAIL>
  requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)
struct Swap<ExecSpace,
            Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
            Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,
            ETI_SPEC_AVAIL> {
  using XVector = Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using YVector = Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  static void swap(ExecSpace const& /*space*/, XVector const& X, YVector const& Y) {
    Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_BLAS,double]");
    HostBlas<double>::swap(X.extent_int(0), X.data(), 1, Y.data(), 1);
    Kokkos::Profiling::popRegion();
  }
};

template <typename ExecSpace, bool ETI_SPEC_AVAIL>
  requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)
struct Swap<ExecSpace,
            Kokkos::View<float*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
            Kokkos::View<float*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,
            ETI_SPEC_AVAIL> {
  using XVector = Kokkos::View<float*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using YVector = Kokkos::View<float*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  static void swap(ExecSpace const& /*space*/, XVector const& X, YVector const& Y) {
    Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_BLAS,float]");
    HostBlas<float>::swap(X.extent_int(0), X.data(), 1, Y.data(), 1);
    Kokkos::Profiling::popRegion();
  }
};

template <typename ExecSpace, bool ETI_SPEC_AVAIL>
  requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)
struct Swap<ExecSpace,
            Kokkos::View<Kokkos::complex<double>*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
            Kokkos::View<Kokkos::complex<double>*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
            true, ETI_SPEC_AVAIL> {
  using XVector = Kokkos::View<Kokkos::complex<double>*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using YVector = Kokkos::View<Kokkos::complex<double>*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  static void swap(ExecSpace const& /*space*/, XVector const& X, YVector const& Y) {
    Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_BLAS,complex<double>]");
    HostBlas<std::complex<double>>::swap(X.extent_int(0), reinterpret_cast<std::complex<double>*>(X.data()), 1,
                                         reinterpret_cast<std::complex<double>*>(Y.data()), 1);
    Kokkos::Profiling::popRegion();
  }
};

template <typename ExecSpace, bool ETI_SPEC_AVAIL>
  requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)
struct Swap<ExecSpace,
            Kokkos::View<Kokkos::complex<float>*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
            Kokkos::View<Kokkos::complex<float>*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
            true, ETI_SPEC_AVAIL> {
  using XVector = Kokkos::View<Kokkos::complex<float>*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using YVector = Kokkos::View<Kokkos::complex<float>*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  static void swap(ExecSpace const& /*space*/, XVector const& X, YVector const& Y) {
    Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_BLAS,complex<float>]");
    HostBlas<std::complex<float>>::swap(X.extent_int(0), reinterpret_cast<std::complex<float>*>(X.data()), 1,
                                        reinterpret_cast<std::complex<float>*>(Y.data()), 1);
    Kokkos::Profiling::popRegion();
  }
};

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_CUBLAS(LAYOUT)                                                            \
  template <bool ETI_SPEC_AVAIL>                                                                                  \
  struct Swap<Kokkos::Cuda, Kokkos::View<double*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              Kokkos::View<double*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,         \
              ETI_SPEC_AVAIL> {                                                                                   \
    using XVector = Kokkos::View<double*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;         \
    using YVector = Kokkos::View<double*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;         \
    static void swap(Kokkos::Cuda const& space, XVector const& X, YVector const& Y) {                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_CUBLAS,double]");                                       \
      swap_print_specialization<Kokkos::Cuda, XVector, YVector>();                                                \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();          \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                   \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasDswap(singleton.handle, X.extent_int(0), X.data(), 1, Y.data(), 1)); \
      Kokkos::Profiling::popRegion();                                                                             \
    }                                                                                                             \
  };

#define KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_CUBLAS(LAYOUT)                                                            \
  template <bool ETI_SPEC_AVAIL>                                                                                  \
  struct Swap<Kokkos::Cuda, Kokkos::View<float*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
              Kokkos::View<float*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,          \
              ETI_SPEC_AVAIL> {                                                                                   \
    using XVector = Kokkos::View<float*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;          \
    using YVector = Kokkos::View<float*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;          \
    static void swap(Kokkos::Cuda const& space, XVector const& X, YVector const& Y) {                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_CUBLAS,float]");                                        \
      swap_print_specialization<Kokkos::Cuda, XVector, YVector>();                                                \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();          \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                   \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSswap(singleton.handle, X.extent_int(0), X.data(), 1, Y.data(), 1)); \
      Kokkos::Profiling::popRegion();                                                                             \
    }                                                                                                             \
  };

#define KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_CUBLAS(LAYOUT)                                                               \
  template <bool ETI_SPEC_AVAIL>                                                                                     \
  struct Swap<Kokkos::Cuda,                                                                                          \
              Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              true, ETI_SPEC_AVAIL> {                                                                                \
    using XVector =                                                                                                  \
        Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    using YVector =                                                                                                  \
        Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    static void swap(Kokkos::Cuda const& space, XVector const& X, YVector const& Y) {                                \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_CUBLAS,complex<double>]");                                 \
      swap_print_specialization<Kokkos::Cuda, XVector, YVector>();                                                   \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();             \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                      \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasZswap(singleton.handle, X.extent_int(0),                                \
                                                   reinterpret_cast<cuDoubleComplex*>(X.data()), 1,                  \
                                                   reinterpret_cast<cuDoubleComplex*>(Y.data()), 1));                \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_CUBLAS(LAYOUT)                                                              \
  template <bool ETI_SPEC_AVAIL>                                                                                    \
  struct Swap<Kokkos::Cuda,                                                                                         \
              Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              true, ETI_SPEC_AVAIL> {                                                                               \
    using XVector =                                                                                                 \
        Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    using YVector =                                                                                                 \
        Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    static void swap(Kokkos::Cuda const& space, XVector const& X, YVector const& Y) {                               \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_CUBLAS,complex<float>]");                                 \
      swap_print_specialization<Kokkos::Cuda, XVector, YVector>();                                                  \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();            \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                     \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasCswap(singleton.handle, X.extent_int(0),                               \
                                                   reinterpret_cast<cuComplex*>(X.data()), 1,                       \
                                                   reinterpret_cast<cuComplex*>(Y.data()), 1));                     \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight)
}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_ROCBLAS(LAYOUT)                                                              \
  template <bool ETI_SPEC_AVAIL>                                                                                     \
  struct Swap<Kokkos::HIP, Kokkos::View<double*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,      \
              Kokkos::View<double*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,             \
              ETI_SPEC_AVAIL> {                                                                                      \
    using XVector = Kokkos::View<double*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;             \
    using YVector = Kokkos::View<double*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;             \
    static void swap(Kokkos::HIP const& space, XVector const& X, YVector const& Y) {                                 \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_ROCBLAS,double]");                                         \
      swap_print_specialization<Kokkos::HIP, XVector, YVector>();                                                    \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();               \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(singleton.handle, space.hip_stream()));                   \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_dswap(singleton.handle, X.extent_int(0), X.data(), 1, Y.data(), 1)); \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_ROCBLAS(LAYOUT)                                                              \
  template <bool ETI_SPEC_AVAIL>                                                                                     \
  struct Swap<Kokkos::HIP, Kokkos::View<float*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,       \
              Kokkos::View<float*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,              \
              ETI_SPEC_AVAIL> {                                                                                      \
    using XVector = Kokkos::View<float*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;              \
    using YVector = Kokkos::View<float*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;              \
    static void swap(Kokkos::HIP const& space, XVector const& X, YVector const& Y) {                                 \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_ROCBLAS,float]");                                          \
      swap_print_specialization<Kokkos::HIP, XVector, YVector>();                                                    \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();               \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(singleton.handle, space.hip_stream()));                   \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_sswap(singleton.handle, X.extent_int(0), X.data(), 1, Y.data(), 1)); \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_ROCBLAS(LAYOUT)                                                             \
  template <bool ETI_SPEC_AVAIL>                                                                                    \
  struct Swap<Kokkos::HIP,                                                                                          \
              Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              true, ETI_SPEC_AVAIL> {                                                                               \
    using XVector =                                                                                                 \
        Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    using YVector =                                                                                                 \
        Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    static void swap(Kokkos::HIP const& space, XVector const& X, YVector const& Y) {                                \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_ROCBLAS,complex_double]");                                \
      swap_print_specialization<Kokkos::HIP, XVector, YVector>();                                                   \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();              \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(singleton.handle, space.hip_stream()));                  \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_zswap(singleton.handle, X.extent_int(0),                            \
                                                      reinterpret_cast<rocblas_double_complex*>(X.data()), 1,       \
                                                      reinterpret_cast<rocblas_double_complex*>(Y.data()), 1));     \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_ROCBLAS(LAYOUT)                                                            \
  template <bool ETI_SPEC_AVAIL>                                                                                   \
  struct Swap<Kokkos::HIP,                                                                                         \
              Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              true, ETI_SPEC_AVAIL> {                                                                              \
    using XVector =                                                                                                \
        Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    using YVector =                                                                                                \
        Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    static void swap(Kokkos::HIP const& space, XVector const& X, YVector const& Y) {                               \
      Kokkos::Profiling::pushRegion("KokkosBlas::swap[TPL_ROCBLAS,complex_float]");                                \
      swap_print_specialization<Kokkos::HIP, XVector, YVector>();                                                  \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();             \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(singleton.handle, space.hip_stream()));                 \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_cswap(singleton.handle, X.extent_int(0),                           \
                                                      reinterpret_cast<rocblas_float_complex*>(X.data()), 1,       \
                                                      reinterpret_cast<rocblas_float_complex*>(Y.data()), 1));     \
      Kokkos::Profiling::popRegion();                                                                              \
    }                                                                                                              \
  };

KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_DSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_SSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_ZSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_CSWAP_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight)
}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCBLAS

#endif
