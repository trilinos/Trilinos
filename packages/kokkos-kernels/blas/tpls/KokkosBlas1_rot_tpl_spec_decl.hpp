// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_ROT_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_ROT_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
template <class ExecutionSpace, class VectorView, class ScalarView>
inline void rot_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
  printf("KokkosBlas::rot<> TPL Blas specialization for < %s, %s, %s >\n", typeid(VectorView).name(),
         typeid(ScalarView).name(), typeid(ExecutionSpace).name());
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
struct Rot<ExecSpace,
           Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
           Kokkos::View<double, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
           Kokkos::View<double, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,
           ETI_SPEC_AVAIL> {
  using VectorView =
      Kokkos::View<double*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using MagnitudeView =
      Kokkos::View<double, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using ScalarView =
      Kokkos::View<double, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  static void rot(ExecSpace const& /*space*/, VectorView const& X, VectorView const& Y, MagnitudeView const& c,
                  ScalarView const& s) {
    Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_BLAS,double]");
    HostBlas<double>::rot(X.extent_int(0), X.data(), 1, Y.data(), 1, c.data(), s.data());
    Kokkos::Profiling::popRegion();
  }
};

template <typename ExecSpace, bool ETI_SPEC_AVAIL>
  requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)
struct Rot<ExecSpace,
           Kokkos::View<float*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
           Kokkos::View<float, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
           Kokkos::View<float, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,
           ETI_SPEC_AVAIL> {
  using VectorView =
      Kokkos::View<float*, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using MagnitudeView =
      Kokkos::View<float, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using ScalarView =
      Kokkos::View<float, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  static void rot(ExecSpace const& /*space*/, VectorView const& X, VectorView const& Y, MagnitudeView const& c,
                  ScalarView const& s) {
    Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_BLAS,float]");
    HostBlas<float>::rot(X.extent_int(0), X.data(), 1, Y.data(), 1, c.data(), s.data());
    Kokkos::Profiling::popRegion();
  }
};

template <typename ExecSpace, bool ETI_SPEC_AVAIL>
  requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)
struct Rot<ExecSpace,
           Kokkos::View<Kokkos::complex<double>*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
           Kokkos::View<double, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
           Kokkos::View<Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HostSpace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
           true, ETI_SPEC_AVAIL> {
  using VectorView = Kokkos::View<Kokkos::complex<double>*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                                  Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using MagnitudeView =
      Kokkos::View<double, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using ScalarView = Kokkos::View<Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HostSpace,
                                  Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  static void rot(ExecSpace const& /*space*/, VectorView const& X, VectorView const& Y, MagnitudeView const& c,
                  ScalarView const& s) {
    Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_BLAS,complex<double>]");
    HostBlas<std::complex<double>>::rot(X.extent_int(0), reinterpret_cast<std::complex<double>*>(X.data()), 1,
                                        reinterpret_cast<std::complex<double>*>(Y.data()), 1, c.data(),
                                        reinterpret_cast<std::complex<double>*>(s.data()));
    Kokkos::Profiling::popRegion();
  }
};

template <typename ExecSpace, bool ETI_SPEC_AVAIL>
  requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)
struct Rot<ExecSpace,
           Kokkos::View<Kokkos::complex<float>*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
           Kokkos::View<float, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
           Kokkos::View<Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HostSpace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
           true, ETI_SPEC_AVAIL> {
  using VectorView = Kokkos::View<Kokkos::complex<float>*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                                  Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using MagnitudeView =
      Kokkos::View<float, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using ScalarView = Kokkos::View<Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::HostSpace,
                                  Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  static void rot(ExecSpace const& /*space*/, VectorView const& X, VectorView const& Y, MagnitudeView const& c,
                  ScalarView const& s) {
    Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_BLAS,complex<float>]");
    HostBlas<std::complex<float>>::rot(X.extent_int(0), reinterpret_cast<std::complex<float>*>(X.data()), 1,
                                       reinterpret_cast<std::complex<float>*>(Y.data()), 1, c.data(),
                                       reinterpret_cast<std::complex<float>*>(s.data()));
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

#define KOKKOSBLAS1_DROT_TPL_SPEC_DECL_CUBLAS(LAYOUT)                                                            \
  template <bool ETI_SPEC_AVAIL>                                                                                 \
  struct Rot<Kokkos::Cuda, Kokkos::View<double*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
             Kokkos::View<double, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                \
             Kokkos::View<double, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,          \
             ETI_SPEC_AVAIL> {                                                                                   \
    using VectorView    = Kokkos::View<double*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;  \
    using MagnitudeView = Kokkos::View<double, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;   \
    using ScalarView    = Kokkos::View<double, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;   \
    static void rot(Kokkos::Cuda const& space, VectorView const& X, VectorView const& Y, MagnitudeView const& c, \
                    ScalarView const& s) {                                                                       \
      Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_CUBLAS,double]");                                       \
      rot_print_specialization<Kokkos::Cuda, VectorView, ScalarView>();                                          \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();         \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                  \
      cublasPointerMode_t pointer_mode;                                                                          \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasGetPointerMode(singleton.handle, &pointer_mode));                   \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));      \
      cublasDrot(singleton.handle, X.extent_int(0), X.data(), 1, Y.data(), 1, c.data(), s.data());               \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, pointer_mode));                    \
      Kokkos::Profiling::popRegion();                                                                            \
    }                                                                                                            \
  };

#define KOKKOSBLAS1_SROT_TPL_SPEC_DECL_CUBLAS(LAYOUT)                                                            \
  template <bool ETI_SPEC_AVAIL>                                                                                 \
  struct Rot<Kokkos::Cuda, Kokkos::View<float*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
             Kokkos::View<float, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                 \
             Kokkos::View<float, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,           \
             ETI_SPEC_AVAIL> {                                                                                   \
    using VectorView    = Kokkos::View<float*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;   \
    using MagnitudeView = Kokkos::View<float, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;    \
    using ScalarView    = Kokkos::View<float, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;    \
    static void rot(Kokkos::Cuda const& space, VectorView const& X, VectorView const& Y, MagnitudeView const& c, \
                    ScalarView const& s) {                                                                       \
      Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_CUBLAS,float]");                                        \
      rot_print_specialization<Kokkos::Cuda, VectorView, ScalarView>();                                          \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();         \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                  \
      cublasPointerMode_t pointer_mode;                                                                          \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasGetPointerMode(singleton.handle, &pointer_mode));                   \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));      \
      cublasSrot(singleton.handle, X.extent_int(0), X.data(), 1, Y.data(), 1, c.data(), s.data());               \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, pointer_mode));                    \
      Kokkos::Profiling::popRegion();                                                                            \
    }                                                                                                            \
  };

#define KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_CUBLAS(LAYOUT)                                                               \
  template <bool ETI_SPEC_AVAIL>                                                                                    \
  struct Rot<Kokkos::Cuda,                                                                                          \
             Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
             Kokkos::View<double, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                   \
             Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
             true, ETI_SPEC_AVAIL> {                                                                                \
    using VectorView =                                                                                              \
        Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;      \
    using MagnitudeView = Kokkos::View<double, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;      \
    using ScalarView =                                                                                              \
        Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    static void rot(Kokkos::Cuda const& space, VectorView const& X, VectorView const& Y, MagnitudeView const& c,    \
                    ScalarView const& s) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_CUBLAS,complex<double>]");                                 \
      rot_print_specialization<Kokkos::Cuda, VectorView, ScalarView>();                                             \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();            \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                     \
      cublasPointerMode_t pointer_mode;                                                                             \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasGetPointerMode(singleton.handle, &pointer_mode));                      \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));         \
      cublasZrot(singleton.handle, X.extent_int(0), reinterpret_cast<cuDoubleComplex*>(X.data()), 1,                \
                 reinterpret_cast<cuDoubleComplex*>(Y.data()), 1, c.data(),                                         \
                 reinterpret_cast<cuDoubleComplex*>(s.data()));                                                     \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, pointer_mode));                       \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_CROT_TPL_SPEC_DECL_CUBLAS(LAYOUT)                                                              \
  template <bool ETI_SPEC_AVAIL>                                                                                   \
  struct Rot<Kokkos::Cuda,                                                                                         \
             Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
             Kokkos::View<float, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                   \
             Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,  \
             true, ETI_SPEC_AVAIL> {                                                                               \
    using VectorView =                                                                                             \
        Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;      \
    using MagnitudeView = Kokkos::View<float, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;      \
    using ScalarView =                                                                                             \
        Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    static void rot(Kokkos::Cuda const& space, VectorView const& X, VectorView const& Y, MagnitudeView const& c,   \
                    ScalarView const& s) {                                                                         \
      Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_CUBLAS,complex<float>]");                                 \
      rot_print_specialization<Kokkos::Cuda, VectorView, ScalarView>();                                            \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();           \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                    \
      cublasPointerMode_t pointer_mode;                                                                            \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasGetPointerMode(singleton.handle, &pointer_mode));                     \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));        \
      cublasCrot(singleton.handle, X.extent_int(0), reinterpret_cast<cuComplex*>(X.data()), 1,                     \
                 reinterpret_cast<cuComplex*>(Y.data()), 1, c.data(), reinterpret_cast<cuComplex*>(s.data()));     \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, pointer_mode));                      \
      Kokkos::Profiling::popRegion();                                                                              \
    }                                                                                                              \
  };

KOKKOSBLAS1_DROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_DROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_SROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_SROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_CROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_CROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight)
}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#endif
