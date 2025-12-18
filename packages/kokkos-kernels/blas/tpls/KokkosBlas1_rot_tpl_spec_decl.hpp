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
         typeid(ScalarView).name(), typeid(ExecutionSpace).name);
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

#define KOKKOSBLAS1_DROT_TPL_SPEC_DECL_BLAS(LAYOUT, EXECSPACE, ETI_SPEC_AVAIL)                                    \
  template <>                                                                                                     \
  struct Rot<EXECSPACE,                                                                                           \
             Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                          \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                               \
             Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                           \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                               \
             Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                           \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                               \
             true, ETI_SPEC_AVAIL> {                                                                              \
    using VectorView    = Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,             \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                  \
    using MagnitudeView = Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,              \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                  \
    using ScalarView    = Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,              \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                  \
    static void rot(EXECSPACE const& /*space*/, VectorView const& X, VectorView const& Y, MagnitudeView const& c, \
                    ScalarView const& s) {                                                                        \
      Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_BLAS,double]");                                          \
      HostBlas<double>::rot(X.extent_int(0), X.data(), 1, Y.data(), 1, c.data(), s.data());                       \
      Kokkos::Profiling::popRegion();                                                                             \
    }                                                                                                             \
  };

#define KOKKOSBLAS1_SROT_TPL_SPEC_DECL_BLAS(LAYOUT, EXECSPACE, ETI_SPEC_AVAIL)                                    \
  template <>                                                                                                     \
  struct Rot<EXECSPACE,                                                                                           \
             Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                           \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                               \
             Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                            \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                               \
             Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                            \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                               \
             true, ETI_SPEC_AVAIL> {                                                                              \
    using VectorView    = Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,              \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                  \
    using MagnitudeView = Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,               \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                  \
    using ScalarView    = Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,               \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                  \
    static void rot(EXECSPACE const& /*space*/, VectorView const& X, VectorView const& Y, MagnitudeView const& c, \
                    ScalarView const& s) {                                                                        \
      Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_BLAS,float]");                                           \
      HostBlas<float>::rot(X.extent_int(0), X.data(), 1, Y.data(), 1, c.data(), s.data());                        \
      Kokkos::Profiling::popRegion();                                                                             \
    }                                                                                                             \
  };

#define KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_BLAS(LAYOUT, EXECSPACE, ETI_SPEC_AVAIL)                                         \
  template <>                                                                                                          \
  struct Rot<EXECSPACE,                                                                                                \
             Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,              \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                                \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,               \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                    \
             true, ETI_SPEC_AVAIL> {                                                                                   \
    using VectorView    = Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>, \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                       \
    using MagnitudeView = Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                   \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                       \
    using ScalarView    = Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,  \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                       \
    static void rot(EXECSPACE const& /*space*/, VectorView const& X, VectorView const& Y, MagnitudeView const& c,      \
                    ScalarView const& s) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_BLAS,complex<double>]");                                      \
      HostBlas<std::complex<double>>::rot(X.extent_int(0), reinterpret_cast<std::complex<double>*>(X.data()), 1,       \
                                          reinterpret_cast<std::complex<double>*>(Y.data()), 1, c.data(),              \
                                          reinterpret_cast<std::complex<double>*>(s.data()));                          \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS1_CROT_TPL_SPEC_DECL_BLAS(LAYOUT, EXECSPACE, ETI_SPEC_AVAIL)                                        \
  template <>                                                                                                         \
  struct Rot<EXECSPACE,                                                                                               \
             Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,              \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
             Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                                \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
             Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,               \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                   \
             true, ETI_SPEC_AVAIL> {                                                                                  \
    using VectorView    = Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>, \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                      \
    using MagnitudeView = Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                   \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                      \
    using ScalarView    = Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,  \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                      \
    static void rot(EXECSPACE const& /*space*/, VectorView const& X, VectorView const& Y, MagnitudeView const& c,     \
                    ScalarView const& s) {                                                                            \
      Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_BLAS,complex<float>]");                                      \
      HostBlas<std::complex<float>>::rot(X.extent_int(0), reinterpret_cast<std::complex<float>*>(X.data()), 1,        \
                                         reinterpret_cast<std::complex<float>*>(Y.data()), 1, c.data(),               \
                                         reinterpret_cast<std::complex<float>*>(s.data()));                           \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSBLAS1_DROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, true)
KOKKOSBLAS1_DROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, false)

KOKKOSBLAS1_SROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, true)
KOKKOSBLAS1_SROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, false)

KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, true)
KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, false)

KOKKOSBLAS1_CROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, true)
KOKKOSBLAS1_CROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::Serial, false)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSBLAS1_DROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, true)
KOKKOSBLAS1_DROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, false)

KOKKOSBLAS1_SROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, true)
KOKKOSBLAS1_SROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, false)

KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, true)
KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, false)

KOKKOSBLAS1_CROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, true)
KOKKOSBLAS1_CROT_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::OpenMP, false)
#endif

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DROT_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                           \
  template <>                                                                                                        \
  struct Rot<                                                                                                        \
      EXECSPACE,                                                                                                     \
      Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      true, ETI_SPEC_AVAIL> {                                                                                        \
    using VectorView =                                                                                               \
        Kokkos::View<double*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using MagnitudeView =                                                                                            \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;  \
    using ScalarView =                                                                                               \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;  \
    static void rot(EXECSPACE const& space, VectorView const& X, VectorView const& Y, MagnitudeView const& c,        \
                    ScalarView const& s) {                                                                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_CUBLAS,double]");                                           \
      rot_print_specialization<EXECSPACE, VectorView, ScalarView>();                                                 \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();             \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                      \
      cublasPointerMode_t pointer_mode;                                                                              \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasGetPointerMode(singleton.handle, &pointer_mode));                       \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));          \
      cublasDrot(singleton.handle, X.extent_int(0), X.data(), 1, Y.data(), 1, c.data(), s.data());                   \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, pointer_mode));                        \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS1_SROT_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                             \
  template <>                                                                                                          \
  struct Rot<                                                                                                          \
      EXECSPACE,                                                                                                       \
      Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,      \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,       \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, \
      ETI_SPEC_AVAIL> {                                                                                                \
    using VectorView =                                                                                                 \
        Kokkos::View<float*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;    \
    using MagnitudeView =                                                                                              \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    using ScalarView =                                                                                                 \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    static void rot(EXECSPACE const& space, VectorView const& X, VectorView const& Y, MagnitudeView const& c,          \
                    ScalarView const& s) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_CUBLAS,float]");                                              \
      rot_print_specialization<EXECSPACE, VectorView, ScalarView>();                                                   \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();               \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                        \
      cublasPointerMode_t pointer_mode;                                                                                \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasGetPointerMode(singleton.handle, &pointer_mode));                         \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));            \
      cublasSrot(singleton.handle, X.extent_int(0), X.data(), 1, Y.data(), 1, c.data(), s.data());                     \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, pointer_mode));                          \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                          \
  template <>                                                                                                       \
  struct Rot<                                                                                                       \
      EXECSPACE,                                                                                                    \
      Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                        \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                        \
      true, ETI_SPEC_AVAIL> {                                                                                       \
    using VectorView = Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,          \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                       \
    using MagnitudeView =                                                                                           \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using ScalarView = Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,           \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                       \
    static void rot(EXECSPACE const& space, VectorView const& X, VectorView const& Y, MagnitudeView const& c,       \
                    ScalarView const& s) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_CUBLAS,complex<double>]");                                 \
      rot_print_specialization<EXECSPACE, VectorView, ScalarView>();                                                \
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

#define KOKKOSBLAS1_CROT_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                         \
  template <>                                                                                                      \
  struct Rot<                                                                                                      \
      EXECSPACE,                                                                                                   \
      Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
      Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                       \
      true, ETI_SPEC_AVAIL> {                                                                                      \
    using VectorView = Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,          \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                      \
    using MagnitudeView =                                                                                          \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>; \
    using ScalarView = Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,           \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                      \
    static void rot(EXECSPACE const& space, VectorView const& X, VectorView const& Y, MagnitudeView const& c,      \
                    ScalarView const& s) {                                                                         \
      Kokkos::Profiling::pushRegion("KokkosBlas::rot[TPL_CUBLAS,complex<float>]");                                 \
      rot_print_specialization<EXECSPACE, VectorView, ScalarView>();                                               \
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

KOKKOSBLAS1_DROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_DROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_DROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_DROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS1_DROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_DROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_SROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_SROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_SROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_SROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS1_SROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_SROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_ZROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_CROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_CROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_CROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_CROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_CROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_CROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS1_CROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_CROT_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#endif
