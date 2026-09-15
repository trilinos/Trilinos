// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_ROTG_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_ROTG_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
template <class Scalar, class ExecutionSpace>
inline void rotg_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
  printf("KokkosBlas1::rotg<> TPL Blas specialization for < %s, %s >\n", typeid(Scalar).name(),
         typeid(ExecutionSpace).name());
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

#define KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_BLAS(LAYOUT)                                                       \
  template <typename ExecSpace, bool ETI_SPEC_AVAIL>                                                       \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                          \
  struct Rotg<ExecSpace, Kokkos::View<double, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              Kokkos::View<double, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,      \
              ETI_SPEC_AVAIL> {                                                                            \
    using SViewType = Kokkos::View<double, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;    \
    using MViewType = Kokkos::View<double, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;    \
    static void rotg(ExecSpace const, SViewType const& a, SViewType const& b, MViewType const& c,          \
                     SViewType const& s) {                                                                 \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_BLAS,double]");                                  \
      HostBlas<double>::rotg(a.data(), b.data(), c.data(), s.data());                                      \
      Kokkos::Profiling::popRegion();                                                                      \
    }                                                                                                      \
  };

#define KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_BLAS(LAYOUT)                                                                   \
  template <typename ExecSpace, bool ETI_SPEC_AVAIL>                                                                   \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                      \
  struct Rotg<ExecSpace, Kokkos::View<float, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,              \
              Kokkos::View<float, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, ETI_SPEC_AVAIL> { \
    using SViewType = Kokkos::View<float, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                 \
    using MViewType = Kokkos::View<float, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                 \
    static void rotg(ExecSpace const, SViewType const& a, SViewType const& b, MViewType const& c,                      \
                     SViewType const& s) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_BLAS,float]");                                               \
      HostBlas<float>::rotg(a.data(), b.data(), c.data(), s.data());                                                   \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_BLAS(LAYOUT)                                                                \
  template <typename ExecSpace, bool ETI_SPEC_AVAIL>                                                                \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                   \
  struct Rotg<                                                                                                      \
      ExecSpace, Kokkos::View<Kokkos::complex<double>, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<double, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, ETI_SPEC_AVAIL> {     \
    using SViewType =                                                                                               \
        Kokkos::View<Kokkos::complex<double>, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;          \
    using MViewType = Kokkos::View<double, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;             \
    static void rotg(ExecSpace const, SViewType const& a, SViewType const& b, MViewType const& c,                   \
                     SViewType const& s) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_BLAS,complex<double>]");                                  \
      HostBlas<std::complex<double>>::rotg(reinterpret_cast<std::complex<double>*>(a.data()),                       \
                                           reinterpret_cast<std::complex<double>*>(b.data()), c.data(),             \
                                           reinterpret_cast<std::complex<double>*>(s.data()));                      \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_BLAS(LAYOUT)                                                                   \
  template <typename ExecSpace, bool ETI_SPEC_AVAIL>                                                                   \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                      \
  struct Rotg<ExecSpace,                                                                                               \
              Kokkos::View<Kokkos::complex<float>, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,        \
              Kokkos::View<float, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, ETI_SPEC_AVAIL> { \
    using SViewType =                                                                                                  \
        Kokkos::View<Kokkos::complex<float>, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;              \
    using MViewType = Kokkos::View<float, LAYOUT, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                 \
    static void rotg(ExecSpace const, SViewType const& a, SViewType const& b, MViewType const& c,                      \
                     SViewType const& s) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_BLAS,complex<float>]");                                      \
      HostBlas<std::complex<float>>::rotg(reinterpret_cast<std::complex<float>*>(a.data()),                            \
                                          reinterpret_cast<std::complex<float>*>(b.data()), c.data(),                  \
                                          reinterpret_cast<std::complex<float>*>(s.data()));                           \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_BLAS(Kokkos::LayoutRight)
}  // namespace Impl
}  // namespace KokkosBlas

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_CUBLAS(LAYOUT)                                                           \
  template <bool ETI_SPEC_AVAIL>                                                                                 \
  struct Rotg<Kokkos::Cuda, Kokkos::View<double, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              Kokkos::View<double, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,         \
              ETI_SPEC_AVAIL> {                                                                                  \
    using SViewType = Kokkos::View<double, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    using MViewType = Kokkos::View<double, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    static void rotg(Kokkos::Cuda const& space, SViewType const& a, SViewType const& b, MViewType const& c,      \
                     SViewType const& s) {                                                                       \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_CUBLAS,double]");                                      \
      rotg_print_specialization<double, Kokkos::Cuda>();                                                         \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();         \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                  \
      cublasPointerMode_t pointer_mode;                                                                          \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasGetPointerMode(singleton.handle, &pointer_mode));                   \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));      \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasDrotg(singleton.handle, a.data(), b.data(), c.data(), s.data()));   \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, pointer_mode));                    \
      Kokkos::Profiling::popRegion();                                                                            \
    }                                                                                                            \
  };

#define KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_CUBLAS(LAYOUT)                                                          \
  template <bool ETI_SPEC_AVAIL>                                                                                \
  struct Rotg<Kokkos::Cuda, Kokkos::View<float, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              Kokkos::View<float, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,         \
              ETI_SPEC_AVAIL> {                                                                                 \
    using SViewType = Kokkos::View<float, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    using MViewType = Kokkos::View<float, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    static void rotg(Kokkos::Cuda const& space, SViewType const& a, SViewType const& b, MViewType const& c,     \
                     SViewType const& s) {                                                                      \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_CUBLAS,float]");                                      \
      rotg_print_specialization<float, Kokkos::Cuda>();                                                         \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();        \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                 \
      cublasPointerMode_t pointer_mode;                                                                         \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasGetPointerMode(singleton.handle, &pointer_mode));                  \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));     \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSrotg(singleton.handle, a.data(), b.data(), c.data(), s.data()));  \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, pointer_mode));                   \
      Kokkos::Profiling::popRegion();                                                                           \
    }                                                                                                           \
  };

#define KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_CUBLAS(LAYOUT)                                                              \
  template <bool ETI_SPEC_AVAIL>                                                                                    \
  struct Rotg<Kokkos::Cuda,                                                                                         \
              Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              Kokkos::View<double, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,            \
              ETI_SPEC_AVAIL> {                                                                                     \
    using SViewType =                                                                                               \
        Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    using MViewType = Kokkos::View<double, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;          \
    static void rotg(Kokkos::Cuda const& space, SViewType const& a, SViewType const& b, MViewType const& c,         \
                     SViewType const& s) {                                                                          \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_CUBLAS,complex<double>]");                                \
      rotg_print_specialization<Kokkos::complex<double>, Kokkos::Cuda>();                                           \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();            \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                     \
      cublasPointerMode_t pointer_mode;                                                                             \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasGetPointerMode(singleton.handle, &pointer_mode));                      \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));         \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasZrotg(singleton.handle, reinterpret_cast<cuDoubleComplex*>(a.data()),  \
                                                   reinterpret_cast<cuDoubleComplex*>(b.data()), c.data(),          \
                                                   reinterpret_cast<cuDoubleComplex*>(s.data())));                  \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, pointer_mode));                       \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_CUBLAS(LAYOUT)                                                             \
  template <bool ETI_SPEC_AVAIL>                                                                                   \
  struct Rotg<Kokkos::Cuda,                                                                                        \
              Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
              Kokkos::View<float, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,            \
              ETI_SPEC_AVAIL> {                                                                                    \
    using SViewType =                                                                                              \
        Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;       \
    using MViewType = Kokkos::View<float, LAYOUT, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;          \
    static void rotg(Kokkos::Cuda const& space, SViewType const& a, SViewType const& b, MViewType const& c,        \
                     SViewType const& s) {                                                                         \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_CUBLAS,complex<float>]");                                \
      rotg_print_specialization<Kokkos::complex<float>, Kokkos::Cuda>();                                           \
      KokkosBlas::Impl::CudaBlasSingleton& singleton = KokkosBlas::Impl::CudaBlasSingleton::singleton();           \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(singleton.handle, space.cuda_stream()));                    \
      cublasPointerMode_t pointer_mode;                                                                            \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasGetPointerMode(singleton.handle, &pointer_mode));                     \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, CUBLAS_POINTER_MODE_DEVICE));        \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasCrotg(singleton.handle, reinterpret_cast<cuComplex*>(a.data()),       \
                                                   reinterpret_cast<cuComplex*>(b.data()), c.data(),               \
                                                   reinterpret_cast<cuComplex*>(s.data())));                       \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(singleton.handle, pointer_mode));                      \
      Kokkos::Profiling::popRegion();                                                                              \
    }                                                                                                              \
  };

KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_ROCBLAS(LAYOUT)                                                           \
  template <bool ETI_SPEC_AVAIL>                                                                                  \
  struct Rotg<Kokkos::HIP, Kokkos::View<double, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
              Kokkos::View<double, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,           \
              ETI_SPEC_AVAIL> {                                                                                   \
    using SViewType = Kokkos::View<double, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;         \
    using MViewType = Kokkos::View<double, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;         \
    static void rotg(Kokkos::HIP const& space, SViewType const& a, SViewType const& b, MViewType const& c,        \
                     SViewType const& s) {                                                                        \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_ROCBLAS,double]");                                      \
      rotg_print_specialization<double, Kokkos::HIP>();                                                           \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();            \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(singleton.handle, space.hip_stream()));                \
      rocblas_pointer_mode pointer_mode;                                                                          \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_get_pointer_mode(singleton.handle, &pointer_mode));               \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(singleton.handle, rocblas_pointer_mode_device)); \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_drotg(singleton.handle, a.data(), b.data(), c.data(), s.data())); \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(singleton.handle, pointer_mode));                \
      Kokkos::Profiling::popRegion();                                                                             \
    }                                                                                                             \
  };

#define KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_ROCBLAS(LAYOUT)                                                           \
  template <bool ETI_SPEC_AVAIL>                                                                                  \
  struct Rotg<Kokkos::HIP, Kokkos::View<float, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,     \
              Kokkos::View<float, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,            \
              ETI_SPEC_AVAIL> {                                                                                   \
    using SViewType = Kokkos::View<float, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;          \
    using MViewType = Kokkos::View<float, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;          \
    static void rotg(Kokkos::HIP const& space, SViewType const& a, SViewType const& b, MViewType const& c,        \
                     SViewType const& s) {                                                                        \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_ROCBLAS,float]");                                       \
      rotg_print_specialization<float, Kokkos::HIP>();                                                            \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();            \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(singleton.handle, space.hip_stream()));                \
      rocblas_pointer_mode pointer_mode;                                                                          \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_get_pointer_mode(singleton.handle, &pointer_mode));               \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(singleton.handle, rocblas_pointer_mode_device)); \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_srotg(singleton.handle, a.data(), b.data(), c.data(), s.data())); \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(singleton.handle, pointer_mode));                \
      Kokkos::Profiling::popRegion();                                                                             \
    }                                                                                                             \
  };

#define KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_ROCBLAS(LAYOUT)                                                              \
  template <bool ETI_SPEC_AVAIL>                                                                                     \
  struct Rotg<Kokkos::HIP,                                                                                           \
              Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,   \
              Kokkos::View<double, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true,              \
              ETI_SPEC_AVAIL> {                                                                                      \
    using SViewType =                                                                                                \
        Kokkos::View<Kokkos::complex<double>, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;         \
    using MViewType = Kokkos::View<double, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;            \
    static void rotg(Kokkos::HIP const& space, SViewType const& a, SViewType const& b, MViewType const& c,           \
                     SViewType const& s) {                                                                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_ROCBLAS,complex<double>]");                                \
      rotg_print_specialization<Kokkos::complex<double>, Kokkos::HIP>();                                             \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();               \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(singleton.handle, space.hip_stream()));                   \
      rocblas_pointer_mode pointer_mode;                                                                             \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_get_pointer_mode(singleton.handle, &pointer_mode));                  \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(singleton.handle, rocblas_pointer_mode_device));    \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_zrotg(singleton.handle,                                              \
                                                      reinterpret_cast<rocblas_double_complex*>(a.data()),           \
                                                      reinterpret_cast<rocblas_double_complex*>(b.data()), c.data(), \
                                                      reinterpret_cast<rocblas_double_complex*>(s.data())));         \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(singleton.handle, pointer_mode));                   \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_ROCBLAS(LAYOUT)                                                                \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct Rotg<                                                                                                         \
      Kokkos::HIP, Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      Kokkos::View<float, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, true, ETI_SPEC_AVAIL> {       \
    using SViewType =                                                                                                  \
        Kokkos::View<Kokkos::complex<float>, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;            \
    using MViewType = Kokkos::View<float, LAYOUT, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;               \
    static void rotg(Kokkos::HIP const& space, SViewType const& a, SViewType const& b, MViewType const& c,             \
                     SViewType const& s) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotg[TPL_ROCBLAS,complex<float>]");                                   \
      rotg_print_specialization<Kokkos::complex<float>, Kokkos::HIP>();                                                \
      KokkosBlas::Impl::RocBlasSingleton& singleton = KokkosBlas::Impl::RocBlasSingleton::singleton();                 \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(singleton.handle, space.hip_stream()));                     \
      rocblas_pointer_mode pointer_mode;                                                                               \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_get_pointer_mode(singleton.handle, &pointer_mode));                    \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(singleton.handle, rocblas_pointer_mode_device));      \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_crotg(singleton.handle,                                                \
                                                      reinterpret_cast<rocblas_float_complex*>(a.data()),              \
                                                      reinterpret_cast<rocblas_float_complex*>(b.data()), c.data(),    \
                                                      reinterpret_cast<rocblas_float_complex*>(s.data())));            \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(singleton.handle, pointer_mode));                     \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_DROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_SROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_ZROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight)

KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft)
KOKKOSBLAS1_CROTG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#endif
