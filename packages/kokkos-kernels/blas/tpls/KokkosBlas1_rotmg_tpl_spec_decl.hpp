// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
template <class Scalar>
inline void rotmg_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
  printf("KokkosBlas1::rotmg<> TPL Blas specialization for < %s >\n", typeid(Scalar).name());
#endif
}
}  // namespace
}  // namespace Impl
}  // namespace KokkosBlas

// Generic Host side BLAS (could be MKL or whatever)
#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS) && !defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL)
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                    \
  template <>                                                                                                          \
  struct Rotmg<                                                                                                        \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      Kokkos::View<SCALAR const, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<SCALAR[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    using DXView =                                                                                                     \
        Kokkos::View<SCALAR, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;  \
    using YView = Kokkos::View<SCALAR const, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                            \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                               \
    using PView = Kokkos::View<SCALAR[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                               \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                               \
    static void rotmg(EXEC_SPACE const& /* space */, DXView& d1, DXView& d2, DXView& x1, YView& y1, PView& param) {    \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotmg[TPL_BLAS,double]");                                             \
      HostBlas<SCALAR>::rotmg(d1.data(), d2.data(), x1.data(), y1.data(), param.data());                               \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(double, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(double, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(float, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(float, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace, false)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(float, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, true)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(float, Kokkos::LayoutRight, Kokkos::Serial, Kokkos::HostSpace, false)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(double, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(double, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(float, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(float, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace, false)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(float, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, true)
KOKKOSBLAS1_ROTMG_TPL_SPEC_DECL_BLAS(float, Kokkos::LayoutRight, Kokkos::OpenMP, Kokkos::HostSpace, false)
#endif

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DROTMG_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                         \
  template <>                                                                                                          \
  struct Rotmg<                                                                                                        \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      Kokkos::View<double const, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<double[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    using DXView =                                                                                                     \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;  \
    using YView = Kokkos::View<double const, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                            \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                               \
    using PView = Kokkos::View<double[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                               \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                               \
                                                                                                                       \
    static void rotmg(EXEC_SPACE const& space, DXView const& d1, DXView const& d2, DXView const& x1, YView const& y1,  \
                      PView const& param) {                                                                            \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotmg[TPL_CUBLAS,double]");                                           \
      rotmg_print_specialization<double>();                                                                            \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                       \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                                \
      cublasPointerMode_t pointer_mode;                                                                                \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasGetPointerMode(s.handle, &pointer_mode));                                 \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(s.handle, CUBLAS_POINTER_MODE_DEVICE));                    \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(                                                                                \
          cublasDrotmg(s.handle, d1.data(), d2.data(), x1.data(), y1.data(), param.data()));                           \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(s.handle, pointer_mode));                                  \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_DROTMG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DROTMG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DROTMG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_DROTMG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)

#define KOKKOSBLAS1_SROTMG_TPL_SPEC_DECL_CUBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                        \
  template <>                                                                                                         \
  struct Rotmg<                                                                                                       \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      Kokkos::View<float const, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<float[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    using DXView =                                                                                                    \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;  \
    using YView = Kokkos::View<float const, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                            \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                              \
    using PView = Kokkos::View<float[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                               \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                              \
                                                                                                                      \
    static void rotmg(EXEC_SPACE const& space, DXView const& d1, DXView const& d2, DXView const& x1, YView const& y1, \
                      PView const& param) {                                                                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotmg[TPL_CUBLAS,float]");                                           \
      rotmg_print_specialization<float>();                                                                            \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                      \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                               \
      cublasPointerMode_t pointer_mode;                                                                               \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasGetPointerMode(s.handle, &pointer_mode));                                \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(s.handle, CUBLAS_POINTER_MODE_DEVICE));                   \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(                                                                               \
          cublasSrotmg(s.handle, d1.data(), d2.data(), x1.data(), y1.data(), param.data()));                          \
      KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(s.handle, pointer_mode));                                 \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

KOKKOSBLAS1_SROTMG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SROTMG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SROTMG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS1_SROTMG_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// rocBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DROTMG_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                        \
  template <>                                                                                                          \
  struct Rotmg<                                                                                                        \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<double, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      Kokkos::View<double const, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<double[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    using DXView =                                                                                                     \
        Kokkos::View<double, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;  \
    using YView = Kokkos::View<double const, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                            \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                               \
    using PView = Kokkos::View<double[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                               \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                               \
                                                                                                                       \
    static void rotmg(EXEC_SPACE const& space, DXView const& d1, DXView const& d2, DXView const& x1, YView const& y1,  \
                      PView const& param) {                                                                            \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotmg[TPL_ROCBLAS,double]");                                          \
      rotmg_print_specialization<double>();                                                                            \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                         \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                             \
      rocblas_pointer_mode pointer_mode;                                                                               \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_get_pointer_mode(s.handle, &pointer_mode));                            \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(s.handle, rocblas_pointer_mode_device));              \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                               \
          rocblas_drotmg(s.handle, d1.data(), d2.data(), x1.data(), y1.data(), param.data()));                         \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(s.handle, pointer_mode));                             \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_DROTMG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_DROTMG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_DROTMG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS1_DROTMG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

#define KOKKOSBLAS1_SROTMG_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                       \
  template <>                                                                                                         \
  struct Rotmg<                                                                                                       \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<float, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      Kokkos::View<float const, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                          \
      Kokkos::View<float[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    using DXView =                                                                                                    \
        Kokkos::View<float, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;  \
    using YView = Kokkos::View<float const, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                            \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                              \
    using PView = Kokkos::View<float[5], LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                               \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                              \
                                                                                                                      \
    static void rotmg(EXEC_SPACE const& space, DXView const& d1, DXView const& d2, DXView const& x1, YView const& y1, \
                      PView const& param) {                                                                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::rotmg[TPL_ROCBLAS,float]");                                          \
      rotmg_print_specialization<float>();                                                                            \
      KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                        \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                            \
      rocblas_pointer_mode pointer_mode;                                                                              \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_get_pointer_mode(s.handle, &pointer_mode));                           \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(s.handle, rocblas_pointer_mode_device));             \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(                                                                              \
          rocblas_srotmg(s.handle, d1.data(), d2.data(), x1.data(), y1.data(), param.data()));                        \
      KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(s.handle, pointer_mode));                            \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

KOKKOSBLAS1_SROTMG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_SROTMG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_SROTMG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)
KOKKOSBLAS1_SROTMG_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutRight, Kokkos::HIP, Kokkos::HIPSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#endif
