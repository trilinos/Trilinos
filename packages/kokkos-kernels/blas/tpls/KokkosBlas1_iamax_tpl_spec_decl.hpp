// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_IAMAX_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_IAMAX_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {
template <class RV, class XV>
inline void iamax_print_specialization() {
#if defined(KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION)
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
  printf("KokkosBlas1::iamax<> TPL cuBLAS specialization for < %s , %s >\n", typeid(RV).name(), typeid(XV).name());
#elif defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)
  printf("KokkosBlas1::iamax<> TPL rocBLAS specialization for < %s , %s >\n", typeid(RV).name(), typeid(XV).name());
#else
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
  printf("KokkosBlas1::iamax<> TPL Blas specialization for < %s , %s >\n", typeid(RV).name(), typeid(XV).name());
#endif
#endif
#endif
}
}  // namespace Impl
}  // namespace KokkosBlas

// Generic Host side BLAS (could be MKL or whatever)
#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_BLAS(SCALAR_TYPE, BASE_SCALAR_TYPE)                                            \
  template <typename ExecSpace, bool ETI_SPEC_AVAIL>                                                                    \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                       \
  struct Iamax<                                                                                                         \
      ExecSpace,                                                                                                        \
      Kokkos::View<unsigned long, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,     \
      Kokkos::View<const SCALAR_TYPE*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1,     \
      true, ETI_SPEC_AVAIL> {                                                                                           \
    typedef Kokkos::View<unsigned long, Kokkos::LayoutLeft, Kokkos::HostSpace,                                          \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                      \
        RV;                                                                                                             \
    typedef Kokkos::View<const SCALAR_TYPE*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >   \
        XV;                                                                                                             \
    typedef typename XV::size_type size_type;                                                                           \
                                                                                                                        \
    static void iamax(const ExecSpace& space, RV& R, const XV& X) {                                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_BLAS," #SCALAR_TYPE "]");                                    \
      const size_type numElems = X.extent(0);                                                                           \
      if (numElems == 0) {                                                                                              \
        R() = 0;                                                                                                        \
        return;                                                                                                         \
      }                                                                                                                 \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                                 \
        iamax_print_specialization<RV, XV>();                                                                           \
        int N         = static_cast<int>(numElems);                                                                     \
        const int XST = X.stride(0);                                                                                    \
        const int LDX = (XST == 0) ? 1 : XST;                                                                           \
        int idx       = HostBlas<BASE_SCALAR_TYPE>::iamax(N, reinterpret_cast<const BASE_SCALAR_TYPE*>(X.data()), LDX); \
        R()           = static_cast<size_type>(idx);                                                                    \
      } else {                                                                                                          \
        Iamax<ExecSpace, RV, XV, 1, false, ETI_SPEC_AVAIL>::iamax(space, R, X);                                         \
      }                                                                                                                 \
      Kokkos::Profiling::popRegion();                                                                                   \
    }                                                                                                                   \
  };

KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_BLAS(double, double)
KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_BLAS(float, float)
KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_BLAS(Kokkos::complex<double>, std::complex<double>)
KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_BLAS(Kokkos::complex<float>, std::complex<float>)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// cuBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

using CUBLAS_DEVICE_TYPE = Kokkos::Cuda;

#define KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_CUBLAS_WRAPPER(SCALAR_TYPE, CUDA_SCALAR_TYPE, CUBLAS_FN, INDEX_TYPE,          \
                                                        EXEC_SPACE, RET_DEVICE_TYPE, CUBLAS_PTR_MODE_1,                \
                                                        CUBLAS_PTR_MODE_2)                                             \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct Iamax<                                                                                                        \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<INDEX_TYPE, Kokkos::LayoutLeft, RET_DEVICE_TYPE, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
      Kokkos::View<const SCALAR_TYPE*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    typedef Kokkos::View<INDEX_TYPE, Kokkos::LayoutLeft, RET_DEVICE_TYPE, Kokkos::MemoryTraits<Kokkos::Unmanaged> >    \
        RV;                                                                                                            \
    typedef Kokkos::View<const SCALAR_TYPE*, Kokkos::LayoutLeft, Kokkos::Cuda,                                         \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        XV;                                                                                                            \
    typedef typename XV::size_type size_type;                                                                          \
                                                                                                                       \
    static void iamax(const EXEC_SPACE& space, RV& R, const XV& X) {                                                   \
      Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_CUBLAS," #SCALAR_TYPE "]");                                 \
      const size_type numElems = X.extent(0);                                                                          \
      if (numElems == 0) {                                                                                             \
        Kokkos::deep_copy(R, 0);                                                                                       \
        return;                                                                                                        \
      }                                                                                                                \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                                \
        iamax_print_specialization<RV, XV>();                                                                          \
        const int N                            = static_cast<int>(numElems);                                           \
        const int XST                          = X.stride(0);                                                          \
        const int LDX                          = (XST == 0) ? 1 : XST;                                                 \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                     \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                              \
        cublasPointerMode_t prevPtrMode;                                                                               \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasGetPointerMode(s.handle, &prevPtrMode));                                \
        if (prevPtrMode == CUBLAS_PTR_MODE_2) {                                                                        \
          KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(s.handle, CUBLAS_PTR_MODE_1));                         \
        }                                                                                                              \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(CUBLAS_FN(s.handle, N, reinterpret_cast<const CUDA_SCALAR_TYPE*>(X.data()),   \
                                                   LDX, reinterpret_cast<int*>(R.data())));                            \
        if (prevPtrMode == CUBLAS_PTR_MODE_2) {                                                                        \
          KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetPointerMode(s.handle, CUBLAS_PTR_MODE_2));                         \
          KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));                                           \
        }                                                                                                              \
      } else {                                                                                                         \
        Iamax<EXEC_SPACE, RV, XV, 1, false, ETI_SPEC_AVAIL>::iamax(space, R, X);                                       \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_CUBLAS(SCALAR_TYPE, CUDA_SCALAR_TYPE, CUBLAS_FN, INDEX_TYPE)                 \
  KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_CUBLAS_WRAPPER(SCALAR_TYPE, CUDA_SCALAR_TYPE, CUBLAS_FN, INDEX_TYPE, Kokkos::Cuda, \
                                                  Kokkos::HostSpace, CUBLAS_POINTER_MODE_HOST,                        \
                                                  CUBLAS_POINTER_MODE_DEVICE)                                         \
  KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_CUBLAS_WRAPPER(SCALAR_TYPE, CUDA_SCALAR_TYPE, CUBLAS_FN, INDEX_TYPE, Kokkos::Cuda, \
                                                  CUBLAS_DEVICE_TYPE, CUBLAS_POINTER_MODE_DEVICE,                     \
                                                  CUBLAS_POINTER_MODE_HOST)

#define KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_CUBLAS(INDEX_TYPE) \
  KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_CUBLAS(double, double, cublasIdamax, INDEX_TYPE)

#define KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_CUBLAS(INDEX_TYPE) \
  KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_CUBLAS(float, float, cublasIsamax, INDEX_TYPE)

#define KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_CUBLAS(INDEX_TYPE) \
  KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<double>, cuDoubleComplex, cublasIzamax, INDEX_TYPE)

#define KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_CUBLAS(INDEX_TYPE) \
  KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<float>, cuComplex, cublasIcamax, INDEX_TYPE)

KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_CUBLAS(unsigned long)

KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_CUBLAS(unsigned long)

KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_CUBLAS(unsigned long)

KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_CUBLAS(unsigned long)

KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_CUBLAS(unsigned int)

KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_CUBLAS(unsigned int)

KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_CUBLAS(unsigned int)

KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_CUBLAS(unsigned int)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// rocBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

using ROCBLAS_DEVICE_TYPE = Kokkos::HIP;

#define KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_ROCBLAS_WRAPPER(SCALAR_TYPE, ROCBLAS_SCALAR_TYPE, ROCBLAS_FN, INDEX_TYPE,    \
                                                         RET_DEVICE_TYPE, ROCBLAS_PTR_MODE_1, ROCBLAS_PTR_MODE_2)     \
  template <bool ETI_SPEC_AVAIL>                                                                                      \
  struct Iamax<                                                                                                       \
      Kokkos::HIP,                                                                                                    \
      Kokkos::View<INDEX_TYPE, Kokkos::LayoutLeft, RET_DEVICE_TYPE, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,        \
      Kokkos::View<const SCALAR_TYPE*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    using execution_space = Kokkos::HIP;                                                                              \
    typedef Kokkos::View<INDEX_TYPE, Kokkos::LayoutLeft, RET_DEVICE_TYPE, Kokkos::MemoryTraits<Kokkos::Unmanaged> >   \
        RV;                                                                                                           \
    typedef Kokkos::View<const SCALAR_TYPE*, Kokkos::LayoutLeft, Kokkos::HIP,                                         \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        XV;                                                                                                           \
    typedef typename XV::size_type size_type;                                                                         \
                                                                                                                      \
    static void iamax(const execution_space& space, RV& R, const XV& X) {                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::iamax[TPL_ROCBLAS," #SCALAR_TYPE "]");                               \
      const size_type numElems = X.extent(0);                                                                         \
      if (numElems == 0) {                                                                                            \
        Kokkos::deep_copy(R, 0);                                                                                      \
        return;                                                                                                       \
      }                                                                                                               \
      if (numElems < static_cast<size_type>(INT_MAX)) {                                                               \
        iamax_print_specialization<RV, XV>();                                                                         \
        const int N                           = static_cast<int>(numElems);                                           \
        const int XST                         = X.stride(0);                                                          \
        const int LDX                         = (XST == 0) ? 1 : XST;                                                 \
        KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                      \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                          \
        rocblas_pointer_mode prevPtrMode;                                                                             \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_get_pointer_mode(s.handle, &prevPtrMode));                          \
        if (prevPtrMode == ROCBLAS_PTR_MODE_2) {                                                                      \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(s.handle, ROCBLAS_PTR_MODE_1));                  \
        }                                                                                                             \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(ROCBLAS_FN(s.handle, N,                                                     \
                                                     reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(X.data()), LDX,     \
                                                     reinterpret_cast<int*>(R.data())));                              \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, NULL));                                        \
        if (prevPtrMode == ROCBLAS_PTR_MODE_2) {                                                                      \
          KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(s.handle, ROCBLAS_PTR_MODE_2));                  \
        }                                                                                                             \
      } else {                                                                                                        \
        Iamax<execution_space, RV, XV, 1, false, ETI_SPEC_AVAIL>::iamax(space, R, X);                                 \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_ROCBLAS(SCALAR_TYPE, ROCBLAS_SCALAR_TYPE, ROCBLAS_FN, INDEX_TYPE)   \
  KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_ROCBLAS_WRAPPER(SCALAR_TYPE, ROCBLAS_SCALAR_TYPE, ROCBLAS_FN, INDEX_TYPE, \
                                                   Kokkos::HostSpace, rocblas_pointer_mode_host,             \
                                                   rocblas_pointer_mode_device)                              \
  KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_ROCBLAS_WRAPPER(SCALAR_TYPE, ROCBLAS_SCALAR_TYPE, ROCBLAS_FN, INDEX_TYPE, \
                                                   ROCBLAS_DEVICE_TYPE, rocblas_pointer_mode_device,         \
                                                   rocblas_pointer_mode_host)

#define KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_ROCBLAS(INDEX_TYPE) \
  KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_ROCBLAS(double, double, rocblas_idamax, INDEX_TYPE)

#define KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_ROCBLAS(INDEX_TYPE) \
  KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_ROCBLAS(float, float, rocblas_isamax, INDEX_TYPE)

#define KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_ROCBLAS(INDEX_TYPE) \
  KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<double>, rocblas_double_complex, rocblas_izamax, INDEX_TYPE)

#define KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_ROCBLAS(INDEX_TYPE) \
  KOKKOSBLAS1_XIAMAX_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<float>, rocblas_float_complex, rocblas_icamax, INDEX_TYPE)

KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_ROCBLAS(unsigned long)

KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_ROCBLAS(unsigned long)

KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_ROCBLAS(unsigned long)

KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_ROCBLAS(unsigned long)

KOKKOSBLAS1_DIAMAX_TPL_SPEC_DECL_ROCBLAS(unsigned int)

KOKKOSBLAS1_SIAMAX_TPL_SPEC_DECL_ROCBLAS(unsigned int)

KOKKOSBLAS1_ZIAMAX_TPL_SPEC_DECL_ROCBLAS(unsigned int)

KOKKOSBLAS1_CIAMAX_TPL_SPEC_DECL_ROCBLAS(unsigned int)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#endif
