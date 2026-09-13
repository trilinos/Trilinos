// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_SCAL_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_SCAL_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
template <class RV, class AS, class XV>
inline void scal_print_specialization() {
#if defined(KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION)
  printf("KokkosBlas1::scal<> TPL Blas specialization for < %s , %s , %s >\n", typeid(RV).name(), typeid(AS).name(),
         typeid(XV).name());
#endif
}
}  // namespace
}  // namespace Impl
}  // namespace KokkosBlas

#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(SCALAR_TYPE, BASE_SCALAR_TYPE)                                           \
  template <typename ExecSpace, bool ETI_SPEC_AVAIL>                                                                  \
    requires(std::is_same_v<typename ExecSpace::memory_space, Kokkos::HostSpace>)                                     \
  struct Scal<                                                                                                        \
      ExecSpace, Kokkos::View<SCALAR_TYPE*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      SCALAR_TYPE,                                                                                                    \
      Kokkos::View<const SCALAR_TYPE*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1,   \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    typedef Kokkos::View<SCALAR_TYPE*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV;   \
    typedef SCALAR_TYPE AS;                                                                                           \
    typedef Kokkos::View<const SCALAR_TYPE*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
        XV;                                                                                                           \
    typedef typename XV::size_type size_type;                                                                         \
                                                                                                                      \
    static void scal(const ExecSpace& space, const RV& R, const AS& alpha, const XV& X) {                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_BLAS," #SCALAR_TYPE "]");                                   \
      const size_type numElems = X.extent(0);                                                                         \
      if ((numElems < static_cast<size_type>(INT_MAX)) && (R.data() == X.data())) {                                   \
        scal_print_specialization<RV, AS, XV>();                                                                      \
        int N                          = numElems;                                                                    \
        int one                        = 1;                                                                           \
        const BASE_SCALAR_TYPE alpha_b = static_cast<BASE_SCALAR_TYPE>(alpha);                                        \
        HostBlas<BASE_SCALAR_TYPE>::scal(N, alpha_b, reinterpret_cast<BASE_SCALAR_TYPE*>(R.data()), one);             \
      } else {                                                                                                        \
        Scal<ExecSpace, RV, AS, XV, 1, false, ETI_SPEC_AVAIL>::scal(space, R, alpha, X);                              \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(double, double)

KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(float, float)

KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(Kokkos::complex<double>, std::complex<double>)

KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(Kokkos::complex<float>, std::complex<float>)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// cuBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(SCALAR_TYPE, CUDA_SCALAR_TYPE, CUBLAS_FN)                               \
  template <bool ETI_SPEC_AVAIL>                                                                                       \
  struct Scal<                                                                                                         \
      Kokkos::Cuda,                                                                                                    \
      Kokkos::View<SCALAR_TYPE*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,          \
      SCALAR_TYPE,                                                                                                     \
      Kokkos::View<const SCALAR_TYPE*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    typedef Kokkos::View<SCALAR_TYPE*, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
    typedef SCALAR_TYPE AS;                                                                                            \
    typedef Kokkos::View<const SCALAR_TYPE*, Kokkos::LayoutLeft, Kokkos::Cuda,                                         \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        XV;                                                                                                            \
    typedef typename XV::size_type size_type;                                                                          \
                                                                                                                       \
    static void scal(const Kokkos::Cuda& space, const RV& R, const AS& alpha, const XV& X) {                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_CUBLAS," #SCALAR_TYPE "]");                                  \
      const size_type numElems = X.extent(0);                                                                          \
      if ((numElems < static_cast<size_type>(INT_MAX)) && (R.data() == X.data())) {                                    \
        scal_print_specialization<RV, AS, XV>();                                                                       \
        const int N                            = static_cast<int>(numElems);                                           \
        constexpr int one                      = 1;                                                                    \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                     \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, space.cuda_stream()));                              \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(CUBLAS_FN(s.handle, N, reinterpret_cast<const CUDA_SCALAR_TYPE*>(&alpha),     \
                                                   reinterpret_cast<CUDA_SCALAR_TYPE*>(R.data()), one));               \
        KOKKOSBLAS_IMPL_CUBLAS_SAFE_CALL(cublasSetStream(s.handle, NULL));                                             \
      } else {                                                                                                         \
        Scal<Kokkos::Cuda, RV, AS, XV, 1, false, ETI_SPEC_AVAIL>::scal(space, R, alpha, X);                            \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(double, double, cublasDscal)

KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(float, float, cublasSscal)

KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<double>, cuDoubleComplex, cublasZscal)

KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<float>, cuComplex, cublasCscal)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// rocBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(SCALAR_TYPE, ROCBLAS_SCALAR_TYPE, ROCBLAS_FN)                         \
  template <bool ETI_SPEC_AVAIL>                                                                                      \
  struct Scal<                                                                                                        \
      Kokkos::HIP,                                                                                                    \
      Kokkos::View<SCALAR_TYPE*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,          \
      SCALAR_TYPE,                                                                                                    \
      Kokkos::View<const SCALAR_TYPE*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, 1, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    typedef Kokkos::View<SCALAR_TYPE*, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
    typedef SCALAR_TYPE AS;                                                                                           \
    typedef Kokkos::View<const SCALAR_TYPE*, Kokkos::LayoutLeft, Kokkos::HIP,                                         \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        XV;                                                                                                           \
    typedef typename XV::size_type size_type;                                                                         \
                                                                                                                      \
    static void scal(const Kokkos::HIP& space, const RV& R, const AS& alpha, const XV& X) {                           \
      Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_ROCBLAS," #SCALAR_TYPE "]");                                \
      const size_type numElems = X.extent(0);                                                                         \
      if ((numElems < static_cast<size_type>(INT_MAX)) && (R.data() == X.data())) {                                   \
        scal_print_specialization<RV, AS, XV>();                                                                      \
        const int N                           = static_cast<int>(numElems);                                           \
        constexpr int one                     = 1;                                                                    \
        KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                      \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_stream(s.handle, space.hip_stream()));                          \
        rocblas_pointer_mode pointer_mode;                                                                            \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_get_pointer_mode(s.handle, &pointer_mode));                         \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(s.handle, rocblas_pointer_mode_host));             \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(ROCBLAS_FN(s.handle, N,                                                     \
                                                     reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(&alpha),            \
                                                     reinterpret_cast<ROCBLAS_SCALAR_TYPE*>(R.data()), one));         \
        KOKKOSBLAS_IMPL_ROCBLAS_SAFE_CALL(rocblas_set_pointer_mode(s.handle, pointer_mode));                          \
      } else {                                                                                                        \
        Scal<Kokkos::HIP, RV, AS, XV, 1, false, ETI_SPEC_AVAIL>::scal(space, R, alpha, X);                            \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(double, double, rocblas_dscal)

KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(float, float, rocblas_sscal)

KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<double>, rocblas_double_complex, rocblas_zscal)

KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<float>, rocblas_float_complex, rocblas_cscal)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#endif
