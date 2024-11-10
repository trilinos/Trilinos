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

#define KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(SCALAR_TYPE, BASE_SCALAR_TYPE, LAYOUT, MEMSPACE, ETI_SPEC_AVAIL) \
  template <class ExecSpace>                                                                                  \
  struct Scal<ExecSpace,                                                                                      \
              Kokkos::View<SCALAR_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                         \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                         \
              SCALAR_TYPE,                                                                                    \
              Kokkos::View<const SCALAR_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                   \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                         \
              1, true, ETI_SPEC_AVAIL> {                                                                      \
    typedef Kokkos::View<SCALAR_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                           \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                            \
        RV;                                                                                                   \
    typedef SCALAR_TYPE AS;                                                                                   \
    typedef Kokkos::View<const SCALAR_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                            \
        XV;                                                                                                   \
    typedef typename XV::size_type size_type;                                                                 \
                                                                                                              \
    static void scal(const ExecSpace& space, const RV& R, const AS& alpha, const XV& X) {                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_BLAS," #SCALAR_TYPE "]");                           \
      const size_type numElems = X.extent(0);                                                                 \
      if ((numElems < static_cast<size_type>(INT_MAX)) && (R.data() == X.data())) {                           \
        scal_print_specialization<RV, AS, XV>();                                                              \
        int N                          = numElems;                                                            \
        int one                        = 1;                                                                   \
        const BASE_SCALAR_TYPE alpha_b = static_cast<BASE_SCALAR_TYPE>(alpha);                                \
        HostBlas<BASE_SCALAR_TYPE>::scal(N, alpha_b, reinterpret_cast<BASE_SCALAR_TYPE*>(R.data()), one);     \
      } else {                                                                                                \
        Scal<ExecSpace, RV, AS, XV, 1, false, ETI_SPEC_AVAIL>::scal(space, R, alpha, X);                      \
      }                                                                                                       \
      Kokkos::Profiling::popRegion();                                                                         \
    }                                                                                                         \
  };

#define KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(double, double, LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(float, float, LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(Kokkos::complex<double>, std::complex<double>, LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(Kokkos::complex<float>, std::complex<float>, LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)

KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// cuBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(SCALAR_TYPE, CUDA_SCALAR_TYPE, CUBLAS_FN, LAYOUT, MEMSPACE,     \
                                               ETI_SPEC_AVAIL)                                                 \
  template <class ExecSpace>                                                                                   \
  struct Scal<ExecSpace,                                                                                       \
              Kokkos::View<SCALAR_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                          \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                          \
              SCALAR_TYPE,                                                                                     \
              Kokkos::View<const SCALAR_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                    \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                          \
              1, true, ETI_SPEC_AVAIL> {                                                                       \
    typedef Kokkos::View<SCALAR_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                            \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                             \
        RV;                                                                                                    \
    typedef SCALAR_TYPE AS;                                                                                    \
    typedef Kokkos::View<const SCALAR_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,                      \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                             \
        XV;                                                                                                    \
    typedef typename XV::size_type size_type;                                                                  \
                                                                                                               \
    static void scal(const ExecSpace& space, const RV& R, const AS& alpha, const XV& X) {                      \
      Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_CUBLAS," #SCALAR_TYPE "]");                          \
      const size_type numElems = X.extent(0);                                                                  \
      if ((numElems < static_cast<size_type>(INT_MAX)) && (R.data() == X.data())) {                            \
        scal_print_specialization<RV, AS, XV>();                                                               \
        const int N                            = static_cast<int>(numElems);                                   \
        constexpr int one                      = 1;                                                            \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();             \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                          \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(CUBLAS_FN(s.handle, N, reinterpret_cast<const CUDA_SCALAR_TYPE*>(&alpha), \
                                               reinterpret_cast<CUDA_SCALAR_TYPE*>(R.data()), one));           \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                         \
      } else {                                                                                                 \
        Scal<ExecSpace, RV, AS, XV, 1, false, ETI_SPEC_AVAIL>::scal(space, R, alpha, X);                       \
      }                                                                                                        \
      Kokkos::Profiling::popRegion();                                                                          \
    }                                                                                                          \
  };

#define KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(double, double, cublasDscal, LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(float, float, cublasSscal, LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                                  \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<double>, cuDoubleComplex, cublasZscal, LAYOUT, MEMSPACE, \
                                         ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)                           \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<float>, cuComplex, cublasCscal, LAYOUT, MEMSPACE, \
                                         ETI_SPEC_AVAIL)

KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// rocBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(SCALAR_TYPE, ROCBLAS_SCALAR_TYPE, ROCBLAS_FN, LAYOUT, EXECSPACE,    \
                                                MEMSPACE, ETI_SPEC_AVAIL)                                           \
  template <>                                                                                                       \
  struct Scal<EXECSPACE,                                                                                            \
              Kokkos::View<SCALAR_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                               \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                               \
              SCALAR_TYPE,                                                                                          \
              Kokkos::View<const SCALAR_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                         \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                               \
              1, true, ETI_SPEC_AVAIL> {                                                                            \
    using execution_space = EXECSPACE;                                                                              \
    typedef Kokkos::View<SCALAR_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                  \
        RV;                                                                                                         \
    typedef SCALAR_TYPE AS;                                                                                         \
    typedef Kokkos::View<const SCALAR_TYPE*, LAYOUT, Kokkos::Device<EXECSPACE, MEMSPACE>,                           \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                  \
        XV;                                                                                                         \
    typedef typename XV::size_type size_type;                                                                       \
                                                                                                                    \
    static void scal(const execution_space& space, const RV& R, const AS& alpha, const XV& X) {                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_ROCBLAS," #SCALAR_TYPE "]");                              \
      const size_type numElems = X.extent(0);                                                                       \
      if ((numElems < static_cast<size_type>(INT_MAX)) && (R.data() == X.data())) {                                 \
        scal_print_specialization<RV, AS, XV>();                                                                    \
        const int N                           = static_cast<int>(numElems);                                         \
        constexpr int one                     = 1;                                                                  \
        KokkosBlas::Impl::RocBlasSingleton& s = KokkosBlas::Impl::RocBlasSingleton::singleton();                    \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_stream(s.handle, space.hip_stream()));                            \
        rocblas_pointer_mode pointer_mode;                                                                          \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_get_pointer_mode(s.handle, &pointer_mode));                           \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(s.handle, rocblas_pointer_mode_host));               \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(ROCBLAS_FN(s.handle, N, reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(&alpha), \
                                                 reinterpret_cast<ROCBLAS_SCALAR_TYPE*>(R.data()), one));           \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(rocblas_set_pointer_mode(s.handle, pointer_mode));                            \
      } else {                                                                                                      \
        Scal<EXECSPACE, RV, AS, XV, 1, false, ETI_SPEC_AVAIL>::scal(space, R, alpha, X);                            \
      }                                                                                                             \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

#define KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(double, double, rocblas_dscal, LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(float, float, rocblas_sscal, LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                      \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<double>, rocblas_double_complex, rocblas_zscal, LAYOUT, \
                                          EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_ROCBLAS(LAYOUT, EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)                    \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::complex<float>, rocblas_float_complex, rocblas_cscal, LAYOUT, \
                                          EXECSPACE, MEMSPACE, ETI_SPEC_AVAIL)

KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)

KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, true)
KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#endif
