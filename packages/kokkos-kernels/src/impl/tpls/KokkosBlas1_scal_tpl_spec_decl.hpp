/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSBLAS1_SCAL_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_SCAL_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
template <class RV, class AS, class XV>
inline void scal_print_specialization() {
#if defined(KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION)
  printf("KokkosBlas1::scal<> TPL Blas specialization for < %s , %s , %s >\n",
         typeid(RV).name(), typeid(AS).name(), typeid(XV).name());
#endif
}
}  // namespace
}  // namespace Impl
}  // namespace KokkosBlas

#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(SCALAR_TYPE, BASE_SCALAR_TYPE,    \
                                             LAYOUT, MEMSPACE, ETI_SPEC_AVAIL) \
  template <class ExecSpace>                                                   \
  struct Scal<                                                                 \
      Kokkos::View<SCALAR_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      SCALAR_TYPE,                                                             \
      Kokkos::View<const SCALAR_TYPE*, LAYOUT,                                 \
                   Kokkos::Device<ExecSpace, MEMSPACE>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      1, true, ETI_SPEC_AVAIL> {                                               \
    typedef Kokkos::View<SCALAR_TYPE*, LAYOUT,                                 \
                         Kokkos::Device<ExecSpace, MEMSPACE>,                  \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        RV;                                                                    \
    typedef SCALAR_TYPE AS;                                                    \
    typedef Kokkos::View<const SCALAR_TYPE*, LAYOUT,                           \
                         Kokkos::Device<ExecSpace, MEMSPACE>,                  \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XV;                                                                    \
    typedef typename XV::size_type size_type;                                  \
                                                                               \
    static void scal(const RV& R, const AS& alpha, const XV& X) {              \
      Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_BLAS," #SCALAR_TYPE  \
                                    "]");                                      \
      const size_type numElems = X.extent(0);                                  \
      if ((numElems < static_cast<size_type>(INT_MAX)) &&                      \
          (R.data() == X.data())) {                                            \
        scal_print_specialization<RV, AS, XV>();                               \
        int N                          = numElems;                             \
        int one                        = 1;                                    \
        const BASE_SCALAR_TYPE alpha_b = static_cast<BASE_SCALAR_TYPE>(alpha); \
        HostBlas<BASE_SCALAR_TYPE>::scal(                                      \
            N, alpha_b, reinterpret_cast<BASE_SCALAR_TYPE*>(R.data()), one);   \
      } else {                                                                 \
        Scal<RV, AS, XV, 1, false, ETI_SPEC_AVAIL>::scal(R, alpha, X);         \
      }                                                                        \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(double, double, LAYOUT, MEMSPACE,       \
                                       ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(float, float, LAYOUT, MEMSPACE,         \
                                       ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(Kokkos::complex<double>,                \
                                       std::complex<double>, LAYOUT, MEMSPACE, \
                                       ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_BLAS(LAYOUT, MEMSPACE, ETI_SPEC_AVAIL) \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_BLAS(Kokkos::complex<float>,                 \
                                       std::complex<float>, LAYOUT, MEMSPACE,  \
                                       ETI_SPEC_AVAIL)

KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace,
                                     true)
KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace,
                                     false)

KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace,
                                     true)
KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace,
                                     false)

KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace,
                                     true)
KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace,
                                     false)

KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace,
                                     true)
KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_BLAS(Kokkos::LayoutLeft, Kokkos::HostSpace,
                                     false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// cuBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(SCALAR_TYPE, CUDA_SCALAR_TYPE, \
                                               CUBLAS_FN, LAYOUT, MEMSPACE,   \
                                               ETI_SPEC_AVAIL)                \
  template <class ExecSpace>                                                  \
  struct Scal<                                                                \
      Kokkos::View<SCALAR_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                 \
      SCALAR_TYPE,                                                            \
      Kokkos::View<const SCALAR_TYPE*, LAYOUT,                                \
                   Kokkos::Device<ExecSpace, MEMSPACE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                 \
      1, true, ETI_SPEC_AVAIL> {                                              \
    typedef Kokkos::View<SCALAR_TYPE*, LAYOUT,                                \
                         Kokkos::Device<ExecSpace, MEMSPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        RV;                                                                   \
    typedef SCALAR_TYPE AS;                                                   \
    typedef Kokkos::View<const SCALAR_TYPE*, LAYOUT,                          \
                         Kokkos::Device<ExecSpace, MEMSPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        XV;                                                                   \
    typedef typename XV::size_type size_type;                                 \
                                                                              \
    static void scal(const RV& R, const AS& alpha, const XV& X) {             \
      Kokkos::Profiling::pushRegion(                                          \
          "KokkosBlas::scal[TPL_CUBLAS," #SCALAR_TYPE "]");                   \
      const size_type numElems = X.extent(0);                                 \
      if ((numElems < static_cast<size_type>(INT_MAX)) &&                     \
          (R.data() == X.data())) {                                           \
        scal_print_specialization<RV, AS, XV>();                              \
        const int N       = static_cast<int>(numElems);                       \
        constexpr int one = 1;                                                \
        KokkosBlas::Impl::CudaBlasSingleton& s =                              \
            KokkosBlas::Impl::CudaBlasSingleton::singleton();                 \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(CUBLAS_FN(                               \
            s.handle, N, reinterpret_cast<const CUDA_SCALAR_TYPE*>(&alpha),   \
            reinterpret_cast<CUDA_SCALAR_TYPE*>(R.data()), one));             \
      } else {                                                                \
        Scal<RV, AS, XV, 1, false, ETI_SPEC_AVAIL>::scal(R, alpha, X);        \
      }                                                                       \
      Kokkos::Profiling::popRegion();                                         \
    }                                                                         \
  };

#define KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS(LAYOUT, MEMSPACE,              \
                                               ETI_SPEC_AVAIL)                \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(double, double, cublasDscal, LAYOUT, \
                                         MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS(LAYOUT, MEMSPACE,            \
                                               ETI_SPEC_AVAIL)              \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(float, float, cublasSscal, LAYOUT, \
                                         MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS(LAYOUT, MEMSPACE,               \
                                               ETI_SPEC_AVAIL)                 \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<double>,              \
                                         cuDoubleComplex, cublasZscal, LAYOUT, \
                                         MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS(LAYOUT, MEMSPACE,            \
                                               ETI_SPEC_AVAIL)              \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::complex<float>, cuComplex, \
                                         cublasCscal, LAYOUT, MEMSPACE,     \
                                         ETI_SPEC_AVAIL)

KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace,
                                       true)
KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace,
                                       false)

KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace,
                                       true)
KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace,
                                       false)

KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace,
                                       true)
KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace,
                                       false)

KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace,
                                       true)
KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaSpace,
                                       false)

KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                                       true)
KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                                       false)

KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                                       true)
KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                                       false)

KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                                       true)
KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                                       false)

KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                                       true)
KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS(Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                                       false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

// rocBLAS
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)
#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(                               \
    SCALAR_TYPE, ROCBLAS_SCALAR_TYPE, ROCBLAS_FN, LAYOUT, MEMSPACE,            \
    ETI_SPEC_AVAIL)                                                            \
  template <class ExecSpace>                                                   \
  struct Scal<                                                                 \
      Kokkos::View<SCALAR_TYPE*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>,  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      SCALAR_TYPE,                                                             \
      Kokkos::View<const SCALAR_TYPE*, LAYOUT,                                 \
                   Kokkos::Device<ExecSpace, MEMSPACE>,                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      1, true, ETI_SPEC_AVAIL> {                                               \
    typedef Kokkos::View<SCALAR_TYPE*, LAYOUT,                                 \
                         Kokkos::Device<ExecSpace, MEMSPACE>,                  \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        RV;                                                                    \
    typedef SCALAR_TYPE AS;                                                    \
    typedef Kokkos::View<const SCALAR_TYPE*, LAYOUT,                           \
                         Kokkos::Device<ExecSpace, MEMSPACE>,                  \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        XV;                                                                    \
    typedef typename XV::size_type size_type;                                  \
                                                                               \
    static void scal(const RV& R, const AS& alpha, const XV& X) {              \
      Kokkos::Profiling::pushRegion(                                           \
          "KokkosBlas::scal[TPL_ROCBLAS," #SCALAR_TYPE "]");                   \
      const size_type numElems = X.extent(0);                                  \
      if ((numElems < static_cast<size_type>(INT_MAX)) &&                      \
          (R.data() == X.data())) {                                            \
        scal_print_specialization<RV, AS, XV>();                               \
        const int N       = static_cast<int>(numElems);                        \
        constexpr int one = 1;                                                 \
        KokkosBlas::Impl::RocBlasSingleton& s =                                \
            KokkosBlas::Impl::RocBlasSingleton::singleton();                   \
        KOKKOS_ROCBLAS_SAFE_CALL_IMPL(ROCBLAS_FN(                              \
            s.handle, N, reinterpret_cast<const ROCBLAS_SCALAR_TYPE*>(&alpha), \
            reinterpret_cast<ROCBLAS_SCALAR_TYPE*>(R.data()), one));           \
      } else {                                                                 \
        Scal<RV, AS, XV, 1, false, ETI_SPEC_AVAIL>::scal(R, alpha, X);         \
      }                                                                        \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_ROCBLAS(LAYOUT, MEMSPACE,        \
                                                ETI_SPEC_AVAIL)          \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(double, double, rocblas_dscal, \
                                          LAYOUT, MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_ROCBLAS(LAYOUT, MEMSPACE,              \
                                                ETI_SPEC_AVAIL)                \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(float, float, rocblas_sscal, LAYOUT, \
                                          MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_ROCBLAS(LAYOUT, MEMSPACE,             \
                                                ETI_SPEC_AVAIL)               \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(                                    \
      Kokkos::complex<double>, rocblas_double_complex, rocblas_zscal, LAYOUT, \
      MEMSPACE, ETI_SPEC_AVAIL)

#define KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_ROCBLAS(LAYOUT, MEMSPACE,           \
                                                ETI_SPEC_AVAIL)             \
  KOKKOSBLAS1_XSCAL_TPL_SPEC_DECL_ROCBLAS(                                  \
      Kokkos::complex<float>, rocblas_float_complex, rocblas_cscal, LAYOUT, \
      MEMSPACE, ETI_SPEC_AVAIL)

KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft,
                                        Kokkos::Experimental::HIPSpace, true)
KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft,
                                        Kokkos::Experimental::HIPSpace, false)

KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft,
                                        Kokkos::Experimental::HIPSpace, true)
KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft,
                                        Kokkos::Experimental::HIPSpace, false)

KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft,
                                        Kokkos::Experimental::HIPSpace, true)
KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft,
                                        Kokkos::Experimental::HIPSpace, false)

KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft,
                                        Kokkos::Experimental::HIPSpace, true)
KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_ROCBLAS(Kokkos::LayoutLeft,
                                        Kokkos::Experimental::HIPSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif

#endif
