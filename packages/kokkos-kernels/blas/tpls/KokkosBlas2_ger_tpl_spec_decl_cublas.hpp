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

#ifndef KOKKOSBLAS2_GER_TPL_SPEC_DECL_CUBLAS_HPP_
#define KOKKOSBLAS2_GER_TPL_SPEC_DECL_CUBLAS_HPP_

#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS2_GER_CUBLAS_DETERMINE_ARGS(LAYOUT)                        \
  bool A_is_ll      = std::is_same<Kokkos::LayoutLeft, LAYOUT>::value;       \
  bool A_is_lr      = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;      \
  const int M       = static_cast<int>(A_is_lr ? A.extent(1) : A.extent(0)); \
  const int N       = static_cast<int>(A_is_lr ? A.extent(0) : A.extent(1)); \
  constexpr int one = 1;                                                     \
  const int LDA     = A_is_lr ? A.stride(0) : A.stride(1);

#define KOKKOSBLAS2_DGER_CUBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                         \
  template <>                                                                                                          \
  struct GER<                                                                                                          \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const double*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<double**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    typedef double SCALAR;                                                                                             \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        XViewType;                                                                                                     \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        YViewType;                                                                                                     \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        AViewType;                                                                                                     \
                                                                                                                       \
    static void ger(const EXEC_SPACE& space, const char /*trans*/[], typename AViewType::const_value_type& alpha,      \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                      \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_CUBLAS,double]");                                             \
      KOKKOSBLAS2_GER_CUBLAS_DETERMINE_ARGS(LAYOUT);                                                                   \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                       \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                    \
      if (A_is_ll) {                                                                                                   \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasDger(s.handle, M, N, &alpha, X.data(), one, Y.data(), one, A.data(), LDA)); \
      } else {                                                                                                         \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasDger(s.handle, M, N, &alpha, Y.data(), one, X.data(), one, A.data(), LDA)); \
      }                                                                                                                \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                                   \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_SGER_CUBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                         \
  template <>                                                                                                          \
  struct GER<                                                                                                          \
      EXEC_SPACE,                                                                                                      \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<float**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >,  \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    typedef float SCALAR;                                                                                              \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        XViewType;                                                                                                     \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        YViewType;                                                                                                     \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                      \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                     \
        AViewType;                                                                                                     \
                                                                                                                       \
    static void ger(const EXEC_SPACE& space, const char /*trans*/[], typename AViewType::const_value_type& alpha,      \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                      \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_CUBLAS,float]");                                              \
      KOKKOSBLAS2_GER_CUBLAS_DETERMINE_ARGS(LAYOUT);                                                                   \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                       \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                    \
      if (A_is_ll) {                                                                                                   \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSger(s.handle, M, N, &alpha, X.data(), one, Y.data(), one, A.data(), LDA)); \
      } else {                                                                                                         \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSger(s.handle, M, N, &alpha, Y.data(), one, X.data(), one, A.data(), LDA)); \
      }                                                                                                                \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                                   \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_ZGER_CUBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                        \
  template <>                                                                                                         \
  struct GER<EXEC_SPACE,                                                                                              \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,              \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,              \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             true, ETI_SPEC_AVAIL> {                                                                                  \
    typedef Kokkos::complex<double> SCALAR;                                                                           \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        XViewType;                                                                                                    \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        YViewType;                                                                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        AViewType;                                                                                                    \
                                                                                                                      \
    static void ger(const EXEC_SPACE& space, const char trans[], typename AViewType::const_value_type& alpha,         \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_CUBLAS,complex<double>]");                                   \
      KOKKOSBLAS2_GER_CUBLAS_DETERMINE_ARGS(LAYOUT);                                                                  \
      bool justTranspose                     = (trans[0] == 'T') || (trans[0] == 't');                                \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                      \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                   \
      if (A_is_ll) {                                                                                                  \
        if (justTranspose) {                                                                                          \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasZgeru(s.handle, M, N, reinterpret_cast<const cuDoubleComplex*>(&alpha),  \
                                                   reinterpret_cast<const cuDoubleComplex*>(X.data()), one,           \
                                                   reinterpret_cast<const cuDoubleComplex*>(Y.data()), one,           \
                                                   reinterpret_cast<cuDoubleComplex*>(A.data()), LDA));               \
        } else {                                                                                                      \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasZgerc(s.handle, M, N, reinterpret_cast<const cuDoubleComplex*>(&alpha),  \
                                                   reinterpret_cast<const cuDoubleComplex*>(X.data()), one,           \
                                                   reinterpret_cast<const cuDoubleComplex*>(Y.data()), one,           \
                                                   reinterpret_cast<cuDoubleComplex*>(A.data()), LDA));               \
        }                                                                                                             \
      } else {                                                                                                        \
        if (justTranspose) {                                                                                          \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasZgeru(s.handle, M, N, reinterpret_cast<const cuDoubleComplex*>(&alpha),  \
                                                   reinterpret_cast<const cuDoubleComplex*>(Y.data()), one,           \
                                                   reinterpret_cast<const cuDoubleComplex*>(X.data()), one,           \
                                                   reinterpret_cast<cuDoubleComplex*>(A.data()), LDA));               \
        } else {                                                                                                      \
          /* cublasZgerc() + ~A_ll => call kokkos-kernels' implementation */                                          \
          GER<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::ger(space, trans, alpha, X, Y, A); \
        }                                                                                                             \
      }                                                                                                               \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                                  \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS2_CGER_CUBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                        \
  template <>                                                                                                         \
  struct GER<EXEC_SPACE,                                                                                              \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,               \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,               \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                    \
                          Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
             true, ETI_SPEC_AVAIL> {                                                                                  \
    typedef Kokkos::complex<float> SCALAR;                                                                            \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        XViewType;                                                                                                    \
    typedef Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        YViewType;                                                                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                     \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >                                                    \
        AViewType;                                                                                                    \
                                                                                                                      \
    static void ger(const EXEC_SPACE& space, const char trans[], typename AViewType::const_value_type& alpha,         \
                    const XViewType& X, const YViewType& Y, const AViewType& A) {                                     \
      Kokkos::Profiling::pushRegion("KokkosBlas::ger[TPL_CUBLAS,complex<float>]");                                    \
      KOKKOSBLAS2_GER_CUBLAS_DETERMINE_ARGS(LAYOUT);                                                                  \
      bool justTranspose                     = (trans[0] == 'T') || (trans[0] == 't');                                \
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                      \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                   \
      if (A_is_ll) {                                                                                                  \
        if (justTranspose) {                                                                                          \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasCgeru(s.handle, M, N, reinterpret_cast<const cuComplex*>(&alpha),        \
                                                   reinterpret_cast<const cuComplex*>(X.data()), one,                 \
                                                   reinterpret_cast<const cuComplex*>(Y.data()), one,                 \
                                                   reinterpret_cast<cuComplex*>(A.data()), LDA));                     \
        } else {                                                                                                      \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasCgerc(s.handle, M, N, reinterpret_cast<const cuComplex*>(&alpha),        \
                                                   reinterpret_cast<const cuComplex*>(X.data()), one,                 \
                                                   reinterpret_cast<const cuComplex*>(Y.data()), one,                 \
                                                   reinterpret_cast<cuComplex*>(A.data()), LDA));                     \
        }                                                                                                             \
      } else {                                                                                                        \
        if (justTranspose) {                                                                                          \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasCgeru(s.handle, M, N, reinterpret_cast<const cuComplex*>(&alpha),        \
                                                   reinterpret_cast<const cuComplex*>(Y.data()), one,                 \
                                                   reinterpret_cast<const cuComplex*>(X.data()), one,                 \
                                                   reinterpret_cast<cuComplex*>(A.data()), LDA));                     \
        } else {                                                                                                      \
          /* cublasCgerc() + ~A_ll => call kokkos-kernels' implementation */                                          \
          GER<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::ger(space, trans, alpha, X, Y, A); \
        }                                                                                                             \
      }                                                                                                               \
      KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                                  \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

KOKKOSBLAS2_DGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_DGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS2_DGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_DGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)

KOKKOSBLAS2_DGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_DGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS2_DGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_DGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS2_SGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_SGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS2_SGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_SGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)

KOKKOSBLAS2_SGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_SGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS2_SGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_SGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS2_ZGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_ZGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS2_ZGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_ZGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)

KOKKOSBLAS2_ZGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_ZGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS2_ZGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_ZGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS2_CGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_CGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS2_CGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_CGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)

KOKKOSBLAS2_CGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_CGER_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS2_CGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_CGER_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif
