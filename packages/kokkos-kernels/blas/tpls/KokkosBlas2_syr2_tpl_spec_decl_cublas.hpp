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

#ifndef KOKKOSBLAS2_SYR2_TPL_SPEC_DECL_CUBLAS_HPP_
#define KOKKOSBLAS2_SYR2_TPL_SPEC_DECL_CUBLAS_HPP_

#include <KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS2_SYR2_CUBLAS_DETERMINE_ARGS(LAYOUT, uploChar)                     \
  bool A_is_ll              = std::is_same<Kokkos::LayoutLeft, LAYOUT>::value;       \
  bool A_is_lr              = std::is_same<Kokkos::LayoutRight, LAYOUT>::value;      \
  const int N               = static_cast<int>(A_is_lr ? A.extent(0) : A.extent(1)); \
  constexpr int one         = 1;                                                     \
  const int LDA             = A_is_lr ? A.stride(0) : A.stride(1);                   \
  cublasFillMode_t fillMode = (uploChar == 'L' || uploChar == 'l') ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;

#define KOKKOSBLAS2_DSYR2_CUBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                        \
  template <>                                                                                                          \
  struct SYR2<                                                                                                         \
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
    static void syr2(const typename AViewType::execution_space& space, const char trans[], const char uplo[],          \
                     typename AViewType::const_value_type& alpha, const XViewType& X, const YViewType& Y,              \
                     const AViewType& A) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr2[TPL_CUBLAS,double]");                                            \
      KOKKOSBLAS2_SYR2_CUBLAS_DETERMINE_ARGS(LAYOUT, uplo[0]);                                                         \
      if (A_is_ll) {                                                                                                   \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                     \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                  \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(                                                                                  \
            cublasDsyr2(s.handle, fillMode, N, &alpha, X.data(), one, Y.data(), one, A.data(), LDA));                  \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                                 \
      } else {                                                                                                         \
        /* cublasDsyr2() + ~A_ll => call kokkos-kernels' implementation */                                             \
        SYR2<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X,   \
                                                                                       Y, A);                          \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_SSYR2_CUBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                       \
  template <>                                                                                                         \
  struct SYR2<                                                                                                        \
      EXEC_SPACE,                                                                                                     \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<const float*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                         \
      Kokkos::View<float**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
      true, ETI_SPEC_AVAIL> {                                                                                         \
    typedef float SCALAR;                                                                                             \
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
    static void syr2(const typename AViewType::execution_space& space, const char trans[], const char uplo[],         \
                     typename AViewType::const_value_type& alpha, const XViewType& X, const YViewType& Y,             \
                     const AViewType& A) {                                                                            \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr2[TPL_CUBLAS,float]");                                            \
      KOKKOSBLAS2_SYR2_CUBLAS_DETERMINE_ARGS(LAYOUT, uplo[0]);                                                        \
      if (A_is_ll) {                                                                                                  \
        KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                    \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                 \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(                                                                                 \
            cublasSsyr2(s.handle, fillMode, N, &alpha, X.data(), one, Y.data(), one, A.data(), LDA));                 \
        KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                                \
      } else {                                                                                                        \
        /* cublasSsyr2() + ~A_ll => call kokkos-kernels' implementation */                                            \
        SYR2<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X,  \
                                                                                       Y, A);                         \
      }                                                                                                               \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

#define KOKKOSBLAS2_ZSYR2_CUBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                        \
  template <>                                                                                                          \
  struct SYR2<EXEC_SPACE,                                                                                              \
              Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
              Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
              Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                   \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
              true, ETI_SPEC_AVAIL> {                                                                                  \
    typedef Kokkos::complex<double> SCALAR;                                                                            \
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
    static void syr2(const typename AViewType::execution_space& space, const char trans[], const char uplo[],          \
                     typename AViewType::const_value_type& alpha, const XViewType& X, const YViewType& Y,              \
                     const AViewType& A) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr2[TPL_CUBLAS,complex<double>]");                                   \
      KOKKOSBLAS2_SYR2_CUBLAS_DETERMINE_ARGS(LAYOUT, uplo[0]);                                                         \
      bool justTranspose = (trans[0] == 'T') || (trans[0] == 't');                                                     \
      if (justTranspose) {                                                                                             \
        if (A_is_ll) {                                                                                                 \
          KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                   \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasZsyr2(s.handle, fillMode, N,                                              \
                                                   reinterpret_cast<const cuDoubleComplex*>(&alpha),                   \
                                                   reinterpret_cast<const cuDoubleComplex*>(X.data()), one,            \
                                                   reinterpret_cast<const cuDoubleComplex*>(Y.data()), one,            \
                                                   reinterpret_cast<cuDoubleComplex*>(A.data()), LDA));                \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                               \
        } else {                                                                                                       \
          /* cublasZsyr2() + ~A_ll => call kokkos-kernels' implementation */                                           \
          SYR2<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X, \
                                                                                         Y, A);                        \
        }                                                                                                              \
      } else {                                                                                                         \
        if (A_is_ll) {                                                                                                 \
          KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                   \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasZher2(s.handle, fillMode, N,                                              \
                                                   reinterpret_cast<const cuDoubleComplex*>(&alpha),                   \
                                                   reinterpret_cast<const cuDoubleComplex*>(X.data()), one,            \
                                                   reinterpret_cast<const cuDoubleComplex*>(Y.data()), one,            \
                                                   reinterpret_cast<cuDoubleComplex*>(A.data()), LDA));                \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                               \
        } else {                                                                                                       \
          /* cublasZher2() + ~A_ll => call kokkos-kernels' implementation */                                           \
          SYR2<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X, \
                                                                                         Y, A);                        \
        }                                                                                                              \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSBLAS2_CSYR2_CUBLAS(LAYOUT, EXEC_SPACE, MEM_SPACE, ETI_SPEC_AVAIL)                                        \
  template <>                                                                                                          \
  struct SYR2<EXEC_SPACE,                                                                                              \
              Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,               \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
              Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,               \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
              Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::Device<EXEC_SPACE, MEM_SPACE>,                    \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                  \
              true, ETI_SPEC_AVAIL> {                                                                                  \
    typedef Kokkos::complex<float> SCALAR;                                                                             \
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
    static void syr2(const typename AViewType::execution_space& space, const char trans[], const char uplo[],          \
                     typename AViewType::const_value_type& alpha, const XViewType& X, const YViewType& Y,              \
                     const AViewType& A) {                                                                             \
      Kokkos::Profiling::pushRegion("KokkosBlas::syr2[TPL_CUBLAS,complex<float>]");                                    \
      KOKKOSBLAS2_SYR2_CUBLAS_DETERMINE_ARGS(LAYOUT, uplo[0]);                                                         \
      bool justTranspose = (trans[0] == 'T') || (trans[0] == 't');                                                     \
      if (justTranspose) {                                                                                             \
        if (A_is_ll) {                                                                                                 \
          KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                   \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasCsyr2(s.handle, fillMode, N, reinterpret_cast<const cuComplex*>(&alpha),  \
                                                   reinterpret_cast<const cuComplex*>(X.data()), one,                  \
                                                   reinterpret_cast<const cuComplex*>(Y.data()), one,                  \
                                                   reinterpret_cast<cuComplex*>(A.data()), LDA));                      \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                               \
        } else {                                                                                                       \
          /* cublasCsyr2() + ~A_ll => call kokkos-kernels' implementation */                                           \
          SYR2<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X, \
                                                                                         Y, A);                        \
        }                                                                                                              \
      } else {                                                                                                         \
        if (A_is_ll) {                                                                                                 \
          KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();                   \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, space.cuda_stream()));                                \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasCher2(s.handle, fillMode, N, reinterpret_cast<const cuComplex*>(&alpha),  \
                                                   reinterpret_cast<const cuComplex*>(X.data()), one,                  \
                                                   reinterpret_cast<const cuComplex*>(Y.data()), one,                  \
                                                   reinterpret_cast<cuComplex*>(A.data()), LDA));                      \
          KOKKOS_CUBLAS_SAFE_CALL_IMPL(cublasSetStream(s.handle, NULL));                                               \
        } else {                                                                                                       \
          /* cublasCher2() + ~A_ll => call kokkos-kernels' implementation */                                           \
          SYR2<EXEC_SPACE, XViewType, YViewType, AViewType, false, ETI_SPEC_AVAIL>::syr2(space, trans, uplo, alpha, X, \
                                                                                         Y, A);                        \
        }                                                                                                              \
      }                                                                                                                \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

KOKKOSBLAS2_DSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_DSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS2_DSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_DSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)

KOKKOSBLAS2_DSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_DSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS2_DSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_DSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS2_SSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_SSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS2_SSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_SSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)

KOKKOSBLAS2_SSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_SSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS2_SSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_SSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS2_ZSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_ZSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS2_ZSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_ZSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)

KOKKOSBLAS2_ZSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_ZSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS2_ZSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_ZSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS2_CSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_CSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace, false)
KOKKOSBLAS2_CSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, true)
KOKKOSBLAS2_CSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaSpace, false)

KOKKOSBLAS2_CSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_CSYR2_CUBLAS(Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS2_CSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS2_CSYR2_CUBLAS(Kokkos::LayoutRight, Kokkos::Cuda, Kokkos::CudaUVMSpace, false)

}  // namespace Impl
}  // namespace KokkosBlas

#endif
