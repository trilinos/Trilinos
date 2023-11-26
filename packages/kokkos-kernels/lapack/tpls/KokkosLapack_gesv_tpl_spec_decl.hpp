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

#ifndef KOKKOSLAPACK_GESV_TPL_SPEC_DECL_HPP_
#define KOKKOSLAPACK_GESV_TPL_SPEC_DECL_HPP_

namespace KokkosLapack {
namespace Impl {
template <class AViewType, class BViewType, class PViewType>
inline void gesv_print_specialization() {
#ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
  printf("KokkosLapack::gesv<> TPL MAGMA specialization for < %s , %s, %s >\n",
         typeid(AViewType).name(), typeid(BViewType).name(),
         typeid(PViewType).name());
#else
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK
  printf("KokkosLapack::gesv<> TPL Lapack specialization for < %s , %s, %s >\n",
         typeid(AViewType).name(), typeid(BViewType).name(),
         typeid(PViewType).name());
#endif
#endif
#endif
}
}  // namespace Impl
}  // namespace KokkosLapack

// Generic Host side LAPACK (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK
#include <KokkosLapack_Host_tpl.hpp>

namespace KokkosLapack {
namespace Impl {

#define KOKKOSLAPACK_DGESV_LAPACK(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)          \
  template <class ExecSpace>                                                  \
  struct GESV<                                                                \
      Kokkos::View<double**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                 \
      Kokkos::View<double**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                 \
      Kokkos::View<int*, LAYOUT,                                              \
                   Kokkos::Device<Kokkos::DefaultHostExecutionSpace,          \
                                  Kokkos::HostSpace>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                 \
      true, ETI_SPEC_AVAIL> {                                                 \
    typedef double SCALAR;                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        AViewType;                                                            \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        BViewType;                                                            \
    typedef Kokkos::View<                                                     \
        int*, LAYOUT,                                                         \
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >                             \
        PViewType;                                                            \
                                                                              \
    static void gesv(const AViewType& A, const BViewType& B,                  \
                     const PViewType& IPIV) {                                 \
      Kokkos::Profiling::pushRegion("KokkosLapack::gesv[TPL_LAPACK,double]"); \
      gesv_print_specialization<AViewType, BViewType, PViewType>();           \
      const bool with_pivot =                                                 \
          !((IPIV.extent(0) == 0) && (IPIV.data() == nullptr));               \
                                                                              \
      const int N    = static_cast<int>(A.extent(1));                         \
      const int AST  = static_cast<int>(A.stride(1));                         \
      const int LDA  = (AST == 0) ? 1 : AST;                                  \
      const int BST  = static_cast<int>(B.stride(1));                         \
      const int LDB  = (BST == 0) ? 1 : BST;                                  \
      const int NRHS = static_cast<int>(B.extent(1));                         \
                                                                              \
      int info = 0;                                                           \
                                                                              \
      if (with_pivot) {                                                       \
        HostLapack<double>::gesv(N, NRHS, A.data(), LDA, IPIV.data(),         \
                                 B.data(), LDB, info);                        \
      }                                                                       \
      Kokkos::Profiling::popRegion();                                         \
    }                                                                         \
  };

#define KOKKOSLAPACK_SGESV_LAPACK(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)           \
  template <class ExecSpace>                                                   \
  struct GESV<                                                                 \
      Kokkos::View<float**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<float**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,      \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      Kokkos::View<int*, LAYOUT,                                               \
                   Kokkos::Device<Kokkos::DefaultHostExecutionSpace,           \
                                  Kokkos::HostSpace>,                          \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                  \
      true, ETI_SPEC_AVAIL> {                                                  \
    typedef float SCALAR;                                                      \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                     \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        AViewType;                                                             \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                     \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                 \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >             \
        BViewType;                                                             \
    typedef Kokkos::View<                                                      \
        int*, LAYOUT,                                                          \
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>,  \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >                              \
        PViewType;                                                             \
                                                                               \
    static void gesv(const AViewType& A, const BViewType& B,                   \
                     const PViewType& IPIV) {                                  \
      Kokkos::Profiling::pushRegion("KokkosLapack::gesv[TPL_LAPACK,float]");   \
      gesv_print_specialization<AViewType, BViewType, PViewType>();            \
      const bool with_pivot =                                                  \
          !((IPIV.extent(0) == 0) && (IPIV.data() == nullptr));                \
                                                                               \
      const int N    = static_cast<int>(A.extent(1));                          \
      const int AST  = static_cast<int>(A.stride(1));                          \
      const int LDA  = (AST == 0) ? 1 : AST;                                   \
      const int BST  = static_cast<int>(B.stride(1));                          \
      const int LDB  = (BST == 0) ? 1 : BST;                                   \
      const int NRHS = static_cast<int>(B.extent(1));                          \
                                                                               \
      int info = 0;                                                            \
                                                                               \
      if (with_pivot) {                                                        \
        HostLapack<float>::gesv(N, NRHS, A.data(), LDA, IPIV.data(), B.data(), \
                                LDB, info);                                    \
      }                                                                        \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#define KOKKOSLAPACK_ZGESV_LAPACK(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)          \
  template <class ExecSpace>                                                  \
  struct GESV<Kokkos::View<Kokkos::complex<double>**, LAYOUT,                 \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              Kokkos::View<Kokkos::complex<double>**, LAYOUT,                 \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              Kokkos::View<int*, LAYOUT,                                      \
                           Kokkos::Device<Kokkos::DefaultHostExecutionSpace,  \
                                          Kokkos::HostSpace>,                 \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              true, ETI_SPEC_AVAIL> {                                         \
    typedef Kokkos::complex<double> SCALAR;                                   \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        AViewType;                                                            \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        BViewType;                                                            \
    typedef Kokkos::View<                                                     \
        int*, LAYOUT,                                                         \
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >                             \
        PViewType;                                                            \
                                                                              \
    static void gesv(const AViewType& A, const BViewType& B,                  \
                     const PViewType& IPIV) {                                 \
      Kokkos::Profiling::pushRegion(                                          \
          "KokkosLapack::gesv[TPL_LAPACK,complex<double>]");                  \
      gesv_print_specialization<AViewType, BViewType, PViewType>();           \
      const bool with_pivot =                                                 \
          !((IPIV.extent(0) == 0) && (IPIV.data() == nullptr));               \
                                                                              \
      const int N    = static_cast<int>(A.extent(1));                         \
      const int AST  = static_cast<int>(A.stride(1));                         \
      const int LDA  = (AST == 0) ? 1 : AST;                                  \
      const int BST  = static_cast<int>(B.stride(1));                         \
      const int LDB  = (BST == 0) ? 1 : BST;                                  \
      const int NRHS = static_cast<int>(B.extent(1));                         \
                                                                              \
      int info = 0;                                                           \
                                                                              \
      if (with_pivot) {                                                       \
        HostLapack<std::complex<double> >::gesv(                              \
            N, NRHS, reinterpret_cast<std::complex<double>*>(A.data()), LDA,  \
            IPIV.data(), reinterpret_cast<std::complex<double>*>(B.data()),   \
            LDB, info);                                                       \
      }                                                                       \
      Kokkos::Profiling::popRegion();                                         \
    }                                                                         \
  };

#define KOKKOSLAPACK_CGESV_LAPACK(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)          \
  template <class ExecSpace>                                                  \
  struct GESV<Kokkos::View<Kokkos::complex<float>**, LAYOUT,                  \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              Kokkos::View<Kokkos::complex<float>**, LAYOUT,                  \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              Kokkos::View<int*, LAYOUT,                                      \
                           Kokkos::Device<Kokkos::DefaultHostExecutionSpace,  \
                                          Kokkos::HostSpace>,                 \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              true, ETI_SPEC_AVAIL> {                                         \
    typedef Kokkos::complex<float> SCALAR;                                    \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        AViewType;                                                            \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        BViewType;                                                            \
    typedef Kokkos::View<                                                     \
        int*, LAYOUT,                                                         \
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >                             \
        PViewType;                                                            \
                                                                              \
    static void gesv(const AViewType& A, const BViewType& B,                  \
                     const PViewType& IPIV) {                                 \
      Kokkos::Profiling::pushRegion(                                          \
          "KokkosLapack::gesv[TPL_LAPACK,complex<float>]");                   \
      gesv_print_specialization<AViewType, BViewType, PViewType>();           \
      const bool with_pivot =                                                 \
          !((IPIV.extent(0) == 0) && (IPIV.data() == nullptr));               \
                                                                              \
      const int N    = static_cast<int>(A.extent(1));                         \
      const int AST  = static_cast<int>(A.stride(1));                         \
      const int LDA  = (AST == 0) ? 1 : AST;                                  \
      const int BST  = static_cast<int>(B.stride(1));                         \
      const int LDB  = (BST == 0) ? 1 : BST;                                  \
      const int NRHS = static_cast<int>(B.extent(1));                         \
                                                                              \
      int info = 0;                                                           \
                                                                              \
      if (with_pivot) {                                                       \
        HostLapack<std::complex<float> >::gesv(                               \
            N, NRHS, reinterpret_cast<std::complex<float>*>(A.data()), LDA,   \
            IPIV.data(), reinterpret_cast<std::complex<float>*>(B.data()),    \
            LDB, info);                                                       \
      }                                                                       \
      Kokkos::Profiling::popRegion();                                         \
    }                                                                         \
  };

KOKKOSLAPACK_DGESV_LAPACK(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSLAPACK_DGESV_LAPACK(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSLAPACK_SGESV_LAPACK(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSLAPACK_SGESV_LAPACK(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSLAPACK_ZGESV_LAPACK(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSLAPACK_ZGESV_LAPACK(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSLAPACK_CGESV_LAPACK(Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSLAPACK_CGESV_LAPACK(Kokkos::LayoutLeft, Kokkos::HostSpace, false)

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_LAPACK

// MAGMA
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#include <KokkosLapack_magma.hpp>

namespace KokkosLapack {
namespace Impl {

#define KOKKOSLAPACK_DGESV_MAGMA(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)           \
  template <class ExecSpace>                                                  \
  struct GESV<                                                                \
      Kokkos::View<double**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                 \
      Kokkos::View<double**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                 \
      Kokkos::View<magma_int_t*, LAYOUT,                                      \
                   Kokkos::Device<Kokkos::DefaultHostExecutionSpace,          \
                                  Kokkos::HostSpace>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                 \
      true, ETI_SPEC_AVAIL> {                                                 \
    typedef double SCALAR;                                                    \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        AViewType;                                                            \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        BViewType;                                                            \
    typedef Kokkos::View<                                                     \
        magma_int_t*, LAYOUT,                                                 \
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >                             \
        PViewType;                                                            \
                                                                              \
    static void gesv(const AViewType& A, const BViewType& B,                  \
                     const PViewType& IPIV) {                                 \
      Kokkos::Profiling::pushRegion("KokkosLapack::gesv[TPL_MAGMA,double]");  \
      gesv_print_specialization<AViewType, BViewType, PViewType>();           \
      const bool with_pivot =                                                 \
          !((IPIV.extent(0) == 0) && (IPIV.data() == nullptr));               \
                                                                              \
      magma_int_t N    = static_cast<magma_int_t>(A.extent(1));               \
      magma_int_t AST  = static_cast<magma_int_t>(A.stride(1));               \
      magma_int_t LDA  = (AST == 0) ? 1 : AST;                                \
      magma_int_t BST  = static_cast<magma_int_t>(B.stride(1));               \
      magma_int_t LDB  = (BST == 0) ? 1 : BST;                                \
      magma_int_t NRHS = static_cast<magma_int_t>(B.extent(1));               \
                                                                              \
      KokkosLapack::Impl::MagmaSingleton& s =                                 \
          KokkosLapack::Impl::MagmaSingleton::singleton();                    \
      magma_int_t info = 0;                                                   \
                                                                              \
      if (with_pivot) {                                                       \
        magma_dgesv_gpu(N, NRHS, reinterpret_cast<magmaDouble_ptr>(A.data()), \
                        LDA, IPIV.data(),                                     \
                        reinterpret_cast<magmaDouble_ptr>(B.data()), LDB,     \
                        &info);                                               \
      } else {                                                                \
        magma_dgesv_nopiv_gpu(                                                \
            N, NRHS, reinterpret_cast<magmaDouble_ptr>(A.data()), LDA,        \
            reinterpret_cast<magmaDouble_ptr>(B.data()), LDB, &info);         \
      }                                                                       \
      Kokkos::Profiling::popRegion();                                         \
    }                                                                         \
  };

#define KOKKOSLAPACK_SGESV_MAGMA(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)           \
  template <class ExecSpace>                                                  \
  struct GESV<                                                                \
      Kokkos::View<float**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                 \
      Kokkos::View<float**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>,     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                 \
      Kokkos::View<magma_int_t*, LAYOUT,                                      \
                   Kokkos::Device<Kokkos::DefaultHostExecutionSpace,          \
                                  Kokkos::HostSpace>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                 \
      true, ETI_SPEC_AVAIL> {                                                 \
    typedef float SCALAR;                                                     \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        AViewType;                                                            \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        BViewType;                                                            \
    typedef Kokkos::View<                                                     \
        magma_int_t*, LAYOUT,                                                 \
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >                             \
        PViewType;                                                            \
                                                                              \
    static void gesv(const AViewType& A, const BViewType& B,                  \
                     const PViewType& IPIV) {                                 \
      Kokkos::Profiling::pushRegion("KokkosLapack::gesv[TPL_MAGMA,float]");   \
      gesv_print_specialization<AViewType, BViewType, PViewType>();           \
      const bool with_pivot =                                                 \
          !((IPIV.extent(0) == 0) && (IPIV.data() == nullptr));               \
                                                                              \
      magma_int_t N    = static_cast<magma_int_t>(A.extent(1));               \
      magma_int_t AST  = static_cast<magma_int_t>(A.stride(1));               \
      magma_int_t LDA  = (AST == 0) ? 1 : AST;                                \
      magma_int_t BST  = static_cast<magma_int_t>(B.stride(1));               \
      magma_int_t LDB  = (BST == 0) ? 1 : BST;                                \
      magma_int_t NRHS = static_cast<magma_int_t>(B.extent(1));               \
                                                                              \
      KokkosLapack::Impl::MagmaSingleton& s =                                 \
          KokkosLapack::Impl::MagmaSingleton::singleton();                    \
      magma_int_t info = 0;                                                   \
                                                                              \
      if (with_pivot) {                                                       \
        magma_sgesv_gpu(N, NRHS, reinterpret_cast<magmaFloat_ptr>(A.data()),  \
                        LDA, IPIV.data(),                                     \
                        reinterpret_cast<magmaFloat_ptr>(B.data()), LDB,      \
                        &info);                                               \
      } else {                                                                \
        magma_sgesv_nopiv_gpu(                                                \
            N, NRHS, reinterpret_cast<magmaFloat_ptr>(A.data()), LDA,         \
            reinterpret_cast<magmaFloat_ptr>(B.data()), LDB, &info);          \
      }                                                                       \
      Kokkos::Profiling::popRegion();                                         \
    }                                                                         \
  };

#define KOKKOSLAPACK_ZGESV_MAGMA(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)           \
  template <class ExecSpace>                                                  \
  struct GESV<Kokkos::View<Kokkos::complex<double>**, LAYOUT,                 \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              Kokkos::View<Kokkos::complex<double>**, LAYOUT,                 \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              Kokkos::View<magma_int_t*, LAYOUT,                              \
                           Kokkos::Device<Kokkos::DefaultHostExecutionSpace,  \
                                          Kokkos::HostSpace>,                 \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              true, ETI_SPEC_AVAIL> {                                         \
    typedef Kokkos::complex<double> SCALAR;                                   \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        AViewType;                                                            \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        BViewType;                                                            \
    typedef Kokkos::View<                                                     \
        magma_int_t*, LAYOUT,                                                 \
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >                             \
        PViewType;                                                            \
                                                                              \
    static void gesv(const AViewType& A, const BViewType& B,                  \
                     const PViewType& IPIV) {                                 \
      Kokkos::Profiling::pushRegion(                                          \
          "KokkosLapack::gesv[TPL_MAGMA,complex<double>]");                   \
      gesv_print_specialization<AViewType, BViewType, PViewType>();           \
      const bool with_pivot =                                                 \
          !((IPIV.extent(0) == 0) && (IPIV.data() == nullptr));               \
                                                                              \
      magma_int_t N    = static_cast<magma_int_t>(A.extent(1));               \
      magma_int_t AST  = static_cast<magma_int_t>(A.stride(1));               \
      magma_int_t LDA  = (AST == 0) ? 1 : AST;                                \
      magma_int_t BST  = static_cast<magma_int_t>(B.stride(1));               \
      magma_int_t LDB  = (BST == 0) ? 1 : BST;                                \
      magma_int_t NRHS = static_cast<magma_int_t>(B.extent(1));               \
                                                                              \
      KokkosLapack::Impl::MagmaSingleton& s =                                 \
          KokkosLapack::Impl::MagmaSingleton::singleton();                    \
      magma_int_t info = 0;                                                   \
                                                                              \
      if (with_pivot) {                                                       \
        magma_zgesv_gpu(                                                      \
            N, NRHS, reinterpret_cast<magmaDoubleComplex_ptr>(A.data()), LDA, \
            IPIV.data(), reinterpret_cast<magmaDoubleComplex_ptr>(B.data()),  \
            LDB, &info);                                                      \
      } else {                                                                \
        magma_zgesv_nopiv_gpu(                                                \
            N, NRHS, reinterpret_cast<magmaDoubleComplex_ptr>(A.data()), LDA, \
            reinterpret_cast<magmaDoubleComplex_ptr>(B.data()), LDB, &info);  \
      }                                                                       \
      Kokkos::Profiling::popRegion();                                         \
    }                                                                         \
  };

#define KOKKOSLAPACK_CGESV_MAGMA(LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL)           \
  template <class ExecSpace>                                                  \
  struct GESV<Kokkos::View<Kokkos::complex<float>**, LAYOUT,                  \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              Kokkos::View<Kokkos::complex<float>**, LAYOUT,                  \
                           Kokkos::Device<ExecSpace, MEM_SPACE>,              \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              Kokkos::View<magma_int_t*, LAYOUT,                              \
                           Kokkos::Device<Kokkos::DefaultHostExecutionSpace,  \
                                          Kokkos::HostSpace>,                 \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged> >,         \
              true, ETI_SPEC_AVAIL> {                                         \
    typedef Kokkos::complex<float> SCALAR;                                    \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        AViewType;                                                            \
    typedef Kokkos::View<SCALAR**, LAYOUT,                                    \
                         Kokkos::Device<ExecSpace, MEM_SPACE>,                \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >            \
        BViewType;                                                            \
    typedef Kokkos::View<                                                     \
        magma_int_t*, LAYOUT,                                                 \
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >                             \
        PViewType;                                                            \
                                                                              \
    static void gesv(const AViewType& A, const BViewType& B,                  \
                     const PViewType& IPIV) {                                 \
      Kokkos::Profiling::pushRegion(                                          \
          "KokkosLapack::gesv[TPL_MAGMA,complex<float>]");                    \
      gesv_print_specialization<AViewType, BViewType, PViewType>();           \
      const bool with_pivot =                                                 \
          !((IPIV.extent(0) == 0) && (IPIV.data() == nullptr));               \
                                                                              \
      magma_int_t N    = static_cast<magma_int_t>(A.extent(1));               \
      magma_int_t AST  = static_cast<magma_int_t>(A.stride(1));               \
      magma_int_t LDA  = (AST == 0) ? 1 : AST;                                \
      magma_int_t BST  = static_cast<magma_int_t>(B.stride(1));               \
      magma_int_t LDB  = (BST == 0) ? 1 : BST;                                \
      magma_int_t NRHS = static_cast<magma_int_t>(B.extent(1));               \
                                                                              \
      KokkosLapack::Impl::MagmaSingleton& s =                                 \
          KokkosLapack::Impl::MagmaSingleton::singleton();                    \
      magma_int_t info = 0;                                                   \
                                                                              \
      if (with_pivot) {                                                       \
        magma_cgesv_gpu(                                                      \
            N, NRHS, reinterpret_cast<magmaFloatComplex_ptr>(A.data()), LDA,  \
            IPIV.data(), reinterpret_cast<magmaFloatComplex_ptr>(B.data()),   \
            LDB, &info);                                                      \
      } else {                                                                \
        magma_cgesv_nopiv_gpu(                                                \
            N, NRHS, reinterpret_cast<magmaFloatComplex_ptr>(A.data()), LDA,  \
            reinterpret_cast<magmaFloatComplex_ptr>(B.data()), LDB, &info);   \
      }                                                                       \
      Kokkos::Profiling::popRegion();                                         \
    }                                                                         \
  };

KOKKOSLAPACK_DGESV_MAGMA(Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSLAPACK_DGESV_MAGMA(Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSLAPACK_SGESV_MAGMA(Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSLAPACK_SGESV_MAGMA(Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSLAPACK_ZGESV_MAGMA(Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSLAPACK_ZGESV_MAGMA(Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSLAPACK_CGESV_MAGMA(Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSLAPACK_CGESV_MAGMA(Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_MAGMA

#endif
