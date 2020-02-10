/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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

#ifndef KOKKOSBLAS_GESV_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS_GESV_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {
  template<class AViewType, class BViewType, class PViewType>
  inline void gesv_print_specialization() {
      #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
        #ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
          printf("KokkosBlas::gesv<> TPL MAGMA specialization for < %s , %s, %s >\n",typeid(AViewType).name(),typeid(BViewType).name(),typeid(PViewType).name());
        #else
          #ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
            printf("KokkosBlas::gesv<> TPL Blas specialization for < %s , %s, %s >\n",typeid(AViewType).name(),typeid(BViewType).name(),typeid(PViewType).name());
          #endif        
        #endif
      #endif
  }
}
}

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#include<KokkosBlas_Host_tpl.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS_DGESV_BLAS( LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GESV< \
     Kokkos::View<double**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<double**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef double SCALAR; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > PViewType; \
  \
  static void \
  gesv (const AViewType& A, \
        const BViewType& B, \
        const PViewType& IPIV) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gesv[TPL_BLAS,double]"); \
    gesv_print_specialization<AViewType,BViewType,PViewType>(); \
    const bool with_pivot = !((IPIV.extent(0) == 0) && (IPIV.data()==nullptr)); \
    \
    const int N    = static_cast<int> (A.extent(1)); \
    const int AST  = static_cast<int> (A.stride(1)); \
    const int LDA  = (AST == 0) ? 1 : AST; \
    const int BST  = static_cast<int> (B.stride(1)); \
    const int LDB  = (BST == 0) ? 1 : BST; \
    const int NRHS = static_cast<int> (B.extent(1)); \
    \
    int info = 0; \
    \
    if(with_pivot) { \
      HostBlas<double>::gesv (N, NRHS, A.data(), LDA, IPIV.data(), B.data(), LDB, info); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS_SGESV_BLAS( LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GESV< \
     Kokkos::View<float**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<float**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef float SCALAR; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > PViewType; \
  \
  static void \
  gesv (const AViewType& A, \
        const BViewType& B, \
        const PViewType& IPIV) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gesv[TPL_BLAS,float]"); \
    gesv_print_specialization<AViewType,BViewType,PViewType>(); \
    const bool with_pivot = !((IPIV.extent(0) == 0) && (IPIV.data()==nullptr)); \
    \
    const int N    = static_cast<int> (A.extent(1)); \
    const int AST  = static_cast<int> (A.stride(1)); \
    const int LDA  = (AST == 0) ? 1 : AST; \
    const int BST  = static_cast<int> (B.stride(1)); \
    const int LDB  = (BST == 0) ? 1 : BST; \
    const int NRHS = static_cast<int> (B.extent(1)); \
    \
    int info = 0; \
    \
    if(with_pivot) { \
      HostBlas<float>::gesv (N, NRHS, A.data(), LDA, IPIV.data(), B.data(), LDB, info); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS_ZGESV_BLAS( LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GESV< \
     Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef Kokkos::complex<double> SCALAR; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > PViewType; \
  \
  static void \
  gesv (const AViewType& A, \
        const BViewType& B, \
        const PViewType& IPIV) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gesv[TPL_BLAS,complex<double>]"); \
    gesv_print_specialization<AViewType,BViewType,PViewType>(); \
    const bool with_pivot = !((IPIV.extent(0) == 0) && (IPIV.data()==nullptr)); \
    \
    const int N    = static_cast<int> (A.extent(1)); \
    const int AST  = static_cast<int> (A.stride(1)); \
    const int LDA  = (AST == 0) ? 1 : AST; \
    const int BST  = static_cast<int> (B.stride(1)); \
    const int LDB  = (BST == 0) ? 1 : BST; \
    const int NRHS = static_cast<int> (B.extent(1)); \
    \
    int info = 0; \
    \
    if(with_pivot) { \
      HostBlas<std::complex<double> >::gesv \
        (N, NRHS, reinterpret_cast<std::complex<double>*>(A.data()), LDA, IPIV.data(), reinterpret_cast<std::complex<double>*>(B.data()), LDB, info); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
}; \

#define KOKKOSBLAS_CGESV_BLAS( LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GESV< \
     Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef Kokkos::complex<float> SCALAR; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<int*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > PViewType; \
  \
  static void \
  gesv (const AViewType& A, \
        const BViewType& B, \
        const PViewType& IPIV) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gesv[TPL_BLAS,complex<float>]"); \
    gesv_print_specialization<AViewType,BViewType,PViewType>(); \
    const bool with_pivot = !((IPIV.extent(0) == 0) && (IPIV.data()==nullptr)); \
    \
    const int N    = static_cast<int> (A.extent(1)); \
    const int AST  = static_cast<int> (A.stride(1)); \
    const int LDA  = (AST == 0) ? 1 : AST; \
    const int BST  = static_cast<int> (B.stride(1)); \
    const int LDB  = (BST == 0) ? 1 : BST; \
    const int NRHS = static_cast<int> (B.extent(1)); \
    \
    int info = 0; \
    \
    if(with_pivot) { \
      HostBlas<std::complex<float> >::gesv \
        (N, NRHS, reinterpret_cast<std::complex<float>*>(A.data()), LDA, IPIV.data(), reinterpret_cast<std::complex<float>*>(B.data()), LDB, info); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

KOKKOSBLAS_DGESV_BLAS( Kokkos::LayoutLeft,  Kokkos::HostSpace, true)
KOKKOSBLAS_DGESV_BLAS( Kokkos::LayoutLeft,  Kokkos::HostSpace, false)

KOKKOSBLAS_SGESV_BLAS( Kokkos::LayoutLeft,  Kokkos::HostSpace, true)
KOKKOSBLAS_SGESV_BLAS( Kokkos::LayoutLeft,  Kokkos::HostSpace, false)

KOKKOSBLAS_ZGESV_BLAS( Kokkos::LayoutLeft,  Kokkos::HostSpace, true)
KOKKOSBLAS_ZGESV_BLAS( Kokkos::LayoutLeft,  Kokkos::HostSpace, false)

KOKKOSBLAS_CGESV_BLAS( Kokkos::LayoutLeft,  Kokkos::HostSpace, true)
KOKKOSBLAS_CGESV_BLAS( Kokkos::LayoutLeft,  Kokkos::HostSpace, false)

}
}
#endif // KOKKOSKERNELS_ENABLE_TPL_BLAS

// MAGMA
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA
#include<KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS_DGESV_MAGMA( LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GESV< \
     Kokkos::View<double**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<double**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<magma_int_t*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef double SCALAR; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<magma_int_t*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > PViewType; \
  \
  static void \
  gesv (const AViewType& A, \
        const BViewType& B, \
        const PViewType& IPIV) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gesv[TPL_MAGMA,double]"); \
    gesv_print_specialization<AViewType,BViewType,PViewType>(); \
    const bool with_pivot = !((IPIV.extent(0) == 0) && (IPIV.data()==nullptr)); \
    \
    magma_int_t N    = static_cast<magma_int_t> (A.extent(1)); \
    magma_int_t AST  = static_cast<magma_int_t> (A.stride(1)); \
    magma_int_t LDA  = (AST == 0) ? 1 : AST; \
    magma_int_t BST  = static_cast<magma_int_t> (B.stride(1)); \
    magma_int_t LDB  = (BST == 0) ? 1 : BST; \
    magma_int_t NRHS = static_cast<magma_int_t> (B.extent(1)); \
    \
    KokkosBlas::Impl::MagmaSingleton & s = KokkosBlas::Impl::MagmaSingleton::singleton(); \
    magma_int_t  info = 0; \
    \
    if(with_pivot) { \
      magma_dgesv_gpu(N,NRHS,reinterpret_cast<magmaDouble_ptr>(A.data()),LDA,IPIV.data(),reinterpret_cast<magmaDouble_ptr>(B.data()),LDB,&info); \
    } \
    else { \
      magma_dgesv_nopiv_gpu(N,NRHS,reinterpret_cast<magmaDouble_ptr>(A.data()),LDA,reinterpret_cast<magmaDouble_ptr>(B.data()),LDB,&info); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS_SGESV_MAGMA( LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GESV< \
     Kokkos::View<float**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<float**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<magma_int_t*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef float SCALAR; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<magma_int_t*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > PViewType; \
  \
  static void \
  gesv (const AViewType& A, \
        const BViewType& B, \
        const PViewType& IPIV) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gesv[TPL_MAGMA,float]"); \
    gesv_print_specialization<AViewType,BViewType,PViewType>(); \
    const bool with_pivot = !((IPIV.extent(0) == 0) && (IPIV.data()==nullptr)); \
    \
    magma_int_t N    = static_cast<magma_int_t> (A.extent(1)); \
    magma_int_t AST  = static_cast<magma_int_t> (A.stride(1)); \
    magma_int_t LDA  = (AST == 0) ? 1 : AST; \
    magma_int_t BST  = static_cast<magma_int_t> (B.stride(1)); \
    magma_int_t LDB  = (BST == 0) ? 1 : BST; \
    magma_int_t NRHS = static_cast<magma_int_t> (B.extent(1)); \
    \
    KokkosBlas::Impl::MagmaSingleton & s = KokkosBlas::Impl::MagmaSingleton::singleton(); \
    magma_int_t  info = 0; \
    \
    if(with_pivot) { \
      magma_sgesv_gpu(N,NRHS,reinterpret_cast<magmaFloat_ptr>(A.data()),LDA,IPIV.data(),reinterpret_cast<magmaFloat_ptr>(B.data()),LDB,&info); \
    } \
    else { \
      magma_sgesv_nopiv_gpu(N,NRHS,reinterpret_cast<magmaFloat_ptr>(A.data()),LDA,reinterpret_cast<magmaFloat_ptr>(B.data()),LDB,&info); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS_ZGESV_MAGMA( LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GESV< \
     Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<Kokkos::complex<double>**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<magma_int_t*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef Kokkos::complex<double> SCALAR; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<magma_int_t*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > PViewType; \
  \
  static void \
  gesv (const AViewType& A, \
        const BViewType& B, \
        const PViewType& IPIV) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gesv[TPL_MAGMA,complex<double>]"); \
    gesv_print_specialization<AViewType,BViewType,PViewType>(); \
    const bool with_pivot = !((IPIV.extent(0) == 0) && (IPIV.data()==nullptr)); \
    \
    magma_int_t N    = static_cast<magma_int_t> (A.extent(1)); \
    magma_int_t AST  = static_cast<magma_int_t> (A.stride(1)); \
    magma_int_t LDA  = (AST == 0) ? 1 : AST; \
    magma_int_t BST  = static_cast<magma_int_t> (B.stride(1)); \
    magma_int_t LDB  = (BST == 0) ? 1 : BST; \
    magma_int_t NRHS = static_cast<magma_int_t> (B.extent(1)); \
    \
    KokkosBlas::Impl::MagmaSingleton & s = KokkosBlas::Impl::MagmaSingleton::singleton(); \
    magma_int_t  info = 0; \
    \
    if(with_pivot) { \
      magma_zgesv_gpu(N,NRHS,reinterpret_cast<magmaDoubleComplex_ptr>(A.data()),LDA,IPIV.data(),reinterpret_cast<magmaDoubleComplex_ptr>(B.data()),LDB,&info); \
    } \
    else { \
      magma_zgesv_nopiv_gpu(N,NRHS,reinterpret_cast<magmaDoubleComplex_ptr>(A.data()),LDA,reinterpret_cast<magmaDoubleComplex_ptr>(B.data()),LDB,&info); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
}; \

#define KOKKOSBLAS_CGESV_MAGMA( LAYOUT, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GESV< \
     Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<Kokkos::complex<float>**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<magma_int_t*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef Kokkos::complex<float> SCALAR; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<magma_int_t*, LAYOUT, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > PViewType; \
 \
  static void \
  gesv (const AViewType& A, \
        const BViewType& B, \
        const PViewType& IPIV) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gesv[TPL_MAGMA,complex<float>]"); \
    gesv_print_specialization<AViewType,BViewType,PViewType>(); \
    const bool with_pivot = !((IPIV.extent(0) == 0) && (IPIV.data()==nullptr)); \
    \
    magma_int_t N    = static_cast<magma_int_t> (A.extent(1)); \
    magma_int_t AST  = static_cast<magma_int_t> (A.stride(1)); \
    magma_int_t LDA  = (AST == 0) ? 1 : AST; \
    magma_int_t BST  = static_cast<magma_int_t> (B.stride(1)); \
    magma_int_t LDB  = (BST == 0) ? 1 : BST; \
    magma_int_t NRHS = static_cast<magma_int_t> (B.extent(1)); \
    \
    KokkosBlas::Impl::MagmaSingleton & s = KokkosBlas::Impl::MagmaSingleton::singleton(); \
    magma_int_t  info = 0; \
    \
    if(with_pivot) { \
      magma_cgesv_gpu(N,NRHS,reinterpret_cast<magmaFloatComplex_ptr>(A.data()),LDA,IPIV.data(),reinterpret_cast<magmaFloatComplex_ptr>(B.data()),LDB,&info); \
    } \
    else { \
      magma_cgesv_nopiv_gpu(N,NRHS,reinterpret_cast<magmaFloatComplex_ptr>(A.data()),LDA,reinterpret_cast<magmaFloatComplex_ptr>(B.data()),LDB,&info); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

KOKKOSBLAS_DGESV_MAGMA( Kokkos::LayoutLeft,  Kokkos::CudaSpace, true)
KOKKOSBLAS_DGESV_MAGMA( Kokkos::LayoutLeft,  Kokkos::CudaSpace, false)

KOKKOSBLAS_SGESV_MAGMA( Kokkos::LayoutLeft,  Kokkos::CudaSpace, true)
KOKKOSBLAS_SGESV_MAGMA( Kokkos::LayoutLeft,  Kokkos::CudaSpace, false)

KOKKOSBLAS_ZGESV_MAGMA( Kokkos::LayoutLeft,  Kokkos::CudaSpace, true)
KOKKOSBLAS_ZGESV_MAGMA( Kokkos::LayoutLeft,  Kokkos::CudaSpace, false)

KOKKOSBLAS_CGESV_MAGMA( Kokkos::LayoutLeft,  Kokkos::CudaSpace, true)
KOKKOSBLAS_CGESV_MAGMA( Kokkos::LayoutLeft,  Kokkos::CudaSpace, false)

}
}
#endif // KOKKOSKERNELS_ENABLE_TPL_MAGMA

#endif
