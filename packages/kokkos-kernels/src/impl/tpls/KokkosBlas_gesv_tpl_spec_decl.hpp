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
     true, ETI_SPEC_AVAIL> { \
  typedef double SCALAR; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
 \
  static void \
  gesv (const char pivot[], \
        const AViewType& A, \
        const BViewType& B) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gesv[TPL_MAGMA,double]"); \
    const bool nopivot_t = (pivot[0]=='N') || (pivot[0]=='n'); \
    const bool pivot_t   = (pivot[0]=='Y') || (pivot[0]=='y'); \
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
    if(pivot_t) { \
      magma_int_t *ipiv = NULL; \
      magma_imalloc_cpu( &ipiv, N ); \
      magma_dgesv_gpu ( N, NRHS, reinterpret_cast<magmaDouble_ptr>(A.data()), LDA, ipiv, reinterpret_cast<magmaDouble_ptr>(B.data()), LDB, &info ); \
      magma_free_cpu( ipiv ); \
    } \
    if(nopivot_t) \
      magma_dgesv_nopiv_gpu( N, NRHS, reinterpret_cast<magmaDouble_ptr>(A.data()), LDA, reinterpret_cast<magmaDouble_ptr>(B.data()), LDB, &info ); \
    \
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
     true, ETI_SPEC_AVAIL> { \
  typedef float SCALAR; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
      \
  static void \
  gesv (const char pivot[], \
        const AViewType& A, \
        const BViewType& B) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gesv[TPL_MAGMA,float]"); \
    const bool nopivot_t = (pivot[0]=='N') || (pivot[0]=='n'); \
    const bool pivot_t   = (pivot[0]=='Y') || (pivot[0]=='y'); \
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
    if(pivot_t) { \
      magma_int_t *ipiv = NULL; \
      magma_imalloc_cpu( &ipiv, N ); \
      magma_sgesv_gpu ( N, NRHS, reinterpret_cast<magmaFloat_ptr>(A.data()), LDA, ipiv, reinterpret_cast<magmaFloat_ptr>(B.data()), LDB, &info ); \
      magma_free_cpu( ipiv ); \
    } \
    if(nopivot_t) \
      magma_sgesv_nopiv_gpu( N, NRHS, reinterpret_cast<magmaFloat_ptr>(A.data()), LDA, reinterpret_cast<magmaFloat_ptr>(B.data()), LDB, &info ); \
    \
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
     true, ETI_SPEC_AVAIL> { \
  typedef Kokkos::complex<double> SCALAR; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
      \
  static void \
  gesv (const char pivot[], \
        const AViewType& A, \
        const BViewType& B) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gesv[TPL_MAGMA,complex<double>]"); \
    const bool nopivot_t = (pivot[0]=='N') || (pivot[0]=='n'); \
    const bool pivot_t   = (pivot[0]=='Y') || (pivot[0]=='y'); \
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
    if(pivot_t) { \
      magma_int_t *ipiv = NULL; \
      magma_imalloc_cpu( &ipiv, N ); \
      magma_zgesv_gpu ( N, NRHS, reinterpret_cast<magmaDoubleComplex_ptr>(A.data()), LDA, ipiv, reinterpret_cast<magmaDoubleComplex_ptr>(B.data()), LDB, &info ); \
      magma_free_cpu( ipiv ); \
    } \
    if(nopivot_t) \
      magma_zgesv_nopiv_gpu( N, NRHS, reinterpret_cast<magmaDoubleComplex_ptr>(A.data()), LDA, reinterpret_cast<magmaDoubleComplex_ptr>(B.data()), LDB, &info ); \
    \
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
     true, ETI_SPEC_AVAIL> { \
  typedef Kokkos::complex<float> SCALAR; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
      \
  static void \
  gesv (const char pivot[], \
        const AViewType& A, \
        const BViewType& B) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gesv[TPL_MAGMA,complex<float>]"); \
    const bool nopivot_t = (pivot[0]=='N') || (pivot[0]=='n'); \
    const bool pivot_t   = (pivot[0]=='Y') || (pivot[0]=='y'); \
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
    if(pivot_t) { \
      magma_int_t *ipiv = NULL; \
      magma_imalloc_cpu( &ipiv, N ); \
      magma_cgesv_gpu ( N, NRHS, reinterpret_cast<magmaFloatComplex_ptr>(A.data()), LDA, ipiv, reinterpret_cast<magmaFloatComplex_ptr>(B.data()), LDB, &info ); \
      magma_free_cpu( ipiv ); \
    } \
    if(nopivot_t) \
      magma_cgesv_nopiv_gpu( N, NRHS, reinterpret_cast<magmaFloatComplex_ptr>(A.data()), LDA, reinterpret_cast<magmaFloatComplex_ptr>(B.data()), LDB, &info ); \
    \
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
