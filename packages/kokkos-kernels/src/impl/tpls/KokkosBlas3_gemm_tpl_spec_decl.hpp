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

#ifndef KOKKOSBLAS3_GEMM_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS3_GEMM_TPL_SPEC_DECL_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
extern "C" void dgemm_( const char* transA, const char* transB,
                        const int* M, const int* N, const int* K,
                        const double* alpha,
                        const double* A, const int* LDA,
                        const double* B, const int* LDB,
                        const double* beta,
                        double* C, const int* LDC);
extern "C" void sgemm_( const char* transA, const char* transB,
                        const int* M, const int* N, const int* K,
                        const float* alpha,
                        const float* A, const int* LDA,
                        const float* B, const int* LDB,
                        const float* beta,
                        float* C, const int* LDC);
extern "C" void zgemm_( const char* transA, const char* transB,
                        const int* M, const int* N, const int* K,
                        const std::complex<double>* alpha,
                        const std::complex<double>* A, const int* LDA,
                        const std::complex<double>* B, const int* LDB,
                        const std::complex<double>* beta,
                        std::complex<double>* C, const int* LDC);
extern "C" void cgemm_( const char* transA, const char* transB,
                        const int* M, const int* N, const int* K,
                        const std::complex<float>* alpha,
                        const std::complex<float>* A, const int* LDA,
                        const std::complex<float>* B, const int* LDB,
                        const std::complex<float>* beta,
                        std::complex<float>* C, const int* LDC);

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS3_DGEMM_BLAS( LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GEMM< \
     Kokkos::View<const double**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const double**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<double**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef double SCALAR; \
  typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > CViewType; \
 \
  static void \
  gemm (const char transA[], \
        const char transB[], \
        typename AViewType::const_value_type& alpha, \
        const AViewType& A, \
        const BViewType& B, \
        typename CViewType::const_value_type& beta, \
        const CViewType& C) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_BLAS,double]"); \
    const bool A_t = (transA[0]!='N') && (transA[0]!='n'); \
    const int M = C.extent(0); \
    const int N = C.extent(1); \
    const int K = A.extent(A_t?0:1); \
    \
    bool A_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTA>::value; \
    bool B_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTB>::value; \
    bool C_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTC>::value; \
    \
    const int AST = A_is_lr?A.stride(0):A.stride(1), LDA = AST == 0 ? 1 : AST; \
    const int BST = B_is_lr?B.stride(0):B.stride(1), LDB = BST == 0 ? 1 : BST; \
    const int CST = C_is_lr?C.stride(0):C.stride(1), LDC = CST == 0 ? 1 : CST; \
    \
    if(!A_is_lr && !B_is_lr && !C_is_lr ) \
      dgemm_(transA,transB,&M,&N,&K,&alpha,A.data(),&LDA,B.data(),&LDB,&beta,C.data(),&LDC); \
    if(A_is_lr && B_is_lr && C_is_lr ) \
      dgemm_(transB,transA,&N,&M,&K,&alpha,B.data(),&LDB,A.data(),&LDA,&beta,C.data(),&LDC); \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS3_SGEMM_BLAS( LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GEMM< \
     Kokkos::View<const float**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const float**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<float**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef float SCALAR; \
  typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > CViewType; \
      \
  static void \
  gemm (const char transA[], \
        const char transB[], \
        typename AViewType::const_value_type& alpha, \
        const AViewType& A, \
        const BViewType& B, \
        typename CViewType::const_value_type& beta, \
        const CViewType& C) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_BLAS,float]"); \
    const bool A_t = (transA[0]!='N') && (transA[0]!='n'); \
    const int M = C.extent(0); \
    const int N = C.extent(1); \
    const int K = A.extent(A_t?0:1); \
    \
    bool A_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTA>::value; \
    bool B_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTB>::value; \
    bool C_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTC>::value; \
    \
    const int AST = A_is_lr?A.stride(0):A.stride(1), LDA = AST == 0 ? 1 : AST; \
    const int BST = B_is_lr?B.stride(0):B.stride(1), LDB = BST == 0 ? 1 : BST; \
    const int CST = C_is_lr?C.stride(0):C.stride(1), LDC = CST == 0 ? 1 : CST; \
    \
    if(!A_is_lr && !B_is_lr && !C_is_lr ) \
      sgemm_(transA,transB,&M,&N,&K,&alpha,A.data(),&LDA,B.data(),&LDB,&beta,C.data(),&LDC); \
    if(A_is_lr && B_is_lr && C_is_lr ) \
      sgemm_(transB,transA,&N,&M,&K,&alpha,B.data(),&LDB,A.data(),&LDA,&beta,C.data(),&LDC); \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS3_ZGEMM_BLAS( LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GEMM< \
     Kokkos::View<const Kokkos::complex<double>**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const Kokkos::complex<double>**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<Kokkos::complex<double>**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef Kokkos::complex<double> SCALAR; \
  typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > CViewType; \
      \
  static void \
  gemm (const char transA[], \
        const char transB[], \
        typename AViewType::const_value_type& alpha, \
        const AViewType& A, \
        const BViewType& B, \
        typename CViewType::const_value_type& beta, \
        const CViewType& C) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_BLAS,complex<double>]"); \
    const bool A_t = (transA[0]!='N') && (transA[0]!='n'); \
    const int M = C.extent(0); \
    const int N = C.extent(1); \
    const int K = A.extent(A_t?0:1); \
    \
    bool A_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTA>::value; \
    bool B_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTB>::value; \
    bool C_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTC>::value; \
    \
    const int AST = A_is_lr?A.stride(0):A.stride(1), LDA = AST == 0 ? 1 : AST; \
    const int BST = B_is_lr?B.stride(0):B.stride(1), LDB = BST == 0 ? 1 : BST; \
    const int CST = C_is_lr?C.stride(0):C.stride(1), LDC = CST == 0 ? 1 : CST; \
    \
    if(!A_is_lr && !B_is_lr && !C_is_lr ) \
      zgemm_(transA,transB,&M,&N,&K, \
        reinterpret_cast<const std::complex<double>*>(&alpha),reinterpret_cast<const std::complex<double>*>(A.data()),&LDA, \
        reinterpret_cast<const std::complex<double>*>(B.data()),&LDB, \
        reinterpret_cast<const std::complex<double>*>(&beta),reinterpret_cast<std::complex<double>*>(C.data()),&LDC); \
    if(A_is_lr && B_is_lr && C_is_lr ) \
      zgemm_(transB,transA,&N,&M,&K, \
        reinterpret_cast<const std::complex<double>*>(&alpha),reinterpret_cast<const std::complex<double>*>(B.data()),&LDB, \
        reinterpret_cast<const std::complex<double>*>(A.data()),&LDA, \
        reinterpret_cast<const std::complex<double>*>(&beta),reinterpret_cast<std::complex<double>*>(C.data()),&LDC); \
    Kokkos::Profiling::popRegion(); \
  } \
}; \

#define KOKKOSBLAS3_CGEMM_BLAS( LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GEMM< \
     Kokkos::View<const Kokkos::complex<float>**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const Kokkos::complex<float>**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<Kokkos::complex<float>**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef Kokkos::complex<float> SCALAR; \
  typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > CViewType; \
      \
  static void \
  gemm (const char transA[], \
        const char transB[], \
        typename AViewType::const_value_type& alpha, \
        const AViewType& A, \
        const BViewType& B, \
        typename CViewType::const_value_type& beta, \
        const CViewType& C) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_BLAS,complex<float>]"); \
    const bool A_t = (transA[0]!='N') && (transA[0]!='n'); \
    const int M = C.extent(0); \
    const int N = C.extent(1); \
    const int K = A.extent(A_t?0:1); \
    \
    bool A_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTA>::value; \
    bool B_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTB>::value; \
    bool C_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTC>::value; \
    \
    const int AST = A_is_lr?A.stride(0):A.stride(1), LDA = AST == 0 ? 1 : AST; \
    const int BST = B_is_lr?B.stride(0):B.stride(1), LDB = BST == 0 ? 1 : BST; \
    const int CST = C_is_lr?C.stride(0):C.stride(1), LDC = CST == 0 ? 1 : CST; \
    \
    if(!A_is_lr && !B_is_lr && !C_is_lr ) \
      cgemm_(transA,transB,&M,&N,&K, \
        reinterpret_cast<const std::complex<float>*>(&alpha),reinterpret_cast<const std::complex<float>*>(A.data()),&LDA, \
        reinterpret_cast<const std::complex<float>*>(B.data()),&LDB, \
        reinterpret_cast<const std::complex<float>*>(&beta),reinterpret_cast<std::complex<float>*>(C.data()),&LDC); \
    if(A_is_lr && B_is_lr && C_is_lr ) \
      cgemm_(transB,transA,&N,&M,&K, \
        reinterpret_cast<const std::complex<float>*>(&alpha),reinterpret_cast<const std::complex<float>*>(B.data()),&LDB, \
        reinterpret_cast<const std::complex<float>*>(A.data()),&LDA, \
        reinterpret_cast<const std::complex<float>*>(&beta),reinterpret_cast<std::complex<float>*>(C.data()),&LDC); \
    Kokkos::Profiling::popRegion(); \
  } \
};

KOKKOSBLAS3_DGEMM_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_DGEMM_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_DGEMM_BLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_DGEMM_BLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS3_SGEMM_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_SGEMM_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_SGEMM_BLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_SGEMM_BLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS3_ZGEMM_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_ZGEMM_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_ZGEMM_BLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_ZGEMM_BLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

KOKKOSBLAS3_CGEMM_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS3_CGEMM_BLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::HostSpace, false)
KOKKOSBLAS3_CGEMM_BLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, true)
KOKKOSBLAS3_CGEMM_BLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::HostSpace, false)

}
}
#endif // KOKKOSKERNELS_ENABLE_TPL_BLAS

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include<KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS3_DGEMM_CUBLAS( LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GEMM< \
     Kokkos::View<const double**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const double**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<double**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef double SCALAR; \
  typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > CViewType; \
 \
  static void \
  gemm (const char transA[], \
        const char transB[], \
        typename AViewType::const_value_type& alpha, \
        const AViewType& A, \
        const BViewType& B, \
        typename CViewType::const_value_type& beta, \
        const CViewType& C) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_BLAS,double]"); \
    const bool A_t = (transA[0]!='N') && (transA[0]!='n'); \
    const int M = static_cast<int> (C.extent(0)); \
    const int N = static_cast<int> (C.extent(1)); \
    const int K = static_cast<int> (A.extent(A_t?0:1)); \
    \
    bool A_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTA>::value; \
    bool B_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTB>::value; \
    bool C_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTC>::value; \
    \
    const int AST = A_is_lr?A.stride(0):A.stride(1), LDA = AST == 0 ? 1 : AST; \
    const int BST = B_is_lr?B.stride(0):B.stride(1), LDB = BST == 0 ? 1 : BST; \
    const int CST = C_is_lr?C.stride(0):C.stride(1), LDC = CST == 0 ? 1 : CST; \
    \
    cublasOperation_t transa, transb; \
    if ((transA[0]=='N')||(transA[0]=='n')) \
      transa = CUBLAS_OP_N; \
    else if ((transA[0]=='T')||(transA[0]=='t')) \
      transa = CUBLAS_OP_T; \
    else \
      transa = CUBLAS_OP_C; \
    if ((transB[0]=='N')||(transB[0]=='n')) \
      transb = CUBLAS_OP_N; \
    else if ((transB[0]=='T')||(transB[0]=='t')) \
      transb = CUBLAS_OP_T; \
    else \
      transb = CUBLAS_OP_C; \
    KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
    if(!A_is_lr && !B_is_lr && !C_is_lr ) \
      cublasDgemm(s.handle, transa, transb, M, N, K, &alpha, A.data(), LDA, B.data(), LDB, &beta, C.data(), LDC); \
    if(A_is_lr && B_is_lr && C_is_lr ) \
      cublasDgemm(s.handle, transb, transa, N, M, K, &alpha, B.data(), LDB, A.data(), LDA, &beta, C.data(), LDC); \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS3_SGEMM_CUBLAS( LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GEMM< \
     Kokkos::View<const float**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const float**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<float**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef float SCALAR; \
  typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > CViewType; \
      \
  static void \
  gemm (const char transA[], \
        const char transB[], \
        typename AViewType::const_value_type& alpha, \
        const AViewType& A, \
        const BViewType& B, \
        typename CViewType::const_value_type& beta, \
        const CViewType& C) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_BLAS,float]"); \
    const bool A_t = (transA[0]!='N') && (transA[0]!='n'); \
    const int M = static_cast<int> (C.extent(0)); \
    const int N = static_cast<int> (C.extent(1)); \
    const int K = static_cast<int> (A.extent(A_t?0:1)); \
    \
    bool A_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTA>::value; \
    bool B_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTB>::value; \
    bool C_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTC>::value; \
    \
    const int AST = A_is_lr?A.stride(0):A.stride(1), LDA = AST == 0 ? 1 : AST; \
    const int BST = B_is_lr?B.stride(0):B.stride(1), LDB = BST == 0 ? 1 : BST; \
    const int CST = C_is_lr?C.stride(0):C.stride(1), LDC = CST == 0 ? 1 : CST; \
    \
    cublasOperation_t transa, transb; \
    if ((transA[0]=='N')||(transA[0]=='n')) \
      transa = CUBLAS_OP_N; \
    else if ((transA[0]=='T')||(transA[0]=='t')) \
      transa = CUBLAS_OP_T; \
    else \
      transa = CUBLAS_OP_C; \
    if ((transB[0]=='N')||(transB[0]=='n')) \
      transb = CUBLAS_OP_N; \
    else if ((transB[0]=='T')||(transB[0]=='t')) \
      transb = CUBLAS_OP_T; \
    else \
      transb = CUBLAS_OP_C; \
    KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
    if(!A_is_lr && !B_is_lr && !C_is_lr ) \
      cublasSgemm(s.handle, transa, transb, M, N, K, &alpha, A.data(), LDA, B.data(), LDB, &beta, C.data(), LDC); \
    if(A_is_lr && B_is_lr && C_is_lr ) \
      cublasSgemm(s.handle, transb, transa, N, M, K, &alpha, B.data(), LDB, A.data(), LDA, &beta, C.data(), LDC); \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS3_ZGEMM_CUBLAS( LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GEMM< \
     Kokkos::View<const Kokkos::complex<double>**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const Kokkos::complex<double>**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<Kokkos::complex<double>**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef Kokkos::complex<double> SCALAR; \
  typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > CViewType; \
      \
  static void \
  gemm (const char transA[], \
        const char transB[], \
        typename AViewType::const_value_type& alpha, \
        const AViewType& A, \
        const BViewType& B, \
        typename CViewType::const_value_type& beta, \
        const CViewType& C) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_BLAS,complex<double>]"); \
    const bool A_t = (transA[0]!='N') && (transA[0]!='n'); \
    const int M = static_cast<int> (C.extent(0)); \
    const int N = static_cast<int> (C.extent(1)); \
    const int K = static_cast<int> (A.extent(A_t?0:1)); \
    \
    bool A_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTA>::value; \
    bool B_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTB>::value; \
    bool C_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTC>::value; \
    \
    const int AST = A_is_lr?A.stride(0):A.stride(1), LDA = AST == 0 ? 1 : AST; \
    const int BST = B_is_lr?B.stride(0):B.stride(1), LDB = BST == 0 ? 1 : BST; \
    const int CST = C_is_lr?C.stride(0):C.stride(1), LDC = CST == 0 ? 1 : CST; \
    \
    cublasOperation_t transa, transb; \
    if ((transA[0]=='N')||(transA[0]=='n')) \
      transa = CUBLAS_OP_N; \
    else if ((transA[0]=='T')||(transA[0]=='t')) \
      transa = CUBLAS_OP_T; \
    else \
      transa = CUBLAS_OP_C; \
    if ((transB[0]=='N')||(transB[0]=='n')) \
      transb = CUBLAS_OP_N; \
    else if ((transB[0]=='T')||(transB[0]=='t')) \
      transb = CUBLAS_OP_T; \
    else \
      transb = CUBLAS_OP_C; \
    KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
    if(!A_is_lr && !B_is_lr && !C_is_lr ) \
      cublasZgemm(s.handle, transa, transb, M, N, K, reinterpret_cast<const cuDoubleComplex*>(&alpha), reinterpret_cast<const cuDoubleComplex*>(A.data()), LDA, reinterpret_cast<const cuDoubleComplex*>(B.data()), LDB, reinterpret_cast<const cuDoubleComplex*>(&beta), reinterpret_cast<cuDoubleComplex*>(C.data()), LDC); \
    if(A_is_lr && B_is_lr && C_is_lr ) \
      cublasZgemm(s.handle, transb, transa, N, M, K, reinterpret_cast<const cuDoubleComplex*>(&alpha), reinterpret_cast<const cuDoubleComplex*>(B.data()), LDB, reinterpret_cast<const cuDoubleComplex*>(A.data()), LDA, reinterpret_cast<const cuDoubleComplex*>(&beta), reinterpret_cast<cuDoubleComplex*>(C.data()), LDC); \
    Kokkos::Profiling::popRegion(); \
  } \
}; \

#define KOKKOSBLAS3_CGEMM_CUBLAS( LAYOUTA, LAYOUTB, LAYOUTC, MEM_SPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct GEMM< \
     Kokkos::View<const Kokkos::complex<float>**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const Kokkos::complex<float>**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<Kokkos::complex<float>**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     true, ETI_SPEC_AVAIL> { \
  typedef Kokkos::complex<float> SCALAR; \
  typedef Kokkos::View<const SCALAR**, LAYOUTA, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > AViewType; \
  typedef Kokkos::View<const SCALAR**, LAYOUTB, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > BViewType; \
  typedef Kokkos::View<SCALAR**, LAYOUTC, Kokkos::Device<ExecSpace, MEM_SPACE>, \
      Kokkos::MemoryTraits<Kokkos::Unmanaged> > CViewType; \
      \
  static void \
  gemm (const char transA[], \
        const char transB[], \
        typename AViewType::const_value_type& alpha, \
        const AViewType& A, \
        const BViewType& B, \
        typename CViewType::const_value_type& beta, \
        const CViewType& C) { \
    \
    Kokkos::Profiling::pushRegion("KokkosBlas::gemm[TPL_BLAS,complex<float>]"); \
    const bool A_t = (transA[0]!='N') && (transA[0]!='n'); \
    const int M = static_cast<int> (C.extent(0)); \
    const int N = static_cast<int> (C.extent(1)); \
    const int K = static_cast<int> (A.extent(A_t?0:1)); \
    \
    bool A_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTA>::value; \
    bool B_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTB>::value; \
    bool C_is_lr = std::is_same<Kokkos::LayoutRight,LAYOUTC>::value; \
    \
    const int AST = A_is_lr?A.stride(0):A.stride(1), LDA = AST == 0 ? 1 : AST; \
    const int BST = B_is_lr?B.stride(0):B.stride(1), LDB = BST == 0 ? 1 : BST; \
    const int CST = C_is_lr?C.stride(0):C.stride(1), LDC = CST == 0 ? 1 : CST; \
    \
    cublasOperation_t transa, transb; \
    if ((transA[0]=='N')||(transA[0]=='n')) \
      transa = CUBLAS_OP_N; \
    else if ((transA[0]=='T')||(transA[0]=='t')) \
      transa = CUBLAS_OP_T; \
    else \
      transa = CUBLAS_OP_C; \
    if ((transB[0]=='N')||(transB[0]=='n')) \
      transb = CUBLAS_OP_N; \
    else if ((transB[0]=='T')||(transB[0]=='t')) \
      transb = CUBLAS_OP_T; \
    else \
      transb = CUBLAS_OP_C; \
    KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
    if(!A_is_lr && !B_is_lr && !C_is_lr ) \
      cublasCgemm(s.handle, transa, transb, M, N, K, reinterpret_cast<const cuComplex*>(&alpha), reinterpret_cast<const cuComplex*>(A.data()), LDA, reinterpret_cast<const cuComplex*>(B.data()), LDB, reinterpret_cast<const cuComplex*>(&beta), reinterpret_cast<cuComplex*>(C.data()), LDC); \
    if(A_is_lr && B_is_lr && C_is_lr ) \
      cublasCgemm(s.handle, transb, transa, N, M, K, reinterpret_cast<const cuComplex*>(&alpha), reinterpret_cast<const cuComplex*>(B.data()), LDB, reinterpret_cast<const cuComplex*>(A.data()), LDA, reinterpret_cast<const cuComplex*>(&beta), reinterpret_cast<cuComplex*>(C.data()), LDC); \
    Kokkos::Profiling::popRegion(); \
  } \
};

KOKKOSBLAS3_DGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_DGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_DGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_DGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_SGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_SGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_SGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_SGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_ZGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_ZGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_ZGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_ZGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_CGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_CGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_CGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_CGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

}
}
#endif // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#endif
