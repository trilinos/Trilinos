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
#include "KokkosBlas_Host_tpl.hpp"

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
    if(!A_is_lr && !B_is_lr && !C_is_lr )                               \
      HostBlas<double>::gemm                                            \
        (transA[0],transB[0],                                                 \
         M,N,K,                                                         \
         alpha,                                                         \
         A.data(),LDA,                                                  \
         B.data(),LDB,                                                  \
         beta,                                                          \
         C.data(),LDC);                                                 \
    if(A_is_lr && B_is_lr && C_is_lr )                                  \
      HostBlas<double>::gemm                                            \
        (transB[0],transA[0],                                                 \
         N,M,K,                                                         \
         alpha,                                                         \
         B.data(),LDB,                                                  \
         A.data(),LDA,                                                  \
         beta,                                                          \
         C.data(),LDC);                                                 \
    Kokkos::Profiling::popRegion();                                     \
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
    if(!A_is_lr && !B_is_lr && !C_is_lr )                              \
      HostBlas<float>::gemm                                             \
        (transA[0],transB[0],                                                 \
         M,N,K,                                                         \
         alpha,                                                         \
         A.data(),LDA,                                                  \
         B.data(),LDB,                                                  \
         beta,                                                          \
         C.data(),LDC);                                                 \
    if(A_is_lr && B_is_lr && C_is_lr )                                  \
      HostBlas<float>::gemm                                             \
        (transB[0],transA[0],                                                 \
         N,M,K,                                                         \
         alpha,                                                         \
         B.data(),LDB,                                                  \
         A.data(),LDA,                                                  \
         beta,                                                          \
         C.data(),LDC);                                                 \
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
    const std::complex<double> alpha_val = alpha, beta_val = beta;     \
    if(!A_is_lr && !B_is_lr && !C_is_lr )                              \
      HostBlas<std::complex<double> >::gemm                             \
        (transA[0],transB[0],                                                 \
         M,N,K,                                                         \
         alpha_val,          \
         reinterpret_cast<const std::complex<double>*>(A.data()),LDA,   \
         reinterpret_cast<const std::complex<double>*>(B.data()),LDB,   \
         beta_val,           \
         reinterpret_cast<      std::complex<double>*>(C.data()),LDC);  \
    if(A_is_lr && B_is_lr && C_is_lr )                                  \
      HostBlas<std::complex<double> >::gemm                              \
        (transB[0],transA[0],                                                 \
         N,M,K,                                                         \
         alpha_val,          \
         reinterpret_cast<const std::complex<double>*>(B.data()),LDB,   \
         reinterpret_cast<const std::complex<double>*>(A.data()),LDA,   \
         beta_val,           \
         reinterpret_cast<      std::complex<double>*>(C.data()),LDC);  \
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
    const std::complex<float> alpha_val = alpha, beta_val = beta;      \
    if(!A_is_lr && !B_is_lr && !C_is_lr )                              \
      HostBlas<std::complex<float> >::gemm                             \
        (transA[0],transB[0],                                                 \
         M,N,K,                                                         \
         alpha_val,          \
         reinterpret_cast<const std::complex<float>*>(A.data()),LDA,   \
         reinterpret_cast<const std::complex<float>*>(B.data()),LDB,   \
         beta_val,           \
         reinterpret_cast<      std::complex<float>*>(C.data()),LDC);  \
    if(A_is_lr && B_is_lr && C_is_lr )                                  \
      HostBlas<std::complex<float> >::gemm                              \
        (transB[0],transA[0],                                                 \
         N,M,K,                                                         \
         alpha_val,          \
         reinterpret_cast<const std::complex<float>*>(B.data()),LDB,   \
         reinterpret_cast<const std::complex<float>*>(A.data()),LDA,   \
         beta_val,           \
         reinterpret_cast<      std::complex<float>*>(C.data()),LDC);  \
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


// DotBasedGEMM implements the optimization for C = beta*C + alpha*A^TB 
// with A and B matrices both being tall and skinny. C matrix is assumably 
// small, so, each entry of C is computed by performing the dot product of 
// respective columns of A and B matrices. Note that the dot products are
// performed on very long vectors, so, each dot product is distributed among
// numDivPerDot teams.     

struct TagZero{};   // The init tag for beta=0 
struct TagInit{};   // The init tag for beta!=0 and beta !=1 
struct TagMult{};   // The multiplication tag for transposed A
struct TagMultCT{};   // The multiplication tag for conjugate-transposed A 
template<class ExecSpace, class AV, class BV, class CV>
struct DotBasedGEMM{

  const AV A;
  const BV B;
  CV C;

  using scalar_A = typename AV::non_const_value_type;
  using size_A = typename AV::size_type;
  using scalar_C = typename CV::non_const_value_type;
  using size_C = typename CV::size_type;
  using AVT = Kokkos::Details::ArithTraits<scalar_A>;
  using CVT = Kokkos::Details::ArithTraits<scalar_C>;

  const scalar_A alpha;
  const scalar_C beta;

  // The following types (especially dotSize) could have simply been int,
  const size_C numCrows;           
  const size_C numCcols;

  size_C numDivPerDot;   // number of teams collectively performing a dot product
  size_C numTeams;       // total number of teams
  
  const size_A dotSize;  // the length of the vectors in the dot products
  size_A chunkSize;      // the local length of each team's share on the dot product  
  

  DotBasedGEMM(const scalar_A& alpha_, const AV& A_, const BV& B_, const scalar_C& beta_, const CV& C_):A(A_),B(B_),C(C_),alpha(alpha_),beta(beta_),numCrows(C.extent(0)),numCcols(C.extent(1)),dotSize(A.extent(0))
  { }

  void run(bool conjugateTranspose) {

    constexpr size_C workPerTeam = 4096;                   // Amount of work per team
    const size_C ndots = numCrows * numCcols;              // Number of dot products
    size_C appxNumTeams = (dotSize * ndots) / workPerTeam; // Estimation for appxNumTeams

    // Adjust appxNumTeams in case it is too small or too large
    if(appxNumTeams < 1)
      appxNumTeams = 1;
    if(appxNumTeams > 1024)
      appxNumTeams = 1024;

    // If there are more dot products than the number of teams,
    // then set the number of teams to be number of dot products
    // and each team will perform only one dot product.
    // We don't want a team to perform more than one dot product.
    if(ndots >= appxNumTeams) {
      numTeams = ndots;
      numDivPerDot = 1;
    }
    // If there are more teams than dot products, each dot product can
    // potentially be performed by multiple teams. First, compute 
    // numDivPerDot as an integer (take the floor, not ceiling), then,
    // compute actual number of teams by using this factor.
    else{
      numDivPerDot = appxNumTeams / ndots;
      numTeams = ndots * numDivPerDot;
    }

    // Determine the local length for the dot product
    chunkSize = dotSize / numDivPerDot;
    if(numDivPerDot > 1)
      chunkSize++;

    // Initialize C matrix if beta != 1
    if(beta == CVT::zero()) {
      Kokkos::MDRangePolicy<TagZero, ExecSpace, Kokkos::Rank<2>> policyInit({0,0}, {numCrows, numCcols});
      Kokkos::parallel_for("Initialize C for Dot Product Based GEMM", policyInit, *this);
    }
    else if(beta != CVT::one()) {
      Kokkos::MDRangePolicy<TagInit, ExecSpace, Kokkos::Rank<2>> policyInit({0,0}, {numCrows, numCcols});
      Kokkos::parallel_for("Initialize C for Dot Product Based GEMM", policyInit, *this);
    }
    
    // Multiply alpha*A^TB and add it to beta*C
    if(conjugateTranspose) {
      Kokkos::TeamPolicy<TagMultCT, ExecSpace> policyMult(numTeams, Kokkos::AUTO);
      Kokkos::parallel_for("Perform Dot Product Based GEMM", policyMult, *this);
    }
    else{
      Kokkos::TeamPolicy<TagMult, ExecSpace> policyMult(numTeams, Kokkos::AUTO);
      Kokkos::parallel_for("Perform Dot Product Based GEMM", policyMult, *this);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const TagZero&, const size_C &rowId, const size_C &colId ) const {
    C(rowId, colId) = CVT::zero(); 
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const TagInit&, const size_C &rowId, const size_C &colId ) const {
    C(rowId, colId) = beta * C(rowId, colId);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const TagMult&, const typename Kokkos::TeamPolicy<>::member_type& teamMember) const {

    const size_C globalRank = teamMember.league_rank();
    const size_C localRank = globalRank % numDivPerDot;
    const size_C i = globalRank / numDivPerDot;
    const size_C rowId = i / numCcols;
    const size_C colId = i % numCcols;
    
    scalar_C result = CVT::zero();
    const size_A baseInd = chunkSize*localRank; 
    Kokkos::parallel_reduce( Kokkos::TeamThreadRange(teamMember, chunkSize), [&]( const size_A k, scalar_C &update ) {
	if(baseInd + k < dotSize)
	  update += alpha * A(baseInd+k, rowId) * B(baseInd+k, colId);
      }, result );

    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () { 
      Kokkos::atomic_add(&C(rowId, colId), result);
      });
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const TagMultCT&, const typename Kokkos::TeamPolicy<>::member_type& teamMember) const {

    const size_C globalRank = teamMember.league_rank();
    const size_C localRank = globalRank % numDivPerDot;
    const size_C i = globalRank / numDivPerDot;
    const size_C rowId = i / numCcols;
    const size_C colId = i % numCcols;
    
    scalar_C result = CVT::zero();
    const size_A baseInd = chunkSize*localRank; 
    Kokkos::parallel_reduce( Kokkos::TeamThreadRange(teamMember, chunkSize), [&]( const size_A k, scalar_C &update ) {
	if(baseInd + k < dotSize)
	  update += alpha * AVT::conj(A(baseInd+k, rowId)) * B(baseInd+k, colId);
      }, result );

    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () { 
      Kokkos::atomic_add(&C(rowId, colId), result);
      });
  }

};

  
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
    \
    constexpr int numDotsLayoutLeftThreshold = 1600; \
    constexpr int numDotsLayoutRightThreshold = 100; \
    if(   (!A_is_lr && transa != CUBLAS_OP_N && transb == CUBLAS_OP_N && M*N < numDotsLayoutLeftThreshold) \
       || ( A_is_lr && transa != CUBLAS_OP_N && transb == CUBLAS_OP_N && M*N < numDotsLayoutRightThreshold)) { \
      DotBasedGEMM<ExecSpace,AViewType,BViewType,CViewType> gemm(alpha,A,B,beta,C); \
      gemm.run(false); \
    } \
    else { \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      if(!A_is_lr && !B_is_lr && !C_is_lr )				\
	cublasDgemm(s.handle, transa, transb, M, N, K, &alpha, A.data(), LDA, B.data(), LDB, &beta, C.data(), LDC); \
      if(A_is_lr && B_is_lr && C_is_lr )				\
	cublasDgemm(s.handle, transb, transa, N, M, K, &alpha, B.data(), LDB, A.data(), LDA, &beta, C.data(), LDC); \
    } \
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
    \
    constexpr int numDotsLayoutLeftThreshold = 1600; \
    constexpr int numDotsLayoutRightThreshold = 100; \
    if(   (!A_is_lr && transa != CUBLAS_OP_N && transb == CUBLAS_OP_N && M*N < numDotsLayoutLeftThreshold) \
       || ( A_is_lr && transa != CUBLAS_OP_N && transb == CUBLAS_OP_N && M*N < numDotsLayoutRightThreshold)) { \
      DotBasedGEMM<ExecSpace,AViewType,BViewType,CViewType> gemm(alpha,A,B,beta,C); \
      gemm.run(false); \
    } \
    else { \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      if(!A_is_lr && !B_is_lr && !C_is_lr ) \
        cublasSgemm(s.handle, transa, transb, M, N, K, &alpha, A.data(), LDA, B.data(), LDB, &beta, C.data(), LDC); \
      if(A_is_lr && B_is_lr && C_is_lr ) \
        cublasSgemm(s.handle, transb, transa, N, M, K, &alpha, B.data(), LDB, A.data(), LDA, &beta, C.data(), LDC); \
    } \
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
    \
    constexpr int numDotsLayoutLeftThreshold = 1600; \
    constexpr int numDotsLayoutRightThreshold = 100; \
    if(   (!A_is_lr && transa != CUBLAS_OP_N && transb == CUBLAS_OP_N && M*N < numDotsLayoutLeftThreshold) \
       || ( A_is_lr && transa != CUBLAS_OP_N && transb == CUBLAS_OP_N && M*N < numDotsLayoutRightThreshold)) { \
      DotBasedGEMM<ExecSpace,AViewType,BViewType,CViewType> gemm(alpha,A,B,beta,C); \
      gemm.run(transa == CUBLAS_OP_C ? true : false);  \
    } \
    else { \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      if(!A_is_lr && !B_is_lr && !C_is_lr ) \
	cublasZgemm(s.handle, transa, transb, M, N, K, reinterpret_cast<const cuDoubleComplex*>(&alpha), reinterpret_cast<const cuDoubleComplex*>(A.data()), LDA, reinterpret_cast<const cuDoubleComplex*>(B.data()), LDB, reinterpret_cast<const cuDoubleComplex*>(&beta), reinterpret_cast<cuDoubleComplex*>(C.data()), LDC); \
      if(A_is_lr && B_is_lr && C_is_lr ) \
	cublasZgemm(s.handle, transb, transa, N, M, K, reinterpret_cast<const cuDoubleComplex*>(&alpha), reinterpret_cast<const cuDoubleComplex*>(B.data()), LDB, reinterpret_cast<const cuDoubleComplex*>(A.data()), LDA, reinterpret_cast<const cuDoubleComplex*>(&beta), reinterpret_cast<cuDoubleComplex*>(C.data()), LDC); \
    } \
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
    \
    constexpr int numDotsLayoutLeftThreshold = 1600; \
    constexpr int numDotsLayoutRightThreshold = 100; \
    if(   (!A_is_lr && transa != CUBLAS_OP_N && transb == CUBLAS_OP_N && M*N < numDotsLayoutLeftThreshold) \
       || ( A_is_lr && transa != CUBLAS_OP_N && transb == CUBLAS_OP_N && M*N < numDotsLayoutRightThreshold)) { \
      DotBasedGEMM<ExecSpace,AViewType,BViewType,CViewType> gemm(alpha,A,B,beta,C); \
      gemm.run(transa == CUBLAS_OP_C ? true : false);  \
    } \
    else { \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      if(!A_is_lr && !B_is_lr && !C_is_lr ) \
	cublasCgemm(s.handle, transa, transb, M, N, K, reinterpret_cast<const cuComplex*>(&alpha), reinterpret_cast<const cuComplex*>(A.data()), LDA, reinterpret_cast<const cuComplex*>(B.data()), LDB, reinterpret_cast<const cuComplex*>(&beta), reinterpret_cast<cuComplex*>(C.data()), LDC); \
      if(A_is_lr && B_is_lr && C_is_lr ) \
	cublasCgemm(s.handle, transb, transa, N, M, K, reinterpret_cast<const cuComplex*>(&alpha), reinterpret_cast<const cuComplex*>(B.data()), LDB, reinterpret_cast<const cuComplex*>(A.data()), LDA, reinterpret_cast<const cuComplex*>(&beta), reinterpret_cast<cuComplex*>(C.data()), LDC); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

KOKKOSBLAS3_DGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_DGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_DGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_DGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_DGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_DGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_DGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_DGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS3_SGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_SGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_SGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_SGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_SGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_SGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_SGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_SGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS3_ZGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_ZGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_ZGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_ZGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_ZGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_ZGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_ZGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_ZGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS3_CGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS3_CGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace, false)
KOKKOSBLAS3_CGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, true)
KOKKOSBLAS3_CGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace, false)

KOKKOSBLAS3_CGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_CGEMM_CUBLAS( Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
KOKKOSBLAS3_CGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS3_CGEMM_CUBLAS( Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace, false)

}
}
#endif // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#endif
