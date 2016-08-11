#ifndef __TACHO_GEMM_NOTRANS_NOTRANS_EXTERNAL_BLAS_HPP__
#define __TACHO_GEMM_NOTRANS_NOTRANS_EXTERNAL_BLAS_HPP__

/// \file Tacho_Gemm_NoTrans_NoTrans_ExternalBlas.hpp
/// \brief BLAS matrix-matrix multiplication 
/// \author Kyungjoo Kim (kyukim@sandia.gov)
#ifdef HAVE_SHYLUTACHO_TEUCHOS
#include "Teuchos_BLAS.hpp"
#endif

namespace Tacho {
  /// BLAS Gemm
  /// =========
  /// Properties:
  /// - Compile with Device (o),
  /// - Callable in KokkosFunctors (o)
  /// - For now, this is for HostSpace only.
  template<>
  template<typename PolicyType,
           typename MemberType,
           typename ScalarType,
           typename DenseExecViewTypeA,
           typename DenseExecViewTypeB,
           typename DenseExecViewTypeC>
  KOKKOS_INLINE_FUNCTION
  int
  Gemm<Trans::NoTranspose,Trans::NoTranspose,
       AlgoGemm::ExternalBlas,Variant::One>
  ::invoke(PolicyType &policy,
           MemberType &member,
           const ScalarType alpha,
           DenseExecViewTypeA &A,
           DenseExecViewTypeB &B,
           const ScalarType beta,
           DenseExecViewTypeC &C) {
    // static_assert( Kokkos::Impl::is_same<
    //                typename DenseMatrixTypeA::space_type,
    //                Kokkos::Cuda
    //                >::value,
    //                "Cuda space is not available for calling external BLAS" );

    // static_assert( Kokkos::Impl::is_same<
    //                typename DenseMatrixTypeA::space_type,
    //                typename DenseMatrixTypeB::space_type
    //                >::value && 
    //                Kokkos::Impl::is_same<
    //                typename DenseMatrixTypeB::space_type,
    //                typename DenseMatrixTypeC::space_type
    //                >::value,
    //                "Space type of input matrices does not match" );
    
    //typedef typename DenseExecViewTypeA::space_type   space_type;
    typedef typename DenseExecViewTypeA::ordinal_type ordinal_type;
    typedef typename DenseExecViewTypeA::value_type   value_type;

#ifndef HAVE_SHYLUTACHO_TEUCHOS    
    TACHO_TEST_FOR_ABORT( true, MSG_NOT_HAVE_PACKAGE("Teuchos") );            
#endif

    if (member.team_size() == 1) {
      // Multithreaded MKL handle its own threads
#ifdef HAVE_SHYLUTACHO_TEUCHOS    
      Teuchos::BLAS<ordinal_type,value_type> blas;
      
      const ordinal_type m = C.NumRows();
      const ordinal_type n = C.NumCols();
      const ordinal_type k = B.NumRows();
      
      if (m > 0 && n > 0 && k > 0)
        blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS,
                  m, n, k,
                  value_type(alpha),
                  A.ValuePtr(), A.BaseObject().ColStride(),
                  B.ValuePtr(), B.BaseObject().ColStride(),
                  value_type(beta),
                  C.ValuePtr(), C.BaseObject().ColStride());
#endif
    } else {
      // Sequential MKL is invoked in the team interface
#ifdef HAVE_SHYLUTACHO_TEUCHOS    
      Teuchos::BLAS<ordinal_type,value_type> blas;
      
      const ordinal_type m = C.NumRows();
      const ordinal_type n = C.NumCols();
      const ordinal_type k = B.NumRows();
      
      if (m > 0 && n > 0 && k > 0) {
        const ordinal_type b = 32;
        const ordinal_type
          pend  = k/b + (k%b > 0),
          k2end = n/b + (n%b > 0),
          k1end = m/b + (m%b > 0);
        const ordinal_type nwork = k1end*k2end;

        for (auto p=0;p<pend;++p) {
          const auto beta_select = (p > 0 ? ScalarType(1.0) : beta);

          const ordinal_type koff = p*b;
          const ordinal_type kk = ((koff + b) > k ? (k - koff) : b);

          // I want parallel chunk size is 1
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, nwork),
                               [&](const ordinal_type iwork) {
                                 const ordinal_type k1 = iwork%k1end;
                                 const ordinal_type k2 = iwork/k1end;
                                 
                                 const ordinal_type moff = k1*b;
                                 const ordinal_type mm = ((moff + b) > m ? (m - moff) : b);

                                 const ordinal_type noff = k2*b;
                                 const ordinal_type nn = ((noff + b) > n ? (n - noff) : b);
                                 
                                 // A part is : k1, p
                                 const auto ldim_a = A.BaseObject().ColStride();
                                 const auto ptr_a  = A.ValuePtr() + moff + koff*ldim_a;
                                 
                                 // B part is : p, k2
                                 const auto ldim_b = B.BaseObject().ColStride();
                                 const auto ptr_b  = B.ValuePtr() + koff + noff*ldim_b;

                                 // C part is : k1, k2
                                 const auto ldim_c = C.BaseObject().ColStride();
                                 const auto ptr_c  = C.ValuePtr() + moff + noff*ldim_c;
                                 
                                 blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS,
                                           mm, nn, kk,
                                           value_type(alpha),
                                           ptr_a, ldim_a, 
                                           ptr_b, ldim_b,
                                           value_type(beta_select),
                                           ptr_c, ldim_c);
                               });
          member.team_barrier();
        }
      }
#endif
    } 

    return 0;
  }
}

#endif
