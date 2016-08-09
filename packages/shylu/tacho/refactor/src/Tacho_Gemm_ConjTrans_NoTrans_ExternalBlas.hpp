#ifndef __TACHO_GEMM_CONJTRANS_NOTRANS_EXTERNAL_BLAS_HPP__
#define __TACHO_GEMM_CONJTRANS_NOTRANS_EXTERNAL_BLAS_HPP__

/// \file Tacho_Gemm_ConjTrans_NoTrans_ExternalBlas.hpp
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
  Gemm<Trans::ConjTranspose,Trans::NoTranspose,
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

    if (member.team_rank() == 0) {
#ifdef HAVE_SHYLUTACHO_TEUCHOS    
      Teuchos::BLAS<ordinal_type,value_type> blas;
      
      const ordinal_type m = C.NumRows();
      const ordinal_type n = C.NumCols();
      const ordinal_type k = B.NumRows();

      if (m > 0 && n > 0 && k > 0)
        blas.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS,
                  m, n, k,
                  alpha,
                  A.ValuePtr(), A.BaseObject().ColStride(),
                  B.ValuePtr(), B.BaseObject().ColStride(),
                  beta,
                  C.ValuePtr(), C.BaseObject().ColStride());
#else
    TACHO_TEST_FOR_ABORT( true, MSG_NOT_HAVE_PACKAGE("Teuchos") );        
#endif
    } 

    return 0;
  }
}

#endif
