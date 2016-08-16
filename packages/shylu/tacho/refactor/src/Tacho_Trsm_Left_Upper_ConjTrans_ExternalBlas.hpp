#ifndef __TACHO_TRSM_LEFT_UPPER_CONJTRANS_EXTERNAL_BLAS_HPP__
#define __TACHO_TRSM_LEFT_UPPER_CONJTRANS_EXTERNAL_BLAS_HPP__


/// \file  Tacho_Trsm_Left_Upper_ConjTrans_ExternalBlas.hpp
/// \brief BLAS triangular solve matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)
#ifdef HAVE_SHYLUTACHO_TEUCHOS
#include "Teuchos_BLAS.hpp"
#endif

namespace Tacho {
  /// BLAS Trsm
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
           typename DenseExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
       AlgoTrsm::ExternalBlas,Variant::One>
  ::invoke(PolicyType &policy,
           MemberType &member,
           const int diagA,
           const ScalarType alpha,
           DenseExecViewTypeA &A,
           DenseExecViewTypeB &B) {
    // static_assert( Kokkos::Impl::is_same<
    //                typename DenseMatrixTypeA::space_type,
    //                typename DenseMatrixTypeB::space_type
    //                >::value,
    //                "Space type of input matrices does not match" );

    if (member.team_rank() == 0) {
#if                                                     \
  defined( HAVE_SHYLUTACHO_TEUCHOS ) &&                 \
  defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      typedef typename DenseExecViewTypeA::ordinal_type ordinal_type;
      typedef typename DenseExecViewTypeA::value_type   value_type;
      
      Teuchos::BLAS<ordinal_type,value_type> blas;

      const ordinal_type m = A.NumRows();
      const ordinal_type n = B.NumCols();

      if (m > 0 && n > 0) 
        blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::CONJ_TRANS,
                  (diagA == Diag::Unit ? Teuchos::UNIT_DIAG : Teuchos::NON_UNIT_DIAG),
                  m, n,
                  alpha,
                  A.ValuePtr(), A.BaseObject().ColStride(),
                  B.ValuePtr(), B.BaseObject().ColStride());
#else
    TACHO_TEST_FOR_ABORT( true, MSG_NOT_HAVE_PACKAGE("Teuchos") );
#endif
    }

    return 0;
  }
}

#endif
