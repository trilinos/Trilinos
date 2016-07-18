#ifndef __TACHO_HERK_UPPER_CONJTRANS_EXTERNAL_BLAS_HPP__
#define __TACHO_HERK_UPPER_CONJTRANS_EXTERNAL_BLAS_HPP__


/// \file  Tacho_Herk_Upper_ConjTrans_ExternalBlas.hpp
/// \brief BLAS hermitian rank-k update
/// \author Kyungjoo Kim (kyukim@sandia.gov)
#ifdef HAVE_SHYLUTACHO_TEUCHOS
#include "Teuchos_BLAS.hpp"
#endif

namespace Tacho {
  /// BLAS Herk
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
           typename DenseExecViewTypeC>
  KOKKOS_INLINE_FUNCTION
  int
  Herk<Uplo::Upper,Trans::ConjTranspose,
       AlgoHerk::ExternalBlas,Variant::One>
  ::invoke(PolicyType &policy,
           const MemberType &member,
           const ScalarType alpha,
           DenseExecViewTypeA &A,
           const ScalarType beta,
           DenseExecViewTypeC &C) {
    // static_assert( Kokkos::Impl::is_same<
    //                typename DenseMatrixTypeA::space_type,
    //                Kokkos::Cuda
    //                >::value,
    //                "Cuda space is not available for calling external BLAS" );

    // static_assert( Kokkos::Impl::is_same<
    //                typename DenseMatrixTypeA::space_type,
    //                typename DenseMatrixTypeC::space_type
    //                >::value,
    //                "Space type of input matrices does not match" );

    //typedef typename DenseExecViewTypeA::space_type   space_type;
    typedef typename DenseExecViewTypeA::ordinal_type ordinal_type;
    typedef typename DenseExecViewTypeA::value_type   value_type;

    if (member.team_rank() == 0) {
#ifdef HAVE_SHYLUTACHO_TEUCHOS
      Teuchos::BLAS<ordinal_type,value_type> blas;

      const ordinal_type n = C.NumRows();
      const ordinal_type k = A.NumRows();

      if (n > 0 && k > 0) 
        blas.HERK(Teuchos::UPPER_TRI, Teuchos::CONJ_TRANS,
                  n, k,
                  alpha,
                  A.ValuePtr(), A.BaseObject().ColStride(),
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
