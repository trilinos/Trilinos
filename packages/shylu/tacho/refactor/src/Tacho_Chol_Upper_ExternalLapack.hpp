#ifndef __TACHO_CHOL_UPPER_EXTERNAL_LAPACK_HPP__
#define __TACHO_CHOL_UPPER_EXTERNAL_LAPACK_HPP__

/// \file  Tacho_Chol_Upper_ExternalLapack.hpp
/// \brief LAPACK Cholesky factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)
#ifdef HAVE_SHYLUTACHO_TEUCHOS
#include "Teuchos_LAPACK.hpp"
#endif

namespace Tacho {
  /// LAPACK Chol
  /// ===========
  /// Properties:
  /// - Compile with Device (o),
  /// - Callable in KokkosFunctors (o)
  /// - For now, this is for HostSpace only.
  template<>
  template<typename PolicyType,
           typename MemberType,
           typename DenseExecViewTypeA>
  KOKKOS_INLINE_FUNCTION
  int
  Chol<Uplo::Upper,
       AlgoChol::ExternalLapack,Variant::One>
  ::invoke(PolicyType &policy,
           const MemberType &member,
           DenseExecViewTypeA &A) {
    // static_assert( Kokkos::Impl::is_same<
    //                typename DenseMatrixTypeA::space_type,
    //                Kokkos::Cuda
    //                >::value,
    //                "Cuda space is not available for calling external BLAS" );

    //typedef typename DenseExecViewTypeA::space_type   space_type;
    typedef typename DenseExecViewTypeA::ordinal_type ordinal_type;
    typedef typename DenseExecViewTypeA::value_type   value_type;

    int r_val = 0;      
    if (member.team_rank() == 0) {
#ifdef HAVE_SHYLUTACHO_TEUCHOS
      Teuchos::LAPACK<ordinal_type,value_type> lapack;

      lapack.POTRF('U',
                   A.NumRows(),
                   A.ValuePtr(), A.BaseObject().ColStride(),
                   &r_val);
#else
    TACHO_TEST_FOR_ABORT( true, MSG_NOT_HAVE_PACKAGE("Teuchos") );
#endif
    }

    return r_val;
  }
}

#endif
