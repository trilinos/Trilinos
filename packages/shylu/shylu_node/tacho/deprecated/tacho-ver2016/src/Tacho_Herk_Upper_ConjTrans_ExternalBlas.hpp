#ifndef __TACHO_HERK_UPPER_CONJTRANS_EXTERNAL_BLAS_HPP__
#define __TACHO_HERK_UPPER_CONJTRANS_EXTERNAL_BLAS_HPP__


/// \file  Tacho_Herk_Upper_ConjTrans_ExternalBlas.hpp
/// \brief BLAS hermitian rank-k update
/// \author Kyungjoo Kim (kyukim@sandia.gov)
#ifdef HAVE_SHYLU_NODETACHO_TEUCHOS
#include "Teuchos_BLAS.hpp"
#endif

#include "Tacho_DenseFlopCount.hpp"

namespace Tacho {

  /// BLAS Herk
  /// =========
  /// Properties:
  /// - Compile with Device (o),
  /// - Callable in KokkosFunctors (o)
  /// - For now, this is for HostSpace only.
  template<>
  template<typename ScalarType,
           typename DenseExecViewTypeA,
           typename DenseExecViewTypeC>
  inline
  Stat
  Herk<Uplo::Upper,Trans::ConjTranspose,
       AlgoHerk::ExternalBlas,Variant::One>
  ::stat(const ScalarType alpha,
         DenseExecViewTypeA &A,
         const ScalarType beta,
         DenseExecViewTypeC &C) {
    Stat r_val;

    typedef typename DenseExecViewTypeA::ordinal_type ordinal_type;
    
    const ordinal_type n = C.NumRows();
    const ordinal_type k = A.NumRows();
    r_val.flop = DenseFlopCount<typename DenseExecViewTypeA::value_type>::Syrk(k, n);
    
    return r_val;
  }

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
           MemberType &member,
           const ScalarType alpha,
           DenseExecViewTypeA &A,
           const ScalarType beta,
           DenseExecViewTypeC &C) {
    // static_assert( std::is_same<
    //                typename DenseMatrixTypeA::space_type,
    //                typename DenseMatrixTypeC::space_type
    //                >::value,
    //                "Space type of input matrices does not match" );
    if (member.team_rank() == 0) {
#if                                                     \
  defined( HAVE_SHYLU_NODETACHO_TEUCHOS ) &&                 \
  defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      typedef typename DenseExecViewTypeA::ordinal_type ordinal_type;
      typedef typename DenseExecViewTypeA::value_type   value_type;
      
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
