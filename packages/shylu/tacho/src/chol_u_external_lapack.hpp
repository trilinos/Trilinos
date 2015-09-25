#pragma once
#ifndef __CHOL_U_EXTERNAL_LAPACK_HPP__
#define __CHOL_U_EXTERNAL_LAPACK_HPP__

/// \file chol_u_external_lapack.hpp
/// \brief BLAS Chloesky factorization.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Teuchos_LAPACK.hpp"

namespace Tacho {

  using namespace std;

  template<>
  template<typename ParallelForType,
           typename DenseExecViewType>
  KOKKOS_INLINE_FUNCTION
  int
  Chol<Uplo::Upper,AlgoChol::ExternalLapack>
  ::invoke(typename DenseExecViewType::policy_type &policy,
           const typename DenseExecViewType::policy_type::member_type &member,
           DenseExecViewType &A) {
    typedef typename DenseExecViewType::ordinal_type ordinal_type;
    typedef typename DenseExecViewType::value_type   value_type;

    if (member.team_rank() == 0) {
      const ordinal_type n = A.NumRows();
      
      int r_val = 0;
      Teuchos::LAPACK<ordinal_type,value_type>::POTRF('U', 
                                                      n, 
                                                      A.ValuePtr(), A.BaseObject->ColStride(),
                                                      &r_val);
    }
    return 0;
  }

}

#endif
