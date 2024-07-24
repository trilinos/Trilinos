// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SILLY_MODIFIED_GRAM_SHMIDT_HPP
#define THYRA_SILLY_MODIFIED_GRAM_SHMIDT_HPP

#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"

namespace Thyra {

/** \brief Silly little implementation of the modified Gram-Schmidt algorithm
 * to compute a QR factorization A=Q*R of a multi-vector A.
 *
 * \param V [in/out] On input, V contains the columns of A to compute the QR
 * factorization.  On output, V contains the columns of Q.
 *
 * \param R [out] On output, contains the upper triangular matrix R.
 *
 * \ingroup Thyra_Op_Vec_examples_cg_grp
 */
template<class Scalar>
void sillyModifiedGramSchmidt(
  const Ptr<MultiVectorBase<Scalar> > &V,
  const Ptr<RCP<MultiVectorBase<Scalar> > > &R_out
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST; using Teuchos::as;
  const int n = V->domain()->dim();
  *R_out = createMembers(V->domain(), V->domain());
  DetachedMultiVectorView<Scalar> R(*(*R_out));
  for (int k = 0; k < n; ++k) {
    R(k,k) = norm(*V->col(k));
    Vt_S(V->col(k).ptr(), ST::one()/R(k,k));
    for (int j = k+1; j < n; ++j) {
      R(k,j) = scalarProd(*V->col(k), *V->col(j));
      Vp_StV(V->col(j).ptr(), -R(k,j), *V->col(k));
    }
  }
} // end sillyModifiedGramSchmidt

} // namespace Thyra

#endif // THYRA_SILLY_MODIFIED_GRAM_SHMIDT_HPP
