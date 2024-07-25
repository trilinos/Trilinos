// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_TPETRA_EUCLIDEAN_SCALAR_PROD_DEF_HPP
#define THYRA_TPETRA_EUCLIDEAN_SCALAR_PROD_DEF_HPP

#include "Thyra_TpetraEuclideanScalarProd_decl.hpp"
#include "Thyra_TpetraMultiVector.hpp"
#include "Thyra_TpetraVector.hpp"


namespace Thyra {


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraEuclideanScalarProd<Scalar,LocalOrdinal,GlobalOrdinal,Node>::scalarProdsImpl(
  const MultiVectorBase<Scalar>& X,
  const MultiVectorBase<Scalar>& Y,
  const ArrayView<Scalar>& scalarProds_out
  ) const
{
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  Teuchos::RCP<const TMV> X_tpetra = this->getConstTpetraMultiVector(Teuchos::rcpFromRef(X));
  Teuchos::RCP<const TMV> Y_tpetra = this->getConstTpetraMultiVector(Teuchos::rcpFromRef(Y));

  if (nonnull(X_tpetra) && nonnull(Y_tpetra)) {
    // Which one do we want transposed?
    // Tpetra transposes the argument of dot.
    // Below is the order from TpetraVectorSpace::scalarProdsImpl,
    // so this would transpose Y. However, Thyra::dots (which calls
    // the RTOp) transposes the first argument, so scalarProdsImpl
    // in EuclideanScalarProd transposes X...
    X_tpetra->dot(*Y_tpetra, scalarProds_out);
  } else {
    EuclideanScalarProd<Scalar>::scalarProdsImpl(X, Y, scalarProds_out);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraEuclideanScalarProd<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getConstTpetraMultiVector(const RCP<const MultiVectorBase<Scalar> >& mv) const
{
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef Thyra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;

  RCP<const TMV> tmv = rcp_dynamic_cast<const TMV>(mv);
  if (nonnull(tmv)) {
    return tmv->getConstTpetraMultiVector();
  }

  RCP<const TV> tv = rcp_dynamic_cast<const TV>(mv);
  if (nonnull(tv)) {
    return tv->getConstTpetraVector();
  }

  return Teuchos::null;
}


} // end namespace Thyra


#endif  // THYRA_EUCLIDEAN_SCALAR_PROD_DEF_HPP
