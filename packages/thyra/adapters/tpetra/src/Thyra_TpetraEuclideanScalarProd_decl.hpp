// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_TPETRA_EUCLIDEAN_SCALAR_PROD_DECL_HPP
#define THYRA_TPETRA_EUCLIDEAN_SCALAR_PROD_DECL_HPP


#include "Thyra_EuclideanScalarProd.hpp"
#include "Tpetra_MultiVector.hpp"


namespace Thyra {

/** \brief Extends concrete implementation of a Euclidean scalar product for
 * specifically Tpetra vectors/multivectors.
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TpetraEuclideanScalarProd : public EuclideanScalarProd<Scalar> {
protected:
  
  /** @name Overridden from EuclideanScalarProd */
  //@{

  /** \brief If X and Y are both Tpetra wrappers, computes the pair-wise
   * scalar products directly with Tpetra calls. Otherwise, this defers to
   * the base class implementaiton, which computes the result with an RTOp.
   */
  virtual void scalarProdsImpl(
    const MultiVectorBase<Scalar>& X,
    const MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar>& scalarProds
    ) const;

  //@}

private:

  /** /brief . */
  Teuchos::RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraMultiVector(const RCP<const MultiVectorBase<Scalar> >& mv) const;

};


/** \brief Nonmember constructor for TpetraEuclideanScalarProd.
 *
 * \relates TpetraEuclideanScalarProd
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
inline
RCP<const TpetraEuclideanScalarProd<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
tpetraEuclideanScalarProd()
{
  return  Teuchos::rcp(new TpetraEuclideanScalarProd<Scalar,LocalOrdinal,GlobalOrdinal,Node>);
}


} // end namespace Thyra


#endif  // THYRA_TPETRA_EUCLIDEAN_SCALAR_PROD_DECL_HPP
