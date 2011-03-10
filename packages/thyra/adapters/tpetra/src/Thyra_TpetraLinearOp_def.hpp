// @HEADER
// ***********************************************************************
// 
//    Thar: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_TPETRA_LINEAR_OP_HPP
#define THYRA_TPETRA_LINEAR_OP_HPP

#include "Thyra_TpetraLinearOp_decl.hpp"
#include "Thyra_TpetraVectorSpace.hpp"

#ifdef HAVE_THYRA_TPETRA_EPETRA
#  include "Thyra_EpetraThyraWrappers.hpp"
#endif

namespace Thyra {


#ifdef HAVE_THYRA_TPETRA_EPETRA

// Utilites


/** \brief Default class returns null. */
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
class GetTpetraEpetraRowMatrixWrapper {
public:
  template<class TpetraMatrixType>
  static
  RCP<Tpetra::EpetraRowMatrix<TpetraMatrixType> >
  get(const RCP<TpetraMatrixType> &tpetraMatrix)
    {
      return Teuchos::null;
    }
};


// NOTE: We could support other ordinal types, but we have to
// specialize the EpetraRowMatrix
template<>
class GetTpetraEpetraRowMatrixWrapper<double, int, int> {
public:
  template<class TpetraMatrixType>
  static
  RCP<Tpetra::EpetraRowMatrix<TpetraMatrixType> >
  get(const RCP<TpetraMatrixType> &tpetraMatrix)
    {
      return Teuchos::rcp(
        new Tpetra::EpetraRowMatrix<TpetraMatrixType>(tpetraMatrix,
          *get_Epetra_Comm(
            *convertTpetraToThyraComm(tpetraMatrix->getRowMap()->getComm())
            )
          )
        );
    }
};


#endif // HAVE_THYRA_TPETRA_EPETRA


// Constructors/initializers


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraLinearOp()
{}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::initialize(
  const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
  const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
  const RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraOperator
  )
{
  initializeImpl(rangeSpace, domainSpace, tpetraOperator);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::constInitialize(
  const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
  const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
  const RCP<const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraOperator
  )
{
  initializeImpl(rangeSpace, domainSpace, tpetraOperator);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getTpetraOperator()
{
  return tpetraOperator_.getNonconstObj();
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getConstTpetraOperator() const
{
  return tpetraOperator_;
}


// Public Overridden functions from LinearOpBase


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Thyra::VectorSpaceBase<Scalar> >
TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::range() const
{
  return rangeSpace_;
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Thyra::VectorSpaceBase<Scalar> >
TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::domain() const
{
  return domainSpace_;
}


// Overridden from EpetraLinearOpBase


#ifdef HAVE_THYRA_TPETRA_EPETRA


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNonconstEpetraOpView(
  const Ptr<RCP<Epetra_Operator> > &epetraOp,
  const Ptr<EOpTransp> &epetraOpTransp,
  const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
  const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport
  )
{
  TEST_FOR_EXCEPT(true);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getEpetraOpView(
  const Ptr<RCP<const Epetra_Operator> > &epetraOp,
  const Ptr<EOpTransp> &epetraOpTransp,
  const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
  const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport
  ) const
{
  using Teuchos::rcp_dynamic_cast;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraRowMatrix_t;
  if (nonnull(tpetraOperator_)) {
    if (is_null(epetraOp_)) {
      epetraOp_ = GetTpetraEpetraRowMatrixWrapper<Scalar,LocalOrdinal,GlobalOrdinal>::get(
        rcp_dynamic_cast<const TpetraRowMatrix_t>(tpetraOperator_.getConstObj(), true));
    }
    *epetraOp = epetraOp_;
    *epetraOpTransp = NOTRANS;
    *epetraOpApplyAs = EPETRA_OP_APPLY_APPLY;
    *epetraOpAdjointSupport = ( tpetraOperator_->hasTransposeApply()
      ? EPETRA_OP_ADJOINT_SUPPORTED : EPETRA_OP_ADJOINT_UNSUPPORTED );
  }
  else {
    *epetraOp = Teuchos::null;
  }
}


#endif // HAVE_THYRA_TPETRA_EPETRA


// Protected Overridden functions from LinearOpBase


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::opSupportedImpl(
  Thyra::EOpTransp M_trans) const
{
  if (is_null(tpetraOperator_))
    return false;
  if (M_trans == NOTRANS)
    return true;
  return tpetraOperator_->hasTransposeApply();
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyImpl(
  const Thyra::EOpTransp M_trans,
  const Thyra::MultiVectorBase<Scalar> &X_in,
  const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> > &Y_inout,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  using Teuchos::rcpFromRef;
  using Teuchos::rcpFromPtr;
  typedef TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node>
    ConverterT;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>
    TpetraMultiVector_t;

  // Get Tpetra::MultiVector objects for X and Y

  const RCP<const TpetraMultiVector_t> tX =
    ConverterT::getConstTpetraMultiVector(rcpFromRef(X_in));

  const RCP<TpetraMultiVector_t> tY =
    ConverterT::getTpetraMultiVector(rcpFromPtr(Y_inout));

  const Teuchos::ETransp tTransp =
    ( M_trans == NOTRANS ? Teuchos::NO_TRANS : Teuchos::CONJ_TRANS );

  // Apply the operator

  tpetraOperator_->apply(*tX, *tY, tTransp, alpha, beta);

}


// private


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template<class TpetraOperator_t>
void TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::initializeImpl(
  const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
  const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
  const RCP<TpetraOperator_t> &tpetraOperator
  )
{
#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT(nonnull(rangeSpace));
  TEUCHOS_ASSERT(nonnull(domainSpace));
  TEUCHOS_ASSERT(nonnull(tpetraOperator));
  // ToDo: Assert that spaces are comparible with tpetraOperator
#endif  
  rangeSpace_ = rangeSpace;
  domainSpace_ = domainSpace;
  tpetraOperator_ = tpetraOperator;
}


} // namespace Thyra


#endif	// THYRA_TPETRA_LINEAR_OP_HPP
