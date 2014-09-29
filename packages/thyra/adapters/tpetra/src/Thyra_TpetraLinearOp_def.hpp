// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_TPETRA_LINEAR_OP_HPP
#define THYRA_TPETRA_LINEAR_OP_HPP

#include "Thyra_TpetraLinearOp_decl.hpp"
#include "Thyra_TpetraVectorSpace.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"

#include "Tpetra_CrsMatrix.hpp"

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


template <class Scalar>
inline
Teuchos::ETransp
convertConjNoTransToTeuchosTransMode()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      Teuchos::ScalarTraits<Scalar>::isComplex,
      Exceptions::OpNotSupported,
      "For complex scalars such as " + Teuchos::TypeNameTraits<Scalar>::name() +
      ", Tpetra does not support conjugation without transposition."
      );
  return Teuchos::NO_TRANS; // For non-complex scalars, CONJ is equivalent to NOTRANS.
}


template <class Scalar>
inline
Teuchos::ETransp
convertToTeuchosTransMode(const Thyra::EOpTransp transp)
{
  switch (transp) {
    case NOTRANS:   return Teuchos::NO_TRANS;
    case CONJ:      return convertConjNoTransToTeuchosTransMode<Scalar>();
    case TRANS:     return Teuchos::TRANS;
    case CONJTRANS: return Teuchos::CONJ_TRANS;
  }

  // Should not escape the switch
  TEUCHOS_TEST_FOR_EXCEPT(true);
}


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
  TEUCHOS_TEST_FOR_EXCEPT(true);
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

  if (M_trans == CONJ) {
    // For non-complex scalars, CONJ is always supported since it is equivalent to NO_TRANS.
    // For complex scalars, Tpetra does not support conjugation without transposition.
    return !Teuchos::ScalarTraits<Scalar>::isComplex;
  }

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

  const Teuchos::ETransp tTransp = convertToTeuchosTransMode<Scalar>(M_trans);

  // Apply the operator

  tpetraOperator_->apply(*tX, *tY, tTransp, alpha, beta);

}

// Protected member functions overridden from ScaledLinearOpBase


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::supportsScaleLeftImpl() const
{
  return true;
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::supportsScaleRightImpl() const
{
  return true;
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void 
TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
scaleLeftImpl(const VectorBase<Scalar> &row_scaling_in)
{
  using Teuchos::rcpFromRef;

  const RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > row_scaling = 
    TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getConstTpetraVector(rcpFromRef(row_scaling_in));

  const RCP<typename Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rowMatrix = 
    Teuchos::rcp_dynamic_cast<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tpetraOperator_.getNonconstObj(),true);

  rowMatrix->leftScale(*row_scaling);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void 
TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
scaleRightImpl(const VectorBase<Scalar> &col_scaling_in)
{
  using Teuchos::rcpFromRef;

  const RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > col_scaling = 
    TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getConstTpetraVector(rcpFromRef(col_scaling_in));

  const RCP<typename Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rowMatrix = 
    Teuchos::rcp_dynamic_cast<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tpetraOperator_.getNonconstObj(),true);

  rowMatrix->rightScale(*col_scaling);
}

// Protected member functions overridden from RowStatLinearOpBase


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
rowStatIsSupportedImpl(
  const RowStatLinearOpBaseUtils::ERowStat rowStat) const
{
  if (is_null(tpetraOperator_))
    return false;

  switch (rowStat) {
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM:
      return true;
    default:
      return false;
  }

  return false; // Will never be called!
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getRowStatImpl(
  const RowStatLinearOpBaseUtils::ERowStat rowStat,
  const Ptr<VectorBase<Scalar> > &rowStatVec_in
  ) const
{
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>
    TpetraVector_t;
  
  if ( (rowStat == RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM) ||
       (rowStat == RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM) ) {
    
    TEUCHOS_ASSERT(nonnull(tpetraOperator_));
    TEUCHOS_ASSERT(nonnull(rowStatVec_in));

    // Currently we only support the case of row sums for a concrete
    // Tpetra::CrsMatrix where (1) the entire row is stored on a
    // single process and (2) that the domain map, the range map and
    // the row map are the SAME.  These checks enforce that.  Later on
    // we hope to add complete support for any mapping to the concrete
    // tpetra matrix types.
    
    const RCP<TpetraVector_t> tRowSumVec =
      TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getTpetraVector(rcpFromPtr(rowStatVec_in));
    
    const RCP<const typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tCrsMatrix = 
      Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tpetraOperator_.getConstObj(),true);
    
    TEUCHOS_ASSERT(tCrsMatrix->getRowMap()->isSameAs(*tCrsMatrix->getDomainMap()));
    TEUCHOS_ASSERT(tCrsMatrix->getRowMap()->isSameAs(*tCrsMatrix->getRangeMap()));
    TEUCHOS_ASSERT(tCrsMatrix->getRowMap()->isSameAs(*tRowSumVec->getMap()));

    size_t numMyRows = tCrsMatrix->getNodeNumRows();
    
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> values;

    for (size_t row=0; row < numMyRows; ++row) {
      Scalar sum = Scalar(0.0);
      tCrsMatrix->getLocalRowView(row,indices,values);

      for (int col = 0; col < values.size(); ++col)
	sum += std::abs(values[col]);
      
      if (rowStat == RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM) {
	if (sum < Teuchos::ScalarTraits<Scalar>::sfmin()) {
	  TEUCHOS_TEST_FOR_EXCEPTION(sum == Teuchos::ScalarTraits<Scalar>::zero(), std::runtime_error,
				     "Error - Thyra::TpetraLinearOp::getRowStatImpl() - Inverse row sum "
				     << "requested for a matrix where one of the rows has a zero row sum!");
	
	  sum = Scalar(1.0) / Teuchos::ScalarTraits<Scalar>::sfmin();
	}
	else 
	  sum = Scalar(1.0) / sum;
      }
      
      tRowSumVec->replaceLocalValue(row,sum);
    }

  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
			       "Error - Thyra::TpetraLinearOp::getRowStatImpl() - Colum sum support not implemented!");
  }
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
