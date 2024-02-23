// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef THYRA_XPETRA_LINEAR_OP_HPP
#define THYRA_XPETRA_LINEAR_OP_HPP

#include "Thyra_XpetraLinearOp_decl.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"

#include "MueLu_XpetraOperator_decl.hpp"
#include "Xpetra_MapExtractor.hpp"

namespace Thyra {

// Constructors/initializers

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>::XpetraLinearOp() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>::initialize(
    const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
    const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
    const RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &xpetraOperator) {
  initializeImpl(rangeSpace, domainSpace, xpetraOperator);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>::constInitialize(
    const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
    const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
    const RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &xpetraOperator) {
  initializeImpl(rangeSpace, domainSpace, xpetraOperator);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getXpetraOperator() {
  return xpetraOperator_.getNonconstObj();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getConstXpetraOperator() const {
  return xpetraOperator_;
}

// Public Overridden functions from LinearOpBase

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Thyra::VectorSpaceBase<Scalar> >
XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>::range() const {
  return rangeSpace_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Thyra::VectorSpaceBase<Scalar> >
XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>::domain() const {
  return domainSpace_;
}

// Protected Overridden functions from LinearOpBase

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>::opSupportedImpl(
    Thyra::EOpTransp M_trans) const {
  if (is_null(xpetraOperator_))
    return false;

  if (M_trans == NOTRANS)
    return true;

  if (M_trans == CONJ) {
    // For non-complex scalars, CONJ is always supported since it is equivalent to NO_TRANS.
    // For complex scalars, Xpetra does not support conjugation without transposition.
    return !Teuchos::ScalarTraits<Scalar>::isComplex;
  }

  return xpetraOperator_->hasTransposeApply();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>::applyImpl(
    const Thyra::EOpTransp M_trans,
    const Thyra::MultiVectorBase<Scalar> &X_in,
    const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> > &Y_inout,
    const Scalar alpha,
    const Scalar beta) const {
  using Teuchos::rcpFromPtr;
  using Teuchos::rcpFromRef;

  TEUCHOS_TEST_FOR_EXCEPTION(getConstXpetraOperator() == Teuchos::null, MueLu::Exceptions::RuntimeError, "XpetraLinearOp::applyImpl: internal Xpetra::Operator is null.");
  RCP<const Teuchos::Comm<int> > comm = getConstXpetraOperator()->getRangeMap()->getComm();

  const RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tX_in =
      Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toXpetra(rcpFromRef(X_in), comm);
  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tY_inout =
      Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toXpetra(rcpFromPtr(Y_inout), comm);
  Teuchos::ETransp transp;
  switch (M_trans) {
    case NOTRANS: transp = Teuchos::NO_TRANS; break;
    case TRANS: transp = Teuchos::TRANS; break;
    case CONJTRANS: transp = Teuchos::CONJ_TRANS; break;
    default: TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::NotImplemented, "Thyra::XpetraLinearOp::apply. Unknown value for M_trans. Only NOTRANS, TRANS and CONJTRANS are supported.");
  }

  xpetraOperator_->apply(*tX_in, *tY_inout, transp, alpha, beta);

  // check whether Y is a product vector
  RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rgMapExtractor = Teuchos::null;
  Teuchos::Ptr<Thyra::ProductMultiVectorBase<Scalar> > prodY_inout =
      Teuchos::ptr_dynamic_cast<Thyra::ProductMultiVectorBase<Scalar> >(Y_inout);
  if (prodY_inout != Teuchos::null) {
    // If Y is a product vector we split up the data from tY and merge them
    // into the product vector. The necessary Xpetra::MapExtractor is extracted
    // from the fine level operator (not this!)

    // get underlying fine level operator (BlockedCrsMatrix)
    // to extract the range MapExtractor
    RCP<MueLu::XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > mueXop =
        Teuchos::rcp_dynamic_cast<MueLu::XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(xpetraOperator_.getNonconstObj());

    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A =
        mueXop->GetHierarchy()->GetLevel(0)->template Get<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >("A");
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(A));

    RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > bA =
        Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(A);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(bA));

    rgMapExtractor = bA->getRangeMapExtractor();
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgMapExtractor));
  }
}

// private

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template <class XpetraOperator_t>
void XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>::initializeImpl(
    const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
    const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
    const RCP<XpetraOperator_t> &xpetraOperator) {
#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT(nonnull(rangeSpace));
  TEUCHOS_ASSERT(nonnull(domainSpace));
  TEUCHOS_ASSERT(nonnull(xpetraOperator));
#endif
  rangeSpace_     = rangeSpace;
  domainSpace_    = domainSpace;
  xpetraOperator_ = xpetraOperator;
}

}  // namespace Thyra

#endif  // THYRA_XPETRA_LINEAR_OP_HPP
