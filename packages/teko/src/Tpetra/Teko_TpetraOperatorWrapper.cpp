// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_TpetraOperatorWrapper.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Tpetra_Vector.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

// #include "Thyra_LinearOperator.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_TpetraLinearOp.hpp"

#include "Teko_TpetraThyraConverter.hpp"
#include "Teuchos_Ptr.hpp"

namespace Teko {
namespace TpetraHelpers {

using namespace Teuchos;
using namespace Thyra;

DefaultMappingStrategy::DefaultMappingStrategy(
    const RCP<const Thyra::LinearOpBase<double> >& thyraOp,
    const Teuchos::Comm<Thyra::Ordinal>& comm) {
  RCP<Teuchos::Comm<Thyra::Ordinal> > newComm = comm.duplicate();

  // extract vector spaces from linear operator
  domainSpace_ = thyraOp->domain();
  rangeSpace_  = thyraOp->range();

  domainMap_ = Teko::TpetraHelpers::thyraVSToTpetraMap(*domainSpace_, newComm);
  rangeMap_  = Teko::TpetraHelpers::thyraVSToTpetraMap(*rangeSpace_, newComm);
}

void DefaultMappingStrategy::copyTpetraIntoThyra(
    const Tpetra::MultiVector<ST, LO, GO, NT>& x,
    const Ptr<Thyra::MultiVectorBase<ST> >& thyraVec) const {
  Teko::TpetraHelpers::blockTpetraToThyra(x, thyraVec);
}

void DefaultMappingStrategy::copyThyraIntoTpetra(
    const RCP<const Thyra::MultiVectorBase<ST> >& thyraVec,
    Tpetra::MultiVector<ST, LO, GO, NT>& v) const {
  Teko::TpetraHelpers::blockThyraToTpetra(thyraVec, v);
}

TpetraOperatorWrapper::TpetraOperatorWrapper() {
  useTranspose_ = false;
  mapStrategy_  = Teuchos::null;
  thyraOp_      = Teuchos::null;
  comm_         = Teuchos::null;
  label_        = Teuchos::null;
}

TpetraOperatorWrapper::TpetraOperatorWrapper(const RCP<const Thyra::LinearOpBase<ST> >& thyraOp) {
  SetOperator(thyraOp);
}

TpetraOperatorWrapper::TpetraOperatorWrapper(const RCP<const Thyra::LinearOpBase<ST> >& thyraOp,
                                             const RCP<const MappingStrategy>& mapStrategy)
    : mapStrategy_(mapStrategy) {
  SetOperator(thyraOp);
}

TpetraOperatorWrapper::TpetraOperatorWrapper(const RCP<const MappingStrategy>& mapStrategy)
    : mapStrategy_(mapStrategy) {
  useTranspose_ = false;
  thyraOp_      = Teuchos::null;
  comm_         = Teuchos::null;
  label_        = Teuchos::null;
}

void TpetraOperatorWrapper::SetOperator(const RCP<const LinearOpBase<ST> >& thyraOp,
                                        bool buildMap) {
  useTranspose_ = false;
  thyraOp_      = thyraOp;
  comm_         = getThyraComm(*thyraOp);
  label_        = thyraOp_->description();
  if (mapStrategy_ == Teuchos::null && buildMap)
    mapStrategy_ = Teuchos::rcp(new DefaultMappingStrategy(thyraOp, *comm_));
}

double TpetraOperatorWrapper::NormInf() const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                             "TpetraOperatorWrapper::NormInf not implemated");
  return 1.0;
}

void TpetraOperatorWrapper::apply(const Tpetra::MultiVector<ST, LO, GO, NT>& X,
                                  Tpetra::MultiVector<ST, LO, GO, NT>& Y,
                                  Teuchos::ETransp /* mode */, ST alpha, ST beta) const {
  if (!useTranspose_) {
    // allocate space for each vector
    RCP<Thyra::MultiVectorBase<ST> > tX;
    RCP<Thyra::MultiVectorBase<ST> > tY;

    tX = Thyra::createMembers(thyraOp_->domain(), X.getNumVectors());
    tY = Thyra::createMembers(thyraOp_->range(), X.getNumVectors());

    Thyra::assign(tX.ptr(), 0.0);
    Thyra::assign(tY.ptr(), 0.0);

    // copy Tpetra X into thyra X
    mapStrategy_->copyTpetraIntoThyra(X, tX.ptr());
    mapStrategy_->copyTpetraIntoThyra(
        Y, tY.ptr());  // if this matrix isn't block square, this probably won't work!

    // perform matrix vector multiplication
    thyraOp_->apply(Thyra::NOTRANS, *tX, tY.ptr(), alpha, beta);

    // copy thyra Y into Tpetra Y
    mapStrategy_->copyThyraIntoTpetra(tY, Y);
  } else {
    TEUCHOS_ASSERT(false);
  }
}

void TpetraOperatorWrapper::applyInverse(const Tpetra::MultiVector<ST, LO, GO, NT>& /* X */,
                                         Tpetra::MultiVector<ST, LO, GO, NT>& /* Y */,
                                         Teuchos::ETransp /* mode */, ST /* alpha */,
                                         ST /* beta */) const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                             "TpetraOperatorWrapper::applyInverse not implemented");
}

RCP<const Teuchos::Comm<Thyra::Ordinal> > TpetraOperatorWrapper::getThyraComm(
    const Thyra::LinearOpBase<ST>& inOp) const {
  RCP<const VectorSpaceBase<ST> > vs = inOp.domain();

  RCP<const SpmdVectorSpaceBase<ST> > spmd;
  RCP<const VectorSpaceBase<ST> > current = vs;
  while (current != Teuchos::null) {
    // try to cast to a product vector space first
    RCP<const ProductVectorSpaceBase<ST> > prod =
        rcp_dynamic_cast<const ProductVectorSpaceBase<ST> >(current);

    // figure out what type it is
    if (prod == Teuchos::null) {
      // hopfully this is a SPMD vector space
      spmd = rcp_dynamic_cast<const SpmdVectorSpaceBase<ST> >(current);

      break;
    } else  // get first convenient vector space
      current = prod->getBlock(0);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(spmd == Teuchos::null, std::runtime_error,
                             "TpetraOperatorWrapper requires std::vector space "
                             "blocks to be SPMD std::vector spaces");

  return spmd->getComm();
}

int TpetraOperatorWrapper::GetBlockRowCount() {
  const RCP<const Thyra::BlockedLinearOpBase<ST> > blkOp =
      Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<ST> >(getThyraOp());

  return blkOp->productRange()->numBlocks();
}

int TpetraOperatorWrapper::GetBlockColCount() {
  const RCP<const Thyra::BlockedLinearOpBase<ST> > blkOp =
      Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<ST> >(getThyraOp());

  return blkOp->productDomain()->numBlocks();
}

Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > TpetraOperatorWrapper::GetBlock(int i,
                                                                                      int j) const {
  const RCP<const Thyra::BlockedLinearOpBase<ST> > blkOp =
      Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<ST> >(getThyraOp());

  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp =
      rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(blkOp->getBlock(i, j));

  return tOp->getConstTpetraOperator();
}

}  // namespace TpetraHelpers
}  // namespace Teko
