// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_BlockedTpetraOperator.hpp"
#include "Teko_TpetraBlockedMappingStrategy.hpp"
#include "Teko_TpetraReorderedMappingStrategy.hpp"

#include "Teuchos_VerboseObject.hpp"

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

#include "MatrixMarket_Tpetra.hpp"

#include "Teko_Utilities.hpp"

namespace Teko {
namespace TpetraHelpers {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

BlockedTpetraOperator::BlockedTpetraOperator(
    const std::vector<std::vector<GO> >& vars,
    const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& content, const std::string& label)
    : Teko::TpetraHelpers::TpetraOperatorWrapper(), label_(label) {
  SetContent(vars, content);
}

void BlockedTpetraOperator::SetContent(
    const std::vector<std::vector<GO> >& vars,
    const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& content) {
  fullContent_    = content;
  blockedMapping_ = rcp(new TpetraBlockedMappingStrategy(vars, fullContent_->getDomainMap(),
                                                         *fullContent_->getDomainMap()->getComm()));
  SetMapStrategy(blockedMapping_);

  // build thyra operator
  BuildBlockedOperator();
}

void BlockedTpetraOperator::BuildBlockedOperator() {
  TEUCHOS_ASSERT(blockedMapping_ != Teuchos::null);

  // get a CRS matrix
  const RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > crsContent =
      rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT> >(fullContent_);

  // ask the strategy to build the Thyra operator for you
  if (blockedOperator_ == Teuchos::null) {
    blockedOperator_ = blockedMapping_->buildBlockedThyraOp(crsContent, label_);
  } else {
    const RCP<Thyra::BlockedLinearOpBase<ST> > blkOp =
        rcp_dynamic_cast<Thyra::BlockedLinearOpBase<ST> >(blockedOperator_, true);
    blockedMapping_->rebuildBlockedThyraOp(crsContent, blkOp);
  }

  // set whatever is returned
  SetOperator(blockedOperator_, false);

  // reorder if neccessary
  if (reorderManager_ != Teuchos::null) Reorder(*reorderManager_);
}

const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > BlockedTpetraOperator::GetBlock(
    int i, int j) const {
  const RCP<const Thyra::BlockedLinearOpBase<ST> > blkOp =
      Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<ST> >(getThyraOp());

  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp =
      rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(blkOp->getBlock(i, j), true);
  return tOp->getConstTpetraOperator();
}

/** Use a reorder manager to block this operator as desired.
 * Multiple calls to the function reorder only the underlying object.
 */
void BlockedTpetraOperator::Reorder(const BlockReorderManager& brm) {
  reorderManager_ = rcp(new BlockReorderManager(brm));

  // build reordered objects
  RCP<const MappingStrategy> reorderMapping =
      rcp(new TpetraReorderedMappingStrategy(*reorderManager_, blockedMapping_));
  RCP<const Thyra::BlockedLinearOpBase<ST> > blockOp =
      rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<ST> >(blockedOperator_);

  RCP<const Thyra::LinearOpBase<ST> > A = buildReorderedLinearOp(*reorderManager_, blockOp);

  // set them as working values
  SetMapStrategy(reorderMapping);
  SetOperator(A, false);
}

//! Remove any reordering on this object
void BlockedTpetraOperator::RemoveReording() {
  SetMapStrategy(blockedMapping_);
  SetOperator(blockedOperator_, false);
  reorderManager_ = Teuchos::null;
}

/** Write out this operator to matrix market files
 */
void BlockedTpetraOperator::WriteBlocks(const std::string& prefix) const {
  RCP<Thyra::PhysicallyBlockedLinearOpBase<ST> > blockOp =
      rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<ST> >(blockedOperator_);

  // get size of blocked block operator
  int rows = Teko::blockRowCount(blockOp);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < rows; j++) {
      // get the row matrix object
      RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp =
          rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(blockOp->getBlock(i, j));
      RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > mat =
          Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT> >(
              tOp->getConstTpetraOperator());

      // write to file
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<ST, LO, GO, NT> >::writeSparseFile(
          formatBlockName(prefix, i, j, rows).c_str(), mat);
    }
  }
}

bool BlockedTpetraOperator::testAgainstFullOperator(int count, ST tol) const {
  Tpetra::Vector<ST, LO, GO, NT> xf(getRangeMap());
  Tpetra::Vector<ST, LO, GO, NT> xs(getRangeMap());
  Tpetra::Vector<ST, LO, GO, NT> y(getDomainMap());

  // test operator many times
  bool result = true;
  ST diffNorm = 0.0, trueNorm = 0.0;
  for (int i = 0; i < count; i++) {
    xf.putScalar(0.0);
    xs.putScalar(0.0);
    y.randomize();

    // apply operator
    apply(y, xs);                // xs = A*y
    fullContent_->apply(y, xf);  // xf = A*y

    // compute norms
    xs.update(-1.0, xf, 1.0);
    diffNorm = xs.norm2();
    trueNorm = xf.norm2();

    // check result
    result &= (diffNorm / trueNorm < tol);
  }

  return result;
}

}  // end namespace TpetraHelpers
}  // end namespace Teko
