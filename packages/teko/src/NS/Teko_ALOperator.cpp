// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * Author: Zhen Wang
 * Email: wangz@ornl.gov
 *        zhen.wang@alum.emory.edu
 */

#include "Teko_ALOperator.hpp"

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

#include "Teko_Utilities.hpp"
#include "Teko_ConfigDefs.hpp"

namespace Teko::NS {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

ALOperator::ALOperator(const std::vector<std::vector<GO> >& vars,
                       const Teuchos::RCP<Tpetra::Operator<ST, LO, GO, NT> >& content,
                       LinearOp pressureMassMatrix, double gamma, const std::string& label)
    : Teko::TpetraHelpers::BlockedTpetraOperator(vars, content, label),
      pressureMassMatrix_(pressureMassMatrix),
      gamma_(gamma) {
  checkDim(vars);
  BuildALOperator();
}

ALOperator::ALOperator(const std::vector<std::vector<GO> >& vars,
                       const Teuchos::RCP<Tpetra::Operator<ST, LO, GO, NT> >& content, double gamma,
                       const std::string& label)
    : Teko::TpetraHelpers::BlockedTpetraOperator(vars, content, label),
      pressureMassMatrix_(Teuchos::null),
      gamma_(gamma) {
  checkDim(vars);
  BuildALOperator();
}

void ALOperator::setPressureMassMatrix(LinearOp pressureMassMatrix) {
  if (pressureMassMatrix != Teuchos::null) pressureMassMatrix_ = pressureMassMatrix;
}

void ALOperator::setGamma(double gamma) {
  TEUCHOS_ASSERT(gamma > 0.0);
  gamma_ = gamma;
}

const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > ALOperator::GetBlock(int i,
                                                                                 int j) const {
  const Teuchos::RCP<Thyra::BlockedLinearOpBase<ST> > blkOp =
      Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<ST> >(blockedOperator_, true);

  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp =
      rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(blkOp->getBlock(i, j), true);

  return tOp->getConstTpetraOperator();
}

void ALOperator::checkDim(const std::vector<std::vector<GO> >& vars) {
  dim_ = static_cast<int>(vars.size()) - 1;
  TEUCHOS_ASSERT(dim_ == 2 || dim_ == 3);
}

void ALOperator::BuildALOperator() {
  TEUCHOS_ASSERT(blockedMapping_ != Teuchos::null);

  // Rebuild the original blocked operator from the stored full matrix
  const Teuchos::RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > crsContent =
      Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT> >(fullContent_, true);

  if (blockedOperator_ == Teuchos::null) {
    blockedOperator_ = blockedMapping_->buildBlockedThyraOp(crsContent, label_);
  } else {
    const Teuchos::RCP<Thyra::BlockedLinearOpBase<ST> > blkOp0 =
        Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<ST> >(blockedOperator_, true);
    blockedMapping_->rebuildBlockedThyraOp(crsContent, blkOp0);
  }

  // Extract original blocks
  const Teuchos::RCP<Thyra::BlockedLinearOpBase<ST> > blkOp =
      Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<ST> >(blockedOperator_, true);

  numBlockRows_ = blkOp->productRange()->numBlocks();

  Teuchos::RCP<const Thyra::LinearOpBase<ST> > blockedOpBlocks[4][4];
  for (int i = 0; i <= dim_; i++) {
    for (int j = 0; j <= dim_; j++) {
      blockedOpBlocks[i][j] = blkOp->getBlock(i, j);
    }
  }

  if (pressureMassMatrix_ != Teuchos::null) {
    invPressureMassMatrix_ = getInvDiagonalOp(pressureMassMatrix_);
  } else {
    std::cout << "Pressure mass matrix is null. Use identity." << std::endl;
    pressureMassMatrix_    = Thyra::identity<ST>(blockedOpBlocks[dim_][0]->range());
    invPressureMassMatrix_ = Thyra::identity<ST>(blockedOpBlocks[dim_][0]->range());
  }

  Teuchos::RCP<Thyra::DefaultBlockedLinearOp<ST> > alOperator = Thyra::defaultBlockedLinearOp<ST>();
  alOperator->beginBlockFill(dim_ + 1, dim_ + 1);

  // Velocity-velocity blocks: A_ij + gamma * B_i^T * W^{-1} * B_j
  for (int i = 0; i < dim_; i++) {
    for (int j = 0; j < dim_; j++) {
      alOperator->setNonconstBlock(
          i, j,
          Teuchos::rcp_const_cast<Thyra::LinearOpBase<ST> >(Thyra::add(
              blockedOpBlocks[i][j],
              Thyra::scale(gamma_, Thyra::multiply(blockedOpBlocks[i][dim_], invPressureMassMatrix_,
                                                   blockedOpBlocks[dim_][j])))));
    }
  }

  // Last row: [ B  -C ]
  for (int j = 0; j <= dim_; j++) {
    alOperator->setNonconstBlock(
        dim_, j, Teuchos::rcp_const_cast<Thyra::LinearOpBase<ST> >(blockedOpBlocks[dim_][j]));
  }

  // Last column: B_i^T - gamma * B_i^T * W^{-1} * C
  for (int i = 0; i < dim_; i++) {
    alOperator->setNonconstBlock(
        i, dim_,
        Teuchos::rcp_const_cast<Thyra::LinearOpBase<ST> >(Thyra::add(
            blockedOpBlocks[i][dim_],
            Thyra::scale(-gamma_, Thyra::multiply(blockedOpBlocks[i][dim_], invPressureMassMatrix_,
                                                  blockedOpBlocks[dim_][dim_])))));
  }

  alOperator->endBlockFill();
  alOperator_ = alOperator;

  // This is the operator used by apply()
  SetOperator(alOperator_, false);

  // Build operator for augmenting RHS:
  // [ I  0  gamma * B^T * W^{-1} ]
  // [ 0  I  gamma * B^T * W^{-1} ]
  // [ 0  0  I                     ]
  Teuchos::RCP<Thyra::DefaultBlockedLinearOp<ST> > alOpRhs = Thyra::defaultBlockedLinearOp<ST>();
  alOpRhs->beginBlockFill(dim_ + 1, dim_ + 1);

  for (int i = 0; i < dim_; i++) {
    alOpRhs->setNonconstBlock(i, i,
                              Teuchos::rcp_const_cast<Thyra::LinearOpBase<ST> >(
                                  Thyra::identity<ST>(blockedOpBlocks[i][i]->range())));
  }
  alOpRhs->setNonconstBlock(dim_, dim_,
                            Teuchos::rcp_const_cast<Thyra::LinearOpBase<ST> >(
                                Thyra::identity<ST>(blockedOpBlocks[dim_][dim_]->range())));

  for (int i = 0; i < dim_; i++) {
    alOpRhs->setNonconstBlock(
        i, dim_,
        Teuchos::rcp_const_cast<Thyra::LinearOpBase<ST> >(Thyra::scale(
            gamma_, Thyra::multiply(blockedOpBlocks[i][dim_], invPressureMassMatrix_))));
  }

  alOpRhs->endBlockFill();
  alOperatorRhs_ = alOpRhs;

  if (reorderManager_ != Teuchos::null) Reorder(*reorderManager_);
}

void ALOperator::augmentRHS(const Tpetra::MultiVector<ST, LO, GO, NT>& b,
                            Tpetra::MultiVector<ST, LO, GO, NT>& bAugmented) {
  Teuchos::RCP<const Teko::TpetraHelpers::MappingStrategy> mapping = this->getMapStrategy();

  Teuchos::RCP<Thyra::MultiVectorBase<ST> > bThyra =
      Thyra::createMembers(thyraOp_->range(), b.getNumVectors());

  Teuchos::RCP<Thyra::MultiVectorBase<ST> > bThyraAugmented =
      Thyra::createMembers(thyraOp_->range(), b.getNumVectors());

  mapping->copyTpetraIntoThyra(b, bThyra.ptr());
  alOperatorRhs_->apply(Thyra::NOTRANS, *bThyra, bThyraAugmented.ptr(), 1.0, 0.0);
  mapping->copyThyraIntoTpetra(bThyraAugmented, bAugmented);
}

}  // end namespace Teko::NS