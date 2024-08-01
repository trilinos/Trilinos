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

#include "Teko_BlockedEpetraOperator.hpp"
#include "Teko_BlockedMappingStrategy.hpp"
#include "Teko_ReorderedMappingStrategy.hpp"

#include "Teuchos_VerboseObject.hpp"

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"

#include "Teko_Utilities.hpp"

#include "Teko_ALOperator.hpp"

namespace Teko {

namespace NS {

ALOperator::ALOperator(const std::vector<std::vector<int> >& vars,
                       const Teuchos::RCP<Epetra_Operator>& content, LinearOp pressureMassMatrix,
                       double gamma, const std::string& label)
    : Teko::Epetra::BlockedEpetraOperator(vars, content, label),
      pressureMassMatrix_(pressureMassMatrix),
      gamma_(gamma) {
  checkDim(vars);
  SetContent(vars, content);
  BuildALOperator();
}

ALOperator::ALOperator(const std::vector<std::vector<int> >& vars,
                       const Teuchos::RCP<Epetra_Operator>& content, double gamma,
                       const std::string& label)
    : Teko::Epetra::BlockedEpetraOperator(vars, content, label),
      pressureMassMatrix_(Teuchos::null),
      gamma_(gamma) {
  checkDim(vars);
  SetContent(vars, content);
  BuildALOperator();
}

void ALOperator::setPressureMassMatrix(LinearOp pressureMassMatrix) {
  if (pressureMassMatrix != Teuchos::null) {
    pressureMassMatrix_ = pressureMassMatrix;
  }
}

void ALOperator::setGamma(double gamma) {
  TEUCHOS_ASSERT(gamma > 0.0);
  gamma_ = gamma;
}

const Teuchos::RCP<const Epetra_Operator> ALOperator::GetBlock(int i, int j) const {
  const Teuchos::RCP<Thyra::BlockedLinearOpBase<double> > blkOp =
      Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(blockedOperator_, true);

  return Thyra::get_Epetra_Operator(*blkOp->getBlock(i, j));
}

void ALOperator::checkDim(const std::vector<std::vector<int> >& vars) {
  dim_ = vars.size() - 1;
  TEUCHOS_ASSERT(dim_ == 2 || dim_ == 3);
}

void ALOperator::BuildALOperator() {
  TEUCHOS_ASSERT(blockedMapping_ != Teuchos::null);

  // Get an Epetra_CrsMatrix.
  const Teuchos::RCP<const Epetra_CrsMatrix> crsContent =
      Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(fullContent_);

  // Ask the strategy to build the Thyra operator for you.
  if (blockedOperator_ == Teuchos::null) {
    blockedOperator_ = blockedMapping_->buildBlockedThyraOp(crsContent, label_);
  } else {
    const Teuchos::RCP<Thyra::BlockedLinearOpBase<double> > blkOp =
        Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(blockedOperator_, true);
    blockedMapping_->rebuildBlockedThyraOp(crsContent, blkOp);
  }

  // Extract blocks.
  const Teuchos::RCP<Thyra::BlockedLinearOpBase<double> > blkOp =
      Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(blockedOperator_, true);
  numBlockRows_ = blkOp->productRange()->numBlocks();
  Teuchos::RCP<const Thyra::LinearOpBase<double> > blockedOpBlocks[4][4];
  for (int i = 0; i <= dim_; i++) {
    for (int j = 0; j <= dim_; j++) {
      blockedOpBlocks[i][j] = blkOp->getBlock(i, j);
    }
  }

  // Pressure mass matrix.
  if (pressureMassMatrix_ != Teuchos::null) {
    invPressureMassMatrix_ = getInvDiagonalOp(pressureMassMatrix_);
  }
  // We need the size of the sub-block to build the identity matrix.
  else {
    std::cout << "Pressure mass matrix is null. Use identity." << std::endl;
    pressureMassMatrix_    = Thyra::identity<double>(blockedOpBlocks[dim_][0]->range());
    invPressureMassMatrix_ = Thyra::identity<double>(blockedOpBlocks[dim_][0]->range());
  }

  // Build the AL operator.
  Teuchos::RCP<Thyra::DefaultBlockedLinearOp<double> > alOperator =
      Thyra::defaultBlockedLinearOp<double>();
  alOperator->beginBlockFill(dim_ + 1, dim_ + 1);

  // Set blocks for the velocity parts and gradient.
  for (int i = 0; i < dim_; i++) {
    for (int j = 0; j < dim_; j++) {
      // build the blocks and place it the right location
      alOperator->setBlock(
          i, j,
          Thyra::add(
              blockedOpBlocks[i][j],
              Thyra::scale(gamma_, Thyra::multiply(blockedOpBlocks[i][dim_], invPressureMassMatrix_,
                                                   blockedOpBlocks[dim_][j]))));
    }  // end for j
  }    // end for i

  // Last row. Divergence and (possible) stabilization matrix.
  for (int j = 0; j <= dim_; j++) {
    alOperator->setBlock(dim_, j, blockedOpBlocks[dim_][j]);
  }

  // Last column. Gradient.
  for (int i = 0; i < dim_; i++) {
    alOperator->setBlock(
        i, dim_,
        Thyra::add(
            blockedOpBlocks[i][dim_],
            Thyra::scale(gamma_, Thyra::multiply(blockedOpBlocks[i][dim_], invPressureMassMatrix_,
                                                 blockedOpBlocks[dim_][dim_]))));
  }

  alOperator->endBlockFill();

  // Set whatever is returned.
  SetOperator(alOperator, false);

  // Set operator for augmenting the right-hand side.
  Teuchos::RCP<Thyra::DefaultBlockedLinearOp<double> > alOpRhs_ =
      Thyra::defaultBlockedLinearOp<double>();
  alOpRhs_->beginBlockFill(dim_ + 1, dim_ + 1);

  // Identity matrices on the main diagonal.
  for (int i = 0; i < dim_; i++) {
    // build the blocks and place it the right location
    alOpRhs_->setBlock(i, i, Thyra::identity<double>(blockedOpBlocks[0][0]->range()));
  }  // end for i
  alOpRhs_->setBlock(dim_, dim_, Thyra::identity<double>(blockedOpBlocks[dim_][dim_]->range()));

  // Last column.
  for (int i = 0; i < dim_; i++) {
    alOpRhs_->setBlock(
        i, dim_,
        Thyra::scale(gamma_, Thyra::multiply(blockedOpBlocks[i][dim_], invPressureMassMatrix_)));
  }

  alOpRhs_->endBlockFill();

  alOperatorRhs_ = alOpRhs_;

  // reorder if necessary
  if (reorderManager_ != Teuchos::null) Reorder(*reorderManager_);
}

void ALOperator::augmentRHS(const Epetra_MultiVector& b, Epetra_MultiVector& bAugmented) {
  Teuchos::RCP<const Teko::Epetra::MappingStrategy> mapping = this->getMapStrategy();
  Teuchos::RCP<Thyra::MultiVectorBase<double> > bThyra =
      Thyra::createMembers(thyraOp_->range(), b.NumVectors());
  Teuchos::RCP<Thyra::MultiVectorBase<double> > bThyraAugmented =
      Thyra::createMembers(thyraOp_->range(), b.NumVectors());
  // std::cout << Teuchos::describe(*bThyra, Teuchos::VERB_EXTREME) << std::endl;
  //  Copy Epetra vector to Thyra vector.
  mapping->copyEpetraIntoThyra(b, bThyra.ptr());
  // Apply operator.
  alOperatorRhs_->apply(Thyra::NOTRANS, *bThyra, bThyraAugmented.ptr(), 1.0, 0.0);
  // Copy Thyra vector to Epetra vector.
  mapping->copyThyraIntoEpetra(bThyraAugmented, bAugmented);
}

}  // end namespace NS

}  // end namespace Teko
