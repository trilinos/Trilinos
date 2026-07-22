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

#include "Teko_ModALPreconditionerFactory.hpp"

#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"

#include "Teko_LU2x2InverseOp.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"
#include "Teko_StaticLSCStrategy.hpp"
#include "Teko_InvLSCStrategy.hpp"
#include "Teko_PresLaplaceLSCStrategy.hpp"

#include "Teuchos_Time.hpp"

namespace Teko {

namespace NS {

ModALPrecondState::ModALPrecondState()
    : pressureMassMatrix_(Teuchos::null), invPressureMassMatrix_(Teuchos::null) {}

ModALPreconditionerFactory::ModALPreconditionerFactory(const Teuchos::RCP<InverseFactory>& factory)
    : invOpsStrategy_(Teuchos::rcp(new InvModALStrategy(factory))), isSymmetric_(true) {}

ModALPreconditionerFactory::ModALPreconditionerFactory(
    const Teuchos::RCP<InverseFactory>& invFactoryA,
    const Teuchos::RCP<InverseFactory>& invFactoryS)
    : invOpsStrategy_(Teuchos::rcp(new InvModALStrategy(invFactoryA, invFactoryS))),
      isSymmetric_(true) {}

ModALPreconditionerFactory::ModALPreconditionerFactory(const Teuchos::RCP<InverseFactory>& factory,
                                                       LinearOp& pressureMassMatrix)
    : invOpsStrategy_(Teuchos::rcp(new InvModALStrategy(factory, pressureMassMatrix))),
      isSymmetric_(true) {}

ModALPreconditionerFactory::ModALPreconditionerFactory(
    const Teuchos::RCP<InverseFactory>& invFactoryA,
    const Teuchos::RCP<InverseFactory>& invFactoryS, LinearOp& pressureMassMatrix)
    : invOpsStrategy_(
          Teuchos::rcp(new InvModALStrategy(invFactoryA, invFactoryS, pressureMassMatrix))),
      isSymmetric_(true) {}

// Construct from a strategy.
ModALPreconditionerFactory::ModALPreconditionerFactory(
    const Teuchos::RCP<InvModALStrategy>& strategy)
    : invOpsStrategy_(strategy), isSymmetric_(true) {}

LinearOp ModALPreconditionerFactory::buildPreconditionerOperator(
    BlockedLinearOp& alOp, BlockPreconditionerState& state) const {
  Teko_DEBUG_SCOPE("ModALPreconditionerFactory::buildPreconditionerOperator()", 10);
  Teko_DEBUG_EXPR(Teuchos::Time timer(""));
  Teko_DEBUG_EXPR(Teuchos::Time totalTimer(""));
  Teko_DEBUG_EXPR(totalTimer.start());

  // Only for 2D or 3D problems.
  int dim = blockRowCount(alOp) - 1;
  TEUCHOS_ASSERT(dim == 2 || dim == 3);

  // Build what is necessary for the state object.
  Teko_DEBUG_EXPR(timer.start(true));
  invOpsStrategy_->buildState(alOp, state);
  Teko_DEBUG_EXPR(timer.stop());
  Teko_DEBUG_MSG("ModALPreconditionerFactory::buildPreconditionerOperator():BuildStateTime = "
                     << timer.totalElapsedTime(),
                 2);

  // Extract inverse operators from strategy
  Teko_DEBUG_EXPR(timer.start(true));
  LinearOp invA11p = invOpsStrategy_->getInvA11p(state);
  LinearOp invA22p = invOpsStrategy_->getInvA22p(state);
  LinearOp invA33p;
  if (dim == 3) {
    invA33p = invOpsStrategy_->getInvA33p(state);
  }

  // The inverse of S can be built from strategy,
  // or just a diagonal matrix.
  ModALPrecondState* modALState = dynamic_cast<ModALPrecondState*>(&state);
  TEUCHOS_ASSERT(modALState != NULL);
  LinearOp invS;
  if (modALState->isStabilized_) {
    invS = invOpsStrategy_->getInvS(state);
  } else {
    invS = scale(modALState->gamma_, modALState->invPressureMassMatrix_);
  }

  Teko_DEBUG_EXPR(timer.stop());
  Teko_DEBUG_MSG(
      "ModALPrecFact::buildPreconditionerOperator(): GetInvTime = " << timer.totalElapsedTime(), 2);

  // Build diagonal operations.
  std::vector<LinearOp> invDiag;
  invDiag.resize(dim + 1);
  invDiag[0] = invA11p;
  invDiag[1] = invA22p;
  if (dim == 2) {
    invDiag[2] = scale(-1.0, invS);
  } else if (dim == 3) {
    invDiag[2] = invA33p;
    invDiag[3] = scale(-1.0, invS);
  }

  // Get the upper triangular matrix.
  BlockedLinearOp U = getUpperTriBlocks(alOp);

  Teko_DEBUG_EXPR(totalTimer.stop());
  Teko_DEBUG_MSG(
      "ModALPrecFact::buildPreconditionerOperator TotalTime = " << totalTimer.totalElapsedTime(),
      2);

  // Create the preconditioner
  return createBlockUpperTriInverseOp(U, invDiag, "Modified AL preconditioner-Upper");
}

}  // end namespace NS

}  // end namespace Teko
