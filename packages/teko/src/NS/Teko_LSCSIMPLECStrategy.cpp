// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_LSCSIMPLECStrategy.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

// Teko includes
#include "Teko_Utilities.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"

using Teuchos::RCP;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

namespace Teko {
namespace NS {

/////////////////////////////////////////////////////////////////////////////
// LSCSIMPLECStrategy Implementation
/////////////////////////////////////////////////////////////////////////////

// constructors
/////////////////////////////////////////////////////////////////////////////
LSCSIMPLECStrategy::LSCSIMPLECStrategy()
    : invFactoryF_(Teuchos::null),
      invFactoryS_(Teuchos::null),
      useFullLDU_(false),
      scaleType_(Diagonal) {}

/////////////////////////////////////////////////////////////////////////////

void LSCSIMPLECStrategy::buildState(BlockedLinearOp& A, BlockPreconditionerState& state) const {
  Teko_DEBUG_SCOPE("LSCSIMPLECStrategy::buildState", 10);

  LSCPrecondState* lscState = dynamic_cast<LSCPrecondState*>(&state);
  TEUCHOS_ASSERT(lscState != 0);

  // if neccessary save state information
  if (not lscState->isInitialized()) {
    Teko_DEBUG_EXPR(Teuchos::Time timer(""));

    // construct operators
    {
      Teko_DEBUG_SCOPE("LSC-SIMPLEC::buildState constructing operators", 1);
      Teko_DEBUG_EXPR(timer.start(true));

      initializeState(A, lscState);

      Teko_DEBUG_EXPR(timer.stop());
      Teko_DEBUG_MSG("LSC-SIMPLEC::buildState BuildOpsTime = " << timer.totalElapsedTime(), 1);
    }

    // Build the inverses
    {
      Teko_DEBUG_SCOPE("LSC-SIMPLEC::buildState calculating inverses", 1);
      Teko_DEBUG_EXPR(timer.start(true));

      computeInverses(A, lscState);

      Teko_DEBUG_EXPR(timer.stop());
      Teko_DEBUG_MSG("LSC-SIMPLEC::buildState BuildInvTime = " << timer.totalElapsedTime(), 1);
    }
  }
}

// functions inherited from LSCStrategy
LinearOp LSCSIMPLECStrategy::getInvBQBt(const BlockedLinearOp& /* A */,
                                        BlockPreconditionerState& state) const {
  return state.getInverse("invBQBtmC");
}

LinearOp LSCSIMPLECStrategy::getInvBHBt(const BlockedLinearOp& /* A */,
                                        BlockPreconditionerState& state) const {
  return state.getInverse("invBQBtmC").getConst();
}

LinearOp LSCSIMPLECStrategy::getInvF(const BlockedLinearOp& /* A */,
                                     BlockPreconditionerState& state) const {
  return state.getInverse("invF");
}

LinearOp LSCSIMPLECStrategy::getInnerStabilization(const BlockedLinearOp& A,
                                                   BlockPreconditionerState& state) const {
  LSCPrecondState* lscState = dynamic_cast<LSCPrecondState*>(&state);
  TEUCHOS_ASSERT(lscState != 0);
  TEUCHOS_ASSERT(lscState->isInitialized())

  const LinearOp C = getBlock(1, 1, A);
  return scale(-1.0, C);
}

LinearOp LSCSIMPLECStrategy::getInvMass(const BlockedLinearOp& /* A */,
                                        BlockPreconditionerState& state) const {
  LSCPrecondState* lscState = dynamic_cast<LSCPrecondState*>(&state);
  TEUCHOS_ASSERT(lscState != 0);
  TEUCHOS_ASSERT(lscState->isInitialized())

  return lscState->invMass_;
}

LinearOp LSCSIMPLECStrategy::getHScaling(const BlockedLinearOp& A,
                                         BlockPreconditionerState& state) const {
  return getInvMass(A, state);
}

//! Initialize the state object using this blocked linear operator
void LSCSIMPLECStrategy::initializeState(const BlockedLinearOp& A, LSCPrecondState* state) const {
  Teko_DEBUG_SCOPE("LSCSIMPLECStrategy::initializeState", 10);

  const LinearOp F  = getBlock(0, 0, A);
  const LinearOp Bt = getBlock(0, 1, A);
  const LinearOp B  = getBlock(1, 0, A);
  const LinearOp C  = getBlock(1, 1, A);

  bool isStabilized = (not isZeroOp(C));

  state->invMass_ = getInvDiagonalOp(F, scaleType_);

  // compute BQBt
  state->BQBt_ = explicitMultiply(B, state->invMass_, Bt, state->BQBt_);
  if (isStabilized) {
    // now build B*Q*Bt-C
    Teko::ModifiableLinearOp BQBtmC = state->getInverse("BQBtmC");
    BQBtmC                          = explicitAdd(state->BQBt_, scale(-1.0, C), BQBtmC);
    state->addInverse("BQBtmC", BQBtmC);
  }
  Teko_DEBUG_MSG("Computed BQBt", 10);

  state->setInitialized(true);
}

/** Compute the inverses required for the LSC Schur complement
 *
 * \note This method assumes that the BQBt and BHBt operators have
 *       been constructed.
 */
void LSCSIMPLECStrategy::computeInverses(const BlockedLinearOp& A, LSCPrecondState* state) const {
  Teko_DEBUG_SCOPE("LSCSIMPLECStrategy::computeInverses", 10);
  Teko_DEBUG_EXPR(Teuchos::Time invTimer(""));

  const LinearOp F = getBlock(0, 0, A);

  /////////////////////////////////////////////////////////

  // (re)build the inverse of F
  Teko_DEBUG_MSG("LSC-SIMPLEC::computeInverses Building inv(F)", 1);
  Teko_DEBUG_EXPR(invTimer.start(true));
  InverseLinearOp invF = state->getInverse("invF");
  if (invF == Teuchos::null) {
    invF = buildInverse(*invFactoryF_, F);
    state->addInverse("invF", invF);
  } else {
    rebuildInverse(*invFactoryF_, F, invF);
  }
  Teko_DEBUG_EXPR(invTimer.stop());
  Teko_DEBUG_MSG("LSC-SIMPLEC::computeInverses GetInvF = " << invTimer.totalElapsedTime(), 1);

  /////////////////////////////////////////////////////////

  // (re)build the inverse of BQBt
  Teko_DEBUG_MSG("LSC-SIMPLEC::computeInverses Building inv(BQBtmC)", 1);
  Teko_DEBUG_EXPR(invTimer.start(true));
  const LinearOp BQBt     = state->getInverse("BQBtmC");
  InverseLinearOp invBQBt = state->getInverse("invBQBtmC");
  if (invBQBt == Teuchos::null) {
    invBQBt = buildInverse(*invFactoryS_, BQBt);
    state->addInverse("invBQBtmC", invBQBt);
  } else {
    rebuildInverse(*invFactoryS_, BQBt, invBQBt);
  }
  Teko_DEBUG_EXPR(invTimer.stop());
  Teko_DEBUG_MSG("LSC-SIMPLEC::computeInverses GetInvBQBt = " << invTimer.totalElapsedTime(), 1);
}

//! Initialize from a parameter list
void LSCSIMPLECStrategy::initializeFromParameterList(const Teuchos::ParameterList& pl,
                                                     const InverseLibrary& invLib) {
  // get string specifying inverse
  std::string invStr = "", invVStr = "", invPStr = "";
  bool useLDU = false;
  scaleType_  = Diagonal;

  // "parse" the parameter list
  if (pl.isParameter("Inverse Type")) invStr = pl.get<std::string>("Inverse Type");
  if (pl.isParameter("Inverse Velocity Type"))
    invVStr = pl.get<std::string>("Inverse Velocity Type");
  if (pl.isParameter("Inverse Pressure Type"))
    invPStr = pl.get<std::string>("Inverse Pressure Type");
  if (pl.isParameter("Use LDU")) useLDU = pl.get<bool>("Use LDU");
  if (pl.isParameter("Scaling Type")) {
    scaleType_ = getDiagonalType(pl.get<std::string>("Scaling Type"));
    TEUCHOS_TEST_FOR_EXCEPT(scaleType_ == NotDiag);
  }

  Teko_DEBUG_MSG_BEGIN(0) DEBUG_STREAM << "LSC Inverse Strategy Parameters: " << std::endl;
  DEBUG_STREAM << "   inv type   = \"" << invStr << "\"" << std::endl;
  DEBUG_STREAM << "   inv v type = \"" << invVStr << "\"" << std::endl;
  DEBUG_STREAM << "   inv p type = \"" << invPStr << "\"" << std::endl;
  DEBUG_STREAM << "   use ldu    = " << useLDU << std::endl;
  DEBUG_STREAM << "   scale type    = " << getDiagonalName(scaleType_) << std::endl;
  DEBUG_STREAM << "LSC  Inverse Strategy Parameter list: " << std::endl;
  pl.print(DEBUG_STREAM);
  Teko_DEBUG_MSG_END()

  // set defaults as needed
#if defined(Teko_ENABLE_Amesos)
      if (invStr == "") invStr = "Amesos";
#elif defined(Teko_ENABLE_Amesos2)
      if (invStr == "") invStr = "Amesos2";
#endif
  if (invVStr == "") invVStr = invStr;
  if (invPStr == "") invPStr = invStr;

  // build velocity inverse factory
  invFactoryF_ = invLib.getInverseFactory(invVStr);
  invFactoryS_ = invFactoryF_;  // by default these are the same
  if (invVStr != invPStr)       // if different, build pressure inverse factory
    invFactoryS_ = invLib.getInverseFactory(invPStr);

  // set other parameters
  setUseFullLDU(useLDU);
}

}  // end namespace NS
}  // end namespace Teko
