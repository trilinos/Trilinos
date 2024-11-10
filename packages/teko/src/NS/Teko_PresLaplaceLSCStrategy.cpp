// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NS/Teko_PresLaplaceLSCStrategy.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

// Teko includes
#include "Teko_Utilities.hpp"
#include "NS/Teko_LSCPreconditionerFactory.hpp"

using Teuchos::RCP;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

namespace Teko {
namespace NS {

/////////////////////////////////////////////////////////////////////////////
// PresLaplaceLSCStrategy Implementation
/////////////////////////////////////////////////////////////////////////////

// constructors
/////////////////////////////////////////////////////////////////////////////
PresLaplaceLSCStrategy::PresLaplaceLSCStrategy()
    : invFactoryV_(Teuchos::null),
      invFactoryP_(Teuchos::null),
      eigSolveParam_(5),
      useFullLDU_(false),
      useMass_(false),
      scaleType_(AbsRowSum) {}

PresLaplaceLSCStrategy::PresLaplaceLSCStrategy(const Teuchos::RCP<InverseFactory>& factory)
    : invFactoryV_(factory),
      invFactoryP_(factory),
      eigSolveParam_(5),
      useFullLDU_(false),
      useMass_(false),
      scaleType_(AbsRowSum) {}

PresLaplaceLSCStrategy::PresLaplaceLSCStrategy(const Teuchos::RCP<InverseFactory>& invFactF,
                                               const Teuchos::RCP<InverseFactory>& invFactS)
    : invFactoryV_(invFactF),
      invFactoryP_(invFactS),
      eigSolveParam_(5),
      useFullLDU_(false),
      useMass_(false),
      scaleType_(AbsRowSum) {}

/////////////////////////////////////////////////////////////////////////////

void PresLaplaceLSCStrategy::buildState(BlockedLinearOp& A, BlockPreconditionerState& state) const {
  Teko_DEBUG_SCOPE("PresLaplaceLSCStrategy::buildState", 10);

  LSCPrecondState* lscState = dynamic_cast<LSCPrecondState*>(&state);
  TEUCHOS_ASSERT(lscState != 0);

  // if neccessary save state information
  if (not lscState->isInitialized()) {
    Teko_DEBUG_EXPR(Teuchos::Time timer(""));

    // construct operators
    {
      Teko_DEBUG_SCOPE("PL-LSC::buildState constructing operators", 1);
      Teko_DEBUG_EXPR(timer.start(true));

      initializeState(A, lscState);

      Teko_DEBUG_EXPR(timer.stop());
      Teko_DEBUG_MSG("PL-LSC::buildState BuildOpsTime = " << timer.totalElapsedTime(), 1);
    }

    // Build the inverses
    {
      Teko_DEBUG_SCOPE("PL-LSC::buildState calculating inverses", 1);
      Teko_DEBUG_EXPR(timer.start(true));

      computeInverses(A, lscState);

      Teko_DEBUG_EXPR(timer.stop());
      Teko_DEBUG_MSG("PL-LSC::buildState BuildInvTime = " << timer.totalElapsedTime(), 1);
    }
  }
}

// functions inherited from LSCStrategy
LinearOp PresLaplaceLSCStrategy::getInvBQBt(const BlockedLinearOp& /* A */,
                                            BlockPreconditionerState& state) const {
  return state.getModifiableOp("invPresLap");
}

LinearOp PresLaplaceLSCStrategy::getInvBHBt(const BlockedLinearOp& /* A */,
                                            BlockPreconditionerState& state) const {
  return state.getModifiableOp("invPresLap");
}

LinearOp PresLaplaceLSCStrategy::getInvF(const BlockedLinearOp& /* A */,
                                         BlockPreconditionerState& state) const {
  return state.getModifiableOp("invF");
}

// LinearOp PresLaplaceLSCStrategy::getInvAlphaD(const BlockedLinearOp & A,BlockPreconditionerState
// & state) const
LinearOp PresLaplaceLSCStrategy::getOuterStabilization(const BlockedLinearOp& /* A */,
                                                       BlockPreconditionerState& state) const {
  LSCPrecondState* lscState = dynamic_cast<LSCPrecondState*>(&state);
  TEUCHOS_ASSERT(lscState != 0);
  TEUCHOS_ASSERT(lscState->isInitialized())

  return lscState->aiD_;
}

LinearOp PresLaplaceLSCStrategy::getInvMass(const BlockedLinearOp& /* A */,
                                            BlockPreconditionerState& state) const {
  LSCPrecondState* lscState = dynamic_cast<LSCPrecondState*>(&state);
  TEUCHOS_ASSERT(lscState != 0);
  TEUCHOS_ASSERT(lscState->isInitialized())

  return lscState->invMass_;
}

LinearOp PresLaplaceLSCStrategy::getHScaling(const BlockedLinearOp& A,
                                             BlockPreconditionerState& state) const {
  return getInvMass(A, state);
}

//! Initialize the state object using this blocked linear operator
void PresLaplaceLSCStrategy::initializeState(const BlockedLinearOp& A,
                                             LSCPrecondState* state) const {
  Teko_DEBUG_SCOPE("PresLaplaceLSCStrategy::initializeState", 10);

  std::string velMassStr = getVelocityMassString();

  const LinearOp F  = getBlock(0, 0, A);
  const LinearOp Bt = getBlock(0, 1, A);
  const LinearOp B  = getBlock(1, 0, A);
  const LinearOp C  = getBlock(1, 1, A);

  LinearOp D = B;
  LinearOp G = Bt;

  bool isStabilized = (not isZeroOp(C));

  // grab operators from state object
  LinearOp massMatrix = state->getLinearOp(velMassStr);

  // The logic follows like this
  //    if there is no mass matrix available --> build from F
  //    if there is a mass matrix and the inverse hasn't yet been built
  //       --> build from the mass matrix
  //    otherwise, there is already an invMass_ matrix that is appropriate
  //       --> use that one
  if (massMatrix == Teuchos::null) {
    Teko_DEBUG_MSG(
        "PL-LSC::initializeState Build Scaling <F> type \"" << getDiagonalName(scaleType_) << "\"",
        1);
    state->invMass_ = getInvDiagonalOp(F, scaleType_);
  } else if (state->invMass_ == Teuchos::null) {
    Teko_DEBUG_MSG("PL-LSC::initializeState Build Scaling <mass> type \""
                       << getDiagonalName(scaleType_) << "\"",
                   1);
    state->invMass_ = getInvDiagonalOp(massMatrix, scaleType_);
  }
  // else "invMass_" should be set and there is no reason to rebuild it

  // if this is a stable discretization...we are done!
  if (not isStabilized) {
    state->aiD_ = Teuchos::null;

    state->setInitialized(true);

    return;
  }

  // compute alpha scaled inv(D): EHSST2007 Eq. 4.29
  // construct B_idF_Bt and save it for refilling later: This could reuse BQBt graph
  LinearOp invDiagF                     = getInvDiagonalOp(F);
  Teko::ModifiableLinearOp& modB_idF_Bt = state->getModifiableOp("BidFBt");
  modB_idF_Bt                           = explicitMultiply(B, invDiagF, Bt, modB_idF_Bt);
  const LinearOp B_idF_Bt               = modB_idF_Bt;

  MultiVector vec_D = getDiagonal(B_idF_Bt);  // this memory could be reused
  update(-1.0, getDiagonal(C), 1.0, vec_D);   // vec_D = diag(B*inv(diag(F))*Bt)-diag(C)
  const LinearOp invD = buildInvDiagonal(vec_D, "inv(D)");

  Teko_DEBUG_MSG("Calculating alpha", 10);
  const LinearOp BidFBtidD = multiply<double>(B_idF_Bt, invD);
  double num = std::fabs(Teko::computeSpectralRad(BidFBtidD, 5e-2, false, eigSolveParam_));
  Teko_DEBUG_MSG("Calculated alpha", 10);
  state->alpha_ = 1.0 / num;
  state->aiD_   = Thyra::scale(state->alpha_, invD);

  Teko_DEBUG_MSG_BEGIN(5) DEBUG_STREAM << "PL-LSC Alpha Parameter = " << state->alpha_ << std::endl;
  Teko_DEBUG_MSG_END()

      state->setInitialized(true);
}

/** Compute the inverses required for the LSC Schur complement
 *
 * \note This method assumes that the BQBt and BHBt operators have
 *       been constructed.
 */
void PresLaplaceLSCStrategy::computeInverses(const BlockedLinearOp& A,
                                             LSCPrecondState* state) const {
  Teko_DEBUG_SCOPE("PresLaplaceLSCStrategy::computeInverses", 10);
  Teko_DEBUG_EXPR(Teuchos::Time invTimer(""));

  std::string presLapStr = getPressureLaplaceString();

  const LinearOp F       = getBlock(0, 0, A);
  const LinearOp presLap = state->getLinearOp(presLapStr);

  /////////////////////////////////////////////////////////

  // (re)build the inverse of F
  Teko_DEBUG_MSG("PL-LSC::computeInverses Building inv(F)", 1);
  Teko_DEBUG_EXPR(invTimer.start(true));
  ModifiableLinearOp& invF = state->getModifiableOp("invF");
  if (invF == Teuchos::null) {
    invF = buildInverse(*invFactoryV_, F);
  } else {
    rebuildInverse(*invFactoryV_, F, invF);
  }
  Teko_DEBUG_EXPR(invTimer.stop());
  Teko_DEBUG_MSG("PL-LSC::computeInverses GetInvF = " << invTimer.totalElapsedTime(), 1);

  /////////////////////////////////////////////////////////

  // (re)build the inverse of P
  Teko_DEBUG_MSG("PL-LSC::computeInverses Building inv(PresLap)", 1);
  Teko_DEBUG_EXPR(invTimer.start(true));
  ModifiableLinearOp& invPresLap = state->getModifiableOp("invPresLap");
  if (invPresLap == Teuchos::null) {
    invPresLap = buildInverse(*invFactoryP_, presLap);
  } else {
    // not need because the pressure laplacian never changes
    // rebuildInverse(*invFactoryP_,presLap,invPresLap);
  }
  Teko_DEBUG_EXPR(invTimer.stop());
  Teko_DEBUG_MSG("PL-LSC::computeInverses GetInvBQBt = " << invTimer.totalElapsedTime(), 1);
}

//! Initialize from a parameter list
void PresLaplaceLSCStrategy::initializeFromParameterList(const Teuchos::ParameterList& pl,
                                                         const InverseLibrary& invLib) {
  // get string specifying inverse
  std::string invStr = "", invVStr = "", invPStr = "";
#if defined(Teko_ENABLE_Amesos)
  invStr = "Amesos";
#elif defined(Teko_ENABLE_Amesos2)
  invStr = "Amesos2";
#endif

  bool useLDU = false;
  scaleType_  = AbsRowSum;

  // "parse" the parameter list
  if (pl.isParameter("Inverse Type")) invStr = pl.get<std::string>("Inverse Type");
  if (pl.isParameter("Inverse Velocity Type"))
    invVStr = pl.get<std::string>("Inverse Velocity Type");
  if (pl.isParameter("Inverse Pressure Type"))
    invPStr = pl.get<std::string>("Inverse Pressure Type");
  if (pl.isParameter("Use LDU")) useLDU = pl.get<bool>("Use LDU");
  if (pl.isParameter("Use Mass Scaling")) useMass_ = pl.get<bool>("Use Mass Scaling");
  if (pl.isParameter("Eigen Solver Iterations"))
    eigSolveParam_ = pl.get<int>("Eigen Solver Iterations");
  if (pl.isParameter("Scaling Type")) {
    scaleType_ = getDiagonalType(pl.get<std::string>("Scaling Type"));
    TEUCHOS_TEST_FOR_EXCEPT(scaleType_ == NotDiag);
  }

  // set defaults as needed
  if (invVStr == "") invVStr = invStr;
  if (invPStr == "") invPStr = invStr;

  Teko_DEBUG_MSG_BEGIN(5) DEBUG_STREAM << "LSC Inverse Strategy Parameters: " << std::endl;
  DEBUG_STREAM << "   inv v type = \"" << invVStr << "\"" << std::endl;
  DEBUG_STREAM << "   inv p type = \"" << invPStr << "\"" << std::endl;
  DEBUG_STREAM << "   use ldu    = " << useLDU << std::endl;
  DEBUG_STREAM << "   use mass    = " << useMass_ << std::endl;
  DEBUG_STREAM << "   scale type    = " << getDiagonalName(scaleType_) << std::endl;
  DEBUG_STREAM << "LSC  Pressure Laplace Strategy Parameter list: " << std::endl;
  pl.print(DEBUG_STREAM);
  Teko_DEBUG_MSG_END()

      // build velocity inverse factory
      invFactoryV_ = invLib.getInverseFactory(invVStr);
  invFactoryP_     = invFactoryV_;  // by default these are the same
  if (invVStr != invPStr)           // if different, build pressure inverse factory
    invFactoryP_ = invLib.getInverseFactory(invPStr);

  // set other parameters
  setUseFullLDU(useLDU);
}

//! For assiting in construction of the preconditioner
Teuchos::RCP<Teuchos::ParameterList> PresLaplaceLSCStrategy::getRequestedParameters() const {
  Teko_DEBUG_SCOPE("PresLaplaceLSCStrategy::getRequestedParameters", 10);
  Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

  // grab parameters from F solver
  RCP<Teuchos::ParameterList> fList = invFactoryV_->getRequestedParameters();
  if (fList != Teuchos::null) {
    Teuchos::ParameterList::ConstIterator itr;
    for (itr = fList->begin(); itr != fList->end(); ++itr) pl->setEntry(itr->first, itr->second);
  }

  // grab parameters from S solver
  RCP<Teuchos::ParameterList> sList = invFactoryP_->getRequestedParameters();
  if (sList != Teuchos::null) {
    Teuchos::ParameterList::ConstIterator itr;
    for (itr = sList->begin(); itr != sList->end(); ++itr) pl->setEntry(itr->first, itr->second);
  }

  // use the mass matrix
  if (useMass_) pl->set(getVelocityMassString(), false, "Velocity mass matrix");
  pl->set(getPressureLaplaceString(), false, "Pressure Laplacian matrix");

  return pl;
}

//! For assiting in construction of the preconditioner
bool PresLaplaceLSCStrategy::updateRequestedParameters(const Teuchos::ParameterList& pl) {
  Teko_DEBUG_SCOPE("PresLaplaceLSCStrategy::updateRequestedParameters", 10);
  bool result = true;

  // update requested parameters in solvers
  result &= invFactoryV_->updateRequestedParameters(pl);
  result &= invFactoryP_->updateRequestedParameters(pl);

  Teuchos::ParameterList hackList(pl);

  // get required operator acknowledgment...user must set these to true
  bool plo = hackList.get<bool>(getPressureLaplaceString(), false);

  bool vmo = true;
  if (useMass_) vmo = hackList.get<bool>(getVelocityMassString(), false);

  if (not plo) {
    Teko_DEBUG_MSG("User must acknowledge the use of the \"" << getPressureLaplaceString() << "\"!",
                   0);
  }
  if (not vmo) {
    Teko_DEBUG_MSG("User must acknowledge the use of the \"" << getVelocityMassString() << "\"!",
                   0);
  }

  result &= (plo & vmo);

  return result;
}

}  // end namespace NS
}  // end namespace Teko
