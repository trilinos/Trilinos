// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_PCDStrategy.hpp"

#include "Teuchos_TimeMonitor.hpp"
#include "Teko_Utilities.hpp"

namespace Teko {
namespace NS {

using Teuchos::TimeMonitor;

Teuchos::RCP<Teuchos::Time> PCDStrategy::initTimer_;
Teuchos::RCP<Teuchos::Time> PCDStrategy::invSTimer_;
Teuchos::RCP<Teuchos::Time> PCDStrategy::invFTimer_;
Teuchos::RCP<Teuchos::Time> PCDStrategy::opsTimer_;

void PCDStrategy::buildTimers() {
  if (initTimer_ == Teuchos::null)
    initTimer_ = TimeMonitor::getNewTimer("PCDStrategy::initializePrec");

  if (invSTimer_ == Teuchos::null)
    invSTimer_ = TimeMonitor::getNewTimer("PCDStrategy::initializePrec invS");

  if (invFTimer_ == Teuchos::null)
    invFTimer_ = TimeMonitor::getNewTimer("PCDStrategy::initializePrec invF");

  if (opsTimer_ == Teuchos::null)
    opsTimer_ = TimeMonitor::getNewTimer("PCDStrategy::initializePrec buildOps");
}

PCDStrategy::PCDStrategy() : massInverseType_(Diagonal), schurCompOrdering_(false) {
  pcdParams_ = Teuchos::rcp(new Teuchos::ParameterList);
  lapParams_ = Teuchos::rcp(new Teuchos::ParameterList);

  lapParams_->set("Name", getPressureLaplaceString());
  pcdParams_->set("Name", getPCDString());

  buildTimers();
}

//! Constructor to set the inverse factories.
PCDStrategy::PCDStrategy(const Teuchos::RCP<InverseFactory>& invFA,
                         const Teuchos::RCP<InverseFactory>& invS)
    : invFactoryF_(invFA),
      invFactoryS_(invS),
      massInverseType_(Diagonal),
      schurCompOrdering_(false) {
  pcdParams_ = Teuchos::rcp(new Teuchos::ParameterList);
  lapParams_ = Teuchos::rcp(new Teuchos::ParameterList);

  lapParams_->set("Name", getPressureLaplaceString());
  pcdParams_->set("Name", getPCDString());

  buildTimers();
}

/** returns the first (approximate) inverse of \f$A_{00}\f$ */
const Teko::LinearOp PCDStrategy::getHatInvA00(const Teko::BlockedLinearOp& A,
                                               BlockPreconditionerState& state) const {
  initializeState(A, state);

  return state.getModifiableOp("invF");
}

/** returns the second (approximate) inverse of \f$A_{00}\f$ */
const Teko::LinearOp PCDStrategy::getTildeInvA00(const Teko::BlockedLinearOp& A,
                                                 BlockPreconditionerState& state) const {
  initializeState(A, state);

  return state.getModifiableOp("invF");
}

/** returns an (approximate) inverse of \f$S = -A_{11} + A_{10} \mbox{diag}(A_{00})^{-1} A_{01}\f$
 */
const Teko::LinearOp PCDStrategy::getInvS(const Teko::BlockedLinearOp& A,
                                          BlockPreconditionerState& state) const {
  initializeState(A, state);

  return state.getLinearOp("invS");
}

void PCDStrategy::initializeState(const Teko::BlockedLinearOp& A,
                                  BlockPreconditionerState& state) const {
  Teko_DEBUG_SCOPE("PCDStrategy::initializeState", 10);
  TEUCHOS_ASSERT(getRequestHandler() != Teuchos::null);

  std::string pcdStr      = getPCDString();
  std::string presLapStr  = getPressureLaplaceString();
  std::string presMassStr = getPressureMassString();

  // no work to be done
  if (state.isInitialized()) return;

  Teuchos::TimeMonitor timer(*initTimer_, true);

  // extract sub blocks
  LinearOp F  = Teko::getBlock(0, 0, A);
  LinearOp Bt = Teko::getBlock(0, 1, A);
  LinearOp B  = Teko::getBlock(1, 0, A);
  LinearOp C  = Teko::getBlock(1, 1, A);

  LinearOp Qp = getRequestHandler()->request<LinearOp>(presMassStr);
  TEUCHOS_ASSERT(Qp != Teuchos::null);

  // build the inverse Laplacian complement
  /////////////////////////////////////////////
  LinearOp iQp;
  if (massInverseType_ == NotDiag) {
    ModifiableLinearOp& invMass = state.getModifiableOp("invMass");
    Teko_DEBUG_SCOPE("Building inv(Mass)", 10);

    if (invMass == Teuchos::null)
      invMass = buildInverse(*invFactoryS_, Qp);
    else
      rebuildInverse(*invFactoryS_, Qp, invMass);

    iQp = invMass;
  } else {
    Teko_DEBUG_MSG(
        "Building inverse mass of type \"" << Teko::getDiagonalName(massInverseType_) << "\"", 10);
    iQp = getInvDiagonalOp(Qp, massInverseType_);
  }

  // build the inverse Laplacian complement
  /////////////////////////////////////////////
  ModifiableLinearOp& invLaplace = state.getModifiableOp("invLaplace");
  {
    Teuchos::TimeMonitor timerInvS(*invSTimer_, true);

    // LinearOp laplace = getRequestHandler()->request<Teko::LinearOp>(presLapStr);
    LinearOp laplace = getRequestHandler()->request<Teko::LinearOp>(RequestMesg(lapParams_));
    TEUCHOS_ASSERT(laplace != Teuchos::null);
    if (invLaplace == Teuchos::null)
      invLaplace = buildInverse(*invFactoryS_, laplace);
    else
      rebuildInverse(*invFactoryS_, laplace, invLaplace);
  }

  // build the inverse Schur complement
  /////////////////////////////////////////////
  {
    Teko_DEBUG_SCOPE("Building S", 10);
    Teuchos::TimeMonitor timerS(*opsTimer_, true);

    // build Schur-complement
    // LinearOp pcd = getRequestHandler()->request<Teko::LinearOp>(pcdStr);
    LinearOp pcd = getRequestHandler()->request<Teko::LinearOp>(RequestMesg(pcdParams_));
    TEUCHOS_ASSERT(pcd != Teuchos::null);
    LinearOp invL = invLaplace;

    LinearOp invS;
    if (schurCompOrdering_ == false)
      invS = multiply(iQp, pcd, invL);
    else
      invS = multiply(invL, pcd, iQp);

    state.addLinearOp("invS", invS);
  }

  // build inverse F
  /////////////////////////////////////////////
  {
    Teko_DEBUG_SCOPE("Building inv(F)", 10);
    Teuchos::TimeMonitor timerInvF(*invFTimer_, true);

    ModifiableLinearOp& invF = state.getModifiableOp("invF");
    if (invF == Teuchos::null)
      invF = buildInverse(*invFactoryF_, F);
    else
      rebuildInverse(*invFactoryF_, F, invF);
  }

  // mark state as initialized
  state.setInitialized(true);
}

/** \brief This function builds the internals of the state from a parameter list.
 *
 * This function builds the internals of the LU 2x2 state
 * from a parameter list. Furthermore, it allows a
 * developer to easily add a factory to the build system.
 *
 * \param[in] settings Parameter list to use as the internal settings
 * \param[in] invLib Inverse library to use for building inverse factory objects
 *
 * \note The default implementation does nothing.
 */
void PCDStrategy::initializeFromParameterList(const Teuchos::ParameterList& pl,
                                              const InverseLibrary& invLib) {
  Teko_DEBUG_SCOPE("PCDStrategy::initializeFromParameterList", 10);

  std::string invStr = "", invFStr = "", invSStr = "";
#if defined(Teko_ENABLE_Amesos)
  invStr = "Amesos";
#elif defined(Teko_ENABLE_Amesos2)
  invStr = "Amesos2";
#endif

  massInverseType_ = Diagonal;

  // "parse" the parameter list
  if (pl.isParameter("Inverse Type")) invStr = pl.get<std::string>("Inverse Type");
  if (pl.isParameter("Inverse F Type")) invFStr = pl.get<std::string>("Inverse F Type");
  if (pl.isParameter("Inverse Laplace Type")) invSStr = pl.get<std::string>("Inverse Laplace Type");
  if (pl.isParameter("Inverse Mass Type")) {
    std::string massInverseStr = pl.get<std::string>("Inverse Mass Type");

    // build inverse types
    massInverseType_ = getDiagonalType(massInverseStr);
  }
  if (pl.isParameter("Flip Schur Complement Ordering"))
    schurCompOrdering_ = pl.get<bool>("Flip Schur Complement Ordering");

  // set defaults as needed
  if (invFStr == "") invFStr = invStr;
  if (invSStr == "") invSStr = invStr;

  // read pressure laplace parameters
  if (pl.isSublist("Pressure Laplace Parameters"))
    lapParams_ =
        Teuchos::rcp(new Teuchos::ParameterList(pl.sublist("Pressure Laplace Parameters")));
  else
    lapParams_ = Teuchos::rcp(new Teuchos::ParameterList);

  // read pcd operator parameters
  if (pl.isSublist("Pressure Convection Diffusion Parameters"))
    pcdParams_ = Teuchos::rcp(
        new Teuchos::ParameterList(pl.sublist("Pressure Convection Diffusion Parameters")));
  else
    pcdParams_ = Teuchos::rcp(new Teuchos::ParameterList);

  // The user should not have already added this parameters
  TEUCHOS_TEST_FOR_EXCEPTION(
      lapParams_->isParameter("Name"), std::logic_error,
      "Teko: Parameter \"Name\" is not allowed in the sublist \"" + lapParams_->name() + "\"");
  TEUCHOS_TEST_FOR_EXCEPTION(
      lapParams_->isParameter("Tag"), std::logic_error,
      "Teko: Parameter \"Tag\" is not allowed in the sublist \"" + lapParams_->name() + "\"");
  TEUCHOS_TEST_FOR_EXCEPTION(
      pcdParams_->isParameter("Name"), std::logic_error,
      "Teko: Parameter \"Name\" is not allowed in the sublist \"" + pcdParams_->name() + "\"");
  TEUCHOS_TEST_FOR_EXCEPTION(
      pcdParams_->isParameter("Tag"), std::logic_error,
      "Teko: Parameter \"Tag\" is not allowed in the sublist \"" + pcdParams_->name() + "\"");

  Teko_DEBUG_MSG_BEGIN(5) DEBUG_STREAM << "PCD Strategy Parameters: " << std::endl;
  DEBUG_STREAM << "   inv type   = \"" << invStr << "\"" << std::endl;
  DEBUG_STREAM << "   inv F type = \"" << invFStr << "\"" << std::endl;
  DEBUG_STREAM << "   inv Laplace type = \"" << invSStr << "\"" << std::endl;
  DEBUG_STREAM << "   inv Mass type = \"" << Teko::getDiagonalName(massInverseType_) << "\""
               << std::endl;
  DEBUG_STREAM << "PCD Strategy Parameter list: " << std::endl;
  pl.print(DEBUG_STREAM);
  Teko_DEBUG_MSG_END()

      // build velocity inverse factory
      invFactoryF_ = invLib.getInverseFactory(invFStr);

  if (invFStr == invSStr)
    invFactoryS_ = invFactoryF_;
  else
    invFactoryS_ = invLib.getInverseFactory(invSStr);

  lapParams_->set("Name", getPressureLaplaceString());
  pcdParams_->set("Name", getPCDString());

  // setup a request for required operators
  getRequestHandler()->preRequest<Teko::LinearOp>(getPressureMassString());
  // getRequestHandler()->preRequest<Teko::LinearOp>(getPCDString());
  // getRequestHandler()->preRequest<Teko::LinearOp>(getPressureLaplaceString());
  getRequestHandler()->preRequest<Teko::LinearOp>(Teko::RequestMesg(lapParams_));
  getRequestHandler()->preRequest<Teko::LinearOp>(Teko::RequestMesg(pcdParams_));
}

//! For assiting in construction of the preconditioner
Teuchos::RCP<Teuchos::ParameterList> PCDStrategy::getRequestedParameters() const {
  TEUCHOS_ASSERT(false);

  return Teuchos::null;
}

//! For assiting in construction of the preconditioner
bool PCDStrategy::updateRequestedParameters(const Teuchos::ParameterList& /* pl */) {
  TEUCHOS_ASSERT(false);

  return true;
}

}  // end namespace NS
}  // end namespace Teko
