// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_TimingsSIMPLEPreconditionerFactory.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"
#ifdef TEKO_HAVE_EPETRA
#include "Teko_DiagonalPreconditionerFactory.hpp"
#endif

#include "Teuchos_Time.hpp"

using Teuchos::RCP;

namespace Teko {
namespace NS {

// Constructor definition
TimingsSIMPLEPreconditionerFactory ::TimingsSIMPLEPreconditionerFactory(
    const RCP<InverseFactory>& inverse, double alpha)
    : invVelFactory_(inverse),
      invPrsFactory_(inverse),
      alpha_(alpha),
      fInverseType_(Diagonal),
      useMass_(false),
      constrTotal_("SIMPLE Constr: Total"),
      subTotal_("SIMPLE Constr: Subs"),
      constrCount_(0) {}

TimingsSIMPLEPreconditionerFactory ::TimingsSIMPLEPreconditionerFactory(
    const RCP<InverseFactory>& invVFact, const RCP<InverseFactory>& invPFact, double alpha)
    : invVelFactory_(invVFact),
      invPrsFactory_(invPFact),
      alpha_(alpha),
      fInverseType_(Diagonal),
      useMass_(false),
      constrTotal_("SIMPLE Constr: Total"),
      subTotal_("SIMPLE Constr: Subs"),
      constrCount_(0) {}

TimingsSIMPLEPreconditionerFactory::TimingsSIMPLEPreconditionerFactory()
    : alpha_(1.0),
      fInverseType_(Diagonal),
      useMass_(false),
      constrTotal_("SIMPLE Constr: Total"),
      subTotal_("SIMPLE Constr: Subs"),
      constrCount_(0) {}

TimingsSIMPLEPreconditionerFactory::~TimingsSIMPLEPreconditionerFactory() {
  if (constrTotal_.totalElapsedTime() > 0.0) {
    Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
    out.setOutputToRootOnly(0);

    out << "==========================================================================="
        << std::endl;
    out << std::endl;
    out << "SIMPLE Construction Count   = " << constrCount_ << std::endl;
    out << "SIMPLE Construction Total   = " << constrTotal_.totalElapsedTime() << std::endl;
    out << "SIMPLE Sub Components Total = " << subTotal_.totalElapsedTime() << std::endl;
    out << std::endl;
    out << "==========================================================================="
        << std::endl;
  }
}

// Use the factory to build the preconditioner (this is where the work goes)
LinearOp TimingsSIMPLEPreconditionerFactory ::buildPreconditionerOperator(
    BlockedLinearOp& blockOp, BlockPreconditionerState& state) const {
  constrTotal_.start();

  Teko_DEBUG_SCOPE("TimingsSIMPLEPreconditionerFactory::buildPreconditionerOperator", 10);
  Teko_DEBUG_EXPR(Teuchos::Time timer(""));

  int rows = blockRowCount(blockOp);
  int cols = blockColCount(blockOp);

  TEUCHOS_ASSERT(rows == 2);  // sanity checks
  TEUCHOS_ASSERT(cols == 2);

  bool buildExplicitSchurComplement = true;

  // extract subblocks
  const LinearOp F  = getBlock(0, 0, blockOp);
  const LinearOp Bt = getBlock(0, 1, blockOp);
  const LinearOp B  = getBlock(1, 0, blockOp);
  const LinearOp C  = getBlock(1, 1, blockOp);

  LinearOp matF = F;
  if (useMass_) {
    TEUCHOS_ASSERT(massMatrix_ != Teuchos::null);
    matF = massMatrix_;
  }

  // get approximation of inv(F) name H
  std::string fApproxStr = "<error>";
  LinearOp H;
  if (fInverseType_ == NotDiag) {
    H          = buildInverse(*customHFactory_, matF);
    fApproxStr = customHFactory_->toString();

    // since H is now implicit, we must build an implicit Schur complement
    buildExplicitSchurComplement = false;
  } else if (fInverseType_ == BlkDiag) {
#ifdef TEKO_HAVE_EPETRA
    // Block diagonal approximation for H
    DiagonalPreconditionerFactory Hfact;
    DiagonalPrecondState Hstate;
    Hfact.initializeFromParameterList(BlkDiagList_);
    H = Hfact.buildPreconditionerOperator(matF, Hstate);

    buildExplicitSchurComplement = true;  // NTS: Do I need this?
                                          // Answer - no, but it is documenting whats going on here.
#else
    throw std::logic_error(
        "TimingsSIMPLEPreconditionerFactory fInverseType_ == "
        "BlkDiag but EPETRA is turned off!");
#endif
  } else {
    // get generic diagonal
    subTotal_.start();
    H = getInvDiagonalOp(matF, fInverseType_);
    subTotal_.stop();
    fApproxStr = getDiagonalName(fInverseType_);
  }

  // adjust H for time scaling if it is a mass matrix
  if (useMass_) {
    RCP<const Teuchos::ParameterList> pl = state.getParameterList();

    if (pl->isParameter("stepsize")) {
      // get the step size
      double stepsize = pl->get<double>("stepsize");

      // scale by stepsize only if it is larger than 0
      if (stepsize > 0.0) H = scale(stepsize, H);
    }
  }

  // build approximate Schur complement: hatS = -C + B*H*Bt
  LinearOp HBt, hatS;

  if (buildExplicitSchurComplement) {
    ModifiableLinearOp& mHBt  = state.getModifiableOp("HBt");
    ModifiableLinearOp& mhatS = state.getModifiableOp("hatS");
    ModifiableLinearOp& BHBt  = state.getModifiableOp("BHBt");

    // build H*Bt
    subTotal_.start();
    mHBt = explicitMultiply(H, Bt, mHBt);
    subTotal_.stop();
    HBt = mHBt;

    // build B*H*Bt
    subTotal_.start();
    BHBt = explicitMultiply(B, HBt, BHBt);
    subTotal_.stop();

    // build C-B*H*Bt
    subTotal_.start();
    mhatS = explicitAdd(C, scale(-1.0, BHBt), mhatS);
    subTotal_.stop();
    hatS = mhatS;
  } else {
    // build an implicit Schur complement
    HBt = multiply(H, Bt);

    hatS = add(C, scale(-1.0, multiply(B, HBt)));
  }

  // time the application of HBt
  if (timed_HBt_ == Teuchos::null) {
    timed_HBt_ = Teuchos::rcp(new DiagnosticLinearOp(getOutputStream(), HBt, "HBt"));
  } else {
    timed_HBt_->setLinearOp(HBt);
  }

  // time the application of B
  if (timed_B_ == Teuchos::null) {
    timed_B_ = Teuchos::rcp(new DiagnosticLinearOp(getOutputStream(), B, "B"));
  } else {
    timed_B_->setLinearOp(B);
  }

  // build the inverse for F
  ModifiableLinearOp& invF = state.getModifiableOp("invF");
  subTotal_.start();
  if (invF == Teuchos::null) {
    invF = buildInverse(*invVelFactory_, F);

    timed_invF_ = Teuchos::rcp(new DiagnosticLinearOp(getOutputStream(), invF, "invF"));
  } else {
    rebuildInverse(*invVelFactory_, F, invF);

    timed_invF_->setLinearOp(invF);
  }
  subTotal_.stop();

  // build the approximate Schur complement
  ModifiableLinearOp& invS = state.getModifiableOp("invS");
  subTotal_.start();
  if (invS == Teuchos::null) {
    invS = buildInverse(*invPrsFactory_, hatS);

    timed_invS_ = Teuchos::rcp(new DiagnosticLinearOp(getOutputStream(), invS, "invS"));
  } else {
    rebuildInverse(*invPrsFactory_, hatS, invS);

    timed_invS_->setLinearOp(invS);
  }
  subTotal_.stop();

  std::vector<LinearOp> invDiag(2);  // vector storing inverses

  // build lower triangular inverse matrix
  BlockedLinearOp L = zeroBlockedOp(blockOp);
  setBlock(1, 0, L, timed_B_);
  endBlockFill(L);

  invDiag[0]    = timed_invF_;
  invDiag[1]    = timed_invS_;
  LinearOp invL = createBlockLowerTriInverseOp(L, invDiag);

  // build upper triangular matrix
  BlockedLinearOp U = zeroBlockedOp(blockOp);
  setBlock(0, 1, U, scale<double>(1.0 / alpha_, timed_HBt_.getConst()));
  endBlockFill(U);

  invDiag[0]    = identity(rangeSpace(invF));
  invDiag[1]    = scale(alpha_, identity(rangeSpace(invS)));
  LinearOp invU = createBlockUpperTriInverseOp(U, invDiag);

  // return implicit product operator
  Teko::LinearOp iU_t_iL = multiply(invU, invL, "SIMPLE_" + fApproxStr);

  // time the application of iU_t_iL
  if (timed_iU_t_iL_ == Teuchos::null)
    timed_iU_t_iL_ = Teuchos::rcp(new DiagnosticLinearOp(getOutputStream(), iU_t_iL, "iU_t_iL"));
  else
    timed_iU_t_iL_->setLinearOp(iU_t_iL);

  constrCount_++;

  constrTotal_.stop();

  return timed_iU_t_iL_;
}

//! Initialize from a parameter list
void TimingsSIMPLEPreconditionerFactory::initializeFromParameterList(
    const Teuchos::ParameterList& pl) {
  RCP<const InverseLibrary> invLib = getInverseLibrary();

  // default conditions
  useMass_        = false;
  customHFactory_ = Teuchos::null;
  fInverseType_   = Diagonal;

  // get string specifying inverse
  std::string invStr = "", invVStr = "", invPStr = "";
  alpha_ = 1.0;

  // "parse" the parameter list
  if (pl.isParameter("Inverse Type")) invStr = pl.get<std::string>("Inverse Type");
  if (pl.isParameter("Inverse Velocity Type"))
    invVStr = pl.get<std::string>("Inverse Velocity Type");
  if (pl.isParameter("Inverse Pressure Type"))
    invPStr = pl.get<std::string>("Inverse Pressure Type");
  if (pl.isParameter("Alpha")) alpha_ = pl.get<double>("Alpha");
  if (pl.isParameter("Explicit Velocity Inverse Type")) {
    std::string fInverseStr = pl.get<std::string>("Explicit Velocity Inverse Type");

    // build inverse types
    fInverseType_ = getDiagonalType(fInverseStr);
    if (fInverseType_ == NotDiag) customHFactory_ = invLib->getInverseFactory(fInverseStr);

    // Grab the sublist if we're using the block diagonal
    if (fInverseType_ == BlkDiag) BlkDiagList_ = pl.sublist("H options");
  }
  if (pl.isParameter("Use Mass Scaling")) useMass_ = pl.get<bool>("Use Mass Scaling");

  Teko_DEBUG_MSG_BEGIN(5) DEBUG_STREAM << "SIMPLE Parameters: " << std::endl;
  DEBUG_STREAM << "   inv type    = \"" << invStr << "\"" << std::endl;
  DEBUG_STREAM << "   inv v type  = \"" << invVStr << "\"" << std::endl;
  DEBUG_STREAM << "   inv p type  = \"" << invPStr << "\"" << std::endl;
  DEBUG_STREAM << "   alpha       = " << alpha_ << std::endl;
  DEBUG_STREAM << "   use mass    = " << useMass_ << std::endl;
  DEBUG_STREAM << "   vel scaling = " << getDiagonalName(fInverseType_) << std::endl;
  DEBUG_STREAM << "SIMPLE Parameter list: " << std::endl;
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

  //  two inverse factory objects
  RCP<InverseFactory> invVFact, invPFact;

  // build velocity inverse factory
  invVFact = invLib->getInverseFactory(invVStr);
  invPFact = invVFact;     // by default these are the same
  if (invVStr != invPStr)  // if different, build pressure inverse factory
    invPFact = invLib->getInverseFactory(invPStr);

  // based on parameter type build a strategy
  invVelFactory_ = invVFact;
  invPrsFactory_ = invPFact;

  if (useMass_) {
    Teuchos::RCP<Teko::RequestHandler> rh = getRequestHandler();
    rh->preRequest<Teko::LinearOp>(Teko::RequestMesg("Velocity Mass Matrix"));
    Teko::LinearOp mass = rh->request<Teko::LinearOp>(Teko::RequestMesg("Velocity Mass Matrix"));
    setMassMatrix(mass);
  }
}

//! For assiting in construction of the preconditioner
Teuchos::RCP<Teuchos::ParameterList> TimingsSIMPLEPreconditionerFactory::getRequestedParameters()
    const {
  Teuchos::RCP<Teuchos::ParameterList> result;
  Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

  // grab parameters from F solver
  RCP<Teuchos::ParameterList> vList = invVelFactory_->getRequestedParameters();
  if (vList != Teuchos::null) {
    Teuchos::ParameterList::ConstIterator itr;
    for (itr = vList->begin(); itr != vList->end(); ++itr) pl->setEntry(itr->first, itr->second);
    result = pl;
  }

  // grab parameters from S solver
  RCP<Teuchos::ParameterList> pList = invPrsFactory_->getRequestedParameters();
  if (pList != Teuchos::null) {
    Teuchos::ParameterList::ConstIterator itr;
    for (itr = pList->begin(); itr != pList->end(); ++itr) pl->setEntry(itr->first, itr->second);
    result = pl;
  }

  // grab parameters from S solver
  if (customHFactory_ != Teuchos::null) {
    RCP<Teuchos::ParameterList> hList = customHFactory_->getRequestedParameters();
    if (hList != Teuchos::null) {
      Teuchos::ParameterList::ConstIterator itr;
      for (itr = hList->begin(); itr != hList->end(); ++itr) pl->setEntry(itr->first, itr->second);
      result = pl;
    }
  }

  return result;
}

//! For assiting in construction of the preconditioner
bool TimingsSIMPLEPreconditionerFactory::updateRequestedParameters(
    const Teuchos::ParameterList& pl) {
  Teko_DEBUG_SCOPE("InvLSCStrategy::updateRequestedParameters", 10);
  bool result = true;

  // update requested parameters in solvers
  result &= invVelFactory_->updateRequestedParameters(pl);
  result &= invPrsFactory_->updateRequestedParameters(pl);
  if (customHFactory_ != Teuchos::null) result &= customHFactory_->updateRequestedParameters(pl);

  return result;
}

}  // end namespace NS
}  // end namespace Teko
