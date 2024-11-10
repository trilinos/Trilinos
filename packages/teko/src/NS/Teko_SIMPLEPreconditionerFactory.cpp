// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_SIMPLEPreconditionerFactory.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"
#include <stdexcept>
#ifdef TEKO_HAVE_EPETRA
#include "Teko_DiagonalPreconditionerFactory.hpp"
#endif

#include "Teuchos_Time.hpp"

using Teuchos::RCP;

namespace Teko {
namespace NS {

// Constructor definition
SIMPLEPreconditionerFactory ::SIMPLEPreconditionerFactory(const RCP<InverseFactory>& inverse,
                                                          double alpha)
    : invVelFactory_(inverse),
      invPrsFactory_(inverse),
      alpha_(alpha),
      fInverseType_(Diagonal),
      useMass_(false) {}

SIMPLEPreconditionerFactory ::SIMPLEPreconditionerFactory(const RCP<InverseFactory>& invVFact,
                                                          const RCP<InverseFactory>& invPFact,
                                                          double alpha)
    : invVelFactory_(invVFact),
      invPrsFactory_(invPFact),
      alpha_(alpha),
      fInverseType_(Diagonal),
      useMass_(false) {}

SIMPLEPreconditionerFactory::SIMPLEPreconditionerFactory()
    : alpha_(1.0), fInverseType_(Diagonal), useMass_(false) {}

// Use the factory to build the preconditioner (this is where the work goes)
LinearOp SIMPLEPreconditionerFactory ::buildPreconditionerOperator(
    BlockedLinearOp& blockOp, BlockPreconditionerState& state) const {
  Teko_DEBUG_SCOPE("SIMPLEPreconditionerFactory::buildPreconditionerOperator", 10);
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

    /*
         // Get a FECrsMarix out of the BDP
         RCP<Epetra_FECrsMatrix> Hcrs=rcp(Hstate.BDP_->CreateFECrsMatrix());
         H=Thyra::epetraLinearOp(Hcrs);
    */

    buildExplicitSchurComplement = true;  // NTS: Do I need this?
                                          // Answer - no, but it is documenting whats going on here.
#else
    throw std::logic_error(
        "SIMPLEPreconditionerFactory fInverseType_ == "
        "BlkDiag but EPETRA is turned off!");
#endif

  } else {
    // get generic diagonal
    H          = getInvDiagonalOp(matF, fInverseType_);
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
    mHBt = explicitMultiply(H, Bt, mHBt);
    HBt  = mHBt;

    // build B*H*Bt
    BHBt = explicitMultiply(B, HBt, BHBt);

    // build C-B*H*Bt
    mhatS = explicitAdd(C, scale(-1.0, BHBt), mhatS);
    hatS  = mhatS;
  } else {
    // build an implicit Schur complement
    HBt = multiply(H, Bt);

    hatS = add(C, scale(-1.0, multiply(B, HBt)));
  }

  Teko::ModifiableLinearOp& precInvF = state.getModifiableOp("precInvF");
  if (precVelFactory_) {
    if (precInvF == Teuchos::null) {
      precInvF = precVelFactory_->buildInverse(F);
      state.addModifiableOp("precInvF", precInvF);
    } else {
      Teko::rebuildInverse(*precVelFactory_, F, precInvF);
    }
  }

  // build the inverse for F
  Teko::ModifiableLinearOp& invF = state.getModifiableOp("invF");
  if (invF == Teuchos::null) {
    if (precInvF.is_null()) {
      invF = Teko::buildInverse(*invVelFactory_, F);
    } else {
      invF = Teko::buildInverse(*invVelFactory_, F, precInvF);
    }
  } else {
    if (precInvF.is_null()) {
      Teko::rebuildInverse(*invVelFactory_, F, invF);
    } else {
      Teko::rebuildInverse(*invVelFactory_, F, precInvF, invF);
    }
  }

  Teko::ModifiableLinearOp& precInvS = state.getModifiableOp("precInvS");
  if (precPrsFactory_) {
    if (precInvS == Teuchos::null) {
      precInvS = precPrsFactory_->buildInverse(hatS);
      state.addModifiableOp("precInvS", precInvS);
    } else {
      Teko::rebuildInverse(*precPrsFactory_, hatS, precInvS);
    }
  }

  // build the approximate Schur complement
  Teko::ModifiableLinearOp& invS = state.getModifiableOp("invS");
  if (invS == Teuchos::null) {
    if (precInvS == Teuchos::null) {
      invS = Teko::buildInverse(*invPrsFactory_, hatS);
    } else {
      invS = Teko::buildInverse(*invPrsFactory_, hatS, precInvS);
    }
  } else {
    if (precInvS == Teuchos::null) {
      Teko::rebuildInverse(*invPrsFactory_, hatS, invS);
    } else {
      Teko::rebuildInverse(*invPrsFactory_, hatS, precInvS, invS);
    }
  }

  std::vector<LinearOp> invDiag(2);  // vector storing inverses

  // build lower triangular inverse matrix
  BlockedLinearOp L = zeroBlockedOp(blockOp);
  setBlock(1, 0, L, B);
  endBlockFill(L);

  invDiag[0]    = invF;
  invDiag[1]    = invS;
  LinearOp invL = createBlockLowerTriInverseOp(L, invDiag);

  // build upper triangular matrix
  BlockedLinearOp U = zeroBlockedOp(blockOp);
  setBlock(0, 1, U, scale(1.0 / alpha_, HBt));
  endBlockFill(U);

  invDiag[0]    = identity(rangeSpace(invF));
  invDiag[1]    = scale(alpha_, identity(rangeSpace(invS)));
  LinearOp invU = createBlockUpperTriInverseOp(U, invDiag);

  // return implicit product operator
  return multiply(invU, invL, "SIMPLE_" + fApproxStr);
}

//! Initialize from a parameter list
void SIMPLEPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList& pl) {
  RCP<const InverseLibrary> invLib = getInverseLibrary();

  // default conditions
  useMass_        = false;
  customHFactory_ = Teuchos::null;
  fInverseType_   = Diagonal;

  // get string specifying inverse
  std::string invStr = "", invVStr = "", invPStr = "", precVStr = "", precPStr = "";
  alpha_ = 1.0;

  // "parse" the parameter list
  if (pl.isParameter("Inverse Type")) invStr = pl.get<std::string>("Inverse Type");
  if (pl.isParameter("Inverse Velocity Type"))
    invVStr = pl.get<std::string>("Inverse Velocity Type");
  if (pl.isParameter("Preconditioner Velocity Type"))
    precVStr = pl.get<std::string>("Preconditioner Velocity Type");
  if (pl.isParameter("Inverse Pressure Type"))
    invPStr = pl.get<std::string>("Inverse Pressure Type");
  if (pl.isParameter("Preconditioner Pressure Type"))
    precPStr = pl.get<std::string>("Preconditioner Pressure Type");
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
  DEBUG_STREAM << "   prec v type  = \"" << precVStr << "\"" << std::endl;
  DEBUG_STREAM << "   inv p type  = \"" << invPStr << "\"" << std::endl;
  DEBUG_STREAM << "   prec p type  = \"" << precPStr << "\"" << std::endl;
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

  RCP<InverseFactory> precVFact, precPFact;
  if (precVStr != "") precVFact = invLib->getInverseFactory(precVStr);

  if (precPStr != "") precPFact = invLib->getInverseFactory(precPStr);

  // based on parameter type build a strategy
  invVelFactory_ = invVFact;
  invPrsFactory_ = invPFact;

  precVelFactory_ = precVFact;
  precPrsFactory_ = precPFact;

  if (useMass_) {
    Teuchos::RCP<Teko::RequestHandler> rh = getRequestHandler();
    rh->preRequest<Teko::LinearOp>(Teko::RequestMesg("Velocity Mass Matrix"));
    Teko::LinearOp mass = rh->request<Teko::LinearOp>(Teko::RequestMesg("Velocity Mass Matrix"));
    setMassMatrix(mass);
  }
}

//! For assiting in construction of the preconditioner
Teuchos::RCP<Teuchos::ParameterList> SIMPLEPreconditionerFactory::getRequestedParameters() const {
  Teuchos::RCP<Teuchos::ParameterList> result;
  Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

  // grab parameters from F solver
  {
    RCP<Teuchos::ParameterList> vList = invVelFactory_->getRequestedParameters();
    if (vList != Teuchos::null) {
      Teuchos::ParameterList::ConstIterator itr;
      for (itr = vList->begin(); itr != vList->end(); ++itr) pl->setEntry(itr->first, itr->second);
      result = pl;
    }
  }

  if (precVelFactory_ != Teuchos::null) {
    RCP<Teuchos::ParameterList> vList = precVelFactory_->getRequestedParameters();
    if (vList != Teuchos::null) {
      Teuchos::ParameterList::ConstIterator itr;
      for (itr = vList->begin(); itr != vList->end(); ++itr) pl->setEntry(itr->first, itr->second);
      result = pl;
    }
  }

  // grab parameters from S solver
  {
    RCP<Teuchos::ParameterList> pList = invPrsFactory_->getRequestedParameters();
    if (pList != Teuchos::null) {
      Teuchos::ParameterList::ConstIterator itr;
      for (itr = pList->begin(); itr != pList->end(); ++itr) pl->setEntry(itr->first, itr->second);
      result = pl;
    }
  }

  if (precPrsFactory_ != Teuchos::null) {
    RCP<Teuchos::ParameterList> pList = precPrsFactory_->getRequestedParameters();
    if (pList != Teuchos::null) {
      Teuchos::ParameterList::ConstIterator itr;
      for (itr = pList->begin(); itr != pList->end(); ++itr) pl->setEntry(itr->first, itr->second);
      result = pl;
    }
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
bool SIMPLEPreconditionerFactory::updateRequestedParameters(const Teuchos::ParameterList& pl) {
  Teko_DEBUG_SCOPE("InvLSCStrategy::updateRequestedParameters", 10);
  bool result = true;

  // update requested parameters in solvers
  result &= invVelFactory_->updateRequestedParameters(pl);
  result &= invPrsFactory_->updateRequestedParameters(pl);
  if (precVelFactory_) result &= precVelFactory_->updateRequestedParameters(pl);
  if (precPrsFactory_) result &= precPrsFactory_->updateRequestedParameters(pl);
  if (customHFactory_ != Teuchos::null) result &= customHFactory_->updateRequestedParameters(pl);

  return result;
}

}  // end namespace NS
}  // end namespace Teko
