// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Teko includes
#include "Teko_DiagnosticPreconditionerFactory.hpp"

#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_DiagnosticLinearOp.hpp"

#include "Teuchos_TimeMonitor.hpp"

namespace Teko {

//! Default constructor, for use with the AutoClone class.
DiagnosticPreconditionerFactory::DiagnosticPreconditionerFactory()
    : outputStream_(Teko::getOutputStream()),
      invFactory_(Teuchos::null),
      diagString_("<label me!>"),
      printResidual_(false) {}

/** Construct a preconditioner factory that applies a specified
 * preconditioner, a fixed number of times.
 */
DiagnosticPreconditionerFactory::DiagnosticPreconditionerFactory(
    const Teuchos::RCP<Teko::InverseFactory>& invFactory, const std::string& label,
    const Teuchos::RCP<std::ostream>& os, bool printResidual)
    : outputStream_(Teko::getOutputStream()),
      invFactory_(invFactory),
      precFactory_(Teuchos::null),
      diagString_(label),
      printResidual_(printResidual) {
  initTimers(diagString_);

  if (os != Teuchos::null) outputStream_ = os;
}

/** Construct a preconditioner factory that applies a specified
 * preconditioned solver, a fixed number of times.
 */
DiagnosticPreconditionerFactory::DiagnosticPreconditionerFactory(
    const Teuchos::RCP<Teko::InverseFactory>& invFactory,
    const Teuchos::RCP<Teko::InverseFactory>& precFactory, const std::string& label,
    const Teuchos::RCP<std::ostream>& os, bool printResidual)
    : outputStream_(Teko::getOutputStream()),
      invFactory_(invFactory),
      precFactory_(precFactory),
      diagString_(label),
      printResidual_(printResidual) {
  initTimers(diagString_);

  if (os != Teuchos::null) outputStream_ = os;
}

double DiagnosticPreconditionerFactory::totalInitialBuildTime() const {
  return buildTimer_->totalElapsedTime() + precBuildTimer_->totalElapsedTime();
}

double DiagnosticPreconditionerFactory::totalRebuildTime() const {
  return rebuildTimer_->totalElapsedTime() + precRebuildTimer_->totalElapsedTime();
}

DiagnosticPreconditionerFactory::~DiagnosticPreconditionerFactory() {
  // check timers for null
  if (buildTimer_ == Teuchos::null || rebuildTimer_ == Teuchos::null) {
    // (*outputStream_) << "DiagnosticPreconditionerFactory \"" << diagString_ << "\": "
    //                  << "Timers not initialized" << std::endl;

    return;
  }

  double initBuildTime = totalInitialBuildTime();
  int initBuilds       = numInitialBuilds();

  double initRebuildTime = totalRebuildTime();
  int initRebuilds       = numRebuilds();

  (*outputStream_) << "DiagnosticPreconditionerFactory \"" << diagString_ << "\":\n";

  // print build string
  (*outputStream_) << "   build elapsed = " << initBuildTime << ", "
                   << "num builds = " << initBuilds << ", ";
  if (initBuilds > 0)
    (*outputStream_) << "timer/app = " << initBuildTime / double(initBuilds) << "\n";
  else
    (*outputStream_) << "timer/app = "
                     << "none"
                     << "\n";

  // print rebuild string
  (*outputStream_) << "   rebuild elapsed = " << initRebuildTime << ", "
                   << "num rebuilds = " << initRebuilds << ", ";
  if (initRebuilds > 0)
    (*outputStream_) << "timer/app = " << initRebuildTime / double(initRebuilds) << "\n";
  else
    (*outputStream_) << "timer/app = "
                     << "none"
                     << "\n";

  // print total string
  (*outputStream_) << "   total elapsed = " << initRebuildTime + initBuildTime << ", "
                   << "num rebuilds = " << initRebuilds + initBuilds << ", ";
  if (initBuilds + initRebuilds > 0)
    (*outputStream_) << "timer/app = "
                     << (initRebuildTime + initBuildTime) / double(initRebuilds + initBuilds)
                     << std::endl;
  else
    (*outputStream_) << "timer/app = "
                     << "none" << std::endl;
}

/** \brief Function that is called to build the preconditioner
 *        for the linear operator that is passed in.
 */
LinearOp DiagnosticPreconditionerFactory::buildPreconditionerOperator(
    LinearOp& lo, PreconditionerState& state) const {
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  TEUCHOS_TEST_FOR_EXCEPTION(
      invFactory_ == Teuchos::null, std::runtime_error,
      "ERROR: Teko::DiagnosticPreconditionerFactory::buildPreconditionerOperator requires that an "
          << "inverse factory has been set. Currently it is null!");

  TEUCHOS_TEST_FOR_EXCEPTION(
      buildTimer_ == Teuchos::null || rebuildTimer_ == Teuchos::null, std::runtime_error,
      "ERROR: Teko::DiagnosticPreconditionerFactory::buildPreconditionerOperator requires that "
          << "the timers be initialized. Currently they are null! (label = \"" << diagString_
          << "\")");

  // build user specified preconditioner
  ModifiableLinearOp& diagOp_ptr      = state.getModifiableOp("diagnosticOp");
  ModifiableLinearOp& diagOp_prec_ptr = state.getModifiableOp("prec_diagnosticOp");

  if (precFactory_ != Teuchos::null) {
    if (diagOp_prec_ptr == Teuchos::null) {
      {
        // start timer on construction, end on destruction
        Teuchos::TimeMonitor monitor(*precBuildTimer_, false);

        diagOp_prec_ptr = precFactory_->buildInverse(lo);
      }

      state.addModifiableOp("prec_diagnosticOp", diagOp_prec_ptr);
    } else {
      Teuchos::TimeMonitor monitor(*precRebuildTimer_, false);
      Teko::rebuildInverse(*precFactory_, lo, diagOp_prec_ptr);
    }
  }

  if (diagOp_ptr == Teuchos::null) {
    ModifiableLinearOp invOp;
    {
      // start timer on construction, end on destruction
      Teuchos::TimeMonitor monitor(*buildTimer_, false);

      if (diagOp_prec_ptr.is_null())
        invOp = Teko::buildInverse(*invFactory_, lo);
      else
        invOp = Teko::buildInverse(*invFactory_, lo, diagOp_prec_ptr);
    }

    // only printing residual requires use of forward operator
    if (printResidual_)
      diagOp_ptr = createDiagnosticLinearOp(outputStream_, lo, invOp, diagString_);
    else
      diagOp_ptr = createDiagnosticLinearOp(outputStream_, invOp, diagString_);
  } else {
    RCP<DiagnosticLinearOp> diagOp = rcp_dynamic_cast<DiagnosticLinearOp>(diagOp_ptr);

    // only printing residual requires use of forward operator
    if (printResidual_) diagOp->setForwardOp(lo);

    ModifiableLinearOp invOp = diagOp->getModifiableOp();
    {
      // start timer on construction, end on destruction
      Teuchos::TimeMonitor monitor(*rebuildTimer_, false);

      if (diagOp_prec_ptr.is_null())
        Teko::rebuildInverse(*invFactory_, lo, invOp);
      else
        Teko::rebuildInverse(*invFactory_, lo, diagOp_prec_ptr, invOp);
    }
  }

  return diagOp_ptr.getConst();
}

/** \brief This function builds the internals of the preconditioner factory
 *        from a parameter list.
 */
void DiagnosticPreconditionerFactory::initializeFromParameterList(
    const Teuchos::ParameterList& settings) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      not settings.isParameter("Inverse Factory"), std::runtime_error,
      "Parameter \"Inverse Factory\" is required by a Teko::DiagnosticPreconditionerFactory");
  TEUCHOS_TEST_FOR_EXCEPTION(
      not settings.isParameter("Descriptive Label"), std::runtime_error,
      "Parameter \"Descriptive Label\" is required by a Teko::DiagnosticPreconditionerFactory");

  // grab library and preconditioner name
  std::string invName  = settings.get<std::string>("Inverse Factory");
  std::string precName = "";
  if (settings.isParameter("Preconditioner Factory"))
    precName = settings.get<std::string>("Preconditioner Factory");
  diagString_ = settings.get<std::string>("Descriptive Label");

  // build preconditioner factory
  Teuchos::RCP<const InverseLibrary> il = getInverseLibrary();
  invFactory_                           = il->getInverseFactory(invName);
  if (precName != "") precFactory_ = il->getInverseFactory(precName);
  TEUCHOS_TEST_FOR_EXCEPTION(invFactory_ == Teuchos::null, std::runtime_error,
                             "ERROR: \"Inverse Factory\" = " << invName << " could not be found");

  if (settings.isParameter("Print Residual")) printResidual_ = settings.get<bool>("Print Residual");

  // build timers to use
  initTimers(diagString_);
}

/** \brief Request the additional parameters this preconditioner factory
 *        needs.
 */
Teuchos::RCP<Teuchos::ParameterList> DiagnosticPreconditionerFactory::getRequestedParameters()
    const {
  TEUCHOS_TEST_FOR_EXCEPTION(
      invFactory_ == Teuchos::null, std::runtime_error,
      "ERROR: Teko::DiagnosticPreconditionerFactory::getRequestedParameters requires that a "
          << "preconditioner factory has been set. Currently it is null!");

  return invFactory_->getRequestedParameters();
}

/** \brief Update this object with the fields from a parameter list.
 */
bool DiagnosticPreconditionerFactory::updateRequestedParameters(const Teuchos::ParameterList& pl) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      invFactory_ == Teuchos::null, std::runtime_error,
      "ERROR: Teko::DiagnosticPreconditionerFactory::updateRequestedParameters requires that a "
          << "preconditioner factory has been set. Currently it is null!");

  bool success = true;
  success &= invFactory_->updateRequestedParameters(pl);
  if (precFactory_) success &= precFactory_->updateRequestedParameters(pl);

  return success;
}

void DiagnosticPreconditionerFactory::initTimers(const std::string& str) {
  buildTimer_       = Teuchos::rcp(new Teuchos::Time(str + " buildTimer"));
  rebuildTimer_     = Teuchos::rcp(new Teuchos::Time(str + " rebuildTimer"));
  precBuildTimer_   = Teuchos::rcp(new Teuchos::Time(str + " precBuildTimer"));
  precRebuildTimer_ = Teuchos::rcp(new Teuchos::Time(str + " precRebuildTimer"));
}

}  // end namespace Teko
