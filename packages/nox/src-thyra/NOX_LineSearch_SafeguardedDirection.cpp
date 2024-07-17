// $Id$
// $Source$

// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_LineSearch_SafeguardedDirection.hpp"

#include "NOX_LineSearch_Utils_Printing.H"
#include "NOX_LineSearch_Utils_Counters.H"
#include "NOX_LineSearch_Utils_Slope.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "NOX_MeritFunction_Generic.H"
#include "NOX_StatusTest_FiniteValue.H"
#include "NOX_GlobalData.H"
#include "NOX_SolverStats.hpp"
#include "NOX_TOpEleWiseMinSwap.hpp"
#include "NOX_Thyra_Vector.H"
#include "Thyra_VectorBase.hpp"
#include <cmath>

NOX::LineSearch::SafeguardedDirection::
SafeguardedDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
             Teuchos::ParameterList& params) :
  paramsPtr_(NULL),
  globalDataPtr_(gd),
  print_(gd->getUtils()),
  counter_(&gd->getNonConstSolverStatistics()->lineSearch)
{
  reset(gd, params);
}

bool NOX::LineSearch::SafeguardedDirection::
reset(const Teuchos::RCP<NOX::GlobalData>& gd,
      Teuchos::ParameterList& params)
{
  globalDataPtr_ = gd;
  print_.reset(gd->getUtils());
  counter_ = &gd->getNonConstSolverStatistics()->lineSearch;
  paramsPtr_ = &params;

  Teuchos::ParameterList& p = params.sublist("Safeguarded Direction");
  p.validateParametersAndSetDefaults(*getValidParameters());

  userLimits_ = p.get<Teuchos::RCP<NOX::Abstract::Vector> >("Update Limit Vector");
  useCounter_ = p.get<bool>("Use Counters");

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(userLimits_), std::runtime_error,
                 "Error: The line search NOX::LineSearch::SafeguardedDirection requires the user to supply the \"Update Limit Vector\" parameter of type NOX::Abstract::Vector for in the parameter list.");

  if (is_null(limitDifference_))
    limitDifference_ = userLimits_->clone(NOX::ShapeCopy);

  if (is_null(scaledUpdate_))
    scaledUpdate_ = userLimits_->clone(NOX::ShapeCopy);

  // Set up counter
  if (useCounter_)
    counter_->reset();

  return true;
}

bool NOX::LineSearch::SafeguardedDirection::compute(Abstract::Group& newGrp,
                            double& step,
                            const Abstract::Vector& dir,
                            const Solver::Generic& /* s */)
{
  printOpeningRemarks();

  if (useCounter_) {
    counter_->incrementNumLineSearches();
    counter_->incrementNumNonTrivialLineSearches();
  }

  // Limit individual entries
  *scaledUpdate_ = dir;
  Teuchos::RCP< ::Thyra::VectorBase<double> > thyraUserLimits = Teuchos::rcp_dynamic_cast<NOX::Thyra::Vector>(userLimits_,true)->getThyraRCPVector();
  Teuchos::RCP< ::Thyra::VectorBase<double> > thyraScaledUpdate = Teuchos::rcp_dynamic_cast<NOX::Thyra::Vector>(scaledUpdate_,true)->getThyraRCPVector();
  ::Thyra::ele_wise_min_swap(*thyraUserLimits,thyraScaledUpdate.ptr());

  // Now compute an equivalent "step" for this limiting (used by
  // status test to make sure no limiting is going on at the converged
  // solution)
  limitDifference_->update(1.0,dir,-1.0,*scaledUpdate_,0.0);
  double dirNorm = dir.norm(NOX::Abstract::Vector::TwoNorm);
  // Protect against divide by zero corner case (for a perfectly converged solution)
  if (dirNorm < 100.0 * Teuchos::ScalarTraits<double>::eps())
    step = 1.0;
  else
    step = 1.0 - limitDifference_->norm(NOX::Abstract::Vector::TwoNorm) / dir.norm(NOX::Abstract::Vector::TwoNorm);

  if (print_.isPrintType(NOX::Utils::Details)) {
    print_.out () << "userLimits_:" << std::endl;
    userLimits_->print(print_.out());
    print_.out () << "scaledUpdate_:" << std::endl;
    scaledUpdate_->print(print_.out());
    print_.out () << "dir:" << std::endl;
    dir.print(print_.out());
    print_.out () << "Old Solution:" << std::endl;
    newGrp.getX().print(print_.out());
  }

  if (print_.isPrintType(NOX::Utils::InnerIteration)) {
    print_.out () << "    approximated equivalent step = " << step << std::endl;
  }

  newGrp.computeX(newGrp,*scaledUpdate_,1.0);

  if (print_.isPrintType(NOX::Utils::Details)) {
   print_.out () << "New Solution:" << std::endl;
    newGrp.getX().print(print_.out());
  }

  if (useCounter_)
    counter_->setValues(*paramsPtr_);

  return true;
}

Teuchos::RCP<const Teuchos::ParameterList>
NOX::LineSearch::SafeguardedDirection::getValidParameters()
{
  if (is_null(validParams_)) {
    validParams_ = Teuchos::parameterList();
    validParams_->set<Teuchos::RCP<NOX::Abstract::Vector> >("Update Limit Vector",Teuchos::null," User defined vector that contains the maximum allowed update for each entry in the direction vector of a nonlinear step.");
    validParams_->set("Use Counters",false);
  }

  return validParams_;
}

void NOX::LineSearch::SafeguardedDirection::printOpeningRemarks() const
{
  if (print_.isPrintType(NOX::Utils::InnerIteration))
  {
    print_.out() << "\n" << NOX::Utils::fill(72) << "\n"
         << "-- SafeguardedDirection Line Search -- \n";
  }
}
