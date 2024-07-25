// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_LineSearch_SafeguardedStep.H"

#include "NOX_LineSearch_Utils_Printing.H"
#include "NOX_LineSearch_Utils_Slope.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_SolverStats.hpp"
#include "NOX_Solver_Generic.H"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "NOX_MeritFunction_Generic.H"
#include "NOX_StatusTest_FiniteValue.H"
#include "NOX_GlobalData.H"
#include <cmath>

NOX::LineSearch::SafeguardedStep::
SafeguardedStep(const Teuchos::RCP<NOX::GlobalData>& gd,
        Teuchos::ParameterList& params) :
  paramsPtr_(NULL),
  globalDataPtr_(gd),
  print_(gd->getUtils()),
  counter_(&gd->getNonConstSolverStatistics()->lineSearch)
{
  reset(gd, params);
}

bool NOX::LineSearch::SafeguardedStep::
reset(const Teuchos::RCP<NOX::GlobalData>& gd,
      Teuchos::ParameterList& params)
{
  globalDataPtr_ = gd;
  print_.reset(gd->getUtils());
  paramsPtr_ = &params;
  counter_ = &gd->getNonConstSolverStatistics()->lineSearch;

  Teuchos::ParameterList& p = params.sublist("Safeguarded Step");
  p.validateParametersAndSetDefaults(*getValidParameters());

  userLimits_ = p.get<Teuchos::RCP<NOX::Abstract::Vector> >("Update Limit Vector");
  lowerStepBound_ = p.get<double>("Step Size Lower Bound");
  upperStepBound_ = p.get<double>("Step Size Upper Bound");
  useCounter_ = p.get<bool>("Use Counters");

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(userLimits_), std::runtime_error,
                 "Error: The line search NOX::LineSearch::SafeguardedStep requires the user to supply the \"Update Limit Vector\" parameter of type NOX::Abstract::Vector for in the parameter list.");

  if (is_null(invLimits_))
    invLimits_ = userLimits_->clone(NOX::ShapeCopy);

  if (is_null(scaledUpdate_))
    scaledUpdate_ = userLimits_->clone(NOX::ShapeCopy);


  // Set up counter
  if (useCounter_)
    counter_->reset();

  return true;
}

bool NOX::LineSearch::SafeguardedStep::compute(Abstract::Group& newGrp,
                           double& step,
                           const Abstract::Vector& dir,
                           const Solver::Generic& /* s */)
{
  printOpeningRemarks();

  if (useCounter_) {
    counter_->incrementNumLineSearches();
    counter_->incrementNumNonTrivialLineSearches();
  }

  invLimits_->reciprocal(*userLimits_);
  *scaledUpdate_ = dir;
  scaledUpdate_->scale(*invLimits_);
  double infNorm = scaledUpdate_->norm(NOX::Abstract::Vector::MaxNorm);

  if (infNorm > 1.0)
    step = 1.0 / infNorm;
  else
    step = 1.0;

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
    print_.out () << "    computed step = " << step << std::endl;
  }

  step = std::max(step,lowerStepBound_);
  step = std::min(step,upperStepBound_);

  if (print_.isPrintType(NOX::Utils::InnerIteration))
    print_.out () << "    bounded step = " << step << std::endl;

  newGrp.computeX(newGrp,dir,step);

  if (print_.isPrintType(NOX::Utils::Details)) {
    // reuse invLimits_, will be reset above
    invLimits_->update(step,dir,0.0);
    print_.out () << "Final Step Scaled Update:" << std::endl;
    invLimits_->print(print_.out());
    print_.out () << "New Solution:" << std::endl;
    newGrp.getX().print(print_.out());
  }

  if (useCounter_)
    counter_->setValues(*paramsPtr_);

  return true;
}

Teuchos::RCP<const Teuchos::ParameterList>
NOX::LineSearch::SafeguardedStep::getValidParameters()
{
  if (is_null(validParams_)) {
    validParams_ = Teuchos::parameterList();
    validParams_->set("Step Size Lower Bound",1.0e-3);
    validParams_->set("Step Size Upper Bound",1.0);
    validParams_->set<Teuchos::RCP<NOX::Abstract::Vector> >("Update Limit Vector",Teuchos::null," User defined vector that contains the maximum allowed update for each entry in the direction vector of a nonlinear step.");
    validParams_->set("Use Counters",false);
  }

  return validParams_;
}

void NOX::LineSearch::SafeguardedStep::printOpeningRemarks() const
{
  if (print_.isPrintType(NOX::Utils::InnerIteration))
  {
    print_.out() << "\n" << NOX::Utils::fill(72) << "\n"
         << "-- SafeguardedStep Line Search -- \n";
  }
}
