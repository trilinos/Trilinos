// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_LineSearch_Polynomial.H"

#include "NOX_LineSearch_Utils_Printing.H"
#include "NOX_SolverStats.hpp"
#include "NOX_LineSearch_Utils_Counters.H"
#include "NOX_LineSearch_Utils_Slope.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_MeritFunction_Generic.H"
#include "NOX_StatusTest_FiniteValue.H"
#include "NOX_GlobalData.H"

NOX::LineSearch::Polynomial::
Polynomial(const Teuchos::RCP<NOX::GlobalData>& gd,
       Teuchos::ParameterList& params) :
  globalDataPtr(gd),
  paramsPtr(NULL),
  print(gd->getUtils()),
  counter(&gd->getNonConstSolverStatistics()->lineSearch),
  slopeUtil(gd)
{
  reset(gd, params);
}

NOX::LineSearch::Polynomial::~Polynomial()
{

}

bool NOX::LineSearch::Polynomial::
reset(const Teuchos::RCP<NOX::GlobalData>& gd,
      Teuchos::ParameterList& params)
{
  globalDataPtr = gd;
  meritFuncPtr = gd->getMeritFunction();
  print.reset(gd->getUtils());
  counter = &gd->getNonConstSolverStatistics()->lineSearch;
  paramsPtr = &params;
  slopeUtil.reset(gd);

  Teuchos::ParameterList& p = params.sublist("Polynomial");

  std::string choice = p.get("Sufficient Decrease Condition", "Armijo-Goldstein");

  if (choice == "Armijo-Goldstein")
    suffDecrCond = ArmijoGoldstein;
  else if (choice == "Ared/Pred")
    suffDecrCond = AredPred;
  else if (choice == "None")
    suffDecrCond = None;
  else
  {
    print.err() << "NOX::LineSearch::Polynomial::reset - Invalid \"Sufficient Decrease Condition\"" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  choice = p.get("Interpolation Type", "Cubic");

  if (choice == "Cubic")
    interpolationType = Cubic;
  else if (choice == "Quadratic")
    interpolationType = Quadratic;
  else if (choice == "Quadratic3")
    interpolationType = Quadratic3;
  else
  {
    print.err() << "NOX::LineSearch::Polynomial::reset - Invalid \"Interpolation Type\"" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  choice = p.get("Recovery Step Type", "Constant");

  if (choice == "Constant")
    recoveryStepType = Constant;
  else if (choice == "Last Computed Step") {
    recoveryStepType = LastComputedStep;
  }
  else {
    print.err() << "NOX::LineSearch::Polynomial::reset - Invalid \"Recovery Step Type\"" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  minStep = p.get("Minimum Step", 1.0e-12);
  defaultStep = p.get("Default Step", 1.0);
  recoveryStep = p.get("Recovery Step", defaultStep);
  maxIters = p.get("Max Iters", 100);
  alpha = p.get("Alpha Factor", 1.0e-4);
  minBoundFactor = p.get("Min Bounds Factor", 0.1);
  maxBoundFactor = p.get("Max Bounds Factor", 0.5);
  doForceInterpolation = p.get("Force Interpolation", false);
  useCounter = p.get("Use Counters", true);
  maxIncreaseIter = p.get("Maximum Iteration for Increase", 0);
  maxRelativeIncrease = p.get("Allowed Relative Increase", 1.e2);

  // Is increase allowed?
  doAllowIncrease = (maxIncreaseIter > 0);

  // Set up counter
  if (useCounter)
    counter->reset();

  return true;
}

bool NOX::LineSearch::Polynomial::compute(Abstract::Group& newGrp,
                      double& step,
                      const Abstract::Vector& dir,
                      const Solver::Generic& s)
{
  printOpeningRemarks();

  int nNonlinearIters = s.getNumIterations();

  if (useCounter)
    counter->incrementNumLineSearches();

  // Get the linear solve tolerance if doing ared/pred for conv criteria
  std::string direction = const_cast<Teuchos::ParameterList&>(s.getList()).
    sublist("Direction").get("Method", "Newton");
  double eta = (suffDecrCond == AredPred) ?
    const_cast<Teuchos::ParameterList&>(s.getList()).
    sublist("Direction").sublist(direction).sublist("Linear Solver").
    get("Tolerance", -1.0) : 0.0;

  // Computations with old group
  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  double oldPhi = meritFuncPtr->computef(oldGrp);    // \phi(0)
  double oldValue = computeValue(oldGrp, oldPhi);
  double oldSlope = meritFuncPtr->computeSlope(dir, oldGrp);

  // Computations with new group
  step = defaultStep;
  updateGrp(newGrp, oldGrp, dir, step);
  double newPhi = meritFuncPtr->computef(newGrp);
  double newValue = computeValue(newGrp, newPhi);

  bool isConverged = false;
  bool isFailed = false;
  int nIters = 1;

  if (oldSlope >= 0.0)
  {
    printBadSlopeWarning(oldSlope);
    isFailed = true;
  }
  else
    isConverged = checkConvergence(newValue, oldValue, oldSlope, step,
                   eta, nIters, nNonlinearIters);

  // Increment the number of newton steps requiring a line search
  if ((useCounter) && (!isConverged))
    counter->incrementNumNonTrivialLineSearches();

  double prevPhi = 0.0;        // \phi(\lambda_{k-1})
  double prevPrevPhi = 0.0;    // \phi(\lambda_{k-2})
  double prevStep = 0.0;    // \lambda_{k-1}
  double prevPrevStep = 0.0;    // \lambda_{k-2}

  while ((!isConverged) && (!isFailed))
  {
    print.printStep(nIters, step, oldValue, newValue,
            "", (suffDecrCond != AredPred));

    if (nIters > maxIters)
    {
      isFailed = true;
      break;
    }

    if (interpolationType == Quadratic3)
    {
      /* 3-Point Quadratic Interpolation */

      prevPrevPhi = prevPhi;
      prevPhi = newPhi;
      prevPrevStep = prevStep;
      prevStep = step;

      if (nIters == 1)
      {
    step = 0.5 * step;
      }
      else
      {
    double c1 = prevStep * prevStep * (prevPrevPhi - oldPhi) -
      prevPrevStep * prevPrevStep * (prevPhi - oldPhi);
    double c2 = prevPrevStep * (prevPhi - oldPhi) -
      prevStep * (prevPrevPhi - oldPhi);

    if (c1 < 0)
      step = -0.5 * c1 / c2;
      }
    }

    else if ((nIters == 1) || (interpolationType == Quadratic))
    {
      /* Quadratic Interpolation */

      prevPhi = newPhi;
      prevStep = step;

      step = - (oldSlope * prevStep * prevStep) /
    (2.0 * (prevPhi - oldPhi - prevStep * oldSlope)) ;

    }

    else
    {
      /*   Cubic Interpolation */

      prevPrevPhi = prevPhi;
      prevPhi = newPhi;
      prevPrevStep = prevStep;
      prevStep = step;

      double term1 = prevPhi - oldPhi - prevStep * oldSlope ;
      double term2 = prevPrevPhi - oldPhi - prevPrevStep * oldSlope ;

      double a = 1.0 / (prevStep - prevPrevStep) *
    (term1 / (prevStep * prevStep) - term2 /
     (prevPrevStep * prevPrevStep)) ;

      double b = 1.0 / (prevStep - prevPrevStep) *
    (-1.0 * term1 * prevPrevStep / (prevStep * prevStep) +
     term2 * prevStep / (prevPrevStep * prevPrevStep)) ;

      double disc = b * b - 3.0 * a * oldSlope;

      if (disc < 0)
      {
    isFailed = true;
    break;
      }

      if (b > 0.0) // Check to prevent round off error (H. Walker)
      {
    step = -oldSlope / (b + sqrt(disc));
      }
      else
      {
    if (fabs(a) < 1.e-12) // check for when a is small
    {
      step = -oldSlope / (2.0 * b);
    }
    else
    {
      step = (-b + sqrt(disc))/ (3.0 * a);
    }
      }
    }

    // Apply bounds
    if (step < minBoundFactor * prevStep)
      step = minBoundFactor * prevStep;
    else if (step > maxBoundFactor * prevStep)
      step = maxBoundFactor * prevStep;

    // Check that step isn't too small
    if (step < minStep)
    {
      isFailed = true;
      break;
    }

    // Update the new group and compute new measures
    updateGrp(newGrp, oldGrp, dir, step);
    newPhi = meritFuncPtr->computef(newGrp);
    newValue = computeValue(newGrp, newPhi);

    nIters ++;

    if (useCounter)
      counter->incrementNumIterations();

    isConverged = checkConvergence(newValue, oldValue, oldSlope, step,
                   eta, nIters, nNonlinearIters);

  } // End while loop


  if (isFailed)
  {
    if (useCounter)
      counter->incrementNumFailedLineSearches();

    if (recoveryStepType == Constant)
      step = recoveryStep;

    if (step == 0.0)
    {
      newGrp = oldGrp;
      newPhi = oldPhi;
      newValue = oldValue;
    }
    else
    {
      updateGrp(newGrp, oldGrp, dir, step);
      newPhi = meritFuncPtr->computef(newGrp);
      newValue = computeValue(newGrp, newPhi);
    }
  }

  std::string message = (isFailed) ? "(USING RECOVERY STEP!)" : "(STEP ACCEPTED!)";
  print.printStep(nIters, step, oldValue, newValue, message, (suffDecrCond != AredPred));

  paramsPtr->set("Adjusted Tolerance", 1.0 - step * (1.0 - eta));

  if (useCounter)
    counter->setValues(*paramsPtr);

  return (!isFailed);
}

bool NOX::LineSearch::Polynomial::
checkConvergence(double newValue, double oldValue,
         double oldSlope,
         double step, double eta,
         int nIters,
         int nNonlinearIters) const
{
  NOX::StatusTest::FiniteValue checkNAN;
  if (checkNAN.finiteNumberTest(newValue) !=0)
    return false;

  if ((nIters == 1) && (doForceInterpolation))
    return false;

  if ((doAllowIncrease) && (nNonlinearIters <= maxIncreaseIter))
  {
    double relativeIncrease = newValue / oldValue;
    if (relativeIncrease < maxRelativeIncrease)
      return true;
  }

  bool returnVal = false;
  switch (suffDecrCond)
  {

  case ArmijoGoldstein:
    returnVal = (newValue <= oldValue + alpha * step * oldSlope);
    break;
  case AredPred:
    {
      double newEta = 1.0 - step * (1.0 - eta);
      returnVal = (newValue <= oldValue * (1.0 - alpha * (1.0 - newEta)));
      break;
    }
  case None:
    returnVal = true;
    break;
  default:

    print.err() << "NOX::LineSearch::Polynomial::isSufficientDecrease - Unknown convergence criteria" << std::endl;
    throw std::runtime_error("NOX Error");

  }
  return returnVal;
}

bool NOX::LineSearch::Polynomial::
updateGrp(NOX::Abstract::Group& newGrp,
      const NOX::Abstract::Group& oldGrp,
      const NOX::Abstract::Vector& dir,
      double step) const
{
  newGrp.computeX(oldGrp, dir, step);

  NOX::Abstract::Group::ReturnType status = newGrp.computeF();
  if (status != NOX::Abstract::Group::Ok)
    return false;

  return true;
}

double NOX::LineSearch::Polynomial::
computeValue(const NOX::Abstract::Group& grp, double phi)
{
  double value = phi;

  if (suffDecrCond == AredPred)
  {
    value = grp.getNormF();
  }

  return value;
}

void NOX::LineSearch::Polynomial::printOpeningRemarks() const
{
  if (print.isPrintType(NOX::Utils::InnerIteration))
  {
    print.out() << "\n" << NOX::Utils::fill(72) << "\n"
        << "-- Polynomial Line Search -- \n";
  }

  if (print.isPrintType(NOX::Utils::Details))
  {
    if (!Teuchos::is_null(meritFuncPtr))
      print.out() << "       Merit Function = " << meritFuncPtr->name()
          << std::endl;
  }
}

void NOX::LineSearch::Polynomial::printBadSlopeWarning(double slope) const
{
    print.out(NOX::Utils::Warning)
      << "WARNING: Computed slope is positive (slope = "
      << slope << ").\n" << "Using recovery step!" << std::endl;
}
