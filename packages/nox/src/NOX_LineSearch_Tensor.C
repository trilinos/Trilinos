// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*  Notes:
**
**  (.) The TensorBasedTest solver is virtually identical to the
**  LineSearchBased solver.  Thus, maybe at some point I should remove
**  this TensorBasedTest solver and convert LineSearchBased solver and the
**  LineSearch classes to work with the tensor linesearches.  Right
**  now I see 2 options for this conversion: Make optional argument in
**  LineSearch::compute to allow for a direction argument or add
**  getDirection method to solver so that linesearch object can
**  compute the curvilinear linesearch.  The latter option might have
**  trouble using the const direction.  Need to investigate...
**     //NOX::Abstract::Vector dir2 = dir.clone(ShapeCopy);
**     //const NOX::Direction::Tensor& direction = s.getDirection();
**
**  (.)  Should change to *sufficient* decrease condition instead of
**  just "fprime<0"
**
**  (.)  Maybe move the test of full step into compute instead of in
**  performLinesearch.  However, this might cause trouble with
**  counters and other things.
**
**  (.)  In the dual linesearch, it is checking both full steps and
**  taking the best of either.  This is different from TENSOLVE.
**
**  (.)  Old comment says:
**  "// Note that for Newton direction, fprime = -2.0*oldf"
**  Is this really true?
*/

#include "NOX_Common.H"

#ifdef WITH_PRERELEASE

#include "NOX_LineSearch_Tensor.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"
#include "NOX_Direction_Tensor.H"
#include "NOX_Solver_TensorBasedTest.H"
#include "NOX_LineSearch_Utils_Counters.H"
#include "NOX_SolverStats.hpp"

NOX::LineSearch::Tensor::
Tensor(const Teuchos::RCP<NOX::GlobalData>& gd,
       Teuchos::ParameterList& params) :
  globalDataPtr(gd),
  paramsPtr(NULL),
  print(gd->getUtils()),
  counter(&gd->getNonConstSolverStatistics()->lineSearch),
  slopeObj(gd)
{
  //  reset(paramsPtr->sublist("Line Search"));
  reset(gd, params);
}

NOX::LineSearch::Tensor::~Tensor()
{
  printf("multsJv = %d   (linesearch)\n", multsJv);
}

bool NOX::LineSearch::Tensor::
reset(const Teuchos::RCP<NOX::GlobalData>& gd,
      Teuchos::ParameterList& lsParams)
{
  globalDataPtr = gd;
  utils = *(gd->getUtils());
  print.reset(gd->getUtils());
  counter = &gd->getNonConstSolverStatistics()->lineSearch;
  paramsPtr = &lsParams;
  slopeObj.reset(gd);

  multsJv = 0;

  // Determine the specific type of tensor linesearch to perform
  std::string choice = lsParams.get("Method", "Curvilinear");

  utils.out() << choice << std::endl;

  if (choice == "Curvilinear")
    lsType = Curvilinear;
  else if (choice == "Dual")
    lsType = Dual;
  else if (choice == "Standard")
    lsType = Standard;
  //  else if (choice == "Full Step")
  //    lsType = FullStep;
  else if (choice == "Newton")
    lsType = Newton;
  else
  {
    if (utils.isPrintType(NOX::Utils::Error))
      utils.err() << "NOX::Direction::Tensor::reset() - The choice of "
       << "\"Line Search\" parameter " << choice
       << " is invalid." << std::endl;
    throw std::runtime_error("NOX Error");
  }
  //  Copy Method into "Submethod" (temporary hack for data scripts)
  lsParams.set("Submethod", choice);

  // Make a reference to the sublist holding the global strategy parameters
  Teuchos::ParameterList& gsParams = lsParams.sublist(choice);

#ifdef CODE_FROM_TENSORBASED
  // Decide what step to use in case of linesearch failure
  choice = gsParams.get("Recovery Step Type", "Constant");
  if (choice == "Constant")
    recoveryStepType = Constant;          // Use value in "Recovery Step"
  else if (choice == "Last Computed Step")
    recoveryStepType = LastComputedStep;  // Use last step from linesearch
  else
  {
    utils.err() << "NOX::Solver::TensorBased::reset() - "
     << "Invalid \"Recovery Step Type\"" << std::endl;
    throw std::runtime_error("NOX Error");
  }
#endif

  // Initialize linesearch parameters for this object
  minStep = gsParams.get("Minimum Step", 1.0e-12);
  defaultStep = gsParams.get("Default Step", 1.0);
  recoveryStep = gsParams.get("Recovery Step", 0.0); // exit on fail
  maxIters = gsParams.get("Max Iters", 40);
  alpha = gsParams.get("Alpha Factor", 1.0e-4);

  choice = gsParams.get("Lambda Selection", "Halving");
  if (choice == "Halving")
    lambdaSelection = Halving;
  else if (choice == "Quadratic")
    lambdaSelection = Quadratic;
  else
  {
    if (utils.isPrintType(NOX::Utils::Error))
      utils.err() << "NOX::Solver::TensorBased::reset() - The choice of "
       << "\"Lambda Selection\" parameter " << choice
       << " is invalid." << std::endl;
    throw std::runtime_error("NOX Error");
  }

  choice = gsParams.get("Sufficient Decrease Condition",
                 "Armijo-Goldstein");
  if (choice == "Armijo-Goldstein")
    suffDecrCond = ArmijoGoldstein;     // This is the only one implemented
  else if (choice == "Ared/Pred")
    suffDecrCond = AredPred;
  else if (choice == "None")
    suffDecrCond = None;
  else
  {
    if (utils.isPrintType(NOX::Utils::Error))
      utils.err() << "NOX::Solver::TensorBased::reset() - The choice of "
       << "\"Sufficient Decrease Condition\" parameter " << choice
       << " is invalid." << std::endl;
    throw std::runtime_error("NOX Error");
  }


#ifdef OLD_CODE
  // Initialize linesearch parameters for this object
  minStep = lsparams.get("Minimum Step", 1.0e-12);
  defaultStep = lsparams.get("Default Step", 1.0);
  recoveryStep = lsparams.get("Recovery Step", 0.0); // exit on fail
  maxIters = lsparams.get("Max Iters", 40);
  alpha = lsparams.get("Alpha Factor", 1.0e-4);
  paramsPtr = &params;

  // Do line search and compute new soln.
  std::string choice = lsparams.get("Submethod", "Curvilinear");

  if (choice == "Curvilinear")
    lsType = Curvilinear;
  else if (choice == "Dual")
    lsType = Dual;
  else if (choice == "Standard")
    lsType = Standard;
  else if (choice == "Newton")
    lsType = Newton;
  else {
    if (utils.isPrintType(NOX::Utils::Warning)) {
      utils.out() << "Warning: NOX::Direction::Tensor::reset() - the choice of "
       << "\"Line Search\" \nparameter is invalid.  Using curvilinear "
       << "line search." << std::endl;
    }
    lsparams.set("Submethod", "Curvilinear");
    lsType = Curvilinear;
  }


  choice = lsparams.get("Lambda Selection", "Halving");
  if (choice == "Halving") {
    lambdaSelection = Halving;
  }
  else if (choice == "Quadratic") {
    lambdaSelection = Quadratic;
  }
  else {
    utils.out() << "Warning: NOX::Solver::TensorBasedTest::init() - the choice of "
     << "\"Lambda Selection\" parameter is invalid." << std::endl;
    lambdaSelection = Halving;
  }


  choice = lsparams.get("Sufficient Decrease Condition",
                 "Armijo-Goldstein");
  if (choice == "Ared/Pred")
    convCriteria = AredPred;
  else if (choice == "None")
    convCriteria = None;
  else
    convCriteria = ArmijoGoldstein;     // bwb - the others aren't implemented
#endif

  counter->reset();

  return true;
}


bool NOX::LineSearch::Tensor::compute(NOX::Abstract::Group& newGrp,
                      double& step,
                      const NOX::Abstract::Vector& dir,
                      const NOX::Solver::Generic& s)
{
  bool ok;
  counter->incrementNumLineSearches();
  isNewtonDirection = false;

  const NOX::Direction::Tensor& direction =
    (dynamic_cast<const Solver::TensorBasedTest*>(&s))->getDirection();

#ifdef TRIAL_CODE
  // New code added
  const Solver::Generic* test = 0;
  test = dynamic_cast<const Solver::TensorBasedTest*>(&s);
  if (test == 0)
    {
      //printf("Not a TensorBasedTest solver...\n");
    }
  else
    {
      //printf("IS a TensorBasedTest solver...\n");
      direction = (dynamic_cast<const Solver::TensorBasedTest*>(&s))->getDirection();
    }
#endif


  if (counter->getNumLineSearches() == 1  ||  lsType == Newton)
    isNewtonDirection = true;

  // Do line search and compute new soln.
  if (lsType != Dual || isNewtonDirection)
    ok = performLinesearch(newGrp, step, dir, s, direction);
  else if (lsType == Dual) {
    const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
    double fprime = slopeObj.computeSlope(dir, oldGrp);
    double tensorf = 0.0;
    double tensorStep = 1.0;
    bool isTensorDescent = false;

    if (fprime < 0) {
      ok = performLinesearch(newGrp, step, dir, s, direction);
      tensorf = 0.5*newGrp.getNormF()*newGrp.getNormF();
      tensorStep = step;
      isTensorDescent = true;
    }

    //const NOX::Abstract::Vector& dir2 = direction.getNewton();
    //ok = performLinesearch(newGrp, step, dir2, s, direction);
    ok = performLinesearch(newGrp, step, direction.getNewton(), s, direction);
    double newValue = 0.5*newGrp.getNormF()*newGrp.getNormF();

    // If backtracking on the tensor step produced a better step, then use it.
    if (isTensorDescent  &&  tensorf < newValue) {
      newGrp.computeX(oldGrp, dir, tensorStep);
      newGrp.computeF();
    }
  }

  return ok;
}


bool NOX::LineSearch::Tensor::performLinesearch(NOX::Abstract::Group& newsoln,
                double& step,
                const NOX::Abstract::Vector& dir,
                const NOX::Solver::Generic& s,
                const NOX::Direction::Tensor& direction)
{
  if (utils.isPrintType(NOX::Utils::InnerIteration)) {
    utils.out() << "\n" << NOX::Utils::fill(72) << "\n";
    utils.out() << "-- Tensor Line Search ("
     << paramsPtr->get("Submethod","Curvilinear")
     << ") -- \n";
  }

  // Local variables
  Teuchos::RCP<NOX::Abstract::Vector> dir2;
  bool isFailed = false;
  bool isAccepted = false;
  bool isFirstPass = true;
  std::string message = "(STEP ACCEPTED!)";

  // Set counters
  int lsIterations = 1;

  // Get the linear solve tolerance if doing ared/pred for conv criteria
  std::string dirString = const_cast<Teuchos::ParameterList&>(s.getList()).
    sublist("Direction").get("Method", "Tensor");
  double eta = (suffDecrCond == AredPred) ?
    const_cast<Teuchos::ParameterList&>(s.getList()).sublist("Direction").
    sublist(dirString).sublist("Linear Solver").get("Tolerance", -1.0) : 0.0;

  // Get Old function value
  const Abstract::Group& oldsoln = s.getPreviousSolutionGroup();
  double oldValue = 0.5*oldsoln.getNormF()*oldsoln.getNormF();

  // Compute directional derivative at old solution
  double fprime = (lsType == Curvilinear) ?
    slopeObj.computeSlope(direction.getNewton(), oldsoln) :
    slopeObj.computeSlope(dir, oldsoln);
  multsJv++;

  // Compute first trial point and its function value
  step = defaultStep;
  newsoln.computeX(oldsoln, dir, step);
  newsoln.computeF();
  double newValue = 0.5*newsoln.getNormF()*newsoln.getNormF();

  // Compute the convergence criteria for the line search
  //  double threshold = oldValue + alpha*step*fprime;
  //  isAccepted = (newValue < threshold);
//   if (fprime >= 0.0)
//   {
//     printBadSlopeWarning(fprime);
//     isFailed = true;
//   }
//   else
    isAccepted = checkConvergence(newValue, oldValue, fprime, step, eta,
                  lsIterations, 0);

  // Update counter and allocate memory for dir2 if a linesearch is needed
  if (!isAccepted) {
    counter->incrementNumNonTrivialLineSearches();
    dir2 = dir.clone(ShapeCopy);
    *dir2 = dir;
  }

  // Iterate until the trial point is accepted....
  while ((!isAccepted) && (!isFailed)) {

    // Check for linesearch failure
    if (lsIterations > maxIters) {
      isFailed = true;
      message = "(FAILED - Max Iters)";
      break;
    }

    print.printStep(lsIterations, step, oldValue, newValue);

    // Is the full tensor step a descent direction?  If not, switch to Newton
    if ((!isNewtonDirection) && (isFirstPass && fprime >= 0)) {
      // dir = oldsoln.getNewton();   // bwb - for when Newton put in group
      *dir2 = direction.getNewton();
      fprime = slopeObj.computeSlope(*dir2, oldsoln);
      multsJv++;

      printf("  Switching to Newton.  New fprime = %e\n", fprime);
    }
    else {
      step = selectLambda(newValue, oldValue, fprime, step);
    }

    isFirstPass = false;

    // Check for linesearch failure
    if (step < minStep) {
      isFailed = true;
      message = "(FAILED - Min Step)";
      break;
    }

    // Update the number of linesearch iterations
    counter->incrementNumIterations();
    lsIterations ++;

    // Compute new trial point and its function value
    if (lsType == Curvilinear) {

      // bwb - oldsoln needed for preconditioner, test when right pre available
      //bool ok = direction.computeCurvilinearStep2(*dir2, newsoln, s, step);
      //bool ok = direction.computeCurvilinearStep(*dir2, newsoln, s, step);
      direction.computeCurvilinearStep(*dir2, oldsoln, s, step);
      newsoln.computeX(oldsoln, *dir2, 1.0);
    }
    else {
      newsoln.computeX(oldsoln, *dir2, step);
    }
    newsoln.computeF();
    newValue = 0.5*newsoln.getNormF()*newsoln.getNormF();

    // Recompute convergence criteria based on new step
    // threshold = oldValue + alpha*step*fprime;
    // isAccepted = (newValue < threshold);
    isAccepted = checkConvergence(newValue, oldValue, fprime, step, eta,
                  lsIterations, 0);
  }

  if (isFailed) {
    counter->incrementNumFailedLineSearches();
    step = recoveryStep;

    if (step != 0.0) {
      // Update the group using Newton direction and recovery step
      newsoln.computeX(oldsoln, direction.getNewton(), step);
      newsoln.computeF();
      newValue = 0.5*newsoln.getNormF()*newsoln.getNormF();

      message = "(USING RECOVERY STEP!)";
    }
  }

  print.printStep(lsIterations, step, oldValue, newValue, message);
  counter->setValues(*paramsPtr);

  dir2 = Teuchos::null;

  if (suffDecrCond == AredPred)
    paramsPtr->set("Adjusted Tolerance", 1.0 - step * (1.0 - eta));

  return (!isFailed);
}


bool NOX::LineSearch::Tensor::checkConvergence(double newValue, double oldValue,
                           double oldSlope,
                           double step, double eta,
                           int nIters,
                           int nNonlinearIters) const
{
  /*
  if ((nIters == 1) && (doForceInterpolation))
    return false;

  if ((doAllowIncrease) && (nNonlinearIters <= maxIncreaseIter))
  {
    double relativeIncrease = newValue / oldValue;
    if (relativeIncrease < maxRelativeIncrease)
      return true;
  }
  */

  switch (suffDecrCond)
  {

  case ArmijoGoldstein:

    return (newValue <= oldValue + alpha * step * oldSlope);
    break;

  case AredPred:
    {
      double newEta = 1.0 - step * (1.0 - eta);
      return (newValue <= oldValue * (1.0 - alpha * (1.0 - newEta)));
      break;
    }
  case None:

    return true;
    break;

  default:

    utils.err() << "NOX::LineSearch::Tensor::checkConvergence - Unknown convergence criteria" << std::endl;
    throw std::runtime_error("NOX Error");

  }
}

double NOX::LineSearch::Tensor::selectLambda(double newf, double oldf,
                         double oldfprime, double lambda)
{
  double lambdaRet;
  double temp;

  if (lambdaSelection == Quadratic) {
    temp = -oldfprime / (2.0*(newf - oldf - oldfprime));
    if (temp < 0.1)
      temp = 0.1;
    lambdaRet = temp *lambda;
  }
  else {
    lambdaRet = 0.5 * lambda;
  }
  return lambdaRet;
}


void NOX::LineSearch::Tensor::printBadSlopeWarning(double slope) const
{
  if (print.isPrintType(NOX::Utils::Warning))
    utils.out() << "WARNING: Computed slope is positive (slope = "
     << slope
     << ").\n" << "Using recovery step!"
     << std::endl;
}

#endif
