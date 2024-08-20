// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Solver_TensorBased.H"    // class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"
#include "NOX_Solver_SolverUtils.H"
#include "NOX_Observer.hpp"
#include "NOX_LineSearch_Utils_Printing.H"  // class data member
#include "NOX_LineSearch_Utils_Counters.H"  // class data member
#include "NOX_LineSearch_Utils_Slope.H"     // class data member
#include "NOX_SolverStats.hpp"

#define CHECK_RESIDUALS
#define DEBUG_LEVEL 0
#define DEVELOPER_CODE
#define USE_INITIAL_GUESS_LOGIC

NOX::Solver::TensorBased::
TensorBased(const Teuchos::RCP<NOX::Abstract::Group>& xGrp,
        const Teuchos::RCP<NOX::StatusTest::Generic>& t,
        const Teuchos::RCP<Teuchos::ParameterList>& p) :
  solnPtr(xGrp),
  oldSolnPtr(xGrp->clone(DeepCopy)), // create via clone
  newtonVecPtr(xGrp->getX().clone(ShapeCopy)), // create via clone
  tensorVecPtr(xGrp->getX().clone(ShapeCopy)), // create via clone
  aVecPtr(xGrp->getX().clone(ShapeCopy)), // create via clone
  sVecPtr(xGrp->getX().clone(ShapeCopy)), // create via clone
  tmpVecPtr(xGrp->getX().clone(ShapeCopy)), // create via clone
  residualVecPtr(xGrp->getX().clone(ShapeCopy)), // create via clone
  testPtr(t),
  paramsPtr(p),
  linearParamsPtr(0),
  beta(0.),
  sTinvJF(0.),
  sTinvJa(0.)
{
  reset(xGrp, t, p);
}

// Protected
void NOX::Solver::TensorBased::init()
{
  // Initialize
  stepSize = 0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;

  // Reset counters
  counter->reset();
  numJvMults = 0;
  numJ2vMults = 0;

  // Print out parameters
  if (utilsPtr->isPrintType(NOX::Utils::Parameters))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utilsPtr->out(),5);
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
  }

}


bool NOX::Solver::TensorBased::
reset(const Teuchos::RCP<NOX::Abstract::Group>& xGrp,
      const Teuchos::RCP<NOX::StatusTest::Generic>& t,
      const Teuchos::RCP<Teuchos::ParameterList>& p)
{
  solnPtr = xGrp;
  testPtr = t;
  paramsPtr = p;

  NOX::Solver::validateSolverOptionsSublist(p->sublist("Solver Options"));
  globalDataPtr = Teuchos::rcp(new NOX::GlobalData(p));
  utilsPtr = globalDataPtr->getUtils();
  print = Teuchos::rcp(new NOX::LineSearch::Utils::Printing(utilsPtr));
  counter = &globalDataPtr->getNonConstSolverStatistics()->lineSearch;
  slopeObj.reset(globalDataPtr);
  observer = NOX::Solver::parseObserver(p->sublist("Solver Options"));

  // *** Reset direction parameters ***
  Teuchos::ParameterList& dirParams = paramsPtr->sublist("Direction");

  // Determine the specific type of direction to compute
  std::string choice = dirParams.get("Method", "Tensor");
  if (choice == "Tensor")
    requestedBaseStep = TensorStep;
  else if (choice == "Newton")
    requestedBaseStep = NewtonStep;
  else
  {
    if (utilsPtr->isPrintType(NOX::Utils::Error))
      utilsPtr->err() << "NOX::Direction::Tensor::reset() - The choice of "
       << "\"Method\" parameter \"" << choice
       << "\" is invalid." << std::endl;
    throw std::runtime_error("NOX Error");
  }

  // Make a reference to the sublist holding the global strategy parameters
  Teuchos::ParameterList& teParams = dirParams.sublist(choice);

  //  Copy Method into "Compute Step" (temporary hack for data scripts)
  dirParams.set("Compute Step", choice);

  // Initialize direction parameters for this object
  doRescue = teParams.get("Rescue Bad Newton Solve", true);

  checkType = parseStatusTestCheckType(paramsPtr->sublist("Solver Options"));

  // Determine whether we should use the Modified Tensor method
  useModifiedMethod = false;
  if (requestedBaseStep == TensorStep)
  {
    useModifiedMethod =
      dirParams.get("Use Modified Bouaricha", true);
    if (useModifiedMethod  &&
    utilsPtr->isPrintType(NOX::Utils::Parameters))
      utilsPtr->out() << "Using Modified Bouaricha method" << std::endl;
  }


  // *** Reset parameters for Line Search ***
  Teuchos::ParameterList& lsParams = paramsPtr->sublist("Line Search");

  // Determine the specific type of tensor linesearch to perform
  choice = lsParams.get("Method", "Curvilinear");

  if (choice == "Curvilinear")
    lsType = Curvilinear;
  else if (choice == "Dual")
    lsType = Dual;
  else if (choice == "Standard")
    lsType = Standard;
  else if (choice == "Full Step")
    lsType = FullStep;
  else if (choice == "Newton")
    lsType = Newton;
  else
  {
    if (utilsPtr->isPrintType(NOX::Utils::Error))
      utilsPtr->err() << "NOX::Direction::Tensor::reset() - The choice of "
       << "\"Line Search\" parameter " << choice
       << " is invalid." << std::endl;
    throw std::runtime_error("NOX Error");
  }
  //  Copy Method into "Submethod" (temporary hack for data scripts)
  lsParams.set("Submethod", choice);

  // Make a reference to the sublist holding the global strategy parameters
  Teuchos::ParameterList& gsParams = lsParams.sublist(choice);

  // Decide what step to use in case of linesearch failure
  choice = gsParams.get("Recovery Step Type", "Constant");
  if (choice == "Constant")
    recoveryStepType = Constant;          // Use value in "Recovery Step"
  else if (choice == "Last Computed Step")
    recoveryStepType = LastComputedStep;  // Use last step from linesearch
  else
  {
    utilsPtr->err() << "NOX::Solver::TensorBased::reset() - "
     << "Invalid \"Recovery Step Type\"" << std::endl;
    throw std::runtime_error("NOX Error");
  }

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
    if (utilsPtr->isPrintType(NOX::Utils::Error))
      utilsPtr->err() << "NOX::Solver::TensorBased::reset() - The choice of "
       << "\"Lambda Selection\" parameter " << choice
       << " is invalid." << std::endl;
    throw std::runtime_error("NOX Error");
  }

  choice = gsParams.get("Sufficient Decrease Condition",
                 "Armijo-Goldstein");
  if (choice == "Armijo-Goldstein")
    convCriteria = ArmijoGoldstein;     // This is the only one implemented
  //else if (choice == "Ared/Pred")
  //  convCriteria = AredPred;
  //else if (choice == "None")
  //  convCriteria = None;
  else
  {
    if (utilsPtr->isPrintType(NOX::Utils::Error))
      utilsPtr->err() << "NOX::Solver::TensorBased::reset() - The choice of "
       << "\"Sufficient Decrease Condition\" parameter " << choice
       << " is invalid." << std::endl;
    throw std::runtime_error("NOX Error");
  }

  init();
  return true;
}

void NOX::Solver::TensorBased::
reset(const NOX::Abstract::Vector& initialGuess,
      const Teuchos::RCP<NOX::StatusTest::Generic>& t)
{
  solnPtr->setX(initialGuess);
  testPtr = t;
  init();
}

void NOX::Solver::TensorBased::
reset(const NOX::Abstract::Vector& initialGuess)
{
  solnPtr->setX(initialGuess);
  init();
}

void NOX::Solver::TensorBased::
reset()
{
  stepSize = 0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
  counter->reset();
  numJvMults = 0;
  numJ2vMults = 0;
}

NOX::Solver::TensorBased::~TensorBased()
{
#ifdef DEVELOPER_CODE
  if (utilsPtr->isPrintType(NOX::Utils::Details))
  {
    utilsPtr->out() << "multsJv = " << numJvMults << "   (linesearch)" << std::endl;
    utilsPtr->out() << "mults2Jv = " << numJ2vMults << std::endl;
  }
#endif
}


NOX::StatusTest::StatusType  NOX::Solver::TensorBased::getStatus() const
{
  return status;
}

NOX::StatusTest::StatusType  NOX::Solver::TensorBased::step()
{
  observer->runPreIterate(*this);

  // On the first step, perform some initl checks
  if (nIter ==0) {
    // Compute F of initial guess
    NOX::Abstract::Group::ReturnType rtype = solnPtr->computeF();
    if (rtype != NOX::Abstract::Group::Ok) {
      utilsPtr->err() << "NOX::Solver::TensorBased::init - "
              << "Unable to compute F" << std::endl;
      throw std::runtime_error("NOX Error");
    }

    // Test the initial guess
    status = testPtr->checkStatus(*this, checkType);
    if ((status == NOX::StatusTest::Converged) &&
    (utilsPtr->isPrintType(NOX::Utils::Warning))) {
      utilsPtr->out() << "Warning: NOX::Solver::TensorBased::init() - "
              << "The solution passed into the solver (either "
              << "through constructor or reset method) "
              << "is already converged!  The solver will not "
              << "attempt to solve this system since status "
              << "is flagged as converged." << std::endl;
    }

    printUpdate();
  }

  // First check status
  if (status != NOX::StatusTest::Unconverged)
  {
    observer->runPostIterate(*this);
    printUpdate();
    return status;
  }

  // Copy pointers into temporary references
  NOX::Abstract::Group& soln = *solnPtr;
  NOX::StatusTest::Generic& test = *testPtr;

  // Compute the direction for the update vector at the current solution.
  bool ok = computeTensorDirection(soln, *this);
  if (!ok)
  {
    if (utilsPtr->isPrintType(NOX::Utils::Error))
      utilsPtr->out() << "NOX::Solver::TensorBased::iterate - "
       << "unable to calculate direction" << std::endl;
    status = NOX::StatusTest::Failed;
    observer->runPostIterate(*this);
    printUpdate();
    return status;
  }

  // Update iteration count.
  nIter ++;

  // Copy current soln to the old soln.
  *oldSolnPtr = soln;

  // Do line search and compute new soln.
  ok = implementGlobalStrategy(soln, stepSize, *this);
  if (!ok)
  {
    if (stepSize == 0.0)
    {
      if (utilsPtr->isPrintType(NOX::Utils::Error))
    utilsPtr->out() << "NOX::Solver::TensorBased::iterate - line search failed"
         << std::endl;
      status = NOX::StatusTest::Failed;
      observer->runPostIterate(*this);
      printUpdate();
      return status;
    }
    else if (utilsPtr->isPrintType(NOX::Utils::Warning))
      utilsPtr->out() << "NOX::Solver::TensorBased::iterate - "
       << "using recovery step for line search" << std::endl;
  }

  // Compute F for new current solution.
  NOX::Abstract::Group::ReturnType rtype = soln.computeF();
  if (rtype != NOX::Abstract::Group::Ok)
  {
    if (utilsPtr->isPrintType(NOX::Utils::Error))
      utilsPtr->out() << "NOX::Solver::TensorBased::iterate - "
       << "unable to compute F" << std::endl;
    status = NOX::StatusTest::Failed;
    observer->runPostIterate(*this);
    printUpdate();
    return status;
  }

  status = test.checkStatus(*this, checkType);

  observer->runPostIterate(*this);

  printUpdate();

  return status;
}


NOX::StatusTest::StatusType  NOX::Solver::TensorBased::solve()
{
  observer->runPreSolve(*this);

  this->reset();

  // Iterate until converged or failed
  while (status == NOX::StatusTest::Unconverged)
  {
    status = step();
  }

  Teuchos::ParameterList& outputParams = paramsPtr->sublist("Output");
  outputParams.set("Nonlinear Iterations", nIter);
  outputParams.set("2-Norm of Residual", solnPtr->getNormF());

  observer->runPostSolve(*this);

  return status;
}

const NOX::Abstract::Group&
NOX::Solver::TensorBased::getSolutionGroup() const
{
  return *solnPtr;
}

const NOX::Abstract::Group&
NOX::Solver::TensorBased::getPreviousSolutionGroup() const
{
  return *oldSolnPtr;
}

int NOX::Solver::TensorBased::getNumIterations() const
{
  return nIter;
}

const Teuchos::ParameterList&
NOX::Solver::TensorBased::getList() const
{
  return *paramsPtr;
}

Teuchos::RCP<const NOX::SolverStats>
NOX::Solver::TensorBased::getSolverStatistics() const
{
  return globalDataPtr->getSolverStatistics();
}

// protected
void NOX::Solver::TensorBased::printUpdate()
{
  double normSoln = 0;
  double normStep = 0;

  // Print the status test parameters at each iteration if requested
  if ((status == NOX::StatusTest::Unconverged) &&
      (utilsPtr->isPrintType(NOX::Utils::OuterIterationStatusTest)))
  {
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration))
  {
    normSoln = solnPtr->getNormF();
    normStep = (nIter > 0) ? tensorVecPtr->norm() : 0;
  }

  // ...But only the print process actually prints the result.
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Nonlinear Solver Step " << nIter << " -- \n";
    utilsPtr->out() << "f = " << utilsPtr->sciformat(normSoln);
    utilsPtr->out() << "  step = " << utilsPtr->sciformat(stepSize);
    utilsPtr->out() << "  dx = " << utilsPtr->sciformat(normStep);
    if (status == NOX::StatusTest::Converged)
      utilsPtr->out() << " (Converged!)";
    if (status == NOX::StatusTest::Failed)
      utilsPtr->out() << " (Failed!)";
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n" << std::endl;
  }

  // Print the final parameter values of the status test
  if ((status != NOX::StatusTest::Unconverged) &&
      (utilsPtr->isPrintType(NOX::Utils::OuterIteration)))
  {
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Final Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
  }
}


bool
NOX::Solver::TensorBased::computeTensorDirection(NOX::Abstract::Group& soln,
                     const NOX::Solver::Generic& solver)
{
  NOX::Abstract::Group::ReturnType dir_status;

  Teuchos::ParameterList& linearParams = paramsPtr->sublist("Direction").
    sublist(paramsPtr->sublist("Direction").
        get("Method","Tensor")).
    sublist("Linear Solver");

  // Compute F at current solution.
  dir_status = soln.computeF();
  if (dir_status != NOX::Abstract::Group::Ok)
    throwError("computeTensorDirection", "Unable to compute F");

  // Compute Jacobian at current solution.
  dir_status = soln.computeJacobian();
  if (dir_status != NOX::Abstract::Group::Ok)
    throwError("computeTensorDirection", "Unable to compute Jacobian");

  // Begin processing for the tensor step, if necessary.
  double sDotS = 0.0;
  int tempVal1 = 0;
  if ((nIter > 0)  &&  (requestedBaseStep == TensorStep))
  {
    // Compute the tensor term s = x_{k-1} - x_k
    *sVecPtr = soln.getX();
    sVecPtr->update(1.0, solver.getPreviousSolutionGroup().getX(), -1.0);
    double normS = sVecPtr->norm();
    sDotS = normS * normS;

    // Form the tensor term a = (F_{k-1} - F_k - J*s) / (s^T s)^2
    soln.applyJacobian(*sVecPtr, *aVecPtr);
    numJvMults++;
    aVecPtr->update(1.0, solver.getPreviousSolutionGroup().getF(), -1.0);
    aVecPtr->update(-1.0, soln.getF(), 1.0);
    if (sDotS != 0)
      aVecPtr->scale(1.0 / (sDotS * sDotS));

    // Save old Newton step as initial guess to second system
    *tmpVecPtr = *newtonVecPtr;
    tmpVecPtr->scale(-1.0);   // Rewrite to avoid this?

    // Compute residual of linear system using initial guess...
    soln.applyJacobian(*tmpVecPtr, *residualVecPtr);
    numJvMults++;
    residualVecPtr->update(1.0, solver.getPreviousSolutionGroup().getF(),-1.0);
    double residualNorm = residualVecPtr->norm();

#if DEBUG_LEVEL > 0
    double tmpVecNorm = tmpVecPtr->norm();
    double residualNormRel = residualNorm /
      solver.getPreviousSolutionGroup().getNormF();
    if (utilsPtr->isPrintType(NOX::Utils::Details))
    {
      utilsPtr->out() << "  Norm of initial guess: " << utilsPtr->sciformat(tmpVecNorm, 6)
       << std::endl;
      utilsPtr->out() << "  initg norm of model residual =   "
       << utilsPtr->sciformat(residualNorm, 6) << " (abs)     "
       << utilsPtr->sciformat(residualNormRel, 6) << " (rel)" << std::endl;
    }
#endif

    // Save some parameters and use them later...
    double tol = linearParams.get("Tolerance", 1e-4);
    double relativeResidual = residualNorm /
      solver.getPreviousSolutionGroup().getNormF();

    // Decide whether to use initial guess...
    bool isInitialGuessGood = false;
#ifdef USE_INITIAL_GUESS_LOGIC
    if (relativeResidual < 1.0)
    {
      if (utilsPtr->isPrintType(NOX::Utils::Details))
    utilsPtr->out() << "  Initial guess is good..." << std::endl;
      isInitialGuessGood = true;
      // RPP - Brett please make sure the line below is correct.
      *tensorVecPtr = *tmpVecPtr;
      double newTol = tol / relativeResidual;
      if (newTol > 0.99)
    newTol = 0.99;  // force at least one iteration
      linearParams.set("Tolerance",  newTol);
      if (utilsPtr->isPrintType(NOX::Utils::Details))
    utilsPtr->out() << "  Setting tolerance to " << utilsPtr->sciformat(newTol,6) << std::endl;
    }
    else
#endif // USE_INITIAL_GUESS_LOGIC
    {
      //utilsPtr->out() << "  Initial guess is BAD... do not use!\n";
      isInitialGuessGood = false;
      *residualVecPtr = solver.getPreviousSolutionGroup().getF();
    }

    // Compute the term inv(J)*Fp....
    tmpVecPtr->init(0.0);
    dir_status = soln.applyJacobianInverse(linearParams, *residualVecPtr,
                       *tmpVecPtr);

    // If it didn't converge, maybe we can recover.
    if (dir_status != NOX::Abstract::Group::Ok)
    {
      if (doRescue == false)
    throwError("computeTensorDirection", "Unable to apply Jacobian inverse");
      else if ((doRescue == true) &&
           (utilsPtr->isPrintType(NOX::Utils::Warning)))
    utilsPtr->out() << "WARNING: NOX::Solver::TensorBased::computeTensorDirection() - "
         << "Linear solve failed to achieve convergence - "
         << "using the step anyway "
         << "since \"Rescue Bad Newton Solve\" is true." << std::endl;
    }

    // Continue processing
#ifdef USE_INITIAL_GUESS_LOGIC
    if (isInitialGuessGood)
    {
      tmpVecPtr->update(1.0, *tensorVecPtr, 1.0);
      linearParams.set("Tolerance",  tol);
    }
#endif

    // Save iteration count for comparison later
    if (linearParams.sublist("Output").
    isParameter("Number of Linear Iterations"))
      tempVal1 = linearParams.sublist("Output").
    get("Number of Linear Iterations",0);

#if DEBUG_LEVEL > 0
    // Compute residual of linear system with initial guess...
    soln.applyJacobian(*tmpVecPtr, *residualVecPtr);
    numJvMults++;
    residualVec.update(-1.0, solver.getPreviousSolutionGroup().getF(),1.0);
    double residualNorm2 = residualVec.norm();
    double residualNorm2Rel = residualNorm2 /
      solver.getPreviousSolutionGroup().getNormF();
    if (utilsPtr->isPrintType(NOX::Utils::Details))
      utilsPtr->out() << " jifp norm of model residual =   "
       << utilsPtr->sciformat(residualNorm2, 6) << " (abs)     "
       << utilsPtr->sciformat(residualNorm2Rel, 6) << " (rel)" << std::endl;
#endif
  }

  // Compute the Newton direction
  dir_status = soln.computeNewton(linearParams);

  // If it didn't converge, maybe we can recover.
  if (dir_status != NOX::Abstract::Group::Ok)
  {
    if (doRescue == false)
      throwError("computeTensorDirection", "Unable to apply Jacobian inverse");
    else if ((doRescue == true) &&
         (utilsPtr->isPrintType(NOX::Utils::Warning)))
      utilsPtr->out() << "WARNING: NOX::Solver::TensorBased::computeTensorDirection() - "
       << "Linear solve failed to achieve convergence - "
       << "using the step anyway "
       << "since \"Rescue Bad Newton Solve\" is true." << std::endl;
  }

  // Set Newton direction
  *newtonVecPtr = soln.getNewton();

  // Update counter
  int tempVal2 = 0;
  if (linearParams.sublist("Output").
      isParameter("Number of Linear Iterations"))
    tempVal2 = linearParams.sublist("Output").
      get("Number of Linear Iterations",0);
  numJ2vMults += (tempVal1 > tempVal2) ? tempVal1 : tempVal2;

#ifdef CHECK_RESIDUALS
  printDirectionInfo("newtonVec", *newtonVecPtr, soln, false);
#endif // CHECK_RESIDUALS

  // Continue processing the tensor step, if necessary
  if ((nIter > 0)  &&  (requestedBaseStep == TensorStep))
  {
    // Form the term inv(J)*a...  (note that a is not multiplied by 2)
    // The next line does not work in some implementations for some reason
    //tmpVec.update(1.0, newtonVec, -1.0, sVec, 1.0);
    tmpVecPtr->update(1.0, *newtonVecPtr, 1.0);
    tmpVecPtr->update(-1.0, *sVecPtr, 1.0);
    if (sDotS != 0.0)
      tmpVecPtr->scale( 1.0 / (sDotS * sDotS));

    // Calculate value of beta
    sTinvJF = -sVecPtr->innerProduct(*newtonVecPtr);
    sTinvJa = sVecPtr->innerProduct(*tmpVecPtr);
    double qval = 0;
    double lambdaBar = 1;
    beta = calculateBeta(sTinvJa, 1.0, sTinvJF, qval, lambdaBar);

    double sVecNorm = sVecPtr->norm();
    double aVecNorm = aVecPtr->norm();
    if (utilsPtr->isPrintType(NOX::Utils::Details))
    {
      utilsPtr->out() << " sTinvJF = " << utilsPtr->sciformat(sTinvJF, 6)
       << "  sTinvJa = " << utilsPtr->sciformat(sTinvJa, 6) << std::endl;
      utilsPtr->out() << " norm(s) = " << utilsPtr->sciformat(sVecNorm, 6)
       << "  norm(a) = " << utilsPtr->sciformat(aVecNorm, 6) << std::endl;
    }

    if (useModifiedMethod)
    {
      double alpha2 = lambdaBar;
      if (utilsPtr->isPrintType(NOX::Utils::Details))
    utilsPtr->out() << " Beta = " << utilsPtr->sciformat(beta, 6)
         << "  Alpha2 = " << utilsPtr->sciformat(alpha2, 6) << std::endl;
      if (alpha2 != 1.0)
      {
    if (utilsPtr->isPrintType(NOX::Utils::Details))
      utilsPtr->out() << "   *** Scaling tensor term a ***" << std::endl;
    aVecPtr->scale(alpha2);
    tmpVecPtr->scale(alpha2);
    sTinvJa *= alpha2;
    beta /= alpha2;
    lambdaBar = 1.0;
    qval = 0;
      }
    }

    // Form the tensor step
    tensorVecPtr->update(1.0, *newtonVecPtr, -beta*beta, *tmpVecPtr, 0.0);

#ifdef CHECK_RESIDUALS
    printDirectionInfo("tensorVec", *tensorVecPtr, soln, true);
#endif // CHECK_RESIDUALS
#if DEBUG_LEVEL > 0
    double sDotT = tensorVecPtr->innerProduct(sVec);
    if (utilsPtr->isPrintType(NOX::Utils::Details))
      utilsPtr->out() << "  Beta = " << utilsPtr->sciformat(beta, 6)
       << "  std = " << utilsPtr->sciformat(sDotT, 6)
       << "  qval = " << utilsPtr->sciformat(qval, 2)
       << "  lambdaBar = " << lambdaBar << std::endl;
#endif
  }
  else
    *tensorVecPtr = *newtonVecPtr;

  return true;
}


double NOX::Solver::TensorBased::calculateBeta(double qa,
                           double qb,
                           double qc,
                           double& qval,
                           double& lambdaBar,
                           double lambda) const
{
  double beta_value = 0.0;
  double discriminant = qb*qb - 4*qa*qc*lambda;

  if (discriminant < 0.0)
  {
    // no real root
    beta_value = -qb / qa / 2.0;
    qval = (qa * beta_value * beta_value) + (qb * beta_value) + (lambda * qc);
    lambdaBar = qb*qb / (4*qa*qc);
#if DEBUG_LEVEL > 0
    if (utilsPtr->isPrintType(NOX::Utils::Details))
      utilsPtr->out() << "  ####  LambdaBar = " << lambdaBar << "  ####\n";
#endif
  }
  else
  {
    qval = 0;
    lambdaBar = 1.0;
    if ( (fabs(qa / qb) < 1e-8)  &&  (fabs(lambda * qc / qb) < 1) )
    {
#if DEBUG_LEVEL > 0
      if (utilsPtr->isPrintType(NOX::Utils::Details))
    utilsPtr->out() << "  qa is relatively small\n";
#endif
      beta_value = -lambda * qc / qb;
    }
    else
    {
      double tmp1 = (-qb + sqrt(discriminant)) / (2*qa);
      double tmp2 = (-qb - sqrt(discriminant)) / (2*qa);
      beta_value = (fabs(tmp1) < fabs(tmp2)) ? tmp1 : tmp2; // bwb - temporary test
#if DEBUG_LEVEL > 1
      if (utilsPtr->isPrintType(NOX::Utils::Details))
    utilsPtr->out() << "  tmp1 = " << utilsPtr->sciformat(tmp1, 6)
         << "  tmp2 = " << utilsPtr->sciformat(tmp2, 6)
         << "  dir0xsc = " << utilsPtr->sciformat(dir0xsc, 6)
         << "  normS = " << utilsPtr->sciformat(normS, 6)
         << std::endl;
#endif
    }
  }
#if DEBUG_LEVEL > 1
  if (utilsPtr->isPrintType(NOX::Utils::Details))
    utilsPtr->out() << "  qa,qb,qc = " << utilsPtr->sciformat(qa, 6)
     << utilsPtr->sciformat(qb, 6)
     << utilsPtr->sciformat(qc, 6)
     << "   beta = " << utilsPtr->sciformat(beta_value, 6)
     << std::endl;
#endif

  return beta_value;
}


bool
NOX::Solver::TensorBased::computeCurvilinearStep(NOX::Abstract::Vector& dir,
                     const NOX::Abstract::Group& /* soln */,
                     const NOX::Solver::Generic& /* s */,
                     double& lambda)
{
  double qval = 0;
  double lambdaBar = 1;
  double beta1 = calculateBeta(sTinvJa, 1, sTinvJF, qval, lambdaBar, lambda);
  double betaFactor = ( (beta == 0.0) ? 0.0 : beta1*beta1 / (beta*beta));

  dir.update(lambda - betaFactor, *newtonVecPtr, betaFactor, *tensorVecPtr, 0.0);

#if DEBUG_LEVEL > 0
  double sDotD = dir.innerProduct(sVec);
  if (utilsPtr->isPrintType(NOX::Utils::Details))
  {
    utilsPtr->out() << "  Beta = " << utilsPtr->sciformat(beta, 6)
     << "  std = " << utilsPtr->sciformat(sDotD, 6)
     << "  qval = " << qval
     << "  lambdaBar = " << lambdaBar
     << std::endl;
    utilsPtr->out() << "  betaFactor = " << utilsPtr->sciformat(betaFactor,6)
     << "  beta1 = " << utilsPtr->sciformat(beta1, 6)
     << std::endl;
  }
#endif

  return true;
}


bool
NOX::Solver::TensorBased::implementGlobalStrategy(NOX::Abstract::Group& newGrp,
                      double& in_stepSize,
                      const NOX::Solver::Generic& s)
{
  bool ok;
  counter->incrementNumLineSearches();
  isNewtonDirection = false;
  NOX::Abstract::Vector& searchDirection = *tensorVecPtr;

  if ((counter->getNumLineSearches() == 1)  ||  (lsType == Newton))
  {
    isNewtonDirection = true;
    searchDirection = *newtonVecPtr;
  }

  // Do line search and compute new soln.
  if ((lsType != Dual) || (isNewtonDirection))
    ok = performLinesearch(newGrp, in_stepSize, searchDirection, s);
  else if (lsType == Dual)
  {
    double fTensor = 0.0;
    double fNew = 0.0;
    double tensorStep = 1.0;
    bool isTensorDescent = false;

    const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
    double fprime = slopeObj.computeSlope(searchDirection, oldGrp);

    // Backtrack along tensor direction if it is descent direction.
    if (fprime < 0)
    {
      ok = performLinesearch(newGrp, in_stepSize, searchDirection, s);
      assert(ok);
      fTensor = 0.5 * newGrp.getNormF() * newGrp.getNormF();
      tensorStep = in_stepSize;
      isTensorDescent = true;
    }

    // Backtrack along the Newton direction.
    ok = performLinesearch(newGrp, in_stepSize, *newtonVecPtr, s);
    fNew = 0.5 * newGrp.getNormF() * newGrp.getNormF();

    // If backtracking on the tensor step produced a better step, then use it.
    if (isTensorDescent  &&  (fTensor <= fNew))
    {
      newGrp.computeX(oldGrp, *tensorVecPtr, tensorStep);
      newGrp.computeF();
    }
  }

  return ok;
}


bool
NOX::Solver::TensorBased::performLinesearch(NOX::Abstract::Group& newSoln,
                        double& in_stepSize,
                        const NOX::Abstract::Vector& lsDir,
                        const NOX::Solver::Generic& s)
{
  if (print->isPrintType(NOX::Utils::InnerIteration))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Tensor Line Search (";
    if (lsType == Curvilinear)
      utilsPtr->out() << "Curvilinear";
    else if (lsType == Standard)
      utilsPtr->out() << "Standard";
    else if (lsType == FullStep)
      utilsPtr->out() << "Full Step";
    else if (lsType == Dual)
      utilsPtr->out() << "Dual";
    utilsPtr->out() << ") -- " << std::endl;
  }

  // Local variables
  bool isFailed = false;
  bool isAcceptable = false;
  bool isFirstPass = true;
  std::string message = "(STEP ACCEPTED!)";

  // Set counters
  int lsIterations = 1;

  // Get Old f
  const Abstract::Group& oldSoln = s.getPreviousSolutionGroup();
  double fOld = 0.5 * oldSoln.getNormF() * oldSoln.getNormF();

  // Compute first trial point and its function value
  in_stepSize = defaultStep;
  newSoln.computeX(oldSoln, lsDir, in_stepSize);
  newSoln.computeF();
  double fNew = 0.5 * newSoln.getNormF() * newSoln.getNormF();

  // Stop here if only using the full step
  if (lsType == FullStep)
  {
      print->printStep(lsIterations, in_stepSize, fOld, fNew, message);
      return (!isFailed);
  }

  // Compute directional derivative
  double fprime;
  if ((lsType == Curvilinear)  &&  !(isNewtonDirection))
    fprime = slopeObj.computeSlope(*newtonVecPtr, oldSoln);
  else
    fprime = slopeObj.computeSlope(lsDir, oldSoln);
  numJvMults++;  // computeSlope() has J*v inside of it

  // Compute the convergence criteria for the line search
  double threshold = fOld + alpha*in_stepSize*fprime;
  isAcceptable = (fNew < threshold);

  // Update counter and temporarily hold direction if a linesearch is needed
  if (!isAcceptable)
  {
    counter->incrementNumNonTrivialLineSearches();
    *tmpVecPtr = lsDir;
  }

  // Iterate until the trial point is accepted....
  while (!isAcceptable)
  {
    // Check for linesearch failure
    if (lsIterations > maxIters)
    {
      isFailed = true;
      message = "(FAILED - Max Iters)";
      break;
    }

    print->printStep(lsIterations, in_stepSize, fOld, fNew);

    // Is the full tensor step a descent direction?  If not, switch to Newton
    if (isFirstPass &&
    (!isNewtonDirection) &&
    (fprime >= 0) &&
    (lsType != Curvilinear) )
    {
      *tmpVecPtr = *newtonVecPtr;
      fprime = slopeObj.computeSlope(*tmpVecPtr, oldSoln);
      numJvMults++;

      if (utilsPtr->isPrintType(NOX::Utils::Details))
    utilsPtr->out() << "  Switching to Newton step.  New fprime = "
         << utilsPtr->sciformat(fprime, 6) << std::endl;
    }
    else
    {
      in_stepSize = selectLambda(fNew, fOld, fprime, in_stepSize);
    }

    isFirstPass = false;

    // Check for linesearch failure
    if (in_stepSize < minStep)
    {
      isFailed = true;
      message = "(FAILED - Min Step)";
      break;
    }

    // Update the number of linesearch iterations
    counter->incrementNumIterations();
    lsIterations ++;

    // Compute new trial point and its function value
    if ((lsType == Curvilinear) && !(isNewtonDirection))
    {
      computeCurvilinearStep(*tmpVecPtr, oldSoln, s, in_stepSize);
      // Note: oldSoln is needed above to get correct preconditioner
      newSoln.computeX(oldSoln, *tmpVecPtr, 1.0);
    }
    else
    {
      newSoln.computeX(oldSoln, *tmpVecPtr, in_stepSize);
    }
    newSoln.computeF();
    fNew = 0.5 * newSoln.getNormF() * newSoln.getNormF();

    // Recompute convergence criteria based on new step
    threshold = fOld + alpha*in_stepSize*fprime;
    isAcceptable = (fNew < threshold);
  }


  if (isFailed)
  {
    counter->incrementNumFailedLineSearches();

    if (recoveryStepType == Constant)
    {
      in_stepSize = recoveryStep;
      if (in_stepSize == 0.0)
      {
    newSoln = oldSoln;
    newSoln.computeF();
    fNew = fOld;
      }
      else
      {
    // Update the group using recovery step
    if ((lsType == Curvilinear) && !(isNewtonDirection))
    {
      computeCurvilinearStep(*tmpVecPtr, oldSoln, s, in_stepSize);
      // Note: oldSoln is needed above to get correct preconditioner
      newSoln.computeX(oldSoln, *tmpVecPtr, 1.0);
    }
    else
    {
      newSoln.computeX(oldSoln, *tmpVecPtr, in_stepSize);
    }
    //newSoln.computeX(oldSoln, lsDir, in_stepSize);
    newSoln.computeF();
    fNew = 0.5 * newSoln.getNormF() * newSoln.getNormF();
    message = "(USING RECOVERY STEP!)";
      }
    }
    else
      message = "(USING LAST STEP!)";
  }

  print->printStep(lsIterations, in_stepSize, fOld, fNew, message);
  counter->setValues(paramsPtr->sublist("Line Search"));

  return (!isFailed);
}


double
NOX::Solver::TensorBased::getNormModelResidual(
                                       const NOX::Abstract::Vector& dir,
                       const NOX::Abstract::Group& soln,
                       bool isTensorModel) const
{

  // Compute residual of Newton model...
  Teuchos::RCP<NOX::Abstract::Vector> residualPtr =
    soln.getF().clone(ShapeCopy);
  soln.applyJacobian(dir, *residualPtr);
  numJvMults++;
  residualPtr->update(1.0, soln.getF(), 1.0);

  // Compute residual of Tensor model, if requested...
  if (isTensorModel)
  {
    double tmp = sVecPtr->innerProduct(dir);
    if (utilsPtr->isPrintType(NOX::Utils::Details))
      utilsPtr->out() << " sc'*dt   = " << utilsPtr->sciformat(tmp, 6) << std::endl;
    residualPtr->update(tmp*tmp, *aVecPtr, 1.0);
  }

  double modelNorm = residualPtr->norm();
  return modelNorm;
}


void
NOX::Solver::TensorBased::printDirectionInfo(std::string dirName,
                    const NOX::Abstract::Vector& dir,
                    const NOX::Abstract::Group& soln,
                    bool isTensorModel) const
{
  double dirNorm = dir.norm();

  double residual = getNormModelResidual(dir, soln, isTensorModel);
  double residualRel = residual / soln.getNormF();

  double fprime = getDirectionalDerivative(dir, soln);
  double fprimeRel = fprime / dirNorm;

  if (utilsPtr->isPrintType(NOX::Utils::Details))
  {
    utilsPtr->out() << " " << dirName << " norm of model residual =   "
     << utilsPtr->sciformat(residual, 6) << " (abs)     "
     << utilsPtr->sciformat(residualRel, 6) << " (rel)" << std::endl;
    utilsPtr->out() << " " << dirName << " directional derivative =  "
     << utilsPtr->sciformat(fprime, 6) << " (abs)    "
     << utilsPtr->sciformat(fprimeRel, 6) << " (rel)" << std::endl;
    utilsPtr->out() << " " << dirName << " norm = "
       << utilsPtr->sciformat(dirNorm, 6) << std::endl;
  }
}


double NOX::Solver::TensorBased::getDirectionalDerivative(
                       const NOX::Abstract::Vector& dir,
                       const NOX::Abstract::Group& soln) const
{
  Teuchos::RCP<NOX::Abstract::Vector> tmpPtr =
    soln.getF().clone(ShapeCopy);
  soln.applyJacobian(dir, *tmpPtr);
  numJvMults++;
  double fprime = tmpPtr->innerProduct(soln.getF());
  return fprime;
}


double NOX::Solver::TensorBased::selectLambda(double fNew, double fOld,
                          double fOldPrime,
                          double lambda)
{
  double lambdaRet;
  double temp;

  if (lambdaSelection == Quadratic)
  {
    temp = -fOldPrime / (2.0*(fNew - fOld - fOldPrime));
    if (temp < 0.1)
      temp = 0.1;
    lambdaRet = temp * lambda;
  }
  else
  {
    lambdaRet = 0.5 * lambda;
  }
  return lambdaRet;
}


void NOX::Solver::TensorBased::throwError(const std::string& functionName,
                      const std::string& errorMsg) const
{
  if (utilsPtr->isPrintType(NOX::Utils::Error))
    utilsPtr->err() << "NOX::Solver::TensorBased::" << functionName
     << " - " << errorMsg << std::endl;
  throw std::runtime_error("NOX Error");
}

