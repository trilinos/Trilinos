// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Solver_TrustRegionBased.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_MeritFunction_Generic.H"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"
#include "NOX_Solver_SolverUtils.H"
#include "NOX_Direction_Generic.H"
#include "NOX_Direction_Factory.H"
#include "NOX_Observer.hpp"
#include "NOX_SolverStats.hpp"
#include <cmath>

using namespace NOX;
using namespace NOX::Solver;

TrustRegionBased::
TrustRegionBased(const Teuchos::RCP<NOX::Abstract::Group>& grp,
		 const Teuchos::RCP<NOX::StatusTest::Generic>& t,
		 const Teuchos::RCP<Teuchos::ParameterList>& p) :
  solnPtr(grp),
  oldSolnPtr(grp->clone(DeepCopy)),
  newtonVecPtr(grp->getX().clone(ShapeCopy)),
  cauchyVecPtr(grp->getX().clone(ShapeCopy)),
  aVecPtr(grp->getX().clone(ShapeCopy)),
  bVecPtr(grp->getX().clone(ShapeCopy)),
  testPtr(t),
  paramsPtr(p),
  useAredPredRatio(false)
{
  NOX::Solver::validateSolverOptionsSublist(p->sublist("Solver Options"));
  globalDataPtr = Teuchos::rcp(new NOX::GlobalData(p));
  utilsPtr = globalDataPtr->getUtils();
  meritFuncPtr = globalDataPtr->getMeritFunction();
  observer = NOX::Solver::parseObserver(p->sublist("Solver Options"));
  init();
}

// Protected
void TrustRegionBased::init()
{
  // Initialize
  nIter = 0;
  dx = 0;
  status = StatusTest::Unconverged;
  radius = 0.;
  
  // Print out initialization information
  if (utilsPtr->isPrintType(NOX::Utils::Parameters)) {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utilsPtr->out(),5);
  }

  // Set default parameter settings using get() if they are not set
  paramsPtr->sublist("Direction").get("Method", "Newton");
  paramsPtr->sublist("Cauchy Direction").
    get("Method", "Steepest Descent");
  paramsPtr->sublist("Cauchy Direction").sublist("Steepest Descent").
    get("Scaling Type", "Quadratic Model Min");

  newtonPtr = NOX::Direction::
    buildDirection(globalDataPtr, paramsPtr->sublist("Direction"));

  cauchyPtr = NOX::Direction::
    buildDirection(globalDataPtr, paramsPtr->sublist("Cauchy Direction"));

  minRadius = paramsPtr->sublist("Trust Region").
    get("Minimum Trust Region Radius", 1.0e-6);
  if (minRadius <= 0)
    invalid("Minimum Trust Region Radius", minRadius);

  maxRadius = paramsPtr->sublist("Trust Region").get("Maximum Trust Region Radius", 1.0e+10);
  if (maxRadius <= minRadius)
    invalid("Maximum Trust Region Radius", maxRadius);

  initRadius = paramsPtr->sublist("Trust Region").
    get("Initial Trust Region Radius", 1.0e0);
  if (initRadius <= minRadius || initRadius >= maxRadius)
    invalid("Initial Trust Region Radius", initRadius);

  minRatio = paramsPtr->sublist("Trust Region").get("Minimum Improvement Ratio", 1.0e-4);
  if (minRatio <= 0)
    invalid("Minimum Improvement Ratio", minRatio);

  contractTriggerRatio = paramsPtr->sublist("Trust Region").get("Contraction Trigger Ratio", 0.1);
  if (contractTriggerRatio < minRatio)
    invalid("Contraction Trigger Ratio", contractTriggerRatio);

  expandTriggerRatio = paramsPtr->sublist("Trust Region").get("Expansion Trigger Ratio", 0.75);
  if (expandTriggerRatio <= contractTriggerRatio)
    invalid("Expansion Trigger Ratio", expandTriggerRatio);

  contractFactor = paramsPtr->sublist("Trust Region").get("Contraction Factor", 0.25);
  if ((contractFactor <= 0) || (contractFactor >= 1))
    invalid("Contraction Factor", contractFactor);

  expandFactor = paramsPtr->sublist("Trust Region").get("Expansion Factor", 4.0);
  if (expandFactor <= 1)
    invalid("Expansion Factor", expandFactor);

  recoveryStep = paramsPtr->sublist("Trust Region").get("Recovery Step", 1.0);
  if (recoveryStep < 0)
    invalid("Recovery Step", recoveryStep);

  checkType = parseStatusTestCheckType(paramsPtr->sublist("Solver Options"));

  // Check for the using Homer Walker's Ared/Pred ratio calculation
  useAredPredRatio =
    paramsPtr->sublist("Trust Region").
    get("Use Ared/Pred Ratio Calculation", false);

}

//PRIVATE
void NOX::Solver::TrustRegionBased::
invalid(const std::string& name, double value) const
{
  utilsPtr->err() << "NOX::Solver::TrustRegionBased::init - "
       << "Invalid \"" << name << "\" (" << value << ")"
       << std::endl;
  throw std::runtime_error("NOX Error");
}

void TrustRegionBased::
reset(const NOX::Abstract::Vector& initialGuess,
      const Teuchos::RCP<NOX::StatusTest::Generic>& t)
{
  // New initial guess and status test
  solnPtr->setX(initialGuess);
  testPtr = t;

  // Initialize
  nIter = 0;
  dx = 0;
  status = StatusTest::Unconverged;
}

void TrustRegionBased::
reset(const NOX::Abstract::Vector& initialGuess)
{
  // New initial guess and status test
  solnPtr->setX(initialGuess);

  // Initialize
  nIter = 0;
  dx = 0;
  status = StatusTest::Unconverged;
}

void TrustRegionBased::
reset()
{
  // Initialize
  nIter = 0;
  dx = 0;
  status = StatusTest::Unconverged;
}

TrustRegionBased::~TrustRegionBased()
{

}


NOX::StatusTest::StatusType TrustRegionBased::getStatus() const
{
  return status;
}

NOX::StatusTest::StatusType TrustRegionBased::step()
{
  observer->runPreIterate(*this);

  if (nIter == 0) {
    globalDataPtr->getNonConstSolverStatistics()->incrementNumNonlinearSolves();

    // Compute F of initital guess
    solnPtr->computeF();
    newF = meritFuncPtr->computef(*solnPtr);

    // Test the initial guess
    status = testPtr->checkStatus(*this, checkType);

    printUpdate();
  }

  // First check status
  if (status != StatusTest::Unconverged)
    return status;

  // Copy pointers into temporary references
  Abstract::Group& soln = *solnPtr;
  StatusTest::Generic& test = *testPtr;

  // Compute Cauchy and Newton points
  bool ok;
  ok = newtonPtr->compute(*newtonVecPtr, soln, *this);
  if (!ok)
  {
    utilsPtr->out() << "NOX::Solver::TrustRegionBased::iterate - unable to calculate Newton direction" << std::endl;
    status = StatusTest::Failed;
    observer->runPostIterate(*this);
    printUpdate();
    return status;
  }

  ok = cauchyPtr->compute(*cauchyVecPtr, soln, *this);
  if (!ok)
  {
    utilsPtr->err() << "NOX::Solver::TrustRegionBased::iterate - unable to calculate Cauchy direction" << std::endl;
    status = StatusTest::Failed;
    observer->runPostIterate(*this);
    printUpdate();
    return status;
  }

  if (nIter == 0)
  {
    radius = newtonVecPtr->norm();
    radius = std::min( initRadius, radius ) ;

    if (radius < minRadius)
      radius = 2 * minRadius;
  }

  // Update iteration count.
  nIter ++;
  globalDataPtr->getNonConstSolverStatistics()->incrementNumNonlinearIterations();

  // Copy current soln to the old soln.
  *oldSolnPtr = *solnPtr;
  // RPP: Can't just copy over oldf.  Scaling could change between iterations
  // so user Merit Functions could be out of sync
  oldF = meritFuncPtr->computef(*oldSolnPtr);

  // Improvement ratio = (oldF - newF) / (mold - mnew)
  double ratio = -1;

  if (utilsPtr->isPrintType(NOX::Utils::InnerIteration))
  {
    utilsPtr->out() << NOX::Utils::fill(72) << std::endl;
    utilsPtr->out() << "-- Trust Region Inner Iteration --" << std::endl;
  }

  // Trust region subproblem loop
  while ((ratio < minRatio) && (radius > minRadius))
  {

    Teuchos::RCP<NOX::Abstract::Vector> dirPtr;
    double step;

    // Trust region step
    double newtonVecNorm = newtonVecPtr->norm();
    double cauchyVecNorm = cauchyVecPtr->norm();

    if (newtonVecNorm <= radius)
    {
      stepType = TrustRegionBased::Newton;
      step = 1.0;
      dirPtr = newtonVecPtr;
    }
    else if (cauchyVecNorm >= radius)
    {
      stepType = TrustRegionBased::Cauchy;
      step = radius / cauchyVecNorm;
      dirPtr = cauchyVecPtr;
    }
    else
    {			// Dogleg computation

      // aVec = newtonVec - cauchyVec
      aVecPtr->update(1.0, *newtonVecPtr, -1.0, *cauchyVecPtr, 0.0);

      // cta = cauchyVec' * aVec
      double cta = cauchyVecPtr->innerProduct(*aVecPtr);
      // ctc = cauchyVec' * cauchyVec
      double ctc = cauchyVecPtr->innerProduct(*cauchyVecPtr);
      // ata = aVec' * aVec
      double ata = aVecPtr->innerProduct(*aVecPtr);

      // sqrt of quadratic equation
      double tmp = (cta * cta) - ((ctc - (radius * radius)) * ata);
      if (tmp < 0) {
	utilsPtr->err() << "NOX::Solver::TrustRegionBased::iterate - invalid computation" << std::endl;
	throw std::runtime_error("NOX Error");
      }

      // final soln to quadratic equation
      double gamma = (sqrt(tmp) - cta) / ata;
      if ((gamma < 0) || (gamma > 1)) {
	utilsPtr->err() << "NOX::Solver::TrustRegionBased::iterate - invalid trust region step" << std::endl;
	throw std::runtime_error("NOX Error");
      }

      // final direction computation
      aVecPtr->update(1.0 - gamma, *cauchyVecPtr, gamma, *newtonVecPtr, 0.0);

      // solution
      stepType = TrustRegionBased::Dogleg;
      dirPtr = aVecPtr;
      step = 1.0;
    }

    // Local reference to use in the remaining computation
    const Abstract::Vector& dir = *dirPtr;

    // Compute new X
    observer->runPreSolutionUpdate(dir,*this);
    soln.computeX(*oldSolnPtr, dir, step);
    observer->runPostSolutionUpdate(*this);

    // Calculate true step length
    dx = step * dir.norm();

    // Compute F for new current solution.
    NOX::Abstract::Group::ReturnType rtype = soln.computeF();
    if (rtype != NOX::Abstract::Group::Ok)
    {
      utilsPtr->err() << "NOX::Solver::TrustRegionBased::iterate - unable to compute F" << std::endl;
      throw std::runtime_error("NOX Error");
    }

    // Compute ratio of actual to predicted reduction
    // If using Homer Walker's Ared/Pred ratio computation,
    // we use F, NOT the merit function, f.
    if (useAredPredRatio) {

      // bVec = F(x) + J d
      rtype = oldSolnPtr->applyJacobian(*dirPtr, *bVecPtr);
      if (rtype != NOX::Abstract::Group::Ok)
      {
	utilsPtr->out() << "NOX::Solver::TrustRegionBased::iterate - "
	     << "unable to compute F" << std::endl;
	throw std::runtime_error("NOX Error");
      }
      bVecPtr->update(1.0, oldSolnPtr->getF(), step);

      // Compute norms
      double oldNormF = oldSolnPtr->getNormF();
      double newNormF = soln.getNormF();
      double normFLinear = bVecPtr->norm();

      ratio = (oldNormF - newNormF) / (oldNormF - normFLinear);

      // Print the ratio values if requested
      if (utilsPtr->isPrintType(NOX::Utils::InnerIteration)) {
	double numerator = oldNormF - newNormF;
	double denominator = oldNormF - normFLinear;
	utilsPtr->out() << "Ratio computation: "
			<< utilsPtr->sciformat(numerator) << "/"
			<< utilsPtr->sciformat(denominator) << "="
			<< ratio << std::endl;
      }

      // Update the merit function (newF used when printing iteration status)
      newF = meritFuncPtr->computef(*solnPtr);

    }
    else {  // Default ratio computation

      newF = meritFuncPtr->computef(*solnPtr);

      if (newF >= oldF)
      {
	ratio = -1;
      }
      else
      {

        // use bVecPtr for scratch space
        bVecPtr->update(step, dir, 0.0);

	double numerator = oldF - newF;
	double denominator = 0.0;

	denominator = fabs(oldF - meritFuncPtr->
			   computeQuadraticModel(*(bVecPtr.get()),*oldSolnPtr));

	ratio = numerator / denominator;
	if (utilsPtr->isPrintType(NOX::Utils::Debug))
	  utilsPtr->out() << "Ratio computation: "
			  << utilsPtr->sciformat(numerator) << "/"
			  << utilsPtr->sciformat(denominator) << "="
			  << utilsPtr->sciformat(ratio) << std::endl;

	// WHY IS THIS CHECK HERE?
	if ((denominator < 1.0e-12) && ((newF / oldF) >= 0.5))
	  ratio = -1;
      }
    }

    if (utilsPtr->isPrintType(Utils::InnerIteration)) {
      utilsPtr->out() << "radius = " << utilsPtr->sciformat(radius, 1);
      utilsPtr->out() << " ratio = " << std::setprecision(2) << std::setw(4) << ratio;
      utilsPtr->out() << " f = " << utilsPtr->sciformat(sqrt(2*newF));
      utilsPtr->out() << " old f = " << utilsPtr->sciformat(sqrt(2*oldF));
      utilsPtr->out() << " ";

      switch(stepType) {
      case TrustRegionBased::Newton:
	utilsPtr->out() << "Newton";
	break;
      case TrustRegionBased::Cauchy:
	utilsPtr->out() << "Cauchy";
	break;
      case TrustRegionBased::Dogleg:
	utilsPtr->out() << "Dogleg";
	break;
      }

      utilsPtr->out() << std::endl;
    }

    // Update trust region
    if (ratio < contractTriggerRatio)
    {
      if (stepType == TrustRegionBased::Newton) {
	radius = newtonVecPtr->norm();
      }
      radius = NOX_MAX(contractFactor * radius, minRadius);
    }
    else if ((ratio > expandTriggerRatio) && (dx == radius))
    {
      radius = NOX_MIN(expandFactor * radius, maxRadius);
    }

  }


  // Evaluate the current status
  if ((radius <= minRadius) && (ratio < minRatio))
  {
    if (utilsPtr->isPrintType(Utils::InnerIteration))
      utilsPtr->out() << "Using recovery step and resetting trust region." << std::endl;
    soln.computeX(*oldSolnPtr, *newtonVecPtr, recoveryStep);
    soln.computeF();
    radius = newtonVecPtr->norm();
    /*if (radius < minRadius)
      radius = 2 * minRadius;*/
  }

  status = test.checkStatus(*this, checkType);

  if (utilsPtr->isPrintType(Utils::InnerIteration))
    utilsPtr->out() << NOX::Utils::fill(72) << std::endl;

  observer->runPostIterate(*this);

  printUpdate();

  return status;
}

NOX::StatusTest::StatusType TrustRegionBased::solve()
{
  observer->runPreSolve(*this);

  this->reset();

  // Iterate until converged or failed
  while (status == StatusTest::Unconverged) {
    status = step();
  }

  Teuchos::ParameterList& outputParams = paramsPtr->sublist("Output");
  outputParams.set("Nonlinear Iterations", nIter);
  outputParams.set("2-Norm of Residual", solnPtr->getNormF());

  observer->runPostSolve(*this);

  return status;
}

const Abstract::Group& TrustRegionBased::getSolutionGroup() const
{
  return *solnPtr;
}

const Abstract::Group& TrustRegionBased::getPreviousSolutionGroup() const
{
  return *oldSolnPtr;
}

int TrustRegionBased::getNumIterations() const
{
  return nIter;
}

const Teuchos::ParameterList& TrustRegionBased::getList() const
{
  return *paramsPtr;
}

Teuchos::RCP<const NOX::SolverStats>
NOX::Solver::TrustRegionBased::getSolverStatistics() const
{
  return globalDataPtr->getSolverStatistics();
}

// protected
void TrustRegionBased::printUpdate()
{
  // Print the status test parameters at each iteration if requested
  if ((status == StatusTest::Unconverged) &&
      (utilsPtr->isPrintType(NOX::Utils::OuterIterationStatusTest))) {
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
  }

  double fmax = solnPtr->getF().norm(Abstract::Vector::MaxNorm);
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration)) {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Newton Trust-Region Step " << nIter << " -- \n";
    utilsPtr->out() << "f = " << utilsPtr->sciformat(sqrt(2*newF));
    utilsPtr->out() << " fmax = " << utilsPtr->sciformat(fmax);
    utilsPtr->out() << "  dx = " << utilsPtr->sciformat(dx);
    utilsPtr->out() << "  radius = " << utilsPtr->sciformat(radius);
    if (status == StatusTest::Converged)
      utilsPtr->out() << " (Converged!)";
    if (status == StatusTest::Failed)
      utilsPtr->out() << " (Failed!)";
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n" << std::endl;
  }

  if ((status != StatusTest::Unconverged) &&
      (utilsPtr->isPrintType(NOX::Utils::OuterIteration))) {
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Final Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
  }
}
