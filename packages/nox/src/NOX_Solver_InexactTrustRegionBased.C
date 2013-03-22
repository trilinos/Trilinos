// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_Common.H"
#include "NOX_Solver_InexactTrustRegionBased.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"
#include "NOX_MeritFunction_Generic.H"
#include "NOX_Solver_PrePostOperator.H"
#include "NOX_Solver_SolverUtils.H"
#include "NOX_Direction_Generic.H"
#include "NOX_Direction_Factory.H"

using namespace NOX;
using namespace NOX::Solver;

//*************************************************************************
//**** Constructor
//*************************************************************************
NOX::Solver::InexactTrustRegionBased::
InexactTrustRegionBased(const Teuchos::RCP<NOX::Abstract::Group>& grp, 
			const Teuchos::RCP<NOX::StatusTest::Generic>& t, 
			const Teuchos::RCP<Teuchos::ParameterList>& p) :
  globalDataPtr(Teuchos::rcp(new NOX::GlobalData(p))),
  utils(globalDataPtr->getUtils()), 
  solnPtr(grp),		// pointer to grp
  oldSolnPtr(grp->clone(DeepCopy)), // create via clone
  newtonVecPtr(grp->getX().clone(ShapeCopy)), // create via clone 
  cauchyVecPtr(grp->getX().clone(ShapeCopy)), // create via clone 
  rCauchyVecPtr(grp->getX().clone(ShapeCopy)), // create via clone 
  residualVecPtr(grp->getX().clone(ShapeCopy)), // create via clone 
  aVecPtr(grp->getX().clone(ShapeCopy)), // create via clone 
  bVecPtr(grp->getX().clone(ShapeCopy)), // create via clone 
  testPtr(t),			// pointer to t
  paramsPtr(p),			// copy p
  inNewtonUtils(globalDataPtr, paramsPtr->sublist("Direction")),
  radius(0.0),
  meritFuncPtr(globalDataPtr->getMeritFunction()),
  useCauchyInNewtonDirection(false),
  writeOutputParamsToList(true),
  useCounters(true),
  numCauchySteps(0),
  numNewtonSteps(0),
  numDoglegSteps(0),
  numTrustRegionInnerIterations(0),
  sumDoglegFracCauchyToNewton(0.0),
  sumDoglegFracNewtonLength(0.0),
  useAredPredRatio(false),
  useDoglegMinimization(false),
  prePostOperator(utils, paramsPtr->sublist("Solver Options"))
{
  init();
}

//*************************************************************************
//**** Destructor
//*************************************************************************
NOX::Solver::InexactTrustRegionBased::~InexactTrustRegionBased() 
{

}

//*************************************************************************
//**** init
//*************************************************************************
// Protected
void NOX::Solver::InexactTrustRegionBased::init()
{
  // Initialize 
  nIter = 0;
  dx = 0.0;
  status = StatusTest::Unconverged;
  if (useCounters)
    resetCounters();

  checkType = parseStatusTestCheckType(paramsPtr->sublist("Solver Options"));

  // Print out initialization information
  if (utils->isPrintType(NOX::Utils::Parameters)) {
    utils->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utils->out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utils->out(),5);
  }

  // Get the trust region method
  std::string methodChoice = 
    paramsPtr->sublist("Trust Region").
    get("Inner Iteration Method", "Inexact Trust Region");
  if (methodChoice == "Standard Trust Region")
    method = Standard;
  else if (methodChoice == "Inexact Trust Region")
    method = Inexact;
  else {
    utils->err() << "NOX::Solver::InexactTrustRegionBased::init - \"" << methodChoice
	 << "\" is an invalid choice for \"Method\" key!"
	 << std::endl;
    throw "NOX Error";
  }

  // Set default parameter settings using get() if they are not set
  // Default directions 
  paramsPtr->sublist("Direction").get("Method", "Newton");
  paramsPtr->sublist("Cauchy Direction")
    .get("Method", "Steepest Descent");
  paramsPtr->sublist("Cauchy Direction").sublist("Steepest Descent")
    .get("Scaling Type", "Quadratic Model Min");

  newtonPtr = NOX::Direction::
    buildDirection(globalDataPtr, paramsPtr->sublist("Direction"));
  cauchyPtr = NOX::Direction::
    buildDirection(globalDataPtr, paramsPtr->sublist("Cauchy Direction"));
  inNewtonUtils.reset(globalDataPtr, paramsPtr->sublist("Direction"));

  minRadius = paramsPtr->sublist("Trust Region")
    .get("Minimum Trust Region Radius", 1.0e-6);
  if (minRadius <= 0.0) 
    invalid("Minimum Trust Region Radius", minRadius);

  maxRadius = paramsPtr->sublist("Trust Region")
    .get("Maximum Trust Region Radius", 1.0e+10);
  if (maxRadius <= minRadius) 
    invalid("Maximum Trust Region Radius", maxRadius);

  minRatio = paramsPtr->sublist("Trust Region")
    .get("Minimum Improvement Ratio", 1.0e-4);
  if (minRatio <= 0.0) 
    invalid("Minimum Improvement Ratio", minRatio);

  contractTriggerRatio = paramsPtr->sublist("Trust Region")
    .get("Contraction Trigger Ratio", 0.1);
  if (contractTriggerRatio < minRatio) 
    invalid("Contraction Trigger Ratio", contractTriggerRatio);

  expandTriggerRatio = paramsPtr->sublist("Trust Region")
    .get("Expansion Trigger Ratio", 0.75);
  if (expandTriggerRatio <= contractTriggerRatio) 
    invalid("Expansion Trigger Ratio", expandTriggerRatio);

  contractFactor = paramsPtr->sublist("Trust Region")
    .get("Contraction Factor", 0.25);
  if ((contractFactor <= 0.0) || (contractFactor >= 1)) 
    invalid("Contraction Factor", contractFactor);

  expandFactor = paramsPtr->sublist("Trust Region")
    .get("Expansion Factor", 4.0);
  if (expandFactor <= 1.0) 
    invalid("Expansion Factor", expandFactor);

  recoveryStep = paramsPtr->sublist("Trust Region")
    .get("Recovery Step", 1.0);
  if (recoveryStep < 0.0) 
    invalid("Recovery Step", recoveryStep);

  useCauchyInNewtonDirection = paramsPtr->sublist("Trust Region")
    .get("Use Cauchy in Newton Direction", false);

  // Check for using Homer Walker's Ared/Pred ratio calculation
  useAredPredRatio = paramsPtr->sublist("Trust Region")
    .get("Use Ared/Pred Ratio Calculation", false);

  // Check for dogleg minimization routine (only vaild for inexact algorithm)
  useDoglegMinimization = paramsPtr->sublist("Trust Region")
    .get("Use Dogleg Segment Minimization", false);

  // Check for statistics tracking
  useCounters = paramsPtr->sublist("Trust Region")
    .get("Use Counters", true);

  // Check for writing statistics to the parameter list
  useCounters = paramsPtr->sublist("Trust Region")
    .get("Write Output Parameters", true);

}

//*************************************************************************
//**** invalid
//*************************************************************************
void NOX::Solver::InexactTrustRegionBased::invalid(const std::string& name, 
						   double value) const
{
  utils->out() << "NOX::Solver::InexactTrustRegionBased::init - " 
       << "Invalid \"" << name << "\" (" << value << ")" 
       << std::endl;
  throw "NOX Error";
}

//*************************************************************************
//**** throwError
//*************************************************************************
void NOX::Solver::InexactTrustRegionBased::throwError(const std::string& method, 
						      const std::string& message) const
{
  utils->out() << "NOX::Solver::InexactTrustRegionBased::" << method << " - " 
       << message << std::endl;
  throw "NOX Error";
}

//*************************************************************************
//**** reset (without reparsing of parameter list)
//*************************************************************************
void NOX::Solver::InexactTrustRegionBased::
reset(const NOX::Abstract::Vector& initialGuess, 
      const Teuchos::RCP<NOX::StatusTest::Generic>& t)
{
  solnPtr->setX(initialGuess);
  testPtr = t;

  // Initialize 
  nIter = 0;
  dx = 0.0;
  status = StatusTest::Unconverged;
  if (useCounters)
    resetCounters();

  // Print out initialization information
  if (utils->isPrintType(NOX::Utils::Parameters)) {
    utils->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utils->out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utils->out(),5);
  }

  // Compute F of initital guess
  solnPtr->computeF();
  newF = meritFuncPtr->computef(*solnPtr);

  // Test the initial guess
  status = testPtr->checkStatus(*this, checkType);

  if (utils->isPrintType(NOX::Utils::Parameters)) {
    utils->out() << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testPtr->print(utils->out(), 5);
    utils->out() <<"\n" << NOX::Utils::fill(72) << "\n";
  }
}

//*************************************************************************
//**** reset (without reparsing of parameter list or status tests)
//*************************************************************************
void NOX::Solver::InexactTrustRegionBased::
reset(const NOX::Abstract::Vector& initialGuess)
{
  solnPtr->setX(initialGuess);

  // Initialize 
  nIter = 0;
  dx = 0.0;
  status = StatusTest::Unconverged;
  if (useCounters)
    resetCounters();

  // Print out initialization information
  if (utils->isPrintType(NOX::Utils::Parameters)) {
    utils->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utils->out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utils->out(),5);
  }

  // Compute F of initital guess
  solnPtr->computeF();
  newF = meritFuncPtr->computef(*solnPtr);

  // Test the initial guess
  status = testPtr->checkStatus(*this, checkType);

  if (utils->isPrintType(NOX::Utils::Parameters)) {
    utils->out() << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testPtr->print(utils->out(), 5);
    utils->out() <<"\n" << NOX::Utils::fill(72) << "\n";
  }
}

//*************************************************************************
//**** getStatus
//*************************************************************************
NOX::StatusTest::StatusType NOX::Solver::InexactTrustRegionBased::getStatus()
{
  return status;
}

//*************************************************************************
//**** step
//*************************************************************************
NOX::StatusTest::StatusType NOX::Solver::InexactTrustRegionBased::step()
{
  prePostOperator.runPreIterate(*this);

  if (nIter == 0) {
    // Compute F of initital guess
    solnPtr->computeF();
    newF = meritFuncPtr->computef(*solnPtr);

    // Test the initial guess
    status = testPtr->checkStatus(*this, checkType);
    
    printUpdate();
  }

  NOX::StatusTest::StatusType status = NOX::StatusTest::Unconverged;

  switch(method) {
  case Standard:
    status = iterateStandard();
    break;

  case Inexact:
    status = iterateInexact();
    break;

  default:
    status = iterateInexact();  // defaults to non-minimized model
    break;
  }
  
  prePostOperator.runPostIterate(*this);

  printUpdate();

  return status;
}

//*************************************************************************
//**** iterateStandard
//*************************************************************************
NOX::StatusTest::StatusType 
NOX::Solver::InexactTrustRegionBased::iterateStandard()
{
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
    utils->out() << "NOX::Solver::InexactTrustRegionBased::iterate - "
	 << "unable to calculate Newton direction" << std::endl;
    status = StatusTest::Failed;
    return status;
  }

  ok = cauchyPtr->compute(*cauchyVecPtr, soln, *this);
  if (!ok) 
  {
    utils->out() << "NOX::Solver::InexactTrustRegionBased::iterate - "
	 << "unable to calculate Cauchy direction" << std::endl;
    status = StatusTest::Failed;
    return status;
  }

  if (nIter == 0) 
  {
    radius = computeNorm(*newtonVecPtr);

    if (radius < minRadius)
      radius = 2 * minRadius;
  }

  // Update iteration count.
  nIter ++;

  // Copy current soln to the old soln.
  *oldSolnPtr = *solnPtr;
  // RPP: Can't just copy over oldf.  Scaling could change between iterations
  // for user defined Merit Functions and throw things out of sync
  //oldF = computeMeritFunction(*oldSolnPtr);
  oldF = meritFuncPtr->computef(*oldSolnPtr);

  // Improvement ratio = (oldF - newF) / (mold - mnew)
  double ratio = -1;

  if (utils->isPrintType(NOX::Utils::InnerIteration)) 
  {
    utils->out() << NOX::Utils::fill(72) << std::endl;
    utils->out() << "-- Trust Region Inner Iteration --" << std::endl;
  }

  // Dogleg variable
  double gamma = 0.0;

  // Trust region subproblem loop
  while ((ratio < minRatio) && (radius > minRadius)) 
  {
    if (useCounters)
      numTrustRegionInnerIterations += 1;

    Teuchos::RCP<NOX::Abstract::Vector> dirPtr;
    double step;

    // Trust region step
    double newtonVecNorm = computeNorm(*newtonVecPtr);
    double cauchyVecNorm = computeNorm(*cauchyVecPtr);
    
    if (newtonVecNorm <= radius) 
    {
      stepType = InexactTrustRegionBased::Newton;
      step = 1.0;
      dirPtr = newtonVecPtr;
    }
    else if (cauchyVecNorm >= radius) 
    {
      stepType = InexactTrustRegionBased::Cauchy;
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
	utils->err() << "NOX::Solver::InexactTrustRegionBased::iterate - invalid computation" << std::endl;
	throw "NOX Error";
      }
      
      // final soln to quadratic equation
      gamma = (sqrt(tmp) - cta) / ata;
      if ((gamma < 0) || (gamma > 1)) {
	utils->err() << "NOX::Solver::InexactTrustRegionBased::iterate - invalid trust region step" << std::endl;
	throw "NOX Error";
      }
      
      // final direction computation
      aVecPtr->update(1.0 - gamma, *cauchyVecPtr, gamma, *newtonVecPtr, 0.0);

      // solution
      stepType = InexactTrustRegionBased::Dogleg;
      dirPtr = aVecPtr;
      step = 1.0;
    }
    
    // Local reference to use in the remaining computation
    const NOX::Abstract::Vector& dir = *dirPtr;

    // Calculate true step length
    dx = step * (computeNorm(dir));
    
    // Compute new X
    soln.computeX(*oldSolnPtr, dir, step);

    // Compute F for new current solution.
    NOX::Abstract::Group::ReturnType rtype = soln.computeF();
    if (rtype != NOX::Abstract::Group::Ok) 
    {
      utils->err() << "NOX::Solver::InexactTrustRegionBased::iterate - "
	   << "unable to compute F" << std::endl;
      throw "NOX Error";
    }

    // Compute ratio of actual to predicted reduction
    // If using Homer Walker's Ared/Pred ratio computation, 
    // we use the residual F, NOT the merit function, f.
    if (useAredPredRatio) {

      // bVec = F(x) + J d
      rtype = oldSolnPtr->applyJacobian(*dirPtr, *bVecPtr);
      if (rtype != NOX::Abstract::Group::Ok) 
      {
	utils->out() << "NOX::Solver::TrustRegionBased::iterate - "
	     << "unable to compute F" << std::endl;
	throw "NOX Error";
      }
      bVecPtr->update(1.0, oldSolnPtr->getF(), 1.0);

      // Compute norms
      double oldNormF = computeNorm(oldSolnPtr->getF());
      double newNormF = computeNorm(solnPtr->getF());
      double normFLinear = computeNorm(*bVecPtr);;
      
      ratio = (oldNormF - newNormF) / (oldNormF - normFLinear);

      // Print the ratio values if requested
      if (utils->isPrintType(NOX::Utils::InnerIteration)) {
      double numerator = oldNormF - newNormF;
      double denominator = oldNormF - normFLinear;
	utils->out() << "Ratio computation: " << utils->sciformat(numerator) 
		     << "/" << utils->sciformat(denominator) << "=" 
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
	
	rtype = oldSolnPtr->applyJacobian(*dirPtr, *bVecPtr);
	if (rtype != NOX::Abstract::Group::Ok) 
	{
	  utils->err() << "NOX::Solver::InexactTrustRegionBased::iterate - "
	       << "unable to compute F" << std::endl;
	  throw "NOX Error";
	}

	double numerator = oldF - newF;
	double denominator = 0.0;

	if (!Teuchos::is_null(meritFuncPtr)) {
	  denominator = fabs(oldF - meritFuncPtr->
			     computeQuadraticModel(dir, *oldSolnPtr));
	}
	else 
	  denominator = fabs(dir.innerProduct(oldSolnPtr->getGradient()) + 
			     0.5 * bVecPtr->innerProduct(*bVecPtr));

	ratio = numerator / denominator;
	if (utils->isPrintType(NOX::Utils::InnerIteration))
	  utils->out() << "Ratio computation: " << utils->sciformat(numerator) << "/" 
	       << utils->sciformat(denominator) << "=" << ratio << std::endl;

	// WHY IS THIS CHECK HERE?
	if ((denominator < 1.0e-12) && ((newF / oldF) >= 0.5))
	  ratio = -1;
      }
    }  // end ratio computation

    if (utils->isPrintType(Utils::InnerIteration)) {
      utils->out() << "radius = " << utils->sciformat(radius, 1);
      utils->out() << " ratio = " << std::setprecision(1) << std::setw(3) << ratio;
      utils->out() << " f = " << utils->sciformat(sqrt(2*newF));
      utils->out() << " oldF = " << utils->sciformat(sqrt(2*oldF));
      utils->out() << " ";

      switch(stepType) {
      case InexactTrustRegionBased::Newton:
	utils->out() << "Newton";
	break;
      case InexactTrustRegionBased::Cauchy:
	utils->out() << "Cauchy";
	break;
      case InexactTrustRegionBased::Dogleg:
	utils->out() << "Dogleg";
	break;
      }

      utils->out() << std::endl;
    }

    // Update trust region
    if (ratio < contractTriggerRatio) 
    {
      if (stepType == InexactTrustRegionBased::Newton) {
	radius = computeNorm(*newtonVecPtr);
      }
      radius = NOX_MAX(contractFactor * radius, minRadius);
    }
    else if ((ratio > expandTriggerRatio) && (dx == radius))
      // RPP Hack
      //else if (ratio > expandTriggerRatio) 
    {
      radius = NOX_MIN(expandFactor * radius, maxRadius);
    }

  } // End of While loop over TR inner iteration

  // Update Counters
  if (useCounters) {
    if (stepType == InexactTrustRegionBased::Newton)
      numNewtonSteps += 1;
    else if (stepType == InexactTrustRegionBased::Cauchy) 
      numCauchySteps += 1;
    else if (stepType == InexactTrustRegionBased::Dogleg) {
      numDoglegSteps += 1;
      sumDoglegFracCauchyToNewton += gamma;
      double tmp = radius/computeNorm(*newtonVecPtr);
      sumDoglegFracNewtonLength += tmp;

      if (utils->isPrintType(Utils::Details)) {
	utils->out() << "    Fraction of Newton Step Length = " << tmp << std::endl;
	utils->out() << "    Fraction Between Cauchy and Newton Direction = " 
	     << gamma << std::endl; 
      }

    }
  }

  // Evaluate the current status
  if ((radius <= minRadius) && (ratio < minRatio)) 
  {
    if (utils->isPrintType(Utils::InnerIteration))
      utils->out() << "Using recovery step and resetting trust region." << std::endl;
    solnPtr->computeX(*oldSolnPtr, *newtonVecPtr, recoveryStep);
    solnPtr->computeF();
    radius = computeNorm(*newtonVecPtr);
    /*if (radius < minRadius)
      radius = 2 * minRadius;*/
  }

  status = test.checkStatus(*this, checkType);
 
  if (utils->isPrintType(Utils::InnerIteration)) 
    utils->out() << NOX::Utils::fill(72) << std::endl;

  // Return status.
  return status;
}

//*************************************************************************
//**** interateInexact
//*************************************************************************
NOX::StatusTest::StatusType 
NOX::Solver::InexactTrustRegionBased::iterateInexact()
{
  // First check the current status
  if (status != StatusTest::Unconverged) 
    return status;

  // Copy pointers into temporary references
  NOX::Abstract::Group& soln = *solnPtr;
  NOX::Abstract::Group& oldSoln = *oldSolnPtr;
  NOX::Abstract::Vector& newtonVec = *newtonVecPtr;
  NOX::Abstract::Vector& cauchyVec = *cauchyVecPtr;
  NOX::Abstract::Vector& rCauchyVec = *rCauchyVecPtr;
  NOX::Abstract::Vector& residualVec = *residualVecPtr;
  NOX::Abstract::Vector& aVec = *aVecPtr;
  NOX::Abstract::Vector& bVec = *bVecPtr;
  NOX::StatusTest::Generic& test = *testPtr;

  // Newton direction eta (do before setting oldSoln = soln).
  eta_last = eta;
  eta = inNewtonUtils.computeForcingTerm(soln, oldSoln, nIter, 
					 *this, eta_last);
  
  // Copy current soln to the old soln.
  oldSoln = soln;

  // Compute Cauchy direction
  bool ok = cauchyPtr->compute(cauchyVec, oldSoln, *this);
  if (!ok) {
    utils->out() << "NOX::Solver::InexactTrustRegionBased::iterate - "
	 << "unable to calculate Cauchy direction" << std::endl;
    status = StatusTest::Failed;
    return status;
  }

  // RPP: Can't just copy over oldf.  Scaling could change between iterations
  // so user Merit Functions could be out of sync
  oldF = meritFuncPtr->computef(oldSoln);
  
  // Compute linear model residual for cauchy direction and tolerance
  oldSoln.applyJacobian(cauchyVec, rCauchyVec);
  rCauchyVec.update(1.0, oldSoln.getF(), 1.0);
  double rCauchy = computeNorm(rCauchyVec);
  double normF = computeNorm(oldSoln.getF());
  double etaCauchy = rCauchy/normF;
  
  // Improvement ratio = (oldF - newF) / (mold - mnew)
  double ratio = -1.0;

  // Initial Radius
  if (nIter == 0) 
      radius = maxRadius;

  if (utils->isPrintType(NOX::Utils::InnerIteration)) 
  {
    utils->out() << NOX::Utils::fill(72) << std::endl;
    utils->out() << "-- Trust Region Inner Iteration --" << std::endl;
  }

  // Set up variables needed during inner iteration.
  NOX::Abstract::Vector* dirPtr = 0;
  NOX::Abstract::Vector& doglegVec = bVec;
  NOX::Abstract::Vector& zVec = aVec;
  double step = 0.0;
  double tau = 0.0;
  double normS = computeNorm(cauchyVec);
  bool computedNewtonDir = false;
  innerIterationStatus = Unconverged;

  // Update iteration count.
  nIter ++;

  // While "s" is not acceptable:
  while (innerIterationStatus == Unconverged) {

    if (useCounters)
      numTrustRegionInnerIterations += 1;

    // Compute a direction and step length
    if (computeNorm(cauchyVec) >= radius) {
      stepType = InexactTrustRegionBased::Cauchy;
      step = radius / computeNorm(cauchyVec);
      dirPtr = &cauchyVec;
    }
    else {
      if (etaCauchy <= eta) {
	stepType = InexactTrustRegionBased::Cauchy;
	step = 1.0;
	dirPtr = &cauchyVec;
      }
      else {  // compute the Newton direction

	if (!computedNewtonDir) {
	  std::string directionMethod = paramsPtr->sublist("Direction").
	    get("Method", "Newton");
	  Teuchos::ParameterList& lsParams = paramsPtr->sublist("Direction").
	    sublist(directionMethod).sublist("Linear Solver");
	  lsParams.set("Tolerance", eta);
	  if (useCauchyInNewtonDirection)
	    newtonVec.update(1.0, cauchyVec, 0.0);
	  else
	    newtonVec.init(0.0);
	  if (!(oldSoln.isJacobian()))
	    oldSoln.computeJacobian();
	  oldSoln.applyJacobianInverse(lsParams, oldSoln.getF(), newtonVec);
	  newtonVec.scale(-1.0);
	  computedNewtonDir = true;
	}

	// Can we take a full Newton step?
	if ((computeNorm(newtonVec) <= radius) && (!useDoglegMinimization)) {
	  stepType = InexactTrustRegionBased::Newton;
	  step = 1.0;
	  dirPtr = &newtonVec;
	}
	else { 

	  // Compute z
	  zVec.update(1.0, newtonVec, -1.0, cauchyVec, 0.0);
	  
	  // Compute <s,z>
	  double sDotZ = 0.0; 
	  sDotZ = cauchyVec.innerProduct(zVec);

	  // Compute tau
	  double normZ = computeNorm(zVec);
	  if (sDotZ <= 0.0)
	    tau = (-1.0 * sDotZ + sqrt(sDotZ * sDotZ + normZ * normZ * 
		  (radius * radius - normS * normS)))/(normZ * normZ);
	  else
	    tau = (radius * radius - normS * normS)/
	      (sDotZ + sqrt(sDotZ * sDotZ + normZ * normZ * 
			    (radius * radius - normS * normS)));

	  // Adjust tau if using dogleg segment minimization
	  if (useDoglegMinimization) {
	    double tauMin = 1.0;
	    // residualVec = r_sd - r_n where r_i = F + J*s_i
	    oldSoln.applyJacobian(newtonVec, residualVec);
	    residualVec.update(1.0, oldSoln.getF(), 1.0);
	    residualVec.update(1.0, rCauchyVec, -1.0);
	    double norm = computeNorm(residualVec);
	    tauMin = rCauchyVec.innerProduct(residualVec) / (norm * norm);

	    tau = NOX_MIN(tauMin, tau);
	  }

	  // Compute dogleg step
	  doglegVec.update(1.0, cauchyVec, tau, zVec, 0.0);
	  stepType = InexactTrustRegionBased::Dogleg;
	  step = 1.0;
	  dirPtr = &doglegVec;
	} 
      }      
    } // Finished computing direction and step length

    // Check if the computed step is valid for the trust region
    //innerIterationStatus = checkStep();
    NOX::Abstract::Vector& dir = *dirPtr;

    // Update X and compute new F
    soln.computeX(oldSoln, dir, step);
    NOX::Abstract::Group::ReturnType rtype = soln.computeF();
    if (rtype != NOX::Abstract::Group::Ok)
      throwError("iterateInexact", "unable to compute F!");
    
    // Compute ratio of actual to predicted reduction
    newF = meritFuncPtr->computef(*solnPtr);

    if (newF >= oldF) 
      ratio = -1;
    else {

     // If using Homer Walker's Ared/Pred ratio computation, 
     // we use F, NOT the merit function, f.
     if (useAredPredRatio) {

       // bVec = F(x) + J d
       rtype = oldSoln.applyJacobian(*dirPtr, bVec);
       if (rtype != NOX::Abstract::Group::Ok) 
       {
	 utils->out() << "NOX::Solver::TrustRegionBased::iterate - "
	      << "unable to compute F" << std::endl;
	 throw "NOX Error";
       } 
       bVec.update(1.0, oldSoln.getF(), 1.0);

       // Compute norms
       double oldNormF = computeNorm(oldSoln.getF());
       double newNormF = computeNorm(soln.getF());
       double normFLinear = computeNorm(bVec);;
       
       ratio = (oldNormF - newNormF) / (oldNormF - normFLinear);
       
       // Print the ratio values if requested
       if (utils->isPrintType(NOX::Utils::InnerIteration)) {
	 double numerator = oldNormF - newNormF;
	 double denominator = oldNormF - normFLinear;
	 utils->out() << "Ratio computation: " << utils->sciformat(numerator) << "/" 
	      << utils->sciformat(denominator) << "=" << ratio << std::endl;
       }

       // Update the merit function (newF used when printing iteration status)
       newF = meritFuncPtr->computef(*solnPtr);

     }
     else {  // Default ratio computation

       rtype = oldSoln.applyJacobian(*dirPtr, zVec);
       if (rtype != NOX::Abstract::Group::Ok)
	 throwError("iterateInexact", "unable to applyJacobian!");
  
       double numerator = oldF - newF;
       double denominator = 0.0;

       denominator = fabs(oldF - meritFuncPtr->
			  computeQuadraticModel(dir,oldSoln));

       ratio = numerator / denominator;
       if (utils->isPrintType(NOX::Utils::InnerIteration))
	 utils->out() << "Ratio computation: " << utils->sciformat(numerator) << "/" 
	      << utils->sciformat(denominator) << "=" << ratio << std::endl;

       // WHY IS THIS CHECK HERE?
       if ((denominator < 1.0e-12) && ((newF / oldF) >= 0.5))
	 ratio = -1;
     }
    } 

    if (utils->isPrintType(Utils::InnerIteration)) {
      utils->out() << "radius = " << utils->sciformat(radius, 1);
      utils->out() << " ratio = " << std::setprecision(1) << std::setw(3) << ratio;
      utils->out() << " f = " << utils->sciformat(sqrt(2*newF));
      utils->out() << " oldF = " << utils->sciformat(sqrt(2*oldF));
      utils->out() << " ";

      switch(stepType) {
      case InexactTrustRegionBased::Newton:
	utils->out() << "Newton";
	break;
      case InexactTrustRegionBased::Cauchy:
	utils->out() << "Cauchy";
	break;
      case InexactTrustRegionBased::Dogleg:
	utils->out() << "Dogleg";
	break;
      }
      
      utils->out() << std::endl;
    }

    // Check the inner iteration status
    if (ratio >= minRatio) {
      innerIterationStatus = Converged;
    }
    else if ((radius <= minRadius) && (ratio < minRatio)) {
      innerIterationStatus = Failed; 
    }

    // Update trust region radius
    if (ratio < contractTriggerRatio) {

      if (stepType == InexactTrustRegionBased::Newton)
	radius = computeNorm(newtonVec);
  
      radius = NOX_MAX(contractFactor * radius, minRadius);
    }
    else if (ratio > expandTriggerRatio)
    {
      radius = NOX_MIN(expandFactor * radius, maxRadius);
    }

  } // End Trust region inner iteration

  // Update Counters
  if (useCounters) {
    if (stepType == InexactTrustRegionBased::Newton)
      numNewtonSteps += 1;
    else if (stepType == InexactTrustRegionBased::Cauchy) 
      numCauchySteps += 1;
    else if (stepType == InexactTrustRegionBased::Dogleg) {
      numDoglegSteps += 1;
      sumDoglegFracCauchyToNewton += tau;
      double tmp = radius/computeNorm(newtonVec);
      sumDoglegFracNewtonLength += tmp;

      if (utils->isPrintType(Utils::Details)) {
	utils->out() << "    Fraction of Newton Step Length = " << tmp << std::endl;
	utils->out() << "    Fraction Between Cauchy and Newton Direction = " 
	     << tau << std::endl; 
      }

    }
  }

  // Adjust the eta when taking a cauchy or dogleg step
  // This is equivalent to accounting for backtracking in inexact Newton
  // line searches.  Have to check theory w/ Homer Walker on this.
  // 3/10/04 This is the correct theory according to Homer!
  if (stepType == InexactTrustRegionBased::Cauchy)
    eta = 1.0 - (computeNorm(*dirPtr) / computeNorm(cauchyVec))*(1.0-etaCauchy);
  else if (stepType == InexactTrustRegionBased::Dogleg)
    eta = etaCauchy - tau * (etaCauchy - eta);

  // If the inner iteration failed, use a recovery step
  if (innerIterationStatus == Failed) {

    if (utils->isPrintType(Utils::InnerIteration))
      utils->out() << "Inner Iteration Failed!\n ";
    if (computedNewtonDir) {
      soln.computeX(oldSoln, cauchyVec, recoveryStep);
      utils->out() << "Using Newton recovery step and resetting trust region!" << std::endl;
    }
    else {
      soln.computeX(oldSoln, cauchyVec, recoveryStep);
      utils->out() << "Using Cauchy recovery step and resetting trust region!" << std::endl;
    }
    soln.computeF();
    radius = computeNorm(newtonVec);
  }

  // Evaluate the current status
  status = test.checkStatus(*this, checkType);
 
  if (utils->isPrintType(Utils::InnerIteration)) 
    utils->out() << NOX::Utils::fill(72) << std::endl;

  return status;    
}

//*************************************************************************
//**** checkStep
//*************************************************************************
NOX::StatusTest::StatusType 
NOX::Solver::InexactTrustRegionBased::
checkStep(const NOX::Abstract::Vector& step,
	  double& radius)
{
  return NOX::StatusTest::Converged;
}

//*************************************************************************
//**** solve
//*************************************************************************
NOX::StatusTest::StatusType NOX::Solver::InexactTrustRegionBased::solve()
{
  prePostOperator.runPreSolve(*this);

  // Iterate until converged or failed
  while (status == StatusTest::Unconverged) {
    status = step();
  }

  if (writeOutputParamsToList) {
    Teuchos::ParameterList& outputParams = paramsPtr->sublist("Output");
    outputParams.set("Nonlinear Iterations", nIter);
    outputParams.set("2-Norm of Residual", solnPtr->getNormF());
    if (useCounters) {
      Teuchos::ParameterList& trOutputParams = paramsPtr->
	sublist("Trust Region").sublist("Output");
      trOutputParams.set("Number of Cauchy Steps", numCauchySteps);
      trOutputParams.set("Number of Newton Steps", numNewtonSteps);
      trOutputParams.set("Number of Dogleg Steps", numDoglegSteps);
      trOutputParams.set("Number of Trust Region Inner Iterations", 
				  numTrustRegionInnerIterations);
      if (numDoglegSteps != 0) {
	trOutputParams.set("Dogleg Steps: Average Fraction of Newton Step Length", (sumDoglegFracNewtonLength/((double)numDoglegSteps)));
	trOutputParams.set("Dogleg Steps: Average Fraction Between Cauchy and Newton Direction", (sumDoglegFracCauchyToNewton/((double)numDoglegSteps)));
      }

    }
  }

  prePostOperator.runPostSolve(*this);

  return status;
}

//*************************************************************************
//**** getSolutionGroup
//*************************************************************************
const Abstract::Group& 
NOX::Solver::InexactTrustRegionBased::getSolutionGroup() const
{
  return *solnPtr;
}

//*************************************************************************
//**** getPreviousSolutionGroup
//*************************************************************************
const Abstract::Group& 
NOX::Solver::InexactTrustRegionBased::getPreviousSolutionGroup() const
{
  return *oldSolnPtr;
}

//*************************************************************************
//**** getNumIterations
//*************************************************************************
int NOX::Solver::InexactTrustRegionBased::getNumIterations() const
{
  return nIter;
}

//*************************************************************************
//**** getList
//*************************************************************************
const Teuchos::ParameterList& 
NOX::Solver::InexactTrustRegionBased::getList() const
{
  return *paramsPtr;
}

//*************************************************************************
//**** printUpdate
//*************************************************************************
// protected
void NOX::Solver::InexactTrustRegionBased::printUpdate() 
{
  // Print the status test parameters at each iteration if requested  
  if ((status == StatusTest::Unconverged) && 
      (utils->isPrintType(NOX::Utils::OuterIterationStatusTest))) {
    utils->out() << NOX::Utils::fill(72) << "\n";
    utils->out() << "-- Status Test Results --\n";    
    testPtr->print(utils->out());
    utils->out() << NOX::Utils::fill(72) << "\n";
  }
  
  double fmax = solnPtr->getF().norm(Abstract::Vector::MaxNorm);
  if (utils->isPrintType(NOX::Utils::OuterIteration)) {
    utils->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utils->out() << "-- Newton Trust-Region Step " << nIter << " -- \n";
    utils->out() << "f = " << utils->sciformat(sqrt(2*newF));
    utils->out() << " fmax = " << utils->sciformat(fmax);
    utils->out() << "  dx = " << utils->sciformat(dx);
    utils->out() << "  radius = " << utils->sciformat(radius);
    if (status == StatusTest::Converged)
      utils->out() << " (Converged!)";
    if (status == StatusTest::Failed)
      utils->out() << " (Failed!)";
    utils->out() << "\n" << NOX::Utils::fill(72) << "\n" << std::endl;
  }
  
  if ((status != StatusTest::Unconverged) && 
      (utils->isPrintType(NOX::Utils::OuterIteration))) {
    utils->out() << NOX::Utils::fill(72) << "\n";
    utils->out() << "-- Final Status Test Results --\n";    
    testPtr->print(utils->out());
    utils->out() << NOX::Utils::fill(72) << "\n";
  }
}

//*************************************************************************
//**** resetCounters
//*************************************************************************
void NOX::Solver::InexactTrustRegionBased::resetCounters()
{
  numCauchySteps = 0;
  numNewtonSteps = 0;
  numDoglegSteps = 0;
  numTrustRegionInnerIterations = 0;
  sumDoglegFracCauchyToNewton = 0.0;
  sumDoglegFracNewtonLength = 0.0;
  return;
}

//*************************************************************************
//**** computeNorm
//*************************************************************************
double NOX::Solver::InexactTrustRegionBased::
computeNorm(const NOX::Abstract::Vector& v)
{
  double norm = 0.0;
  norm = v.norm();
  return norm;
}
