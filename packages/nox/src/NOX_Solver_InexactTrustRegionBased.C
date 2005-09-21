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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifdef WITH_PRERELEASE

#include "NOX_Solver_InexactTrustRegionBased.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Parameter_List.H"
#include "NOX_Parameter_UserNorm.H"
#include "NOX_Parameter_MeritFunction.H"
#include "NOX_Solver_PrePostOperator.H"

using namespace NOX;
using namespace NOX::Solver;

//*************************************************************************
//**** Constructor
//*************************************************************************
NOX::Solver::InexactTrustRegionBased::
InexactTrustRegionBased(const Teuchos::RefCountPtr<Abstract::Group>& grp, 
			const Teuchos::RefCountPtr<StatusTest::Generic>& t, 
			const Teuchos::RefCountPtr<Parameter::List>& p) :
  solnPtr(grp),		// pointer to grp
  oldSolnPtr(grp->clone(DeepCopy)), // create via clone
  oldSoln(*oldSolnPtr),		// reference to just-created pointer
  newtonVecPtr(grp->getX().clone(ShapeCopy)), // create via clone 
  newtonVec(*newtonVecPtr),	// reference to just-created pointer
  cauchyVecPtr(grp->getX().clone(ShapeCopy)), // create via clone 
  cauchyVec(*cauchyVecPtr),	// reference to just-created pointer
  rCauchyVecPtr(grp->getX().clone(ShapeCopy)), // create via clone 
  rCauchyVec(*rCauchyVecPtr),		// reference to just-created pointer
  residualVecPtr(grp->getX().clone(ShapeCopy)), // create via clone 
  residualVec(*rCauchyVecPtr),		// reference to just-created pointer
  aVecPtr(grp->getX().clone(ShapeCopy)), // create via clone 
  aVec(*aVecPtr),		// reference to just-created pointer
  bVecPtr(grp->getX().clone(ShapeCopy)), // create via clone 
  bVec(*bVecPtr),		// reference to just-created pointer
  testPtr(t),			// pointer to t
  paramsPtr(p),			// copy p
  utils(paramsPtr->sublist("Printing")), // inititalize utils
  inNewtonUtils(utils, paramsPtr->sublist("Direction")),
  newton(utils),		// initialize direction
  cauchy(utils),		// initialize direction
  radius(0.0),
  userNormPtr(0),
  userMeritFuncPtr(0),
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
  delete oldSolnPtr;
  delete cauchyVecPtr;
  delete newtonVecPtr;
  delete rCauchyVecPtr;
  delete residualVecPtr;
  delete aVecPtr;
  delete bVecPtr;
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

  // Print out initialization information
  if (utils.isPrintType(NOX::Utils::Parameters)) {
    utils.out() << "\n" << NOX::Utils::fill(72) << "\n";
    utils.out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utils.out(),5);
  }

  // Get the trust region method
  string methodChoice = 
    paramsPtr->sublist("Trust Region").
    getParameter("Inner Iteration Method", "Inexact Trust Region");
  if (methodChoice == "Standard Trust Region")
    method = Standard;
  else if (methodChoice == "Inexact Trust Region")
    method = Inexact;
  else {
    utils.err() << "NOX::Solver::InexactTrustRegionBased::init - \"" << methodChoice
	 << "\" is an invalid choice for \"Method\" key!"
	 << endl;
    throw "NOX Error";
  }

  // Set default parameter settings using getParameter() if they are not set
  // Default directions 
  paramsPtr->sublist("Direction").getParameter("Method", "Newton");
  paramsPtr->sublist("Cauchy Direction")
    .getParameter("Method", "Steepest Descent");
  paramsPtr->sublist("Cauchy Direction").sublist("Steepest Descent")
    .getParameter("Scaling Type", "Quadratic Model Min");

  newton.reset(paramsPtr->sublist("Direction"));
  cauchy.reset(paramsPtr->sublist("Cauchy Direction"));

  minRadius = paramsPtr->sublist("Trust Region")
    .getParameter("Minimum Trust Region Radius", 1.0e-6);
  if (minRadius <= 0.0) 
    invalid("Minimum Trust Region Radius", minRadius);

  maxRadius = paramsPtr->sublist("Trust Region")
    .getParameter("Maximum Trust Region Radius", 1.0e+10);
  if (maxRadius <= minRadius) 
    invalid("Maximum Trust Region Radius", maxRadius);

  minRatio = paramsPtr->sublist("Trust Region")
    .getParameter("Minimum Improvement Ratio", 1.0e-4);
  if (minRatio <= 0.0) 
    invalid("Minimum Improvement Ratio", minRatio);

  contractTriggerRatio = paramsPtr->sublist("Trust Region")
    .getParameter("Contraction Trigger Ratio", 0.1);
  if (contractTriggerRatio < minRatio) 
    invalid("Contraction Trigger Ratio", contractTriggerRatio);

  expandTriggerRatio = paramsPtr->sublist("Trust Region")
    .getParameter("Expansion Trigger Ratio", 0.75);
  if (expandTriggerRatio <= contractTriggerRatio) 
    invalid("Expansion Trigger Ratio", expandTriggerRatio);

  contractFactor = paramsPtr->sublist("Trust Region")
    .getParameter("Contraction Factor", 0.25);
  if ((contractFactor <= 0.0) || (contractFactor >= 1)) 
    invalid("Contraction Factor", contractFactor);

  expandFactor = paramsPtr->sublist("Trust Region")
    .getParameter("Expansion Factor", 4.0);
  if (expandFactor <= 1.0) 
    invalid("Expansion Factor", expandFactor);

  recoveryStep = paramsPtr->sublist("Trust Region")
    .getParameter("Recovery Step", 1.0);
  if (recoveryStep < 0.0) 
    invalid("Recovery Step", recoveryStep);

  useCauchyInNewtonDirection = paramsPtr->sublist("Trust Region")
    .getParameter("Use Cauchy in Newton Direction", false);

  // Check for the using Homer Walker's Ared/Pred ratio calculation
  useAredPredRatio = paramsPtr->sublist("Trust Region")
    .getParameter("Use Ared/Pred Ratio Calculation", false);

  // Check for dogleg minimization routine (only vaild for inexact algorithm)
  useDoglegMinimization = paramsPtr->sublist("Trust Region")
    .getParameter("Use Dogleg Segment Minimization", false);

  // Check for statistics tracking
  useCounters = paramsPtr->sublist("Trust Region")
    .getParameter("Use Counters", true);

  // Check for writing statistics to the parameter list
  useCounters = paramsPtr->sublist("Trust Region")
    .getParameter("Write Output Parameters", true);

  // Check for a user defined Norm
  if (paramsPtr->sublist("Trust Region").
      isParameterArbitrary("User Defined Norm")) {
    const NOX::Parameter::UserNorm& un = 
      dynamic_cast<const NOX::Parameter::UserNorm&>(paramsPtr->
      sublist("Trust Region").getArbitraryParameter("User Defined Norm"));
    userNormPtr = const_cast<NOX::Parameter::UserNorm*>(&un);

    // RPP: Hack!!  Need to compute the row sum scaling vector 
    // so norms are up-to-date (fix this in Epetra objects later).
    // This only is needed for user defined norms at this point.
    solnPtr->computeF();
    solnPtr->computeJacobian();
    solnPtr->computeNewton(paramsPtr->sublist("Direction").
			   sublist("Newton").sublist("Linear Solver"));
    
  }

  // Check for a user defined Merit Function
  if (paramsPtr->sublist("Trust Region").
      isParameterArbitrary("User Defined Merit Function")) {
    const NOX::Parameter::MeritFunction& mf = 
      dynamic_cast<const NOX::Parameter::MeritFunction&>
      (paramsPtr->sublist("Trust Region").
       getArbitraryParameter("User Defined Merit Function"));
    userMeritFuncPtr = const_cast<NOX::Parameter::MeritFunction*>(&mf);
  }

  // Compute F of initital guess
  solnPtr->computeF();
  if (userMeritFuncPtr != 0) {
    newF = userMeritFuncPtr->computef(*solnPtr);
  }
  else 
    newF = 0.5 * solnPtr->getNormF() * solnPtr->getNormF();

  // Test the initial guess
  status = testPtr->checkStatus(*this);

  if (utils.isPrintType(NOX::Utils::Parameters)) {
    utils.out() << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testPtr->print(utils.out(), 5);
    utils.out() <<"\n" << NOX::Utils::fill(72) << "\n";
  }

}

//*************************************************************************
//**** invalid
//*************************************************************************
void NOX::Solver::InexactTrustRegionBased::invalid(const string& name, 
						   double value) const
{
  utils.out() << "NOX::Solver::InexactTrustRegionBased::init - " 
       << "Invalid \"" << name << "\" (" << value << ")" 
       << endl;
  throw "NOX Error";
}

//*************************************************************************
//**** throwError
//*************************************************************************
void NOX::Solver::InexactTrustRegionBased::throwError(const string& method, 
						      const string& message) const
{
  utils.out() << "NOX::Solver::InexactTrustRegionBased::" << method << " - " 
       << message << endl;
  throw "NOX Error";
}

//*************************************************************************
//**** reset
//*************************************************************************
bool NOX::Solver::InexactTrustRegionBased::
reset(const Teuchos::RefCountPtr<Abstract::Group>& grp, 
      const Teuchos::RefCountPtr<StatusTest::Generic>& t, 
      const Teuchos::RefCountPtr<Parameter::List>& p) 
{
  solnPtr = grp;
  testPtr = t;
  paramsPtr = p;			
  utils.reset(paramsPtr->sublist("Printing"));
  prePostOperator.reset(utils, paramsPtr->sublist("Solver Options"));
  init();
  return true;
}

//*************************************************************************
//**** reset (without reparsing of parameter list)
//*************************************************************************
bool NOX::Solver::InexactTrustRegionBased::
reset(const Teuchos::RefCountPtr<Abstract::Group>& grp, 
      const Teuchos::RefCountPtr<StatusTest::Generic>& t)
{
  solnPtr = grp;
  testPtr = t;

  // Initialize 
  nIter = 0;
  dx = 0.0;
  status = StatusTest::Unconverged;
  if (useCounters)
    resetCounters();

  // Print out initialization information
  if (utils.isPrintType(NOX::Utils::Parameters)) {
    utils.out() << "\n" << NOX::Utils::fill(72) << "\n";
    utils.out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utils.out(),5);
  }

  // Compute F of initital guess
  solnPtr->computeF();
  if (userMeritFuncPtr != 0) {
    newF = userMeritFuncPtr->computef(*solnPtr);
  }
  else 
    newF = 0.5 * solnPtr->getNormF() * solnPtr->getNormF();

  // Test the initial guess
  status = testPtr->checkStatus(*this);

  if (utils.isPrintType(NOX::Utils::Parameters)) {
    utils.out() << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testPtr->print(utils.out(), 5);
    utils.out() <<"\n" << NOX::Utils::fill(72) << "\n";
  }
  return true;
}

//*************************************************************************
//**** getStatus
//*************************************************************************
NOX::StatusTest::StatusType NOX::Solver::InexactTrustRegionBased::getStatus()
{
  return status;
}

//*************************************************************************
//**** iterate
//*************************************************************************
NOX::StatusTest::StatusType NOX::Solver::InexactTrustRegionBased::iterate()
{

  prePostOperator.runPreIterate(*this);

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
  ok = newton.compute(newtonVec, soln, *this);
  if (!ok) 
  {
    utils.out() << "NOX::Solver::InexactTrustRegionBased::iterate - "
	 << "unable to calculate Newton direction" << endl;
    status = StatusTest::Failed;
    return status;
  }

  ok = cauchy.compute(cauchyVec, soln, *this);
  if (!ok) 
  {
    utils.out() << "NOX::Solver::InexactTrustRegionBased::iterate - "
	 << "unable to calculate Cauchy direction" << endl;
    status = StatusTest::Failed;
    return status;
  }

  if (nIter == 0) 
  {
    radius = computeNorm(newtonVec);

    if (radius < minRadius)
      radius = 2 * minRadius;
  }

  // Update iteration count.
  nIter ++;

  // Copy current soln to the old soln.
  oldSoln = soln;
  // RPP: Can't just copy over oldf.  Scaling could change between iterations
  // for user defined Merit Functions and throw things out of sync
  oldF = computeMeritFunction(oldSoln);

  // Improvement ratio = (oldF - newF) / (mold - mnew)
  double ratio = -1;

  if (utils.isPrintType(NOX::Utils::InnerIteration)) 
  {
    utils.out() << NOX::Utils::fill(72) << endl;
    utils.out() << "-- Trust Region Inner Iteration --" << endl;
  }

  // Dogleg variable
  double gamma = 0.0;

  // Trust region subproblem loop
  while ((ratio < minRatio) && (radius > minRadius)) 
  {
    if (useCounters)
      numTrustRegionInnerIterations += 1;

    Abstract::Vector* dirPtr;
    double step;

    // Trust region step
    double newtonVecNorm = computeNorm(newtonVec);
    double cauchyVecNorm = computeNorm(cauchyVec);
    
    if (newtonVecNorm <= radius) 
    {
      stepType = InexactTrustRegionBased::Newton;
      step = 1.0;
      dirPtr = &newtonVec;
    }
    else if (cauchyVecNorm >= radius) 
    {
      stepType = InexactTrustRegionBased::Cauchy;
      step = radius / cauchyVecNorm;
      dirPtr = &cauchyVec;
    }
    else 
    {			// Dogleg computation

      // aVec = newtonVec - cauchyVec
      aVec.update(1.0, newtonVec, -1.0, cauchyVec, 0.0);

      // cta = cauchyVec' * aVec
      double cta = 0.0;
      // ctc = cauchyVec' * cauchyVec
      double ctc = 0.0;
      // ata = aVec' * aVec
      double ata = 0.0;

      if (userNormPtr != 0) {
	cta = userNormPtr->dot(cauchyVec, aVec);
	ctc = userNormPtr->dot(cauchyVec, cauchyVec);
	ata = userNormPtr->dot(aVec, aVec);
      }
      else {
	cta = cauchyVec.dot(aVec);
	ctc = cauchyVec.dot(cauchyVec);
	ata = aVec.dot(aVec);
      }

      // sqrt of quadratic equation
      double tmp = (cta * cta) - ((ctc - (radius * radius)) * ata);
      if (tmp < 0) {
	utils.err() << "NOX::Solver::InexactTrustRegionBased::iterate - invalid computation" << endl;
	throw "NOX Error";
      }
      
      // final soln to quadratic equation
      gamma = (sqrt(tmp) - cta) / ata;
      if ((gamma < 0) || (gamma > 1)) {
	utils.err() << "NOX::Solver::InexactTrustRegionBased::iterate - invalid trust region step" << endl;
	throw "NOX Error";
      }
      
      // final direction computation
      aVec.update(1.0 - gamma, cauchyVec, gamma, newtonVec, 0.0);

      // solution
      stepType = InexactTrustRegionBased::Dogleg;
      dirPtr = &aVec;
      step = 1.0;
    }
    
    // Local reference to use in the remaining computation
    const Abstract::Vector& dir = *dirPtr;

    // Calculate true step length
    dx = step * (computeNorm(dir));
    
    // Compute new X
    soln.computeX(oldSoln, dir, step);

    // Compute F for new current solution.
    NOX::Abstract::Group::ReturnType rtype = soln.computeF();
    if (rtype != NOX::Abstract::Group::Ok) 
    {
      utils.err() << "NOX::Solver::InexactTrustRegionBased::iterate - "
	   << "unable to compute F" << endl;
      throw "NOX Error";
    }

    // Compute ratio of actual to predicted reduction
    // If using Homer Walker's Ared/Pred ratio computation, 
    // we use F, NOT the merit function, f.
    if (useAredPredRatio) {

      // bVec = F(x) + J d
      rtype = oldSoln.applyJacobian(*dirPtr, bVec);
      if (rtype != NOX::Abstract::Group::Ok) 
      {
	utils.out() << "NOX::Solver::TrustRegionBased::iterate - "
	     << "unable to compute F" << endl;
	throw "NOX Error";
      }
      bVec.update(1.0, oldSoln.getF(), 1.0);

      // Compute norms
      double oldNormF = computeNorm(oldSoln.getF());
      double newNormF = computeNorm(soln.getF());
      double normFLinear = computeNorm(bVec);;
      
      ratio = (oldNormF - newNormF) / (oldNormF - normFLinear);

      // Print the ratio values if requested
      if (utils.isPrintType(NOX::Utils::InnerIteration)) {
      double numerator = oldNormF - newNormF;
      double denominator = oldNormF - normFLinear;
	utils.out() << "Ratio computation: " << utils.sciformat(numerator) << "/" 
	     << utils.sciformat(denominator) << "=" << ratio << endl;
      }

      // Update the merit function (newF used when printing iteration status)
      newF = computeMeritFunction(*solnPtr);

    }
    else {  // Default ratio computation

      newF = computeMeritFunction(*solnPtr);

      if (newF >= oldF) 
      {
	ratio = -1;
      }
      else 
      {
	
	rtype = oldSoln.applyJacobian(*dirPtr, bVec);
	if (rtype != NOX::Abstract::Group::Ok) 
	{
	  utils.err() << "NOX::Solver::InexactTrustRegionBased::iterate - "
	       << "unable to compute F" << endl;
	  throw "NOX Error";
	}

	double numerator = oldF - newF;
	double denominator = 0.0;

	if (userMeritFuncPtr != 0) {
	  denominator = fabs(oldF - userMeritFuncPtr->
			     computeQuadraticModel(dir,oldSoln));
	}
	else 
	  denominator = fabs(dir.dot(oldSoln.getGradient()) + 
			     0.5 * bVec.dot(bVec));

	ratio = numerator / denominator;
	if (utils.isPrintType(NOX::Utils::InnerIteration))
	  utils.out() << "Ratio computation: " << utils.sciformat(numerator) << "/" 
	       << utils.sciformat(denominator) << "=" << ratio << endl;

	// WHY IS THIS CHECK HERE?
	if ((denominator < 1.0e-12) && ((newF / oldF) >= 0.5))
	  ratio = -1;
      }
    }  // end ratio computation

    if (utils.isPrintType(Utils::InnerIteration)) {
      utils.out() << "radius = " << utils.sciformat(radius, 1);
      utils.out() << " ratio = " << setprecision(1) << setw(3) << ratio;
      utils.out() << " f = " << utils.sciformat(sqrt(2*newF));
      utils.out() << " oldF = " << utils.sciformat(sqrt(2*oldF));
      utils.out() << " ";

      switch(stepType) {
      case InexactTrustRegionBased::Newton:
	utils.out() << "Newton";
	break;
      case InexactTrustRegionBased::Cauchy:
	utils.out() << "Cauchy";
	break;
      case InexactTrustRegionBased::Dogleg:
	utils.out() << "Dogleg";
	break;
      }

      utils.out() << endl;
    }

    // Update trust region
    if (ratio < contractTriggerRatio) 
    {
      if (stepType == InexactTrustRegionBased::Newton) {
	radius = computeNorm(newtonVec);
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
      double tmp = radius/computeNorm(newtonVec);
      sumDoglegFracNewtonLength += tmp;

      if (utils.isPrintType(Utils::Details)) {
	utils.out() << "    Fraction of Newton Step Length = " << tmp << endl;
	utils.out() << "    Fraction Between Cauchy and Newton Direction = " 
	     << gamma << endl; 
      }

    }
  }

  // Evaluate the current status
  if ((radius <= minRadius) && (ratio < minRatio)) 
  {
    if (utils.isPrintType(Utils::InnerIteration))
      utils.out() << "Using recovery step and resetting trust region." << endl;
    soln.computeX(oldSoln, newtonVec, recoveryStep);
    soln.computeF();
    radius = computeNorm(newtonVec);
    /*if (radius < minRadius)
      radius = 2 * minRadius;*/
  }

  status = test.checkStatus(*this);
 
  if (utils.isPrintType(Utils::InnerIteration)) 
    utils.out() << NOX::Utils::fill(72) << endl;

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
  Abstract::Group& soln = *solnPtr;
  StatusTest::Generic& test = *testPtr;

  // Newton direction eta (do before setting oldSoln = soln).
  eta_last = eta;
  eta = inNewtonUtils.computeForcingTerm(soln, oldSoln, nIter, *this, eta_last);
  
  // Copy current soln to the old soln.
  oldSoln = soln;

  // Compute Cauchy direction
  bool ok = cauchy.compute(cauchyVec, oldSoln, *this);
  if (!ok) {
    utils.out() << "NOX::Solver::InexactTrustRegionBased::iterate - "
	 << "unable to calculate Cauchy direction" << endl;
    status = StatusTest::Failed;
    return status;
  }

  // RPP: Can't just copy over oldf.  Scaling could change between iterations
  // so user Merit Functions could be out of sync
  oldF = computeMeritFunction(oldSoln);
  
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

  if (utils.isPrintType(NOX::Utils::InnerIteration)) 
  {
    utils.out() << NOX::Utils::fill(72) << endl;
    utils.out() << "-- Trust Region Inner Iteration --" << endl;
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
	  string directionMethod = paramsPtr->sublist("Direction").
	    getParameter("Method", "Newton");
	  NOX::Parameter::List& lsParams = paramsPtr->sublist("Direction").
	    sublist(directionMethod).sublist("Linear Solver");
	  lsParams.setParameter("Tolerance", eta);
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
	  if (userNormPtr != 0) 
	    sDotZ = userNormPtr->dot(cauchyVec,zVec);
	  else
	    sDotZ = cauchyVec.dot(zVec);

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
	    if (userNormPtr != 0) 
	      tauMin = userNormPtr->dot(rCauchyVec, residualVec) / (norm*norm);
	    else
	      tauMin = rCauchyVec.dot(residualVec) / (norm * norm);

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
    newF = computeMeritFunction(*solnPtr);

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
	 utils.out() << "NOX::Solver::TrustRegionBased::iterate - "
	      << "unable to compute F" << endl;
	 throw "NOX Error";
       } 
       bVec.update(1.0, oldSoln.getF(), 1.0);

       // Compute norms
       double oldNormF = computeNorm(oldSoln.getF());
       double newNormF = computeNorm(soln.getF());
       double normFLinear = computeNorm(bVec);;
       
       ratio = (oldNormF - newNormF) / (oldNormF - normFLinear);
       
       // Print the ratio values if requested
       if (utils.isPrintType(NOX::Utils::InnerIteration)) {
	 double numerator = oldNormF - newNormF;
	 double denominator = oldNormF - normFLinear;
	 utils.out() << "Ratio computation: " << utils.sciformat(numerator) << "/" 
	      << utils.sciformat(denominator) << "=" << ratio << endl;
       }

       // Update the merit function (newF used when printing iteration status)
       newF = computeMeritFunction(*solnPtr);

     }
     else {  // Default ratio computation

       rtype = oldSoln.applyJacobian(*dirPtr, zVec);
       if (rtype != NOX::Abstract::Group::Ok)
	 throwError("iterateInexact", "unable to applyJacobian!");
  
       double numerator = oldF - newF;
       double denominator = 0.0;

       if (userMeritFuncPtr != 0) {
	 denominator = fabs(oldF - userMeritFuncPtr->
			    computeQuadraticModel(dir,oldSoln));
       }
       else 
	 denominator = fabs(dir.dot(oldSoln.getGradient()) + 
			    0.5 * zVec.dot(zVec));

       ratio = numerator / denominator;
       if (utils.isPrintType(NOX::Utils::InnerIteration))
	 utils.out() << "Ratio computation: " << utils.sciformat(numerator) << "/" 
	      << utils.sciformat(denominator) << "=" << ratio << endl;

       // WHY IS THIS CHECK HERE?
       if ((denominator < 1.0e-12) && ((newF / oldF) >= 0.5))
	 ratio = -1;
     }
    } 

    if (utils.isPrintType(Utils::InnerIteration)) {
      utils.out() << "radius = " << utils.sciformat(radius, 1);
      utils.out() << " ratio = " << setprecision(1) << setw(3) << ratio;
      utils.out() << " f = " << utils.sciformat(sqrt(2*newF));
      utils.out() << " oldF = " << utils.sciformat(sqrt(2*oldF));
      utils.out() << " ";

      switch(stepType) {
      case InexactTrustRegionBased::Newton:
	utils.out() << "Newton";
	break;
      case InexactTrustRegionBased::Cauchy:
	utils.out() << "Cauchy";
	break;
      case InexactTrustRegionBased::Dogleg:
	utils.out() << "Dogleg";
	break;
      }
      
      utils.out() << endl;
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

      if (utils.isPrintType(Utils::Details)) {
	utils.out() << "    Fraction of Newton Step Length = " << tmp << endl;
	utils.out() << "    Fraction Between Cauchy and Newton Direction = " 
	     << tau << endl; 
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

    if (utils.isPrintType(Utils::InnerIteration))
      utils.out() << "Inner Iteration Failed!\n ";
    if (computedNewtonDir) {
      soln.computeX(oldSoln, cauchyVec, recoveryStep);
      utils.out() << "Using Newton recovery step and resetting trust region!" << endl;
    }
    else {
      soln.computeX(oldSoln, cauchyVec, recoveryStep);
      utils.out() << "Using Cauchy recovery step and resetting trust region!" << endl;
    }
    soln.computeF();
    radius = computeNorm(newtonVec);
  }

  // Evaluate the current status
  status = test.checkStatus(*this);
 
  if (utils.isPrintType(Utils::InnerIteration)) 
    utils.out() << NOX::Utils::fill(72) << endl;

  return status;    
}

//*************************************************************************
//**** checkStep
//*************************************************************************
NOX::StatusTest::StatusType 
NOX::Solver::InexactTrustRegionBased::checkStep(const NOX::Abstract::Vector& step,
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

  printUpdate();

  // Iterate until converged or failed
  while (status == StatusTest::Unconverged) {
    status = iterate();
    printUpdate();
  }

  if (writeOutputParamsToList) {
    Parameter::List& outputParams = paramsPtr->sublist("Output");
    outputParams.setParameter("Nonlinear Iterations", nIter);
    outputParams.setParameter("2-Norm of Residual", solnPtr->getNormF());
    if (useCounters) {
      Parameter::List& trOutputParams = paramsPtr->
	sublist("Trust Region").sublist("Output");
      trOutputParams.setParameter("Number of Cauchy Steps", numCauchySteps);
      trOutputParams.setParameter("Number of Newton Steps", numNewtonSteps);
      trOutputParams.setParameter("Number of Dogleg Steps", numDoglegSteps);
      trOutputParams.setParameter("Number of Trust Region Inner Iterations", 
				  numTrustRegionInnerIterations);
      if (numDoglegSteps != 0) {
	trOutputParams.setParameter("Dogleg Steps: Average Fraction of Newton Step Length", (sumDoglegFracNewtonLength/((double)numDoglegSteps)));
	trOutputParams.setParameter("Dogleg Steps: Average Fraction Between Cauchy and Newton Direction", (sumDoglegFracCauchyToNewton/((double)numDoglegSteps)));
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
  return oldSoln;
}

//*************************************************************************
//**** getNumIterations
//*************************************************************************
int NOX::Solver::InexactTrustRegionBased::getNumIterations() const
{
  return nIter;
}

//*************************************************************************
//**** getParameterList
//*************************************************************************
const Parameter::List& 
NOX::Solver::InexactTrustRegionBased::getParameterList() const
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
      (utils.isPrintType(NOX::Utils::OuterIterationStatusTest))) {
    utils.out() << NOX::Utils::fill(72) << "\n";
    utils.out() << "-- Status Test Results --\n";    
    testPtr->print(utils.out());
    utils.out() << NOX::Utils::fill(72) << "\n";
  }
  
  double fmax = solnPtr->getF().norm(Abstract::Vector::MaxNorm);
  if (utils.isPrintType(NOX::Utils::OuterIteration)) {
    utils.out() << "\n" << NOX::Utils::fill(72) << "\n";
    utils.out() << "-- Newton Trust-Region Step " << nIter << " -- \n";
    utils.out() << "f = " << utils.sciformat(sqrt(2*newF));
    utils.out() << " fmax = " << utils.sciformat(fmax);
    utils.out() << "  dx = " << utils.sciformat(dx);
    utils.out() << "  radius = " << utils.sciformat(radius);
    if (status == StatusTest::Converged)
      utils.out() << " (Converged!)";
    if (status == StatusTest::Failed)
      utils.out() << " (Failed!)";
    utils.out() << "\n" << NOX::Utils::fill(72) << "\n" << endl;
  }
  
  if ((status != StatusTest::Unconverged) && 
      (utils.isPrintType(NOX::Utils::OuterIteration))) {
    utils.out() << NOX::Utils::fill(72) << "\n";
    utils.out() << "-- Final Status Test Results --\n";    
    testPtr->print(utils.out());
    utils.out() << NOX::Utils::fill(72) << "\n";
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
  if (userNormPtr != 0)
    norm = userNormPtr->norm(v);
  else
    norm = v.norm();

  return norm;
}

//*************************************************************************
//**** computeMeritFunction
//*************************************************************************
double NOX::Solver::InexactTrustRegionBased::
computeMeritFunction(NOX::Abstract::Group& grp)
{
  double f = 0.0;

  if (userMeritFuncPtr != 0)
    f = userMeritFuncPtr->computef(grp);
  else 
    f = 0.5 * grp.getNormF() * grp.getNormF();

  return f;
}

#endif // WITH_PRERELEASE
