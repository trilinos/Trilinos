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

#include "NOX_Solver_TrustRegionBased.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Parameter_List.H"
#include "NOX_Parameter_UserNorm.H"
#include "NOX_Parameter_MeritFunction.H"
#include "NOX_Parameter_PrePostOperator.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::Solver;

/* Some compilers (in particular the SGI and ASCI Red - TFLOP) 
 * fail to find the max and min function.  Therfore we redefine them 
 * here. 
 */ 
#ifdef max
#undef max
#endif

#define max(a,b) ((a)>(b)) ? (a) : (b);

#ifdef min
#undef min
#endif

#define min(a,b) ((a)<(b)) ? (a) : (b);

TrustRegionBased::TrustRegionBased(Abstract::Group& grp, StatusTest::Generic& t, Parameter::List& p) :
  solnPtr(&grp),		// pointer to grp
  oldSolnPtr(grp.clone(DeepCopy)), // create via clone
  oldSoln(*oldSolnPtr),		// reference to just-created pointer
  newtonVecPtr(grp.getX().clone(ShapeCopy)), // create via clone 
  newtonVec(*newtonVecPtr),	// reference to just-created pointer
  cauchyVecPtr(grp.getX().clone(ShapeCopy)), // create via clone 
  cauchyVec(*cauchyVecPtr),	// reference to just-created pointer
  aVecPtr(grp.getX().clone(ShapeCopy)), // create via clone 
  aVec(*aVecPtr),		// reference to just-created pointer
  bVecPtr(grp.getX().clone(ShapeCopy)), // create via clone 
  bVec(*bVecPtr),		// reference to just-created pointer
  testPtr(&t),			// pointer to t
  paramsPtr(&p),			// copy p
  utils(paramsPtr->sublist("Printing")), // inititalize utils
  newton(utils),		// initialize direction
  cauchy(utils),       		// initialize direction
  userNormPtr(0),
  userMeritFuncPtr(0),
  useAredPredRatio(false),
  prePostOperatorPtr(0),
  havePrePostOperator(false)
{
  init();
}

// Protected
void TrustRegionBased::init()
{
  // Initialize 
  nIter = 0;
  dx = 0;
  status = StatusTest::Unconverged;
  havePrePostOperator = false;

  // Print out initialization information
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << "\n" << NOX::Utils::fill(72) << "\n";
    cout << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(cout,5);
  }

  // Set default parameter settings using getParameter() if they are not set
  paramsPtr->sublist("Direction").getParameter("Method", "Newton");
  paramsPtr->sublist("Cauchy Direction").getParameter("Method", "Steepest Descent");
  paramsPtr->sublist("Cauchy Direction").sublist("Steepest Descent").getParameter("Scaling Type", "Quadratic Model Min");

  newton.reset(paramsPtr->sublist("Direction"));
  cauchy.reset(paramsPtr->sublist("Cauchy Direction"));

  minRadius = paramsPtr->sublist("Trust Region").getParameter("Minimum Trust Region Radius", 1.0e-6);
  if (minRadius <= 0) 
    invalid("Minimum Trust Region Radius", minRadius);

  maxRadius = paramsPtr->sublist("Trust Region").getParameter("Maximum Trust Region Radius", 1.0e+10);
  if (maxRadius <= minRadius) 
    invalid("Maximum Trust Region Radius", maxRadius);

  minRatio = paramsPtr->sublist("Trust Region").getParameter("Minimum Improvement Ratio", 1.0e-4);
  if (minRatio <= 0) 
    invalid("Minimum Improvement Ratio", minRatio);

  contractTriggerRatio = paramsPtr->sublist("Trust Region").getParameter("Contraction Trigger Ratio", 0.1);
  if (contractTriggerRatio < minRatio) 
    invalid("Contraction Trigger Ratio", contractTriggerRatio);

  expandTriggerRatio = paramsPtr->sublist("Trust Region").getParameter("Expansion Trigger Ratio", 0.75);
  if (expandTriggerRatio <= contractTriggerRatio) 
    invalid("Expansion Trigger Ratio", expandTriggerRatio);

  contractFactor = paramsPtr->sublist("Trust Region").getParameter("Contraction Factor", 0.25);
  if ((contractFactor <= 0) || (contractFactor >= 1)) 
    invalid("Contraction Factor", contractFactor);

  expandFactor = paramsPtr->sublist("Trust Region").getParameter("Expansion Factor", 4.0);
  if (expandFactor <= 1) 
    invalid("Expansion Factor", expandFactor);

  recoveryStep = paramsPtr->sublist("Trust Region").getParameter("Recovery Step", 1.0);
  if (recoveryStep < 0) 
    invalid("Recovery Step", recoveryStep);


  // Check for a user defined Norm
  if (paramsPtr->sublist("Trust Region").
      isParameterArbitrary("User Defined Norm")) {
    const NOX::Parameter::UserNorm& un = 
      dynamic_cast<const NOX::Parameter::UserNorm&>(paramsPtr->
      sublist("Trust Region").getArbitraryParameter("User Defined Norm"));
    userNormPtr = const_cast<NOX::Parameter::UserNorm*>(&un);

    /*
    // RPP: Hack!!  Need to get the scaling vectors computed.
    solnPtr->computeF();
    solnPtr->computeJacobian();
    solnPtr->computeNewton(paramsPtr->sublist("Direction").sublist("Newton").sublist("Linear Solver"));
    */
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

  // Check for a user defined Pre/Post Operator
  NOX::Parameter::List& p = paramsPtr->sublist("Solver Options");
  havePrePostOperator = false;
  prePostOperatorPtr = 0;
  if (p.isParameter("User Defined Pre/Post Operator")) {
    if (p.isParameterArbitrary("User Defined Pre/Post Operator")) {
      prePostOperatorPtr = dynamic_cast<NOX::Parameter::PrePostOperator*>
	(p.getArbitraryParameter("User Defined Pre/Post Operator").clone());
      if (prePostOperatorPtr != 0)
	havePrePostOperator = true;
      else
	if (utils.isPrintProcessAndType(NOX::Utils::Warning))
	  cout << "Warning: NOX::Solver::LineSearchBased::init() - " 
	       << "\"User Defined Pre/Post Operator\" not derived from " 
	       << "NOX::Parameter::PrePostOperator class!\n" 
	       << "Ignoring this flag!"<< endl;
    }
    else {
      cout << "ERROR: NOX::Solver::LineSearchBased::init() - the parameter "
	   << "\"User Defined Pre/Post Operator\" must be derived from an"
	   << "arbitrary parameter!" << endl;
      throw "NOX Error";
    }
  }

  // Check for the using Homer Walker's Ared/Pred ratio calculation
  useAredPredRatio = 
    paramsPtr->sublist("Trust Region").getParameter("Use Ared/Pred Ratio Calculation", false);

  // Compute F of initital guess
  solnPtr->computeF();
  if (userMeritFuncPtr != 0) {
    newF = userMeritFuncPtr->computef(*solnPtr);
  }
  else 
    newF = 0.5 * solnPtr->getNormF() * solnPtr->getNormF();

  // Test the initial guess
  status = testPtr->checkStatus(*this);

  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testPtr->print(cout, 5);
    cout <<"\n" << NOX::Utils::fill(72) << "\n";
  }

}

//PRIVATE
void NOX::Solver::TrustRegionBased::invalid(const string& name, double value) const
{
  cerr << "NOX::Solver::TrustRegionBased::init - " 
       << "Invalid \"" << name << "\" (" << value << ")" 
       << endl;
  throw "NOX Error";
}

bool TrustRegionBased::reset(Abstract::Group& grp, StatusTest::Generic& t, 
			     Parameter::List& p) 
{
  solnPtr = &grp;
  testPtr = &t;
  paramsPtr = &p;			
  utils.reset(paramsPtr->sublist("Printing"));
  init();
  return true;
}

bool TrustRegionBased::reset(Abstract::Group& grp, StatusTest::Generic& t)
{
  // New initial guess and status test
  solnPtr = &grp;
  testPtr = &t;

  // Initialize 
  nIter = 0;
  dx = 0;
  status = StatusTest::Unconverged;

  // Print out initialization information
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << "\n" << NOX::Utils::fill(72) << "\n";
    cout << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(cout,5);
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

  if (utils.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testPtr->print(cout, 5);
    cout <<"\n" << NOX::Utils::fill(72) << "\n";
  }
  return true;
}

TrustRegionBased::~TrustRegionBased() 
{
  delete prePostOperatorPtr;
  delete oldSolnPtr;
}


NOX::StatusTest::StatusType TrustRegionBased::getStatus()
{
  return status;
}

NOX::StatusTest::StatusType TrustRegionBased::iterate()
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPreIterate(*this);

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
    cout << "NOX::Solver::TrustRegionBased::iterate - unable to calculate Newton direction" << endl;
    status = StatusTest::Failed;
    if (havePrePostOperator)
      prePostOperatorPtr->runPostIterate(*this);
    return status;
  }

  ok = cauchy.compute(cauchyVec, soln, *this);
  if (!ok) 
  {
    cerr << "NOX::Solver::TrustRegionBased::iterate - unable to calculate Cauchy direction" << endl;
    status = StatusTest::Failed;
    if (havePrePostOperator)
      prePostOperatorPtr->runPostIterate(*this);
    return status;
  }

  if (nIter == 0) 
  {
    if (userNormPtr != 0) {
      radius = userNormPtr->norm(newtonVec);
    }
    else 
      radius = newtonVec.norm();
    
    if (radius < minRadius)
      radius = 2 * minRadius;
  }

  // Update iteration count.
  nIter ++;

  // Copy current soln to the old soln.
  oldSoln = soln;
  // RPP: Can't just copy over oldf.  Scaling could change between iterations
  // so user Merit Functions could be out of sync
  if (userMeritFuncPtr != 0) {
    oldF = userMeritFuncPtr->computef(oldSoln);
  }
  else 
    oldF = newF;

  // Improvement ratio = (oldF - newF) / (mold - mnew)
  double ratio = -1;

  if (utils.isPrintProcessAndType(NOX::Utils::InnerIteration)) 
  {
    cout << NOX::Utils::fill(72) << endl;
    cout << "-- Trust Region Inner Iteration --" << endl;
  }

  // Trust region subproblem loop
  while ((ratio < minRatio) && (radius > minRadius)) 
  {

    Abstract::Vector* dirPtr;
    double step;

    // Trust region step
    double newtonVecNorm = 0.0;
    double cauchyVecNorm = 0.0;
    if (userNormPtr != 0) {
      newtonVecNorm = userNormPtr->norm(newtonVec);
      cauchyVecNorm = userNormPtr->norm(cauchyVec);
    }
    else { 
      newtonVecNorm = newtonVec.norm();
      cauchyVecNorm = cauchyVec.norm();
    }

    if (newtonVecNorm <= radius) 
    {
      stepType = TrustRegionBased::Newton;
      step = 1.0;
      dirPtr = &newtonVec;
    }
    else if (cauchyVecNorm >= radius) 
    {
      stepType = TrustRegionBased::Cauchy;
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
	cerr << "NOX::Solver::TrustRegionBased::iterate - invalid computation" << endl;
	throw "NOX Error";
      }
      
      // final soln to quadratic equation
      double gamma = (sqrt(tmp) - cta) / ata;
      if ((gamma < 0) || (gamma > 1)) {
	cerr << "NOX::Solver::TrustRegionBased::iterate - invalid trust region step" << endl;
	throw "NOX Error";
      }
      
      // final direction computation
      aVec.update(1.0 - gamma, cauchyVec, gamma, newtonVec, 0.0);

      // solution
      stepType = TrustRegionBased::Dogleg;
      dirPtr = &aVec;
      step = 1.0;
    }
    
    // Local reference to use in the remaining computation
    const Abstract::Vector& dir = *dirPtr;

    // Calculate true step length
    if (userNormPtr != 0) {
      dx = step * userNormPtr->norm(dir);
    }
    else
      dx = step * dir.norm();

    // Compute new X
    soln.computeX(oldSoln, dir, step);

    // Compute F for new current solution.
    NOX::Abstract::Group::ReturnType rtype = soln.computeF();
    if (rtype != NOX::Abstract::Group::Ok) 
    {
      cerr << "NOX::Solver::TrustRegionBased::iterate - unable to compute F" << endl;
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
	cout << "NOX::Solver::TrustRegionBased::iterate - "
	     << "unable to compute F" << endl;
	throw "NOX Error";
      }
      bVec.update(1.0, oldSoln.getF(), 1.0);

      // Compute norms
      double oldNormF = 0.0;
      double newNormF = 0.0;
      double normFLinear = 0.0;
      if (userNormPtr != 0) {
	oldNormF = userNormPtr->norm(oldSoln.getF());
	newNormF = userNormPtr->norm(soln.getF());
	normFLinear = userNormPtr->norm(bVec);
      }
      else {
	oldNormF = oldSoln.getNormF();
	newNormF = soln.getNormF();
	normFLinear = bVec.norm();
      }

      ratio = (oldNormF - newNormF) / (oldNormF - normFLinear);

      // Print the ratio values if requested
      if (utils.isPrintProcessAndType(NOX::Utils::InnerIteration)) {
      double numerator = oldNormF - newNormF;
      double denominator = oldNormF - normFLinear;
	cout << "Ratio computation: " << utils.sciformat(numerator) << "/" 
	     << utils.sciformat(denominator) << "=" << ratio << endl;
      }

      // Update the merit function (newF used when printing iteration status)
      if (userMeritFuncPtr != 0) {
	newF = userMeritFuncPtr->computef(*solnPtr);
      }
      else 
	newF = 0.5 * solnPtr->getNormF() * solnPtr->getNormF();

    }
    else {  // Default ratio computation

      if (userMeritFuncPtr != 0) {
	newF = userMeritFuncPtr->computef(*solnPtr);
      }
      else 
	newF = 0.5 * solnPtr->getNormF() * solnPtr->getNormF();
      
      if (newF >= oldF) 
      {
	ratio = -1;
      }
      else 
      {
	  
	rtype = oldSoln.applyJacobian(*dirPtr, bVec);
	if (rtype != NOX::Abstract::Group::Ok) 
	{
	  cerr << "NOX::Solver::TrustRegionBased::iterate - unable to compute F" << endl;
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
	if (utils.isPrintProcessAndType(NOX::Utils::InnerIteration))
	  cout << "Ratio computation: " << utils.sciformat(numerator) << "/" 
	       << utils.sciformat(denominator) << "=" << ratio << endl;
	
	// WHY IS THIS CHECK HERE?
	if ((denominator < 1.0e-12) && ((newF / oldF) >= 0.5))
	  ratio = -1;
      }
    }

    if (utils.isPrintProcessAndType(Utils::InnerIteration)) {
      cout << "radius = " << utils.sciformat(radius, 1);
      cout << " ratio = " << setprecision(1) << setw(3) << ratio;
      cout << " f = " << utils.sciformat(sqrt(2*newF));
      cout << " oldF = " << utils.sciformat(sqrt(2*oldF));
      cout << " ";

      switch(stepType) {
      case TrustRegionBased::Newton:
	cout << "Newton";
	break;
      case TrustRegionBased::Cauchy:
	cout << "Cauchy";
	break;
      case TrustRegionBased::Dogleg:
	cout << "Dogleg";
	break;
      }

      cout << endl;
    }

    // Update trust region
    if (ratio < contractTriggerRatio) 
    {
      if (stepType == TrustRegionBased::Newton) {
	if (userNormPtr != 0) {
	  radius = userNormPtr->norm(newtonVec);
	}
	else 
	  radius = newtonVec.norm();
      }
      radius = max(contractFactor * radius, minRadius);
    }
    else if ((ratio > expandTriggerRatio) && (dx == radius)) 
    {
      radius = min(expandFactor * radius, maxRadius);
    }

  }


  // Evaluate the current status
  if ((radius <= minRadius) && (ratio < minRatio)) 
  {
    if (utils.isPrintProcessAndType(Utils::InnerIteration))
      cout << "Using recovery step and resetting trust region." << endl;
    soln.computeX(oldSoln, newtonVec, recoveryStep);
    soln.computeF();
    if (userNormPtr != 0) {
      radius = userNormPtr->norm(newtonVec);
    }
    else
      radius = newtonVec.norm();
    /*if (radius < minRadius)
      radius = 2 * minRadius;*/
  }

  status = test.checkStatus(*this);
 
  if (utils.isPrintProcessAndType(Utils::InnerIteration)) 
    cout << NOX::Utils::fill(72) << endl;

  if (havePrePostOperator)
    prePostOperatorPtr->runPostIterate(*this);

  // Return status.
  return status;
}

NOX::StatusTest::StatusType TrustRegionBased::solve()
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPreSolve(*this);

  printUpdate();

  // Iterate until converged or failed
  while (status == StatusTest::Unconverged) {
    status = iterate();
    printUpdate();
  }

  Parameter::List& outputParams = paramsPtr->sublist("Output");
  outputParams.setParameter("Nonlinear Iterations", nIter);
  outputParams.setParameter("2-Norm of Residual", solnPtr->getNormF());

  if (havePrePostOperator)
    prePostOperatorPtr->runPostSolve(*this);

  return status;
}

const Abstract::Group& TrustRegionBased::getSolutionGroup() const
{
  return *solnPtr;
}

const Abstract::Group& TrustRegionBased::getPreviousSolutionGroup() const
{
  return oldSoln;
}

int TrustRegionBased::getNumIterations() const
{
  return nIter;
}

const Parameter::List& TrustRegionBased::getParameterList() const
{
  return *paramsPtr;
}

// protected
void TrustRegionBased::printUpdate() 
{
  // Print the status test parameters at each iteration if requested  
  if ((status == StatusTest::Unconverged) && 
      (utils.isPrintProcessAndType(NOX::Utils::OuterIterationStatusTest))) {
    cout << NOX::Utils::fill(72) << "\n";
    cout << "-- Status Test Results --\n";    
    testPtr->print(cout);
    cout << NOX::Utils::fill(72) << "\n";
  }
  
  double fmax = solnPtr->getF().norm(Abstract::Vector::MaxNorm);
  if (utils.isPrintProcessAndType(NOX::Utils::OuterIteration)) {
    cout << "\n" << NOX::Utils::fill(72) << "\n";
    cout << "-- Newton Trust-Region Step " << nIter << " -- \n";
    cout << "f = " << utils.sciformat(sqrt(2*newF));
    cout << " fmax = " << utils.sciformat(fmax);
    cout << "  dx = " << utils.sciformat(dx);
    cout << "  radius = " << utils.sciformat(radius);
    if (status == StatusTest::Converged)
      cout << " (Converged!)";
    if (status == StatusTest::Failed)
      cout << " (Failed!)";
    cout << "\n" << NOX::Utils::fill(72) << "\n" << endl;
  }
  
  if ((status != StatusTest::Unconverged) && 
      (utils.isPrintProcessAndType(NOX::Utils::OuterIteration))) {
    cout << NOX::Utils::fill(72) << "\n";
    cout << "-- Final Status Test Results --\n";    
    testPtr->print(cout);
    cout << NOX::Utils::fill(72) << "\n";
  }
}

