// $Id$ 
// $Source$ 

#ifdef WITH_PRERELEASE
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
//
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "NOX_Solver_TensorBased.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Parameter_List.H"
#include "NOX_Parameter_PrePostOperator.H"
#include "NOX_Utils.H"

#include "NOX_LineSearch_Utils_Printing.H"  // class data member
#include "NOX_LineSearch_Utils_Counters.H"  // class data member
#include "NOX_LineSearch_Utils_Slope.H"     // class data member

//   include "stdio.h"  // for printf()


#define CHECK_RESIDUALS
#define DEBUG_LEVEL 0

NOX::Solver::TensorBased::TensorBased(NOX::Abstract::Group& xGrp,
				      NOX::StatusTest::Generic& t,
				      NOX::Parameter::List& p) :
  solnPtr(&xGrp),		// pointer to xGrp
  oldSolnPtr(xGrp.clone(DeepCopy)), // create via clone
  oldSoln(*oldSolnPtr),		// reference to just-created pointer
  newtonVecPtr(xGrp.getX().clone(ShapeCopy)), // create via clone 
  newtonVec(*newtonVecPtr),	// reference to just-created pointer
  tensorVecPtr(xGrp.getX().clone(ShapeCopy)), // create via clone 
  tensorVec(*tensorVecPtr),	// reference to just-created pointer
  acVecPtr(xGrp.getX().clone(ShapeCopy)), // create via clone 
  acVec(*acVecPtr),		// reference to just-created pointer
  scVecPtr(xGrp.getX().clone(ShapeCopy)), // create via clone 
  scVec(*scVecPtr),		// reference to just-created pointer
  tmpVecPtr(xGrp.getX().clone(ShapeCopy)), // create via clone 
  tmpVec(*tmpVecPtr),		// reference to just-created pointer
  residualVecPtr(xGrp.getX().clone(ShapeCopy)), // create via clone 
  testPtr(&t),			// pointer to t
  paramsPtr(&p),		// copy p
  lsParams(paramsPtr->sublist("Line Search")),  // reference to list
  dirParams(paramsPtr->sublist("Direction")),   // reference to list
  localParams(paramsPtr->sublist("Direction").sublist("Tensor").
	      sublist("Linear Solver")),   // reference to list
  utils(paramsPtr->sublist("Printing")),               // initialize utils
  print(utils),
  prePostOperatorPtr(0),
  havePrePostOperator(false)
{
  //init();
  reset(xGrp, t, p);
}

// Protected
void NOX::Solver::TensorBased::init()
{
  // Initialize
  step = 0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;

  // Reset counters 
  counter.reset();
  numJvMults = 0;
  numJ2vMults = 0;
  
  // Check for a user defined Pre/Post Operator
  NOX::Parameter::List& p = paramsPtr->sublist("Solver Options");
  havePrePostOperator = false;
  prePostOperatorPtr = 0;
  if (p.isParameter("User Defined Pre/Post Operator"))
  {
    if (p.isParameterArbitrary("User Defined Pre/Post Operator"))
    {
      prePostOperatorPtr = dynamic_cast<NOX::Parameter::PrePostOperator*>
	(p.getArbitraryParameter("User Defined Pre/Post Operator").clone());
      if (prePostOperatorPtr != 0)
	havePrePostOperator = true;
      else
	if (utils.isPrintProcessAndType(NOX::Utils::Warning))
	  cout << "Warning: NOX::Solver::TensorBased::init() - " 
	       << "\"User Defined Pre/Post Operator\" not derived from " 
	       << "NOX::Parameter::PrePostOperator class!\n" 
	       << "Ignoring this flag!"<< endl;
    }
    else
    {
      cout << "ERROR: NOX::Solver::TensorBased::init() - the parameter "
	   << "\"User Defined Pre/Post Operator\" must be derived from an"
	   << "arbitrary parameter!" << endl;
      throw "NOX Error";
    }
  }
  
  // Print out parameters
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters))
  {
    cout << "\n" << NOX::Utils::fill(72) << "\n";
    cout << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(cout,5);
    cout << "\n" << NOX::Utils::fill(72) << "\n";
  }

  // Compute F of initial guess
  NOX::Abstract::Group::ReturnType rtype = solnPtr->computeF();
  if (rtype != NOX::Abstract::Group::Ok)
  {
    cout << "NOX::Solver::TensorBased::init - Unable to compute F" << endl;
    throw "NOX Error";
  }

  // Test the initial guess
  status = testPtr->checkStatus(*this);
  if ((status == NOX::StatusTest::Converged) &&
      (utils.isPrintProcessAndType(NOX::Utils::Warning)))
  {
    cout << "Warning: NOX::Solver::TensorBased::init() - The solution passed "
	 << "into the solver (either through constructor or reset method) "
	 << "is already converged!  The solver will not "
	 << "attempt to solve this system since status is flagged as "
	 << "converged." << endl;
  }

  // Print out status tests
  if (utils.isPrintProcessAndType(NOX::Utils::Parameters))
  {
    cout << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testPtr->print(cout, 5);
    cout <<"\n" << NOX::Utils::fill(72) << "\n";
  }

}


bool NOX::Solver::TensorBased::reset(NOX::Abstract::Group& xGrp,
				     NOX::StatusTest::Generic& t,
				     NOX::Parameter::List& p)
{
  solnPtr = &xGrp;
  testPtr = &t;
  paramsPtr = &p;

  utils.reset(paramsPtr->sublist("Printing"));

  // *** Reset direction parameters ***
  dirParams = paramsPtr->sublist("Direction");

  // bwb: Do we just want to use Method instead of "Compute Step"?
  string choice = dirParams.getParameter("Compute Step", "Tensor");
  if (choice == "Tensor")
  {
    requestedBaseStep = TensorStep;
  }
  else if (choice == "Newton")
  {
    requestedBaseStep = NewtonStep;
    if (utils.isPrintProcessAndType(NOX::Utils::Parameters))
      cout << "\n\n    **** Newton step is requested ***** \n\n";
  }
  else
  {
    if (utils.isPrintProcessAndType(NOX::Utils::Warning))
      cout << "Warning: NOX::Direction::Tensor::reset() - The choice of "
	   << "\"Compute Step\" \nparameter \"" << choice
	   << "\" is invalid.  Using \"Tensor\" instead." << endl;
    requestedBaseStep = TensorStep;
    dirParams.setParameter("Compute Step", "Tensor");
  }

  bool useModifiedMethod =
    dirParams.getParameter("Use Modified Bouaricha", true);
  if (useModifiedMethod  &&  utils.isPrintProcessAndType(NOX::Utils::Parameters))
    cout << "Using ALPHA scaling" << endl;
  
  //NOX::Parameter::List& teParams = dirParams.sublist("Tensor");
  //doRescue = teParams.getParameter("Rescue Bad Newton Solve", true);

  
  // *** Reset parameters for Line Search ***
  lsParams = paramsPtr->sublist("Line Search");

  NOX::Parameter::List& lsparams =
    lsParams.sublist(lsParams.getParameter("Method", "Tensor"));
  
  // Initialize linesearch parameters for this object
  minStep = lsparams.getParameter("Minimum Step", 1.0e-12);
  defaultStep = lsparams.getParameter("Default Step", 1.0);
  recoveryStep = lsparams.getParameter("Recovery Step", 0.0); // force exit on failure
  maxIters = lsparams.getParameter("Max Iters", 40);
  alpha = lsparams.getParameter("Alpha Factor", 1.0e-4);

  // Determine the specific type of tensor linesearch to perform
  choice = lsparams.getParameter("Submethod", "Curvilinear");

  if (choice == "Curvilinear")
    lsType = Curvilinear;
  else if (choice == "Dual")
    lsType = Dual;
  else if (choice == "Standard")
    lsType = Standard;
  else if (choice == "Newton")
    lsType = Newton;
  else
  {
    if (print.isPrintProcessAndType(NOX::Utils::Warning)) {
      cout << "Warning: NOX::Direction::Tensor::reset() - the choice of "
	   << "\"Line Search\" \nparameter " << choice
	   << " is invalid.  Using \"Curvilinear\" instead." << endl;
    }
    lsparams.setParameter("Submethod", "Curvilinear");
    lsType = Curvilinear;
  }

  choice = lsparams.getParameter("Lambda Selection", "Halving");
  if (choice == "Halving")
    lambdaSelection = Halving;
  else if (choice == "Quadratic") 
    lambdaSelection = Quadratic;
  else
  {
    if (utils.isPrintProcessAndType(NOX::Utils::Warning))
      cout << "Warning: NOX::Solver::TensorBased::reset() - the choice of "
	   << "\"Lambda Selection\" \nparameter " << choice
	   << " is invalid.  Using halving instead." << endl;
    lambdaSelection = Halving;
  }

  choice = lsparams.getParameter("Sufficient Decrease Condition",
				 "Armijo-Goldstein");
  if (choice == "Ared/Pred") 
    convCriteria = AredPred;
  else if (choice == "None")
    convCriteria = None;
  else 
    convCriteria = ArmijoGoldstein;     // bwb - the others aren't implemented

  init();
  return true;
}

bool NOX::Solver::TensorBased::reset(NOX::Abstract::Group& xGrp,
				     NOX::StatusTest::Generic& t)
{
  solnPtr = &xGrp;
  testPtr = &t;
  init();
  return true;
}

NOX::Solver::TensorBased::~TensorBased() 
{
  if (utils.isPrintProcessAndType(NOX::Utils::Details))
  {
    cout << "multsJv = " << numJvMults << "   (linesearch)" << endl;
    cout << "mults2Jv = " << numJ2vMults << endl;
  }
  delete prePostOperatorPtr;
  delete oldSolnPtr;
  delete newtonVecPtr;
  delete tensorVecPtr;
  delete acVecPtr;
  delete scVecPtr;
  delete tmpVecPtr;
  delete residualVecPtr;
}


NOX::StatusTest::StatusType  NOX::Solver::TensorBased::getStatus()
{
  return status;
}

NOX::StatusTest::StatusType  NOX::Solver::TensorBased::iterate()
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPreIterate(*this);

  // First check status
  if (status != NOX::StatusTest::Unconverged)
  { 
    if (havePrePostOperator)
      prePostOperatorPtr->runPostIterate(*this);
    return status;
  }

  // Copy pointers into temporary references
  NOX::Abstract::Group& soln = *solnPtr;
  NOX::StatusTest::Generic& test = *testPtr;

  // Compute the direction for the update vector at the current solution.
  bool ok = computeTensorDirection(soln, *this);
  if (!ok)
  {
    cout << "NOX::Solver::TensorBased::iterate - "
	 << "unable to calculate direction" << endl;
    status = NOX::StatusTest::Failed;
    if (havePrePostOperator)
      prePostOperatorPtr->runPostIterate(*this);
    return status;
  }

  // Update iteration count.
  nIter ++;

  // Copy current soln to the old soln.
  oldSoln = soln;

  // Do line search and compute new soln.
  //ok = lineSearch.compute(soln, step, dir, *this);
  ok = implementGlobalStrategy(soln, step, *this);
  if (!ok)
  {
    if (step == 0)
    {
      cout << "NOX::Solver::TensorBased::iterate - line search failed"
	   << endl;
      status = NOX::StatusTest::Failed;
      if (havePrePostOperator)
	prePostOperatorPtr->runPostIterate(*this);
      return status;
    }
    else if (utils.isPrintProcessAndType(NOX::Utils::Warning))
      cout << "NOX::Solver::TensorBased::iterate - "
	   << "using recovery step for line search" << endl;
  }

  // Compute F for new current solution.
  NOX::Abstract::Group::ReturnType rtype = soln.computeF();
  if (rtype != NOX::Abstract::Group::Ok)
  {
    cout << "NOX::Solver::LineSearchBased::iterate - "
	 << "unable to compute F" << endl;
    status = NOX::StatusTest::Failed;
    if (havePrePostOperator)
      prePostOperatorPtr->runPostIterate(*this);
    return status;
  }

  status = test.checkStatus(*this);
 
  if (havePrePostOperator)
    prePostOperatorPtr->runPostIterate(*this);

  return status;
}


NOX::StatusTest::StatusType  NOX::Solver::TensorBased::solve()
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPreSolve(*this);

  printUpdate();

  // Iterate until converged or failed
  while (status == NOX::StatusTest::Unconverged)
  {
    status = iterate();
    printUpdate();
  }

  NOX::Parameter::List& outputParams = paramsPtr->sublist("Output");
  outputParams.setParameter("Nonlinear Iterations", nIter);
  outputParams.setParameter("2-Norm of Residual", solnPtr->getNormF());

  if (havePrePostOperator)
    prePostOperatorPtr->runPostSolve(*this);

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
  return oldSoln;
}

int NOX::Solver::TensorBased::getNumIterations() const
{
  return nIter;
}

const NOX::Parameter::List&
NOX::Solver::TensorBased::getParameterList() const
{
  return *paramsPtr;
}

// protected
void NOX::Solver::TensorBased::printUpdate() 
{
  double normSoln = 0;
  double normStep = 0;

  // Print the status test parameters at each iteration if requested  
  if ((status == NOX::StatusTest::Unconverged) &&
      (utils.isPrintProcessAndType(NOX::Utils::OuterIterationStatusTest)))
  {
    cout << NOX::Utils::fill(72) << "\n";
    cout << "-- Status Test Results --\n";    
    testPtr->print(cout);
    cout << NOX::Utils::fill(72) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (utils.isPrintType(NOX::Utils::InnerIteration))
  {
    normSoln = solnPtr->getNormF();
    normStep = (nIter > 0) ? tensorVec.norm() : 0;
  }

  // ...But only the print process actually prints the result.
  if (utils.isPrintProcessAndType(NOX::Utils::OuterIteration))
  {
    cout << "\n" << NOX::Utils::fill(72) << "\n";
    cout << "-- Nonlinear Solver Step " << nIter << " -- \n";
    cout << "f = " << utils.sciformat(normSoln);
    cout << "  step = " << utils.sciformat(step);
    cout << "  dx = " << utils.sciformat(normStep);
    if (status == NOX::StatusTest::Converged)
      cout << " (Converged!)";
    if (status == NOX::StatusTest::Failed)
      cout << " (Failed!)";
    cout << "\n" << Utils::fill(72) << "\n" << endl;
  }

  // Print the final parameter values of the status test
  if ((status != NOX::StatusTest::Unconverged) && 
      (utils.isPrintProcessAndType(NOX::Utils::OuterIteration)))
  {
    cout << NOX::Utils::fill(72) << "\n";
    cout << "-- Final Status Test Results --\n";    
    testPtr->print(cout);
    cout << NOX::Utils::fill(72) << "\n";
  }
}


bool
NOX::Solver::TensorBased::computeTensorDirection(NOX::Abstract::Group& soln,
					 const NOX::Solver::Generic& solver)
{
  NOX::Abstract::Group::ReturnType status;
  
  // Compute F at current solution.
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok) 
    throwError("computeTensorDirection", "Unable to compute F");

  // Compute Jacobian at current solution.
  status = soln.computeJacobian();
  if (status != NOX::Abstract::Group::Ok) 
    throwError("computeTensorDirection", "Unable to compute Jacobian");
  
  // Begin processing for the tensor step, if necessary.
  double sDotS = 0.0;
  int tempVal1 = 0;
  if ((nIter > 0)  &&  (requestedBaseStep == TensorStep))
  {
    // Save old Newton step as initial guess to second system  (not necessary)
    tmpVec = newtonVec;
    tmpVecPtr->scale(-1.0);   // could probably rewrite to avoid this...

    // Compute the tensor term s = x_{k-1} - x_k
    scVec = soln.getX();
    scVec.update(1.0, solver.getPreviousSolutionGroup().getX(), -1.0);
    double normS = scVec.norm();
    sDotS = normS * normS;

    // Form the tensor term a = (F_{k-1} - F_k - J*s) / (s^T s)^2
    soln.applyJacobian(*scVecPtr, *acVecPtr);
    numJvMults++;
    acVec.update(1.0, solver.getPreviousSolutionGroup().getF(), -1.0);
    acVec.update(-1.0, soln.getF(), 1.0);
    if (sDotS != 0)
      acVec.scale(1 / (sDotS * sDotS));
    
#undef OLD_WAY
#ifdef OLD_WAY
    // Compute inv(J)*a
    //tmpVec.init(0.0);
    //printf("\n\n\nNorm of tmpVec = %8e\n\n\n", tmpVec.norm()); 
    status = soln.applyJacobianInverse(localParams, acVec, tmpVec);
    if (status != NOX::Abstract::Group::Ok)
      throwError("computeTensorDirection", "Unable to apply Jacobian inverse");
    //printf("\n\n\nNorm of tmpVec = %8e\n\n\n", tmpVec.norm()); 
#endif // OLD_WAY

    // Compute residual of linear system using initial guess...
    soln.applyJacobian(tmpVec, *residualVecPtr);
    numJvMults++;
    residualVecPtr->update(1.0, solver.getPreviousSolutionGroup().getF(),-1.0);
    double residualNorm = residualVecPtr->norm();

#if DEBUG_LEVEL > 0
    cout << " Norm of initial guess: " << utils.sciformat(tmpVec.norm(), 6)
	 << endl;
    cout << " initg norm of model residual =   "
	 << utils.sciformat(residualNorm, 6) << " (abs)     "
	 << utils.sciformat(residualNorm /
			    solver.getPreviousSolutionGroup().getNormF(), 6)
	 << " (rel)" << endl;
#endif

    // Save some parameters and use them later...
    double tol = localParams.getParameter("Tolerance", 1e-4);
    double relativeResidual = residualNorm /
      solver.getPreviousSolutionGroup().getNormF();

    // Decide whether to use initial guess...
    bool isInitialGuessGood = false;
#ifdef USE_INITIAL_GUESS_LOGIC    
    if (relativeResidual < 1.0)
    {
      cout << "Initial guess is good..." << endl;
      isInitialGuessGood = true;
      tensorVec = tmpVec;
      double newTol = tol / relativeResidual;
      if (newTol > 0.99)
	newTol = 0.99;  // force at least one iteration
      localParams.setParameter("Tolerance",  newTol);
      cout << "Setting tolerance to " << utils.sciformat(newTol,6) << endl;
    }
    else
#endif // USE_INITIAL_GUESS_LOGIC    
    {
      //printf("Initial guess is BAD... do not use!\n");
      isInitialGuessGood = false;
      *residualVecPtr = solver.getPreviousSolutionGroup().getF();
    }
    
    // Compute the term inv(J)*Fp....
    //tmpVec.init(0.0);
    status = soln.applyJacobianInverse(localParams, *residualVecPtr, tmpVec);
    if (status != NOX::Abstract::Group::Ok)
      throwError("computeTensorDirection", "Unable to apply Jacobian inverse");
    if (isInitialGuessGood) 
      tmpVec.update(1.0, tensorVec, 1.0);
    localParams.setParameter("Tolerance",  tol);
    tempVal1 = localParams.sublist("Output").
    getParameter("Number of Linear Iterations", 0);
    
#if DEBUG_LEVEL > 0
    // Compute residual of linear system with initial guess...
    soln.applyJacobian(tmpVec, *residualVecPtr);
    numJvMults++;
    residualVecPtr->update(-1.0, solver.getPreviousSolutionGroup().getF(),1.0);
    double residualNorm2 = residualVecPtr->norm();
    cout << " jifp norm of model residual =   "
	 << utils.sciformat(residualNorm2, 6) << " (abs)     "
	 << utils.sciformat(residualNorm2 /
			    solver.getPreviousSolutionGroup().getNormF(), 6)
	 << " (rel)" << endl;

    //printf(" jifp norm of model residual = %14.6e (abs)   %14.6e (rel)\n",
    //   residualNorm2,
    //   residualNorm2 / solver.getPreviousSolutionGroup().getNormF());
#endif
  }

  // Compute the Newton direction
  status = soln.computeNewton(localParams);
  if (status != NOX::Abstract::Group::Ok)
    throwError("computeTensorDirection", "Unable to compute Newton step");
  newtonVec = soln.getNewton();
  int tempVal2 = localParams.sublist("Output").
    getParameter("Number of Linear Iterations", 0);

  numJ2vMults += (tempVal1 > tempVal2) ? tempVal1 : tempVal2;
  
#ifdef CHECK_RESIDUALS
  printDirectionInfo("newtonVec", newtonVec, soln, false);
#endif // CHECK_RESIDUALS

  // Continue processing the tensor step, if necessary
  if ((nIter > 0)  &&  (requestedBaseStep == TensorStep))
  {

    // Form the term inv(J)*a...  (note that a is not multiplied by 2)
    tmpVecPtr->update(1.0, newtonVec, -1.0, scVec, 1.0);
    if (sDotS != 0)
      tmpVecPtr->scale(1/(sDotS * sDotS));

    // Calculate value of beta
    sctjf = -scVec.dot(newtonVec);
    sctja = scVec.dot(tmpVec);
    double qval = 0;
    double lambdaBar = 1;
    beta = calculateBeta(sctja, 1.0, sctjf, qval, lambdaBar);

    cout << " sctjf = " << utils.sciformat(sctjf, 6)
	 << "  sctja = " << utils.sciformat(sctja, 6) << endl;
    cout << " norm(s) = " << utils.sciformat(scVec.norm(), 6)
	 << "  norm(a) = " << utils.sciformat(acVecPtr->norm(), 6) << endl;
 
    //printf(" sctjf = %e  sctja = %e\n", sctjf, sctja);
    //printf(" norm(s) = %e  norm(a) = %e\n", scVec.norm(), acVecPtr->norm());
    
    if (useModifiedMethod)
    {
      double alpha2 = lambdaBar;
      cout << " Beta = " << utils.sciformat(beta, 6)
	   << "  Alpha2 = " << utils.sciformat(alpha2, 6) << endl;
      //printf(" Beta = %e   Alpha2 = %e\n", beta, alpha2);
      if (alpha2 != 1.0)
      {
	cout << "   *** Scaling tensor term a ***" << endl;
	//printf("   *** Scaling tensor term a ***\n");
	acVecPtr->scale(alpha2);
	tmpVec.scale(alpha2);
	sctja *= alpha2;
	beta /= alpha2;
	lambdaBar = 1.0;
	qval = 0;
      }
    }
    
    // Form the tensor step
    tensorVec.update(1.0, newtonVec, -beta*beta, tmpVec, 0.0);
    
#ifdef CHECK_RESIDUALS
    printDirectionInfo("tensorVec", tensorVec, soln, true);
#endif // CHECK_RESIDUALS
#if DEBUG_LEVEL > 0
    cout << "Beta = " << utils.sciformat(beta, 6)
	 << "  std = " << utils.sciformat(tensorVec.dot(scVec), 6)
         << "  qval = " << utils.sciformat(qval, 2)
         << "  lambdaBar = " << lambdaBar << endl;
    //printf("Beta = %e  std = %e  qval = %.2f   lambdaBar = %f\n",
    //   beta, tensorVec.dot(scVec), qval, lambdaBar);
#endif
  }
  else
    tensorVec = newtonVec;
  
  return true;
}


double NOX::Solver::TensorBased::calculateBeta(double qa,
					       double qb,
					       double qc,
					       double& qval,
					       double& lambdaBar,
					       double lambda) const
{
  double beta = 0.0;
  double discriminant = qb*qb - 4*qa*qc*lambda;

  if (discriminant < 0.0)
  {
    // no real root
    beta = -qb / qa / 2.0;
    qval = (qa * beta * beta) + (qb * beta) + (lambda * qc);
    lambdaBar = qb*qb / (4*qa*qc);
#if DEBUG_LEVEL > 0
    cout << "  ####  LambdaBar = " << lambdaBar << "  ####\n";
#endif
  }
  else
  {
    qval = 0;
    lambdaBar = 1.0;
    if ( (fabs(qa / qb) < 1e-8)  &&  (fabs(lambda * qc / qb) < 1) )
    {
#if DEBUG_LEVEL > 0
      cout << "qa is relatively small\n";
#endif 
      beta = -lambda * qc / qb;
    }
    else
    {
      double tmp1 = (-qb + sqrt(discriminant)) / (2*qa);
      double tmp2 = (-qb - sqrt(discriminant)) / (2*qa);
      beta = (fabs(tmp1) < fabs(tmp2)) ? tmp1 : tmp2; // bwb - temporary test
      //beta = (fabs(dir0xsc + normS*tmp1) < fabs(dir0xsc + normS*tmp2)) ? tmp1 : tmp2;
#if DEBUG_LEVEL > 1
      cout << "  tmp1 = " << utils.sciformat(tmp1, 6)
	   << "  tmp2 = " << utils.sciformat(tmp2, 6)
	   << endl;
      //printf("  tmp1 = %e  tmp2 = %e  dir0xsc = %e  normS = %e\n",
      //     tmp1, tmp2, dir0xsc, normS);
#endif
    }
  }
#if DEBUG_LEVEL > 1
  cout << "  qa,qb,qc = " << utils.sciformat(qa, 6)
       << utils.sciformat(qb, 6)
       << utils.sciformat(qc, 6)
       << "   beta = " << utils.sciformat(beta, 6)
       << endl;
  //printf("  qa,qb,qc = %e  %e  %e   beta = %e\n", qa, qb, qc, beta);
#endif

  return beta;
}


bool
NOX::Solver::TensorBased::computeCurvilinearStep(NOX::Abstract::Vector& dir,
					 const NOX::Abstract::Group& soln,
					 const NOX::Solver::Generic& s,
					 double& lambda)
{
  //dir = soln.getNewton();
  //dir.scale(step);

  double qval = 0;
  double lambdaBar = 1;
  double beta1 = calculateBeta(sctja, 1, sctjf, qval, lambdaBar, lambda);
  double betaFactor = ( (beta == 0.0) ? 0.0 : beta1*beta1 / (beta*beta));
  
  dir.update(lambda - betaFactor, newtonVec, betaFactor, tensorVec, 0.0);

#if DEBUG_LEVEL > 0
  cout << "  beta = " << utils.sciformat(beta, 6)
       << "  std = " << utils.sciformat(dir.dot(scVec), 6)
       << "  qval = " << qval
       << "  lambda = " << lambda
       << endl;
  cout << "betaFactor = " << utils.sciformat(betaFactor,6)
       << "  beta1 = " << utils.sciformat(beta1, 6)
       << endl;
  //printf("Beta = %e  std = %e  qval = %.2f   lambdaBar = %f\n",
  // beta, dir.dot(scVec), qval, lambdaBar);
  //printf("betaFactor = %e  beta1 = %e\n", betaFactor, beta1);
#endif
  
  return true;
}


bool
NOX::Solver::TensorBased::implementGlobalStrategy(NOX::Abstract::Group& newGrp,
					  double& step, 
					  const NOX::Solver::Generic& s)
{
  bool ok;
  counter.incrementNumLineSearches();
  isNewtonDirection = false;
  NOX::Abstract::Vector& searchDirection = tensorVec;

  if ((counter.getNumLineSearches() == 1)  ||  (lsType == Newton))
  {
    isNewtonDirection = true;
    searchDirection = newtonVec;
  }

  // Do line search and compute new soln.
  if ((lsType != Dual) || (isNewtonDirection))
    ok = performLinesearch(newGrp, step, searchDirection, s);
  else if (lsType == Dual)
  {
    double fTensor;
    double tensorStep;
    bool isTensorDescent = false;

    const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
    double fprime = slopeObj.computeSlope(searchDirection, oldGrp);

    // Backtrack along tensor direction if it is descent direction.
    if (fprime < 0)
    {
      ok = performLinesearch(newGrp, step, searchDirection, s);
      fTensor = 0.5 * newGrp.getNormF() * newGrp.getNormF();
      tensorStep = step;
      isTensorDescent = true;
    }

    // Backtrack along the Newton direction.
    ok = performLinesearch(newGrp, step, newtonVec, s);
    double fNew = 0.5 * newGrp.getNormF() * newGrp.getNormF();

    // If backtracking on the tensor step produced a better step, then use it.
    if (isTensorDescent  &&  (fTensor < fNew))
    {
      newGrp.computeX(oldGrp, tensorVec, tensorStep);
      newGrp.computeF();    
    }
  }
  
  return ok;
}


bool
NOX::Solver::TensorBased::performLinesearch(NOX::Abstract::Group& newsoln,
					    double& step,
					    const NOX::Abstract::Vector& dir,
					    const NOX::Solver::Generic& s)
{
  if (print.isPrintProcessAndType(NOX::Utils::InnerIteration))
  {
    cout << "\n" << NOX::Utils::fill(72) << "\n";
    cout << "-- Tensor Line Search ("
	 << paramsPtr->sublist("Line Search").sublist("Tensor").
      getParameter("Submethod","Curvilinear")
	 << ") -- \n";
  }

  // Local variables
  bool isFailed = false;
  bool isAcceptable = false;
  bool isFirstPass = true;
  string message = "(STEP ACCEPTED!)";

  // Set counters
  int lsIterations = 1;

  // Get Old f
  const Abstract::Group& oldSoln = s.getPreviousSolutionGroup();
  double fOld = 0.5 * oldSoln.getNormF() * oldSoln.getNormF();  

  // Compute first trial point and its function value
  step = defaultStep;
  newsoln.computeX(oldSoln, dir, step);
  newsoln.computeF();    
  double fNew = 0.5 * newsoln.getNormF() * newsoln.getNormF();  

  // Compute directional derivative
  double fprime;
  if ((lsType == Curvilinear)  &&  !(isNewtonDirection)) 
    fprime = slopeObj.computeSlope(newtonVec, oldSoln);
  else 
    fprime = slopeObj.computeSlope(dir, oldSoln);
  numJvMults++;  // computeSlope() has J*v inside of it
  
  // Compute the convergence criteria for the line search 
  double threshold = fOld + alpha*step*fprime;
  isAcceptable = (fNew < threshold);

  // Update counter and allocate memory for dir2 if a linesearch is needed
  if (!isAcceptable)
  { 
    counter.incrementNumNonTrivialLineSearches();
    tmpVec = dir;
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

    print.printStep(lsIterations, step, fOld, fNew);

    // Is the full tensor step a descent direction?  If not, switch to Newton
    if ((!isNewtonDirection) && (isFirstPass && fprime >= 0))
    {
      tmpVec = newtonVec;
      fprime = slopeObj.computeSlope(tmpVec, oldSoln);
      numJvMults++;
      
      cout << "  Switching to Newton.  New fprime = "
	   << utils.sciformat(fprime, 6) << endl;
      //printf("  Switching to Newton.  New fprime = %e\n", fprime);
    }
    else
    {
      step = selectLambda(fNew, fOld, fprime, step);
    }
    
    isFirstPass = false;

    // Check for linesearch failure
    if (step < minStep)
    {
      isFailed = true;
      message = "(FAILED - Min Step)";
      break;
    }

    // Update the number of linesearch iterations
    counter.incrementNumIterations();
    lsIterations ++;

    // Compute new trial point and its function value
    if ((lsType == Curvilinear) && !(isNewtonDirection))
    {
      bool ok = computeCurvilinearStep(tmpVec, oldSoln, s, step);
      // Note: oldSoln is needed above to get correct preconditioner 
      newsoln.computeX(oldSoln, tmpVec, 1.0);
    }
    else
    {
      newsoln.computeX(oldSoln, tmpVec, step);
    }
    newsoln.computeF();    
    fNew = 0.5 * newsoln.getNormF() * newsoln.getNormF();

    // Recompute convergence criteria based on new step
    threshold = fOld + alpha*step*fprime;
    isAcceptable = (fNew < threshold);
  }


  if (isFailed)
  {
    counter.incrementNumFailedLineSearches();
    step = recoveryStep;

    if (step != 0.0)
    {
      // Update the group using Newton direction and recovery step
      newsoln.computeX(oldSoln, newtonVec, step);
      newsoln.computeF();    
      fNew = 0.5 * newsoln.getNormF() * newsoln.getNormF();

      message = "(USING RECOVERY STEP!)";
    }
  }
  
  print.printStep(lsIterations, step, fOld, fNew, message);
  counter.setValues(paramsPtr->sublist("Line Search"));

  return (!isFailed);
}


double
NOX::Solver::TensorBased::getNormModelResidual(
                                       const NOX::Abstract::Vector& dir,
				       const NOX::Abstract::Group& soln,
				       bool isTensorModel) const
{
  NOX::Abstract::Vector* residualPtr = NULL;

  // Compute residual of Newton model...
  residualPtr = soln.getF().clone(ShapeCopy);
  soln.applyJacobian(dir, *residualPtr);
  numJvMults++;
  residualPtr->update(1.0, soln.getF(), 1.0);

  // Compute residual of Tensor model, if requested...
  if (isTensorModel)
  {
    double beta = scVecPtr->dot(dir);
    cout << " sc'*dt   = " << utils.sciformat(beta, 6) << endl;
    cout << " norm(dt) = " << utils.sciformat(dir.norm(), 6) << endl;
    residualPtr->update(beta*beta, *acVecPtr, 1.0);
  }

#ifdef LEAVE_OUT
  if (precondition == Left)
  {
    NOX::Abstract::Vector* tmpPtr = soln.getF().clone(ShapeCopy);
    *tmpPtr = *residualPtr;
    applyPreconditioner(false, soln, *localParamsPtr, *tmpPtr, *residualPtr,
			"compute");
    delete tmpPtr;
  }
#endif
  
  double modelNorm = residualPtr->norm();
  delete residualPtr;
  return modelNorm;
}


void
NOX::Solver::TensorBased::printDirectionInfo(char* dirName,
					const NOX::Abstract::Vector& dir,
					const NOX::Abstract::Group& soln,
					bool isTensorModel) const
{
  double residual = getNormModelResidual(dir, soln, isTensorModel);
 
  cout << " " << dirName << " norm of model residual =   "
       << utils.sciformat(residual, 6) << " (abs)     "
       << utils.sciformat(residual / soln.getNormF(), 6) << " (rel)" << endl;
  double fprime = getDirectionalDerivative(dir, soln);
  cout << " " << dirName << " directional derivative =  "
       << utils.sciformat(fprime, 6) << " (abs)    "
       << utils.sciformat(fprime / dir.norm(), 6) << " (rel)" << endl;
  cout << " " << dirName << " norm = "
       << utils.sciformat(dir.norm(), 6) << endl;
}


double NOX::Solver::TensorBased::getDirectionalDerivative(
				       const NOX::Abstract::Vector& dir,
				       const NOX::Abstract::Group& soln) const
{
  NOX::Abstract::Vector* tmpPtr = soln.getF().clone(ShapeCopy);
  soln.applyJacobian(dir,*tmpPtr);
  numJvMults++;
  double fprime = tmpPtr->dot(soln.getF());
  delete tmpPtr;
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


void NOX::Solver::TensorBased::throwError(const string& functionName,
					  const string& errorMsg) const
{
  if (utils.isPrintProcessAndType(NOX::Utils::Error))
    cerr << "NOX::Solver::TensorBased::" << functionName
	 << " - " << errorMsg << endl;
  throw "NOX Error";
}

#endif  // WITH_PRERELEASE
