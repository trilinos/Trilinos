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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "NOX_Solver_TrustRegion.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Parameter_List.H"
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

TrustRegion::TrustRegion(Abstract::Group& grp, Status::Test& t, const Parameter::List& p) :
  solnPtr(&grp),		// pointer to grp
  oldSolnPtr(grp.clone(DeepCopy)), // create via clone
  oldSoln(*oldSolnPtr),		// reference to just-created pointer
  newtonVecPtr(grp.getX().clone(CopyShape)), // create via clone 
  newtonVec(*newtonVecPtr),	// reference to just-created pointer
  cauchyVecPtr(grp.getX().clone(CopyShape)), // create via clone 
  cauchyVec(*cauchyVecPtr),	// reference to just-created pointer
  aVecPtr(grp.getX().clone(CopyShape)), // create via clone 
  aVec(*aVecPtr),		// reference to just-created pointer
  bVecPtr(grp.getX().clone(CopyShape)), // create via clone 
  bVec(*bVecPtr),		// reference to just-created pointer
  testPtr(&t),			// pointer to t
  iparams(p),			// copy p
  newton(),			// initialize direction
  cauchy()			// initialize direction
{
  init();
}

// Protected
void TrustRegion::init()
{
  // Initialize 
  niter = 0;
  dx = 0;
  status = Status::Unconverged;

  // Set up utilities (i.e., set print processor, etc)
  Utils::setUtils(iparams);
  
  // Print out initialization information
  if (Utils::doPrint(Utils::Parameters)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    iparams.print(cout,5);
  }

  // Compute RHS of initital guess
  solnPtr->computeRHS();
  fnew = 0.5 * solnPtr->getNormRHS() * solnPtr->getNormRHS();

  // Get parameter settings
  if (!iparams.sublist("Direction").isParameter("Method"))
    iparams.sublist("Direction").setParameter("Method", "Newton");

  if (!iparams.sublist("Cauchy Direction").isParameter("Method"))
    iparams.sublist("Cauchy Direction").setParameter("Method", "Steepest Descent");

  if (!iparams.sublist("Cauchy Direction").isParameter("Scaling Type"))
    iparams.sublist("Cauchy Direction").setParameter("Scaling Type", "Quadratic Model Min");

  newton.reset(iparams.sublist("Direction"));
  cauchy.reset(iparams.sublist("Cauchy Direction"));

  initRadius = iparams.getParameter("Initial Radius", 0.0);
  radius = initRadius;

  /*
  if (radius <= 0) {
    cerr << "NOX::Solver::TrustRegion::init - Invalid \"Initial Radius\" (" 
	 << radius << ")" << endl;
    throw "NOX Error";
  }
  */

  minRadius = iparams.getParameter("Minimum Trust Region Radius", 1.0e-6);

  if (minRadius <= 0) {
    cerr << "NOX::Solver::TrustRegion::init - Invalid \"Minimum Trust Region Radius\" (" 
	 << minRadius << ")" << endl;
    throw "NOX Error";
  }

  maxRadius = iparams.getParameter("Maximum Trust Region Radius", 1.0e+10);

  if (maxRadius <= minRadius) {
    cerr << "NOX::Solver::TrustRegion::init - Invalid \"Maximum Trust Region Radius\" (" 
	 << maxRadius << ")" << endl;
    throw "NOX Error";
  }

  minRatio = iparams.getParameter("Minimum Improvement Ratio", 1.0e-4);

  if (minRatio <= 0) {
    cerr << "NOX::Solver::TrustRegion::init - Invalid \"Minimum Improvement Ratio\" (" 
	 << minRatio << ")" << endl;
    throw "NOX Error";
  }

  contractTriggerRatio = iparams.getParameter("Contraction Trigger Ratio", 0.1);

  if (contractTriggerRatio < minRatio) {
    cerr << "NOX::Solver::TrustRegion::init - Invalid \"Contraction Trigger Ratio\" (" 
	 << contractTriggerRatio << ")" << endl;
    throw "NOX Error";
  }


  expandTriggerRatio = iparams.getParameter("Expansion Trigger Ratio", 0.75);

  if (expandTriggerRatio <= contractTriggerRatio) {
    cerr << "NOX::Solver::TrustRegion::init - Invalid \"Expansion Trigger Ratio\" (" 
	 << expandTriggerRatio << ")" << endl;
    throw "NOX Error";
  }

  contractFactor = iparams.getParameter("Contraction Factor", 0.25);

  if ((contractFactor <= 0) || (contractFactor >= 1)) {
    cerr << "NOX::Solver::TrustRegion::init - Invalid \"Contraction Factor\" (" 
	 << contractFactor << ")" << endl;
    throw "NOX Error";
  }

  expandFactor = iparams.getParameter("Expansion Factor", 4.0);

  if (expandFactor <= 1) {
    cerr << "NOX::Solver::TrustRegion::init - Invalid \"Expansion Factor\" (" 
	 << expandFactor << ")" << endl;
    throw "NOX Error";
  }

  recoveryStep = iparams.getParameter("Recovery Step", 1.0);

  if (recoveryStep < 0) {
    cerr << "NOX::Solver::TrustRegion::init - Invalid \"Recovery Step\" (" 
	 << recoveryStep << ")" << endl;
    throw "NOX Error";
  }


  // Test the initial guess
  status = testPtr->operator()(*this);

  if (Utils::doPrint(Utils::Parameters)) {
    cout << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testPtr->print(cout, 5);
    cout <<"\n" << Utils::fill(72) << "\n";
  }

}

bool TrustRegion::reset(Abstract::Group& grp, Status::Test& t, const Parameter::List& p) 
{
  solnPtr = &grp;
  testPtr = &t;
  iparams = p;			
  init();
  return true;
}

TrustRegion::~TrustRegion() 
{
  delete oldSolnPtr;
}


Status::StatusType TrustRegion::getStatus()
{
  return status;
}

Status::StatusType TrustRegion::iterate()
{
  // First check status
  if (status != Status::Unconverged) 
    return status;

  // Copy pointers into temporary references
  Abstract::Group& soln = *solnPtr;
  Status::Test& test = *testPtr;

  // Compute Cauchy and Newton points
  bool ok;
  ok = newton(newtonVec, soln, *this);
  if (!ok) {
    cout << "NOX::Solver::TrustRegion::iterate - unable to calculate Newton direction" << endl;
    status = Status::Failed;
    return status;
  }

  ok = cauchy(cauchyVec, soln, *this);
  if (!ok) {
    cerr << "NOX::Solver::TrustRegion::iterate - unable to calculate Cauchy direction" << endl;
    status = Status::Failed;
    return status;
  }

  if (niter == 0) {
    radius = newtonVec.norm();
    if (radius < minRadius)
      radius = 2 * minRadius;
  }

  // Copy current soln to the old soln.
  oldSoln = soln;
  fold = fnew;

  //! Improvement ratio = (fold - fnew) / (mold - mnew)
  double ratio = -1;

  if (Utils::doPrint(Utils::InnerIteration)) {
    cout << Utils::fill(72) << endl;
    cout << "-- Trust Region Inner Iteration --" << endl;
  }


  // Trust region subproblem loop
  while ((ratio < minRatio) && (radius > minRadius)) {

    Abstract::Vector* dirPtr;
    double step;

    // Trust region step
    if (newtonVec.norm() <= radius) {
      stepType = TrustRegion::Newton;
      step = 1.0;
      dirPtr = &newtonVec;
    }
    else if (cauchyVec.norm() >= radius) {
      stepType = TrustRegion::Cauchy;
      step = radius / cauchyVec.norm();
      dirPtr = &cauchyVec;
    }
    else {			// Dogleg computation
      

      // aVec = newtonVec - cauchyVec
      aVec.update(1.0, newtonVec, -1.0, cauchyVec, 0.0);

      // cta = cauchyVec' * aVec
      double cta = cauchyVec.dot(aVec);	
      // ctc = cauchyVec' * cauchyVec
      double ctc = cauchyVec.dot(cauchyVec);
      // ata = aVec' * aVec
      double ata = aVec.dot(aVec);

      // sqrt of quadratic equation
      double tmp = (cta * cta) - ((ctc - (radius * radius)) * ata);
      if (tmp < 0) {
	cerr << "NOX::Solver::TrustRegion::iterate - invalid computation" << endl;
	throw "NOX Error";
      }
      
      // final soln to quadratic equation
      double gamma = (sqrt(tmp) - cta) / ata;
      if ((gamma < 0) || (gamma > 1)) {
	cerr << "NOX::Solver::TrustRegion::iterate - invalid trust region step" << endl;
	throw "NOX Error";
      }
      
      // final direction computation
      aVec.update(1.0 - gamma, cauchyVec, gamma, newtonVec, 0.0);

      // solution
      stepType = TrustRegion::Dogleg;
      dirPtr = &aVec;
      step = 1.0;
    }
    
    // Local reference to use in the remaining computation
    const Abstract::Vector& dir = *dirPtr;

    // Calculate true step length
    dx = step * dir.norm();
    
    // Compute new X
    soln.computeX(oldSoln, dir, step);

    // Compute RHS for new current solution.
    ok = soln.computeRHS();
    if (!ok) {
      cerr << "NOX::Solver::TrustRegion::iterate - unable to compute RHS" << endl;
      throw "NOX Error";
    }

    // Compute ratio of actual to predicted reduction
    fnew = 0.5 * solnPtr->getNormRHS() * solnPtr->getNormRHS();
    if (fnew >= fold) {
      ratio = -1;
    }
    else {

      ok = oldSoln.applyJacobian(*dirPtr, bVec);
      if (!ok) {
	cerr << "NOX::Solver::TrustRegion::iterate - unable to compute RHS" << endl;
	throw "NOX Error";
      }
      double numerator = fold - fnew;
      double denominator = abs(dir.dot(oldSoln.getGrad()) + 0.5 * bVec.dot(bVec));
      ratio = numerator / denominator;
      cout << "Ratio computation: " << Utils::sci(numerator) << "/" 
	   << Utils::sci(denominator) << "=" << ratio << endl;
      if ((denominator < 1.0e-12) && ((fnew / fold) >= 0.5))
	ratio = -1;
    }


    if (Utils::doPrint(Utils::InnerIteration)) {
      cout << "radius = " << Utils::sci(radius, 1);
      cout << " ratio = " << setprecision(1) << setw(3) << ratio;
      cout << " f = " << Utils::sci(sqrt(2*fnew));
      cout << " fold = " << Utils::sci(sqrt(2*fold));
      cout << " ";

      switch(stepType) {
      case TrustRegion::Newton:
	cout << "Newton";
	break;
      case TrustRegion::Cauchy:
	cout << "Cauchy";
	break;
      case TrustRegion::Dogleg:
	cout << "Dogleg";
	break;
      }

      cout << endl;
    }

    // Update trust region
    if (ratio < contractTriggerRatio) {
      if (stepType == TrustRegion::Newton)
	radius = newtonVec.norm();
      radius = max(contractFactor * radius, minRadius);
    }
    else if ((ratio > expandTriggerRatio) && (dx == radius)) {
      radius = min(expandFactor * radius, maxRadius);
    }

  }

  // Update iteration count.
  niter ++;


  // Evaluate the current status
  if ((radius <= minRadius) && (ratio < minRatio)) {
    if (Utils::doPrint(Utils::InnerIteration))
      cout << "Using recovery step and resetting trust region." << endl;
    soln.computeX(oldSoln, newtonVec, recoveryStep);
    soln.computeRHS();
    radius = newtonVec.norm();
    /*if (radius < minRadius)
      radius = 2 * minRadius;*/
  }

  status = test(*this);
 
  if (Utils::doPrint(Utils::InnerIteration)) 
    cout << Utils::fill(72) << endl;

  // Return status.
  return status;
}

Status::StatusType TrustRegion::solve()
{
  printUpdate();

  // Iterate until converged or failed
  while (status == Status::Unconverged) {
    status = iterate();
    printUpdate();
  }

  return status;
}

const Abstract::Group& TrustRegion::getSolutionGroup() const
{
  return *solnPtr;
}

const Abstract::Group& TrustRegion::getPreviousSolutionGroup() const
{
  return oldSoln;
}

int TrustRegion::getNumIterations() const
{
  return niter;
}

const Parameter::List& TrustRegion::getOutputParameters() const
{
  oparams.setParameter("Nonlinear Iterations", niter);
  oparams.setParameter("2-Norm of Residual", solnPtr->getNormRHS());
  return oparams;
}

// protected
void TrustRegion::printUpdate() 
{
  double fmax = solnPtr->getRHS().norm(Abstract::Vector::INF);
  if (Utils::doPrint(Utils::OuterIteration)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "-- Newton Trust-Region Step " << niter << " -- \n";
    cout << "f = " << Utils::sci(sqrt(2*fnew));
    cout << " fmax = " << Utils::sci(fmax);
    cout << "  dx = " << Utils::sci(dx);
    cout << "  radius = " << Utils::sci(radius);
    if (status == Status::Converged)
      cout << " (Converged!)";
    if (status == Status::Failed)
      cout << " (Failed!)";
    cout << "\n" << Utils::fill(72) << "\n" << endl;
  }
  
  if ((status != Status::Unconverged) && 
      (Utils::doPrint(Utils::OuterIteration))) {
    cout << Utils::fill(72) << "\n";
    cout << "-- Final Status Test Results --\n";    
    testPtr->print(cout);
    cout << Utils::fill(72) << "\n";
  }
}

