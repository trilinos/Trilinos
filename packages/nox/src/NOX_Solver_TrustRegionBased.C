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
#include "Teuchos_ParameterList.hpp"
#include "NOX_MeritFunction_Generic.H"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"

using namespace NOX;
using namespace NOX::Solver;

TrustRegionBased::
TrustRegionBased(const Teuchos::RefCountPtr<Abstract::Group>& grp,
		 const Teuchos::RefCountPtr<StatusTest::Generic>& t,
		 const Teuchos::RefCountPtr<Teuchos::ParameterList>& p) :
  globalDataPtr(Teuchos::rcp(new NOX::GlobalData(p))),
  utilsPtr(globalDataPtr->getUtils()), 
  solnPtr(grp),		
  oldSolnPtr(grp->clone(DeepCopy)), 
  newtonVecPtr(grp->getX().clone(ShapeCopy)), 
  cauchyVecPtr(grp->getX().clone(ShapeCopy)), 
  aVecPtr(grp->getX().clone(ShapeCopy)), 
  bVecPtr(grp->getX().clone(ShapeCopy)), 
  testPtr(t),			
  paramsPtr(p),			
  newton(globalDataPtr),
  cauchy(globalDataPtr),
  meritFuncPtr(globalDataPtr->getMeritFunction()),
  useAredPredRatio(false),
  prePostOperator(utilsPtr, paramsPtr->sublist("Solver Options"))
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

  newton.reset(globalDataPtr, paramsPtr->sublist("Direction"));
  cauchy.reset(globalDataPtr, paramsPtr->sublist("Cauchy Direction"));

  minRadius = paramsPtr->sublist("Trust Region").
    get("Minimum Trust Region Radius", 1.0e-6);
  if (minRadius <= 0) 
    invalid("Minimum Trust Region Radius", minRadius);

  maxRadius = paramsPtr->sublist("Trust Region").get("Maximum Trust Region Radius", 1.0e+10);
  if (maxRadius <= minRadius) 
    invalid("Maximum Trust Region Radius", maxRadius);

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

  // Get the checktype
  checkType = (NOX::StatusTest::CheckType) paramsPtr->
    sublist("Solver Options").get("Status Test Check Type", 
					   NOX::StatusTest::Minimal);

  // Check for the using Homer Walker's Ared/Pred ratio calculation
  useAredPredRatio = 
    paramsPtr->sublist("Trust Region").
    get("Use Ared/Pred Ratio Calculation", false);

}

//PRIVATE
void NOX::Solver::TrustRegionBased::
invalid(const string& name, double value) const
{
  utilsPtr->err() << "NOX::Solver::TrustRegionBased::init - " 
       << "Invalid \"" << name << "\" (" << value << ")" 
       << endl;
  throw "NOX Error";
}

bool TrustRegionBased::
reset(const Teuchos::RefCountPtr<Abstract::Group>& grp, 
      const Teuchos::RefCountPtr<StatusTest::Generic>& t, 
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& p) 
{
  globalDataPtr = Teuchos::rcp(new NOX::GlobalData(p));
  utilsPtr = globalDataPtr->getUtils(); 
  solnPtr = grp;
  testPtr = t;
  paramsPtr = p;			
  prePostOperator.reset(utilsPtr, paramsPtr->sublist("Solver Options"));
  init();
  return true;
}

bool TrustRegionBased::
reset(const Teuchos::RefCountPtr<Abstract::Group>& grp,
      const Teuchos::RefCountPtr<StatusTest::Generic>& t)
{
  // New initial guess and status test
  solnPtr = grp;
  testPtr = t;

  // Initialize 
  nIter = 0;
  dx = 0;
  status = StatusTest::Unconverged;

  // Print out initialization information
  if (utilsPtr->isPrintType(NOX::Utils::Parameters)) {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utilsPtr->out(),5);
  }

  return true;
}

bool TrustRegionBased::
reset(const Teuchos::RefCountPtr<Abstract::Group>& grp)
{
  // New initial guess and status test
  solnPtr = grp;

  // Initialize 
  nIter = 0;
  dx = 0;
  status = StatusTest::Unconverged;

  // Print out initialization information
  if (utilsPtr->isPrintType(NOX::Utils::Parameters)) {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utilsPtr->out(),5);
  }

  return true;
}

TrustRegionBased::~TrustRegionBased() 
{
  
}


NOX::StatusTest::StatusType TrustRegionBased::getStatus()
{
  return status;
}

NOX::StatusTest::StatusType TrustRegionBased::step()
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

  // First check status
  if (status != StatusTest::Unconverged) 
    return status;

  // Copy pointers into temporary references
  Abstract::Group& soln = *solnPtr;
  StatusTest::Generic& test = *testPtr;

  // Compute Cauchy and Newton points
  bool ok;
  ok = newton.compute(*newtonVecPtr, soln, *this);
  if (!ok) 
  {
    utilsPtr->out() << "NOX::Solver::TrustRegionBased::iterate - unable to calculate Newton direction" << endl;
    status = StatusTest::Failed;
    prePostOperator.runPostIterate(*this);
    return status;
  }

  ok = cauchy.compute(*cauchyVecPtr, soln, *this);
  if (!ok) 
  {
    utilsPtr->err() << "NOX::Solver::TrustRegionBased::iterate - unable to calculate Cauchy direction" << endl;
    status = StatusTest::Failed;
    prePostOperator.runPostIterate(*this);
    return status;
  }

  if (nIter == 0) 
  {
    radius = newtonVecPtr->norm();
    
    if (radius < minRadius)
      radius = 2 * minRadius;
  }

  // Update iteration count.
  nIter ++;

  // Copy current soln to the old soln.
  *oldSolnPtr = *solnPtr;
  // RPP: Can't just copy over oldf.  Scaling could change between iterations
  // so user Merit Functions could be out of sync
  oldF = meritFuncPtr->computef(*oldSolnPtr);

  // Improvement ratio = (oldF - newF) / (mold - mnew)
  double ratio = -1;

  if (utilsPtr->isPrintType(NOX::Utils::InnerIteration)) 
  {
    utilsPtr->out() << NOX::Utils::fill(72) << endl;
    utilsPtr->out() << "-- Trust Region Inner Iteration --" << endl;
  }

  // Trust region subproblem loop
  while ((ratio < minRatio) && (radius > minRadius)) 
  {

    Teuchos::RefCountPtr<NOX::Abstract::Vector> dirPtr;
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
	utilsPtr->err() << "NOX::Solver::TrustRegionBased::iterate - invalid computation" << endl;
	throw "NOX Error";
      }
      
      // final soln to quadratic equation
      double gamma = (sqrt(tmp) - cta) / ata;
      if ((gamma < 0) || (gamma > 1)) {
	utilsPtr->err() << "NOX::Solver::TrustRegionBased::iterate - invalid trust region step" << endl;
	throw "NOX Error";
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

    // Calculate true step length
    dx = step * dir.norm();

    // Compute new X
    soln.computeX(*oldSolnPtr, dir, step);

    // Compute F for new current solution.
    NOX::Abstract::Group::ReturnType rtype = soln.computeF();
    if (rtype != NOX::Abstract::Group::Ok) 
    {
      utilsPtr->err() << "NOX::Solver::TrustRegionBased::iterate - unable to compute F" << endl;
      throw "NOX Error";
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
	     << "unable to compute F" << endl;
	throw "NOX Error";
      }
      bVecPtr->update(1.0, oldSolnPtr->getF(), 1.0);

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
			<< ratio << endl;
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
	  utilsPtr->err() << "NOX::Solver::TrustRegionBased::iterate - unable to compute F" << endl;
	  throw "NOX Error";
	}
	double numerator = oldF - newF;
	double denominator = 0.0;
	  
	denominator = fabs(oldF - meritFuncPtr->
			   computeQuadraticModel(dir,*oldSolnPtr));
	
	ratio = numerator / denominator;
	if (utilsPtr->isPrintType(NOX::Utils::Debug))
	  utilsPtr->out() << "Ratio computation: " 
			  << utilsPtr->sciformat(numerator) << "/" 
			  << utilsPtr->sciformat(denominator) << "=" 
			  << utilsPtr->sciformat(ratio) << endl;
	
	// WHY IS THIS CHECK HERE?
	if ((denominator < 1.0e-12) && ((newF / oldF) >= 0.5))
	  ratio = -1;
      }
    }

    if (utilsPtr->isPrintType(Utils::InnerIteration)) {
      utilsPtr->out() << "radius = " << utilsPtr->sciformat(radius, 1);
      utilsPtr->out() << " ratio = " << setprecision(2) << setw(4) << ratio;
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

      utilsPtr->out() << endl;
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
      utilsPtr->out() << "Using recovery step and resetting trust region." << endl;
    soln.computeX(*oldSolnPtr, *newtonVecPtr, recoveryStep);
    soln.computeF();
    radius = newtonVecPtr->norm();
    /*if (radius < minRadius)
      radius = 2 * minRadius;*/
  }

  status = test.checkStatus(*this, checkType);
 
  if (utilsPtr->isPrintType(Utils::InnerIteration)) 
    utilsPtr->out() << NOX::Utils::fill(72) << endl;

  prePostOperator.runPostIterate(*this);

  return status;
}

NOX::StatusTest::StatusType TrustRegionBased::solve()
{
  prePostOperator.runPreSolve(*this);

  // Iterate until converged or failed
  while (status == StatusTest::Unconverged) {
    status = step();
    printUpdate();
  }

  Teuchos::ParameterList& outputParams = paramsPtr->sublist("Output");
  outputParams.set("Nonlinear Iterations", nIter);
  outputParams.set("2-Norm of Residual", solnPtr->getNormF());

  prePostOperator.runPostSolve(*this);

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
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n" << endl;
  }
  
  if ((status != StatusTest::Unconverged) && 
      (utilsPtr->isPrintType(NOX::Utils::OuterIteration))) {
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Final Status Test Results --\n";    
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
  }
}

