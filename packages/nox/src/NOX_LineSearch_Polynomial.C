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

#include "NOX_LineSearch_Polynomial.H"

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

NOX::LineSearch::Polynomial::Polynomial(const NOX::Utils& u, Parameter::List& params) :
  Common(u, params),
  paramsPtr(NULL)
{
  reset(params);
}

NOX::LineSearch::Polynomial::~Polynomial()
{

}

bool NOX::LineSearch::Polynomial::reset(Parameter::List& params)
{ 
  NOX::Parameter::List& p = params.sublist("Polynomial");
  

  string choice = p.getParameter("Convergence Criteria", "Armijo-Goldstein");

  if (choice == "Ared/Pred") 
    convCriteria = AredPred;
  else if (choice == "None")
    convCriteria = None;
  else {
    convCriteria = ArmijoGoldstein;
  }

  choice = p.getParameter("Interpolation Type", "Cubic");
  if (choice == "Quadratic")
    interpolationType = Quadratic;
  else if (choice == "Cubic") 
    interpolationType = Cubic;
  else {
    cerr << "NOX::LineSearch::Polynomial::reset - Invalid \"Interpolation Type\"" << endl;
    throw "NOX Error";
  }

  minStep = p.getParameter("Minimum Step", 1.0e-12);
  defaultStep = p.getParameter("Default Step", 1.0);
  recoveryStep = p.getParameter("Recovery Step", defaultStep);
  maxIters = p.getParameter("Max Iters", 100);
  alpha = p.getParameter("Alpha Factor", 1.0e-4);
  minBoundFactor = p.getParameter("Min Bounds Factor", 0.1);
  maxBoundFactor = p.getParameter("Max Bounds Factor", 0.9);
  doForceInterpolation = p.getParameter("Force Interpolation", false);
  paramsPtr = &params;
  totalNumLineSearchCalls = 0;
  totalNumNonTrivialLineSearches = 0;
  totalNumFailedLineSearches = 0;
  totalNumIterations = 0;

  allowIncrease = params.isParameter("Allowed Relative Increase");
  if(allowIncrease) {
    relIncrease = params.getParameter("Allowed Relative Increase", 1.e2);
    numAllowed = params.getParameter("Maximum Increase Steps", maxIters);
    if(relIncrease<=0) {
      cerr << "NOX::LineSearch::NM_Polynomial::reset - Invalid \"Allowed Relative Increase\"" << endl;
      throw "NOX Error";
    }
  }

  return true;
}

bool NOX::LineSearch::Polynomial::compute(Abstract::Group& newGrp, double& step, 
			 const Abstract::Vector& dir,
			 const Solver::Generic& s) 
{

  if (utils.isPrintProcessAndType(Utils::InnerIteration)) 
  {
   cout << "\n" << Utils::fill(72) << "\n" << "-- Polynomial Line Search -- \n";
  }

  int nIters = 1;
  totalNumLineSearchCalls += 1;

  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  double oldf = 0.5 * oldGrp.getNormF() * oldGrp.getNormF();  

  step = defaultStep;
  newGrp.computeX(oldGrp, dir, step);
  newGrp.computeF();    
  double newf = 0.5 * newGrp.getNormF() * newGrp.getNormF();  

  

  // Get the linear solve tolerance if doing ared/pred for conv criteria
  double eta_original = 0.0;
  double eta = 0.0;
  if (convCriteria == AredPred) 
  {
    const NOX::Parameter::List& p = s.getParameterList();
    eta_original = p.sublist("Direction").sublist("Newton").sublist("Linear Solver").getParameter("Tolerance", -1.0);
    eta = eta_original;
  }

  bool isConverged = false;
  bool isFailed = false;
  double slope = computeSlope(dir, oldGrp);

  if (slope >= 0)
    isFailed = true;
  else if (doForceInterpolation)
    isConverged = false;
  else {
    isConverged = isSufficientDecrease(newf, oldf, step, slope, eta);
    if(allowIncrease)
      isConverged = (isConverged || 
                     isIncreaseAllowed(newf, oldf, totalNumLineSearchCalls) );
  }

  // Increment the number of newton steps requiring a line search
  if (!isConverged)
    totalNumNonTrivialLineSearches += 1;

  double prevf = 0;
  double previousStep = 0;
  double tempStep;
  bool isFirstPass = true;
  while ((!isConverged) && (!isFailed)) {

    if (nIters > maxIters) {
      isFailed = true;
      break;
    }
	
    printStep(nIters, step, oldf, newf);
    
    totalNumIterations += 1;
    nIters ++;
    
    if ((isFirstPass) || (interpolationType == Quadratic)) 
    {
      
      /* Quadratic Interpolation */
      
      tempStep = - (slope * step * step) / (2.0 * (newf - oldf - step * slope)) ;
      isFirstPass = false;
    }

    else 
    {
    
      /*   Cubic Interpolation */
      
      double term1 = newf - oldf - step * slope ;
      double term2 = prevf - oldf - previousStep * slope ;
      
      double a = 1.0 / (step - previousStep) * 
	(term1 / (step * step) - term2 / (previousStep * previousStep)) ;
      
      double b = 1.0 / (step - previousStep) *
	(-1.0 * term1 * previousStep / (step * step) +
	 term2 * step / (previousStep * previousStep)) ;
      
      double disc = b * b - 3.0 * a * slope ;
      
      if (disc < 0) 
      {
	isFailed = true;
	break;
      }
      
      if (fabs(a) < 1.e-12) 
      {
	tempStep = -slope / (2.0 * b);
      }
      else 
      {
	tempStep = (-b + sqrt(disc))/ (3.0 * a);
      }
      

      if (tempStep > 0.5 * step) 
	tempStep = 0.5 * step ;
    }
    
    previousStep = step ;
    prevf = newf ; 

    if (tempStep < minBoundFactor * step) 
      step *= minBoundFactor;
    else if ((nIters > 2) && (tempStep > maxBoundFactor * step))
      step *= maxBoundFactor;
    else 
      step = tempStep;


    if (step < minStep) 
    {
      isFailed = true;
      break;
    }
    
    newGrp.computeX(oldGrp, dir, step);
    newGrp.computeF();    
    newf = 0.5 * newGrp.getNormF() * newGrp.getNormF(); 
    
    eta = 1.0 - step * (1.0 - eta);
    isConverged = isSufficientDecrease(newf, oldf, step, slope, eta);
    if(allowIncrease)
      isConverged = (isConverged || 
                     isIncreaseAllowed(newf, oldf, totalNumLineSearchCalls) );
    
  } // End while loop 

  string message = "(STEP ACCEPTED!)";

  if (isFailed) {

    totalNumFailedLineSearches += 1;
    step = recoveryStep;
    newGrp.computeX(oldGrp, dir, step);
    newGrp.computeF();    
    newf = 0.5 * newGrp.getNormF() * newGrp.getNormF(); 

    eta = 1.0 - step * (1.0 - eta_original);
    paramsPtr->setParameter("Adjusted Tolerance", eta);
    
    message = "(USING RECOVERY STEP!)";
  }

  printStep(nIters, step, oldf, newf, message);
  paramsPtr->setParameter("Adjusted Tolerance", eta);
  setOutputParameters();
  return (!isFailed);
}

bool NOX::LineSearch::Polynomial::isSufficientDecrease(double newf, double oldf, double step, double slope, double eta) const
{
  double rhs = 0.0;
  if (convCriteria == ArmijoGoldstein) 
  {
    rhs = oldf + alpha * step * slope;
  }
  else if (convCriteria == AredPred) 
  {
    rhs = oldf * (1.0 - alpha * (1.0 - eta));
  }
  else if (convCriteria == None)
    return true;
  else 
  {
    cerr << "NOX::LineSearch::Polynomial::isSufficientDecrease - Invalid convergence criteria" << endl;
    throw "NOX Error";
  }

  return (newf <= rhs);
}

bool NOX::LineSearch::Polynomial::setOutputParameters() {
  NOX::Parameter::List& outputList = paramsPtr->sublist("Output");
  outputList.setParameter("Total Number of Line Search Calls", totalNumLineSearchCalls);
  outputList.setParameter("Total Number of Non-trivial Line Searches", totalNumNonTrivialLineSearches);
  outputList.setParameter("Total Number of Failed Line Searches", totalNumFailedLineSearches);
  outputList.setParameter("Total Number of Line Search Inner Iterations", totalNumIterations);
  return true;
}

