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

#ifdef WITH_PRERELEASE

#include "NOX_LineSearch_PolynomialWalker.H"

#include "NOX_LineSearch_Utils_Printing.H"
#include "NOX_LineSearch_Utils_Slope.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Parameter_List.H"

NOX::LineSearch::PolynomialWalker::PolynomialWalker(const NOX::Utils& u, Parameter::List& params) :
  paramsPtr(0),
  print(u)
{
  reset(params);
}

NOX::LineSearch::PolynomialWalker::~PolynomialWalker()
{

}

bool NOX::LineSearch::PolynomialWalker::reset(Parameter::List& params)
{ 
  NOX::Parameter::List& p = params.sublist("PolynomialWalker");
  
  string choice = p.getParameter("Sufficient Decrease Condition", "Armijo-Goldstein");
  if (choice == "Ared/Pred") 
    convCriteria = AredPred;
  else if (choice == "None")
    convCriteria = None;
  else 
    convCriteria = ArmijoGoldstein;

  choice = p.getParameter("Interpolation Type", "Cubic");
  if (choice == "Quadratic")
    interpolationType = Quadratic;
  else if (choice == "Cubic") 
    interpolationType = Cubic;
  else 
  {
    cerr << "NOX::LineSearch::PolynomialWalker::reset - Invalid \"Interpolation Type\"" << endl;
    throw "NOX Error";
  }

  minStep = p.getParameter("Minimum Step", 1.0e-12);
  defaultStep = p.getParameter("Default Step", 1.0);
  recoveryStep = p.getParameter("Recovery Step", defaultStep);
  maxIters = p.getParameter("Max Iters", 100);
  alpha = p.getParameter("Alpha Factor", 1.0e-4);
  minBoundFactor = p.getParameter("Min Bounds Factor", 0.1);
  maxBoundFactor = p.getParameter("Max Bounds Factor", 0.5);
  doForceInterpolation = p.getParameter("Force Interpolation", false);
  paramsPtr = &params;

  useScaledNorms = p.getParameter("Use Scaled Norms", false);

  choice = p.getParameter("Merit Function", "Norm-F Squared");
  if (choice == "Norm-F")
    meritFunctionType = NormF;
  else if (choice == "Norm-F Squared") 
    meritFunctionType = NormFSquared;
  else 
  {
    cerr << "NOX::LineSearch::PolynomialWalker::reset - Invalid \"Merit Function\"" << endl;
    throw "NOX Error";
  }

  allowIncrease = p.isParameter("Allowed Relative Increase");
  if(allowIncrease) 
  {
    relIncrease = p.getParameter("Allowed Relative Increase", 1.e2);
    numAllowed = p.getParameter("Maximum Increase Steps", maxIters);
    if(relIncrease <= 0) 
    {
      cerr << "NOX::LineSearch::NM_PolynomialWalker::reset - Invalid \"Allowed Relative Increase\"" << endl;
      throw "NOX Error";
    }
  }

  counter.reset();

  return true;
}

bool NOX::LineSearch::PolynomialWalker::compute(Abstract::Group& newGrp, double& step, 
			 const Abstract::Vector& dir,
			 const Solver::Generic& s) 
{

  if (print.isPrintProcessAndType(NOX::Utils::InnerIteration)) 
  {
    cout << "\n" << NOX::Utils::fill(72) << "\n" << "-- PolynomialWalker Line Search -- \n";
  }

  // Set counter to 1
  int nIters = 1;
  counter.incrementNumLineSearches();

  // Compute the scale factor for norms
  if (useScaledNorms) 
    factor = 1.0/(sqrt((double) (dir.length())));
  else
    factor = 1.0;

  // Get Old f(0)
  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  oldf_1 = oldGrp.getNormF() * factor;
  oldf_2 = 0.5 * oldGrp.getNormF() * oldGrp.getNormF() * factor * factor;
  oldf_interp = 0.0;
  if (meritFunctionType == NormF)
    oldf_interp = oldf_1;
  else 
    oldf_interp = oldf_2;

  // Get the slope f'(0)
  slope_2 = slopeObj.computeSlope(dir, oldGrp);
  slope_1 = slope_2 / oldf_1;
  slope_interp = 0.0;
  if (meritFunctionType == NormF)
    slope_interp = slope_1;
  else 
    slope_interp = slope_2;

  // Get New f(step)
  step = defaultStep;
  computeNewF(newGrp, oldGrp, dir, step);

  // Get the linear solve tolerance if doing ared/pred for conv criteria
  double eta_original = 0.0;
  double eta = 0.0;
  if (convCriteria == AredPred) 
  {
    const NOX::Parameter::List& p = s.getParameterList();
    eta_original = p.sublist("Direction").sublist("Newton")
                    .sublist("Linear Solver").getParameter("Tolerance", -1.0);
    eta = eta_original;
  }

  bool isConverged = false;
  bool isFailed = false;

  if (slope_interp >= 0.0)
    isFailed = true;
  else if (doForceInterpolation)
    isConverged = false;
  else 
  {
    isConverged = isSufficientDecrease(step, eta);
    if(allowIncrease) 
    {
      isConverged = ( isConverged || 
		      isIncreaseAllowed(counter.getNumLineSearches()) );
    }
  }

  // Increment the number of newton steps requiring a line search
  if (!isConverged)
    counter.incrementNumNonTrivialLineSearches();

  double prevf = 0.0;
  double previousStep = 0.0;
  double tempStep = 0.0;
  bool isFirstPass = true;
  while ((!isConverged) && (!isFailed)) {

    if (nIters > maxIters) 
    {
      isFailed = true;
      break;
    }
	
    if (convCriteria == AredPred)
      print.printStep(nIters, step, oldf_1, newf_1, "", false);
    else
      print.printStep(nIters, step, oldf_2, newf_2);
    
    counter.incrementNumIterations();
    nIters ++;
    
    if ((isFirstPass) || (interpolationType == Quadratic)) 
    {
      
      /* Quadratic Interpolation */
      
      tempStep = - (slope_interp * step * step) / (2.0 * (newf_interp - oldf_interp - step * slope_interp)) ;
      isFirstPass = false;

    }

    else 
    {
    
      /*   Cubic Interpolation */
      
      double term1 = newf_interp - oldf_interp - step * slope_interp ;
      double term2 = prevf - oldf_interp - previousStep * slope_interp ;
      
      double a = 1.0 / (step - previousStep) * 
	(term1 / (step * step) - term2 / (previousStep * previousStep)) ;
      
      double b = 1.0 / (step - previousStep) *
	(-1.0 * term1 * previousStep / (step * step) +
	 term2 * step / (previousStep * previousStep)) ;
      
      double disc = b * b - 3.0 * a * slope_interp ;
      
      if (disc < 0) 
      {
	isFailed = true;
	break;
      }
      
      // The folowing has been removed at the request of Homwer Walker
      /*
      if (fabs(a) < 1.e-12) 
      {
	tempStep = -slope / (2.0 * b);
      }
      else 
      {
	tempStep = (-b + sqrt(disc))/ (3.0 * a);
      }
      */
      // Homer suggests this test be used instead:
      if (b < 0)
        tempStep = (-b + sqrt(disc))/ (3.0 * a);
      else
        tempStep = -slope_interp/(b + sqrt(disc));
      
      /* Also removed at the request of Homer Walker
      if (tempStep > 0.5 * step) 
	tempStep = 0.5 * step ;
      */

    }
    
    previousStep = step ;
    prevf = newf_interp ; 

    if (tempStep < minBoundFactor * step) 
      step *= minBoundFactor;
    // Changed the following line at Homer Walker's request.
    //else if ((nIters > 2) && (tempStep > maxBoundFactor * step))
    else if (tempStep > maxBoundFactor * step)
      step *= maxBoundFactor;
    else 
      step = tempStep;


    if (step < minStep) 
    {
      isFailed = true;
      break;
    }
    
    computeNewF(newGrp, oldGrp, dir, step);
    
    eta = 1.0 - step * (1.0 - eta_original);
    isConverged = isSufficientDecrease(step, eta);
    if(allowIncrease)
      isConverged = (isConverged || 
                     isIncreaseAllowed(counter.getNumLineSearches()));
    
  } // End while loop 

  string message = "(STEP ACCEPTED!)";

  if (isFailed) {

    counter.incrementNumFailedLineSearches();
    step = recoveryStep;
    computeNewF(newGrp, oldGrp, dir, step);

    eta = 1.0 - step * (1.0 - eta_original);
    paramsPtr->setParameter("Adjusted Tolerance", eta);
    
    message = "(USING RECOVERY STEP!)";
  }

  if (convCriteria == AredPred)
    print.printStep(nIters, step, oldf_1, newf_1, message, false);
  else
    print.printStep(nIters, step, oldf_2, newf_2, message);
  paramsPtr->setParameter("Adjusted Tolerance", eta);
  counter.setValues(*paramsPtr);
  return (!isFailed);
}

bool NOX::LineSearch::PolynomialWalker::isSufficientDecrease(double step, 
							     double eta) const
{
  double rhs = 0.0;
  if (convCriteria == ArmijoGoldstein) 
  {
    rhs = oldf_interp + alpha * step * slope_interp;
    return (newf_interp <= rhs);
  }
  else if (convCriteria == AredPred) 
  {
    rhs = oldf_1 * (1.0 - alpha * (1.0 - eta));
    return (newf_1 <= rhs);
  }
  else if (convCriteria == None)
  {
    return true;
  }

  cerr << "NOX::LineSearch::PolynomialWalker::isSufficientDecrease - Invalid convergence criteria" << endl;
  throw "NOX Error";

}

bool NOX::LineSearch::PolynomialWalker::computeNewF(
                                          NOX::Abstract::Group& newGrp, 
					  const NOX::Abstract::Group& oldGrp, 
					  const NOX::Abstract::Vector& dir, 
					  double step)
{
  newGrp.computeX(oldGrp, dir, step);

  NOX::Abstract::Group::ReturnType status = newGrp.computeF();
  if (status != NOX::Abstract::Group::Ok)
    return false;

  newf_1 = newGrp.getNormF() * factor;
  newf_2 = 0.5 * newGrp.getNormF() * newGrp.getNormF() * factor * factor;  
  newf_interp = 0.0;
  if (meritFunctionType == NormF)
    newf_interp = newf_1;
  else 
    newf_interp = newf_2;

  return true;
}

bool NOX::LineSearch::PolynomialWalker::isIncreaseAllowed(int nOuterIters) const
{
  double increase = 0.0;
  if (convCriteria == AredPred)  
    increase = sqrt(newf_1 / oldf_1);
  else
    increase = sqrt(newf_interp / oldf_interp);

  return ( (increase <= relIncrease) && (nOuterIters <= numAllowed) );
}

#endif
