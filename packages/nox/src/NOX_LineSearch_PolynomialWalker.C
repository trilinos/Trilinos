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
#include "NOX_Parameter_UserNorm.H"
#include "NOX_Parameter_MeritFunction.H"

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
  paramsPtr = &params;

  NOX::Parameter::List& p = params.sublist("PolynomialWalker");
  
  string choice = p.getParameter("Sufficient Decrease Condition", 
				 "Armijo-Goldstein");
  if (choice == "Ared/Pred") 
    suffDecrCond = AredPred;
  else if (choice == "None")
    suffDecrCond = None;
  else 
    suffDecrCond = ArmijoGoldstein;   // default

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

  // Determine if a User Defined Norm is present
  userNormPtr = 0;
  if (p.isParameterArbitrary("User Defined Norm")){
    userNormPtr = const_cast<NOX::Parameter::UserNorm*>(dynamic_cast<const NOX::Parameter::UserNorm*>(&(p.getArbitraryParameter("User Defined Norm"))));
  }

  // Determine the merit function to use
  meritFuncPtr = 0;
  if (p.isParameterArbitrary("Merit Function")) {
    meritFuncPtr = const_cast<NOX::Parameter::MeritFunction*>(dynamic_cast<const NOX::Parameter::MeritFunction*>(&(p.getArbitraryParameter("Merit Function"))));
  }
  else {
    choice = p.getParameter("Merit Function", "Norm-F Squared");
    if (choice == "Norm-F Squared") {
      // Do nothing, this is the default
    }
    else {
      cerr << "NOX::LineSearch::PolynomialWalker::reset - Invalid "
	   << "\"Merit Function\"" << endl;
      throw "NOX Error";
    }
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

bool NOX::LineSearch::PolynomialWalker::compute(Abstract::Group& newGrp, 
						double& step, 
			                        const Abstract::Vector& dir,
						const Solver::Generic& s) 
{

  if (print.isPrintProcessAndType(NOX::Utils::InnerIteration)) 
  {
    cout << "\n" << NOX::Utils::fill(72) << "\n" << "-- PolynomialWalker Line Search -- \n";
  }

  // Print out Details of line search
  if (print.isPrintProcessAndType(NOX::Utils::Details)) {

    if (userNormPtr != 0)
      cout << "       Norms = Using a user defined norm" << endl;
    else 
      cout << "       Norms = L-2" << endl;
  
    if (meritFuncPtr != 0) 
      cout << "       Merit Function = User Defined" << endl;
    else 
      cout << "       Merit Function = 0.5 * || F || * || F ||" << endl;
  }

  // Set counter to 1
  int nIters = 1;
  counter.incrementNumLineSearches();


  // Get Old Group
  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();

  // Compute ||F(0)|| only for Ared/Pred condition
  if (suffDecrCond == AredPred) {
    if (userNormPtr != 0)
      normFOld = userNormPtr->norm(oldGrp.getF());
    else
      normFOld = oldGrp.getNormF();
  }

  // Compute merit function f(0)
  if (meritFuncPtr != 0) {
    meritFuncOld = meritFuncPtr->computef(oldGrp);
  }
  else {
    meritFuncOld = 0.5 * oldGrp.getNormF() * oldGrp.getNormF();
  }

  // Compute the slope f'(0)
  if (meritFuncPtr != 0) {
    slopeOld = meritFuncPtr->computeSlope(dir, oldGrp);
  }
  else {
    slopeOld = slopeObj.computeSlope(dir, oldGrp);
  }
  
  // Compute New f(step) and if doing Ared/Pred also do ||F(step)||
  step = defaultStep;
  computeNewF(newGrp, oldGrp, dir, step);

  //cout << "Old f = " << normFOld << endl;
  //cout << "Slope Old f = " << slope_1 << endl;
  //cout << "New f = " << normFNew << endl;

  // Get the linear solve tolerance if doing ared/pred for conv criteria
  double eta_original = 0.0;
  double eta = 0.0;
  if (suffDecrCond == AredPred) 
  {
    const NOX::Parameter::List& p = s.getParameterList();
    eta_original = p.sublist("Direction").sublist("Newton")
                    .sublist("Linear Solver").getParameter("Tolerance", -1.0);
    eta = eta_original;
  }

  bool isConverged = false;
  bool isFailed = false;

  if (slopeOld >= 0.0)
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
	
    if (suffDecrCond == AredPred)
      print.printStep(nIters, step, normFOld, normFNew, "", false);
    else
      print.printStep(nIters, step, meritFuncOld, meritFuncNew);
    
    counter.incrementNumIterations();
    nIters ++;
    
    if ((isFirstPass) || (interpolationType == Quadratic)) 
    {
      
      /* Quadratic Interpolation */
      
      tempStep = - (slopeOld * step * step) / (2.0 * (meritFuncNew - meritFuncOld - step * slopeOld)) ;
      isFirstPass = false;

    }

    else 
    {
    
      /*   Cubic Interpolation */
      
      double term1 = meritFuncNew - meritFuncOld - step * slopeOld ;
      double term2 = prevf - meritFuncOld - previousStep * slopeOld ;
      
      double a = 1.0 / (step - previousStep) * 
	(term1 / (step * step) - term2 / (previousStep * previousStep)) ;
      
      double b = 1.0 / (step - previousStep) *
	(-1.0 * term1 * previousStep / (step * step) +
	 term2 * step / (previousStep * previousStep)) ;
      
      double disc = b * b - 3.0 * a * slopeOld ;
      
      if (disc < 0) 
      {
	isFailed = true;
	break;
      }
      
      // The folowing has been removed at the request of Homwer Walker
      /*
      if (fabs(a) < 1.e-12) 
      {
	tempStep = -slopeOld / (2.0 * b);
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
        tempStep = -slopeOld/(b + sqrt(disc));
      
      /* Also removed at the request of Homer Walker
      if (tempStep > 0.5 * step) 
	tempStep = 0.5 * step ;
      */

    }
    
    previousStep = step ;
    prevf = meritFuncNew ; 

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

    message = "(USING RECOVERY STEP!)";
  }

  if (suffDecrCond == AredPred)
    print.printStep(nIters, step, normFOld, normFNew, message, false);
  else
    print.printStep(nIters, step, meritFuncOld, meritFuncNew, message);

  counter.setValues(*paramsPtr);
  return (!isFailed);
}

bool NOX::LineSearch::PolynomialWalker::isSufficientDecrease(double step, 
							     double eta) const
{
  double rhs = 0.0;
  if (suffDecrCond == ArmijoGoldstein) 
  {
    rhs = meritFuncOld + alpha * step * slopeOld;
    return (meritFuncNew <= rhs);
  }
  else if (suffDecrCond == AredPred) 
  {
    rhs = normFOld * (1.0 - alpha * (1.0 - eta));
    return (normFNew <= rhs);
  }
  else if (suffDecrCond == None)
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

  // Compute ||F(0)|| ONLY for Ared/Pred condition
  if (suffDecrCond == AredPred) {
    if (userNormPtr != 0)
      normFNew = userNormPtr->norm(newGrp.getF());
    else
      normFNew = newGrp.getNormF();
  }

  // Compute merit function f(0)
  if (meritFuncPtr != 0) {
    meritFuncNew = meritFuncPtr->computef(newGrp);
  }
  else {
    meritFuncNew = 0.5 * newGrp.getNormF() * newGrp.getNormF();
  }

  return true;
}

bool NOX::LineSearch::PolynomialWalker::isIncreaseAllowed(int nOuterIters) const
{
  double increase = 0.0;
  if (suffDecrCond == AredPred)  
    increase = sqrt(normFNew / normFOld);
  else
    increase = sqrt(meritFuncNew / meritFuncOld);

  return ( (increase <= relIncrease) && (nOuterIters <= numAllowed) );
}

#endif
