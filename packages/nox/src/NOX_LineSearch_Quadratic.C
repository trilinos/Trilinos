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

#include "NOX_LineSearch_Quadratic.H"   // class definition

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::LineSearch;

Quadratic::Quadratic(Parameter::List& params) :
  inputList(0)
{
  reset(params);
}

Quadratic::~Quadratic()
{

}

bool Quadratic::reset(Parameter::List& params)
{ 
  if (params.isParameter("Convergence Criteria")) {
    string choice = params.getParameter("Convergence Criteria", "AredPred");
    if (choice == "AredPred") {
      convCriteria = AredPred;
    }
    else if (choice == "Armijo-Goldstein") {
      convCriteria = ArmijoGoldstein;
    }
    else {
      cout << "ERROR: NOX::LineSearch::Quadratic::reset() - the choice of "
	   << "\"Convergence Criteria\" is invalid." << endl;
      throw "NOX Error";
    }
  }
  else {
    // default to "Ared/Pred"
    convCriteria = AredPred;
  } 

  minStep = params.getParameter("Minimum Step", 1.0e-12);
  defaultStep = params.getParameter("Default Step", 1.0);
  recoveryStep = params.getParameter("Recovery Step", defaultStep);
  maxIters = params.getParameter("Max Iters", 100);
  alpha = params.getParameter("Alpha Factor", 1.0e-4);
  boundFactor = params.getParameter("Bounds Factor", 0.1);
  inputList = &params;
  totalNumIterations = 0;
  totalNumFailedLineSearches = 0;
  return true;
}

bool Quadratic::compute(Abstract::Group& newgrp, double& step, 
			const Abstract::Vector& dir,
			const Solver::Generic& s) 
{

  const Abstract::Group& oldgrp = s.getPreviousSolutionGroup();
  double oldf = 0.5*oldgrp.getNormF()*oldgrp.getNormF(); 

  // General computation of directional derivative used in curvature condition
  // Note that for Newton direction, oldfprime = -2.0*oldf
  Abstract::Vector* tmpvecptr = oldgrp.getX().clone(ShapeCopy);
  oldgrp.applyJacobian(dir,*tmpvecptr);
  double oldfprime = tmpvecptr->dot(oldgrp.getF());
  delete tmpvecptr;

  double newf, prevf;
  double tempStep, previousStep;
  bool isfailed = false;

  step = defaultStep;
  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeF();    

  // f = 0.5 * 2-Norm(F) * 2-Norm(F)
  newf = 0.5*newgrp.getNormF()*newgrp.getNormF();  

  int niters = 1;

  if (Utils::doPrint(Utils::InnerIteration)) {
   cout << "\n" << Utils::fill(72) << "\n" << "-- Quadratic Line Search -- \n";
  }
  
  // Get the Linear solver tolerance used in the Newton 
  const NOX::Parameter::List& const_p = s.getParameterList();
  NOX::Parameter::List& p = const_cast <NOX::Parameter::List&>(const_p);
  double eta_original = p.sublist("Linear Solver").getParameter("Tolerance", 0.0);
  double eta = eta_original;

  // Compute the convergence criteria for the line search 
  double convergence = 0.0;
  if (convCriteria == ArmijoGoldstein) {
    convergence = oldf + alpha*step*oldfprime;
  }
  else if (convCriteria == AredPred) {
    convergence = oldf*(1.0-alpha*(1-eta));
  }
  
  // Iterate until convergence criteria is satisfied
  while (newf >= convergence) {  

    if (Utils::doPrint(Utils::InnerIteration)) {
      cout << setw(3) << niters << ":";
      cout << " step = " << Utils::sci(step);
      cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << " newf = " << Utils::sci(sqrt(2.*newf));
      cout << endl;
    }

    // Compute a new step length
    tempStep = -oldfprime/(2.0*(newf - oldf - oldfprime));

    
    //   Enforce bounds on minimum step size
    if(tempStep < boundFactor) 
      tempStep = boundFactor;

    // Safeguard while loop termination by adjusting eta
    // the direction also needs to use this eta in the next computation
    // if using "Type 1" or Type 2" comutation.
    eta = 1.0-tempStep*(1.0-eta);
    tempStep *= step;
    previousStep = step;
    prevf = newf; 
    step = tempStep;

    // update the iteration count 
    totalNumIterations += 1;
    
    // Check for linesearch failure: if true, exit the method
    if ((step < minStep) || (niters > maxIters))
    {
      totalNumFailedLineSearches += 1;
      step = recoveryStep;

      // If we are not using a constant tolerance in the linear solver 
      // (i.e. we are using either "Type 1" or "Type 2") then we must adjust
      // the old tolerance to be used in the next computation of eta in the 
      // direction object based on these linesearch results
      string forcingTermMethod = 
	p.sublist("Direction").getParameter("Forcing Term Method", "");
      if ((forcingTermMethod == "Type 1") ||
	  (forcingTermMethod == "Type 2")) {
	eta = 1.0-step*(1.0-eta_original);
	p.sublist("Linear Solver").setParameter("Tolerance", eta);
      }

      newgrp.computeX(oldgrp, dir, step);
      newgrp.computeF();    
      newf = 0.5*newgrp.getNormF()*newgrp.getNormF();
      cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
      cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << Utils::fill(1,' ') << "newf = " << Utils::sci(sqrt(2.*newf));
      cout << " (USING RECOVERY STEP!)" << endl;
      cout << Utils::fill(72) << "\n" << endl;
      isfailed = true;
      setOutputParameters();
      return(!isfailed);
    }
    
    niters ++;
    newgrp.computeX(oldgrp, dir, step);
    newgrp.computeF();    
    newf = 0.5*newgrp.getNormF()*newgrp.getNormF();

    // Compute the convergence criteria for the line search 
    double convergence = 0.0;
    if (convCriteria == ArmijoGoldstein) {
      convergence = oldf + alpha*step*oldfprime;
    }
    else if (convCriteria == AredPred) {
      convergence = oldf*(1.0-alpha*(1-eta));
    }

  } // end while loop

  // If we are not using a constant tolerance in the linear solver 
  // (i.e. we are using either "Type 1" or "Type 2") then we must adjust
  // the old tolerance to be used in the next computation of eta in the 
  // direction object based on these linesearch results
  string forcingTermMethod = 
    p.sublist("Direction").getParameter("Forcing Term Method", "");
  if ((forcingTermMethod == "Type 1") ||
      (forcingTermMethod == "Type 2")) {
    eta = 1.0-step*(1.0-eta_original);
    p.sublist("Linear Solver").setParameter("Tolerance", eta);
  }

  if (Utils::doPrint(Utils::InnerIteration)) {
      cout << setw(3) << niters << ":";
      cout << " step = " << Utils::sci(step);
      cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << " newf = " << Utils::sci(sqrt(2.*newf));
      cout << " (STEP ACCEPTED!)" << endl;
      cout << Utils::fill(72) << "\n" << endl;
  }

  setOutputParameters();
  return (!isfailed);
}

bool Quadratic::setOutputParameters() {
  NOX::Parameter::List& outputList = inputList->sublist("Output");
  outputList.setParameter("Total Number of Line Search Iterations", totalNumIterations);
  outputList.setParameter("Total Number of Failed Line Searches", totalNumFailedLineSearches);
}
