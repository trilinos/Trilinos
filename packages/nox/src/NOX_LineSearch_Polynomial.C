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

using namespace NOX;
using namespace NOX::LineSearch;

Polynomial::Polynomial(Parameter::List& params) :
  inputList(0)
{
  reset(params);
}

Polynomial::~Polynomial()
{

}

bool Polynomial::reset(Parameter::List& params)
{ 
  if (params.isParameterString("Convergence Criteria")) {
    string choice = params.getParameter("Convergence Criteria", "Ared/Pred");
    if (choice == "Ared/Pred") {
      convCriteria = AredPred;
    }
    else if (choice == "Armijo-Goldstein") {
      convCriteria = ArmijoGoldstein;
    }
    else {
      cout << "ERROR: NOX::LineSearch::Quadratic::reset() - the choice of "
	   << "\"Convergence Criteria\" parameter is invalid." << endl;
      throw "NOX Error";
    }
  }
  else {
    // default to "Armijo-Goldstein" if "Convergence Criteria" is not set
    convCriteria = ArmijoGoldstein;
  } 

  minstep = params.getParameter("Minimum Step", 1.0e-12);
  defaultstep = params.getParameter("Default Step", 1.0);
  recoverystep = params.getParameter("Recovery Step", defaultstep);
  maxiters = params.getParameter("Max Iters", 100);
  alpha = params.getParameter("Alpha Factor", 1.0e-4);
  minBoundFactor = params.getParameter("Min Bounds Factor", 0.1);
  maxBoundFactor = params.getParameter("Max Bounds Factor", 0.9);
  inputList = &params;
  totalNumLineSearchCalls = 0;
  totalNumNonTrivialLineSearches = 0;
  totalNumFailedLineSearches = 0;
  totalNumIterations = 0;
  return true;
}

bool Polynomial::compute(Abstract::Group& newgrp, double& step, 
			 const Abstract::Vector& dir,
			 const Solver::Generic& s) 
{
  totalNumLineSearchCalls += 1;

  const Abstract::Group& oldgrp = s.getPreviousSolutionGroup();
  
  double oldf = 0.5*oldgrp.getNormF()*oldgrp.getNormF();  
                            // Redefined f(), RH

  // General computation of directional derivative used in curvature condition
  // Note that for Newton direction, oldfprime = -2.0*oldf
  Abstract::Vector* tmpvecptr = oldgrp.getX().clone(ShapeCopy);
  oldgrp.applyJacobian(dir,*tmpvecptr);
  double oldfprime = tmpvecptr->dot(oldgrp.getF());
  delete tmpvecptr;

  double newf, prevf;
  double tempStep, previousStep;
  double a,b,term1,term2,disc ;
  bool isfailed = false;
  bool firstPass = true ;

  step = defaultstep;
  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeF();    
  newf = 0.5*newgrp.getNormF()*newgrp.getNormF();  
                            // Redefined f(), RH

  int niters = 1;

  if (Utils::doPrint(Utils::InnerIteration)) {
   cout << "\n" << Utils::fill(72) << "\n" << "-- Polynomial Line Search -- \n";
  }

  // Get the linear solve tolerance if doing ared/pred for conv criteria
  const NOX::Parameter::List& p = s.getParameterList();
  double eta_original = 0.0;
  double eta = 0.0;
  if (convCriteria == AredPred) {
    eta_original = p.sublist("Direction").sublist("Linear Solver").getParameter("Tolerance", -1.0);
    eta = eta_original;
  }

  // Compute the convergence criteria for the line search 
  double convergence = 0.0;
  if (convCriteria == ArmijoGoldstein) {
    convergence = oldf + alpha*step*oldfprime;
  }
  else if (convCriteria == AredPred) {
    convergence = oldf*(1.0-alpha*(1.0-eta));
  }
  
  // Increment the number of newton steps requiring a line search
  if (newf >= convergence)
    totalNumNonTrivialLineSearches += 1;

  while (newf >= convergence) {

    if (Utils::doPrint(Utils::InnerIteration)) {
      cout << setw(3) << niters << ":";
      cout << " step = " << Utils::sci(step);
      cout << " oldf = " << Utils::sci(sqrt(2.*oldf));
      cout << " newf = " << Utils::sci(sqrt(2.*newf));
      cout << endl;
    }

    // update the iteration count 
    totalNumIterations += 1;
    
    if( firstPass == true) {

      /*   First try quadratic  */
      
      tempStep = -oldfprime/(2.0*(newf - oldf - oldfprime)) ;
      firstPass = false;
    }
    else {

      /*   Do cubic as many times as needed */
      
      term1 = newf - oldf - step*oldfprime ;
      term2 = prevf - oldf - previousStep*oldfprime ;
      a = 1.0/(step-previousStep)*( term1/step/step -
                                    term2/previousStep/previousStep) ;
      b = 1.0/(step-previousStep)*( -term1*previousStep/step/step +
                                    term2*step/previousStep/previousStep) ;
      disc = b*b - 3.0*a*oldfprime ;
      if(disc < 0) {
	step = recoverystep;
	newgrp.computeX(oldgrp, dir, step);
	newgrp.computeF();    
	newf = 0.5*newgrp.getNormF()*newgrp.getNormF(); 
	if (Utils::doPrint(Utils::InnerIteration)) { 
	  cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
	  cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(sqrt(2.*oldf));
	  cout << Utils::fill(1,' ') << "newf = " << Utils::sci(sqrt(2.*newf));
	  cout << " (USING RECOVERY STEP!)" << endl;
	  cout << Utils::fill(72) << "\n" << endl;
	}
	isfailed = true;
	setOutputParameters();
	return(!isfailed);
      }
      if( fabs(a) < 1.e-12)
	tempStep = -oldfprime/2.0/b ;
      else
	tempStep = (-b + sqrt(disc))/3.0/a ;
      
      if(tempStep > 0.5*step) tempStep = 0.5*step ;
      
    }
    
    previousStep = step ;
    prevf = newf ; 
    
    if (tempStep < minBoundFactor*step) 
      step *= minBoundFactor;
    else if (tempStep > maxBoundFactor*step) 
      step *= maxBoundFactor;
    else step = tempStep;
    
    // Safeguard while loop termination for by adjusting eta during "Ared/Pred".
    // The "Direction" also needs to use this eta in the next computation
    // if using inexact Newton with "Type 1" or Type 2" eta computation.
    eta = 1.0-tempStep*(1.0-eta);
    
    if ((step < minstep) || (niters > maxiters))
      {
	totalNumFailedLineSearches += 1;
	step = recoverystep;
	
	// For directions that require an iterative linear solve:
	// If we are not using a constant tolerance in the linear solver 
	// (i.e. we are using Homer Walker's "Type 1" or "Type 2" criteria) 
	// then we must adjust the old tolerance to be used in the next 
	// computation of eta in the direction object based on these 
	// linesearch results
	eta = 1.0-step*(1.0-eta_original);
	inputList->setParameter("Adjusted Tolerance", eta);
	
	// Update the group, and exit compute()
	newgrp.computeX(oldgrp, dir, step);
	newgrp.computeF();    
	newf = 0.5*newgrp.getNormF()*newgrp.getNormF(); 
	if (Utils::doPrint(Utils::InnerIteration)) {
	  cout << Utils::fill(5,' ') << "step = " << Utils::sci(step);
	  cout << Utils::fill(1,' ') << "oldf = " << Utils::sci(sqrt(2.*oldf));
	  cout << Utils::fill(1,' ') << "newf = " << Utils::sci(sqrt(2.*newf));
	  cout << " (USING RECOVERY STEP!)" << endl;
	  cout << Utils::fill(72) << "\n" << endl;
	}
	isfailed = true;
	setOutputParameters();
	return(!isfailed);
      }
    
    niters ++;
    newgrp.computeX(oldgrp, dir, step);
    newgrp.computeF();    
    newf = 0.5*newgrp.getNormF()*newgrp.getNormF(); 
    
    // Compute the convergence criteria for the line search 
    if (convCriteria == AredPred) {
      convergence = oldf*(1.0-alpha*(1.0-eta));
    } 
    else {
      // defaults to Armijo-Goldstein
      convergence = oldf + alpha*step*oldfprime;
    }
    
  } // End while loop 
  
  // For directions that require an iterative linear solve:
  // If we are not using a constant tolerance in the linear solver 
  // (i.e. we are using Homer Walker's "Type 1" or "Type 2" criteria) 
  // then we must adjust the old tolerance to be used in the next 
  // computation of eta in the direction object based on these 
  // linesearch results
  inputList->setParameter("Adjusted Tolerance", eta);
  
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

bool Polynomial::setOutputParameters() {
  NOX::Parameter::List& outputList = inputList->sublist("Output");
  outputList.setParameter("Total Number of Line Search Calls", totalNumLineSearchCalls);
  outputList.setParameter("Total Number of Non-trivial Line Searches", totalNumNonTrivialLineSearches);
  outputList.setParameter("Total Number of Failed Line Searches", totalNumFailedLineSearches);
  outputList.setParameter("Total Number of Line Search Inner Iterations", totalNumIterations);
  return true;
}
