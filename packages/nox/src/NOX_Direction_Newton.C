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

#include "NOX_Direction_Newton.H" // class definition
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::Direction;

Newton::Newton(Parameter::List& p) :
  predrhs(0),
  stepdir(0)
  
{
  reset(p);
}

Newton::~Newton()
{
  delete predrhs;
  delete stepdir;
}

bool Newton::reset(Parameter::List& p)
{
  paramsPtr = &p;
  if (!paramsPtr->sublist("Linear Solver").isParameter("Tolerance"))
    paramsPtr->sublist("Linear Solver").setParameter("Tolerance", 1.0e-10);
  return true;
}

bool Newton::compute(Abstract::Vector& dir, 
		     Abstract::Group& soln, 
		     const Solver::Generic& solver)
{
  // Compute F at current solution
  bool ok = soln.computeF();
  double normF = soln.getNormF();

  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Newton::compute - Unable to compute F." << endl;
    return false;
  }

  // Reset the linear solver tolerance
  ok = resetForcingTerm(soln, solver.getPreviousSolutionGroup(), solver.getNumIterations(), solver.getParameterList());
  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Newton::compute - Unable to set Forcing term." << endl;
    return false;
  }

  // Compute Jacobian at current solution.
  ok = soln.computeJacobian();

  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Newton::compute - Unable to compute Jacobian." << endl;
    return false;
  }
  
  // Compute the Newton direction
  ok = soln.computeNewton(paramsPtr->sublist("Linear Solver"));

  // It didn't work, but maybe it's ok anyway...
  if (!ok) {

    double accuracy = soln.getNormNewtonSolveResidual();

    if (accuracy < 0) {
      cerr << "NOX::Direction::Newton::compute " 
	   << "- getNormNewtonSolveResidual returned a negative value" << endl;
    }
      
    // Check if there is any improvement in the relative residual
    if (accuracy < normF) {
      ok = true;
      double tolerance = paramsPtr->sublist("Linear Solver").getParameter("Tolerance", 1.0e-10);
      if (Utils::doPrint(Utils::Warning)) 
	cout << "WARNING: NOX::Direction::Newton::compute - Newton solve failure.\n" 
	     << "Desired accuracy is " << Utils::sci(tolerance) << ".\n"
	     << "Using solution with accuracy of " << Utils::sci(accuracy) << "." << endl;
    }
  }

  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Newton::compute - Unable to compute Newton direction." << endl;
    return false;
  }
  
  // Set search direction.
  dir = soln.getNewton();

  return ok;
}


// protected
bool Newton::resetForcingTerm(const Abstract::Group& soln, 
			      const Abstract::Group& oldsoln, 
			      int niter,
			      const Parameter::List& solverParams)
{
  // Reset the forcing term at the beginning on a nonlinear iteration,
  // based on the last iteration.

  if ((!paramsPtr->isParameter("Forcing Term Method")) ||
      (paramsPtr->isParameterEqual("Forcing Term Method", "None")))
    return true;

  // Get forcing term parameters.
  const string method = paramsPtr->getParameter("Forcing Term Method", "Constant");
  const double eta_min = paramsPtr->getParameter("Forcing Term Minimum Tolerance", 1.0e-4);
  const double eta_max = paramsPtr->getParameter("Forcing Term Maximum Tolerance", 0.9);

  // Get linear solver current tolerance.

  double eta_km1 = 0.0;
  if (((method == "Type 1") || (method == "Type 2")) 
	&& (solverParams.sublist("Line Search").isParameterDouble("Adjusted Tolerance"))) {
    
    // Tolerance may have been adjusted in a line search algorithm   
    eta_km1 = solverParams.sublist("Line Search")
      .getParameter("Adjusted Tolerance", 0.0);
    
  }
  else {
    // Default to the old tolerance
    eta_km1 = paramsPtr->sublist("Linear Solver")
      .getParameter("Tolerance", 0.0);
    
  }

  // New forcing term.
  double eta_k;

  const string indent = "       ";

  if (Utils::doPrint(Utils::Details)) {
    cout << indent << "CALCULATING FORCING TERM" << endl;
    cout << indent << "Method: " << method << endl;
  }

  if (method == "Constant") {    

    if (eta_km1 != 0.0) 
      eta_k = eta_km1;
    else
      eta_k = eta_min;

  }        

  else if (method == "Type 1") {
    
    if (niter == 0) {

      eta_k = paramsPtr->getParameter("Forcing Term Initial Tolerance", 0.01);

    }
    else {

      // Return norm of predicted F

      // do NOT use the following line!! This does NOT account for 
      // line search step length taken.
      //const double normpredf = oldsoln.getNormNewtonSolveResidual();
      
      // Create a new vector to be the predicted RHS
      if (predrhs == NULL) {
	predrhs = oldsoln.getF().clone(ShapeCopy);
      }
      if (stepdir == NULL) {
	stepdir = oldsoln.getF().clone(ShapeCopy);
      }
      
      // stepdir = X - oldX (i.e., the step times the direction)
      stepdir->update(1.0, soln.getX(), -1.0, oldsoln.getX(), 0);
      
      // Compute predrhs = Jacobian * step * dir
      oldsoln.applyJacobian(*stepdir, *predrhs);
      
      // Compute predrhs = RHSVector + predrhs (this is the predicted RHS)
      predrhs->update(1.0, oldsoln.getF(), 1.0);
      
      // Return norm of predicted RHS
      const double normpredf = predrhs->norm();
      
      if (normpredf < 0) {
	cerr << "NOX::Direction::Newton::resetForcingTerm " 
	     << "- getNormNewtonSolveResidual returned a negative value" << endl;
      }

      // Get other norms
      const double normf = soln.getNormF();
      const double normoldf = oldsoln.getNormF();
      
      // Compute forcing term
      eta_k = fabs(normf - normpredf) / normoldf;
      
      // Some output
      if (Utils::doPrint(Utils::Details)) {
	cout << indent << "Residual Norm k-1 =             " << normoldf << "\n";
	cout << indent << "Residual Norm Linear Model k =  " << normpredf << "\n";
	cout << indent << "Residual Norm k =               " << normf << "\n";
	cout << indent << "Calculated eta_k (pre-bounds) = " << eta_k << endl;
      }
      
      // Impose safeguard and constraints ...
      const double alpha = (1.0 + sqrt(5.0)) / 2.0;
      const double eta_km1_alpha = pow(eta_km1, alpha);
      eta_k = max(eta_k, eta_km1_alpha);
      eta_k = max(eta_k, eta_min);
      eta_k = min(eta_max, eta_k);
    }
  }
    
  else if (method == "Type 2") {  
    
    if (niter == 0) {

      eta_k = paramsPtr->getParameter("Forcing Term Initial Tolerance", 0.01);

    }
    else {

      const double normf = soln.getNormF();
      const double normoldf = oldsoln.getNormF();
      const double alpha = paramsPtr->getParameter("Forcing Term Alpha", 1.5);
      const double gamma = paramsPtr->getParameter("Forcing Term Gamma", 0.9);
      const double residual_ratio = normf / normoldf;
      
      eta_k = gamma * pow(residual_ratio, alpha);
      
      // Some output
      if (Utils::doPrint(Utils::Details)) {
	cout << indent << "Residual Norm k-1 =             " << normoldf << "\n";
	cout << indent << "Residual Norm k =               " << normf << "\n";
	cout << indent << "Calculated eta_k (pre-bounds) = " << eta_k << endl;
      }
      
      // Impose safeguard and constraints ... 
      const double eta_k_alpha = gamma * pow(eta_km1, alpha);
      eta_k = max(eta_k, eta_k_alpha);
      eta_k = max(eta_k, eta_min);
      eta_k = min(eta_max, eta_k);
    }

  }

  else {

    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Newton::resetForcingTerm "
	   << "- invalid forcing term method (" << method << ")" << endl;

    return false;
  }

  // Reset linear solver tolerance
  paramsPtr->sublist("Linear Solver").setParameter("Tolerance", eta_k);

  if (Utils::doPrint(Utils::Details)) 
    cout << indent << "Forcing Term: " << eta_k << endl;

  return true;
}





