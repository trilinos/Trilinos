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

Newton::Newton(Parameter::List& p) 
{
  predrhs = NULL;
  stepdir = NULL;
  reset(p);
}

Newton::~Newton()
{
  delete predrhs;
}

bool Newton::reset(Parameter::List& p)
{
  paramsptr = &p;
  if (!paramsptr->sublist("Linear Solver").isParameter("Tolerance"))
    paramsptr->sublist("Linear Solver").setParameter("Tolerance", 1.0e-10);
  return true;
}

bool Newton::operator()(Abstract::Vector& dir, 
			Abstract::Group& soln, 
			const Solver::Generic& solver)
{
  // Compute RHS at current solution
  bool ok = soln.computeRHS();

  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Newton::operator() - Unable to compute RHS." << endl;
    return false;
  }

  // Reset the linear solver tolerance
  ok = resetForcingTerm(soln, solver.getPreviousSolutionGroup(), solver.getNumIterations());
  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Newton::operator() - Unable to set Forcing term." << endl;
    return false;
  }

  // Compute Jacobian at current solution.
  ok = soln.computeJacobian();

  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Newton::operator() - Unable to compute Jacobian." << endl;
    return false;
  }
  
  // Compute the Newton direction
  ok = soln.computeNewton(paramsptr->sublist("Linear Solver"));


  // It didn't work, but maybe it's ok anyway...
  if (!ok) {

    if (predrhs == NULL) {
      predrhs = soln.getRHS().clone(CopyShape);
    }

    soln.applyJacobian(soln.getNewton(), *predrhs);    
    predrhs->update(-1.0, soln.getRHS(), 1.0);
    double accuracy = predrhs->norm();
    if (accuracy < 1) {
      ok = true;
      if (Utils::doPrint(Utils::Warning)) 
	cout << "WARNING: NOX::Direction::Newton::operator() - Newton solve failure.\n" 
	     << "Desired accuracy is " 
	     << Utils::sci(paramsptr->sublist("Linear Solver").getParameter("Tolerance", 1.0e-10)) << ".\n"
	     << "Using solution with accuracy of " << Utils::sci(accuracy) << "." << endl;
    }
  }

  if (!ok) {
    if (Utils::doPrint(Utils::Warning))
      cout << "NOX::Direction::Newton::operator() - Unable to compute Newton direction." << endl;
    return false;
  }
  
  // Set search direction.
  dir = soln.getNewton();

  return ok;
}


// protected
bool Newton::resetForcingTerm(const Abstract::Group& soln, const Abstract::Group& oldsoln, int niter)
{
  // Reset the forcing term at the beginning on a nonlinear iteration,
  // based on the last iteration.

  if ((!paramsptr->isParameter("Forcing Term Method")) ||
      (paramsptr->isParameterEqual("Forcing Term Method", "None")))
    return true;

  // Get forcing term parameters.
  const string method = paramsptr->getParameter("Forcing Term Method", "");
  const double eta_min = paramsptr->getParameter("Forcing Term Minimum Tolerance", 1.0e-6);
  const double eta_max = paramsptr->getParameter("Forcing Term Maximum Tolerance", 0.01);

  // Get linear solver parameter list and current tolerance.
  const double eta_km1 = paramsptr->sublist("Linear Solver").getParameter("Tolerance", 0.0);

  // New forcing term.
  double eta_k;

  const string indent = "       ";

  if (Utils::doPrint(Utils::Details)) {
    cout << indent << "CALCULATING FORCING TERM" << endl;
    cout << indent << "Method: " << method << endl;
  }

  if (method == "Constant") {    

    eta_k = eta_min;

  }        

  else if (method == "Type 1") {
    
    if (niter == 0) {

      eta_k = paramsptr->getParameter("Forcing Term Initial Tolerance", 0.01);

    }
    else {


      // Create a new vector to be the predicted RHS
      if (predrhs == NULL) {
	predrhs = oldsoln.getRHS().clone(CopyShape);
      }
      if (stepdir == NULL) {
	stepdir = oldsoln.getRHS().clone(CopyShape);
      }
      
      // stepdir = X - oldX (i.e., the step times the direction)
      stepdir->update(1.0, soln.getX(), -1.0, oldsoln.getX(), 0);
      
      // Compute predrhs = Jacobian * step * dir
      oldsoln.applyJacobian(*stepdir, *predrhs);
      
      // Compute predrhs = RHSVector + predrhs (this is the predicted RHS)
      predrhs->update(1.0, oldsoln.getRHS(), 1.0);
      
      // Return norm of predicted RHS
      const double normpredrhs = predrhs->norm();
      
      // Get other norms
      const double normrhs = soln.getNormRHS();
      const double normoldrhs = oldsoln.getNormRHS();
      
      // Compute forcing term
      eta_k = fabs(normrhs - normpredrhs) / normoldrhs;
      
      // Some output
      if (Utils::doPrint(Utils::Details)) {
	cout << indent << "Residual Norm k-1 =             " << normoldrhs << "\n";
	cout << indent << "Residual Norm Linear Model k =  " << normpredrhs << "\n";
	cout << indent << "Residual Norm k =               " << normrhs << "\n";
	cout << indent << "Calculated eta_k (pre-bounds) = " << eta_k << endl;
      }
      
      // Impose safeguard and constraints ...
      const double alpha = (1.0 + sqrt(5.0)) / 2.0;
      const double eta_k_alpha = pow(eta_km1, alpha);
      eta_k = max(eta_k, eta_k_alpha);
      eta_k = max(eta_k, eta_min);
      eta_k = min(eta_max, eta_k);
    }
  }
    
  else if (method == "Type 2") {  
    
    if (niter == 0) {

      eta_k = paramsptr->getParameter("Forcing Term Initial Tolerance", 0.01);

    }
    else {

      const double normrhs = soln.getNormRHS();
      const double normoldrhs = oldsoln.getNormRHS();
      const double alpha = paramsptr->getParameter("Forcing Term Alpha", 1.5);
      const double gamma = paramsptr->getParameter("Forcing Term Gamma", 0.9);
      const double residual_ratio = normrhs / normoldrhs;
      
      eta_k = gamma * pow(residual_ratio, alpha);
      
      // Some output
      if (Utils::doPrint(Utils::Details)) {
	cout << indent << "Residual Norm k-1 =             " << normoldrhs << "\n";
	cout << indent << "Residual Norm k =               " << normrhs << "\n";
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
      cout << "NOX::Direction::Newton::resetForcingTerm - invalid forcing term method (" << method << ")" << endl;

    return false;
  }

  // Reset linear solver tolerance
  paramsptr->sublist("Linear Solver").setParameter("Tolerance", eta_k);

  if (Utils::doPrint(Utils::Details)) 
    cout << indent << "Forcing Term: " << eta_k << endl;

  return true;
}





