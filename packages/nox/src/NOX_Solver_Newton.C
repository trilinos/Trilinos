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

#include "NOX_Solver_Newton.H"	// class definition

using namespace NOX;
using namespace NOX::Solver;

Newton::Newton(Abstract::Group& xgrp, Status::Test& t, Parameter::List& p) :
  solnptr(&xgrp),		// pointer to xgrp
  oldsolnptr(xgrp.clone(DeepCopy)), // create via clone
  oldsoln(*oldsolnptr),		// reference to just-created pointer
  dirptr(xgrp.getX().clone(CopyShape)), // create via clone 
  dir(*dirptr),			// reference to just-created pointer
  testptr(&t),			// pointer to t
  iparamsptr(&p),		// pointer to p
  oparams(),			// empty list
  linesearch(p.sublist("Line Search")), // initialize line search
  step(0.0),			// initialize to zero
  niter(0),			// initialize to zero
  status(Status::Unconverged)	// initialize convergence status
{
  // Print out initialization information
  if (Utils::doPrint(Utils::Parameters)) {

    cout << "\n" << Utils::fill(72) << "\n";
    cout << "\n-- Parameters Passed to Nonlinear Solver on Print Processor --\n\n";
    iparamsptr->print(cout,5);
    cout << "\n-- Status Tests Passed to Nonlinear Solver on Print Processor --\n\n";
    testptr->print(cout, 5);
    cout <<"\n" << Utils::fill(72) << "\n";

  }

  // Compute RHS of initital guess
  solnptr->computeRHS();
}

Newton::~Newton() 
{
  delete oldsolnptr;
  delete dirptr;
}

void Newton::resetInputParameters(Parameter::List& p)
{
}

Status::StatusType Newton::getStatus()
{
  status = testptr->operator()(*this);
  return status;
}

Status::StatusType Newton::iterate()
{
  // Copy pointers into temporary references
  Abstract::Group& soln = *solnptr;
  Status::Test& test = *testptr;

  // Reset forcing term.
  resetForcingTerm();

  // Compute Jacobian at current solution.
  soln.computeJacobian();

  // Compute Newton direction for current solution.
  /* NOTE FROM TAMMY: Need to check the return status! */
  soln.computeNewton(iparamsptr->sublist("Linear Solver"));

  // Set search direction.
  dir = soln.getNewton();

  // Copy current soln to the old soln.
  oldsoln = soln;

  // Do line search and compute new soln.
  /* NOTE FROM TAMMY: Need to check the return status! */
  linesearch(soln, step, oldsoln, dir);

  // Compute RHS for new current solution.
  soln.computeRHS();

  // Update iteration count.
  niter ++;

  // Evaluate the current status.
  status = test(*this);
 
  // Return status.
  return status;
}

Status::StatusType Newton::solve()
{
  status = testptr->operator()(*this);
  printUpdate();

  // Iterate until converged or failed
  while (status == Status::Unconverged) {
    status = iterate();
    printUpdate();
  }

  return status;
}

const Abstract::Group& Newton::getSolutionGroup() const
{
  return *solnptr;
}

const Abstract::Group& Newton::getPreviousSolutionGroup() const
{
  return oldsoln;
}

int Newton::getNumIterations() const
{
  return niter;
}

const Parameter::List& Newton::getOutputParameters() const
{
  oparams.setParameter("Nonlinear Iterations", niter);
  oparams.setParameter("2-Norm of Residual", solnptr->getNormRHS());
  return oparams;
}

// protected
void Newton::printUpdate() 
{
  double norm_k;
  double norm_newton;

  // All processes participate in the computation of these norms...
  if (Utils::doAllPrint(Utils::OuterIteration)) {
    norm_k = solnptr->getNormRHS();
    norm_newton = (niter > 0) ? oldsoln.getNewton().norm() : 0;
  }

  // ...But only the print process actually prints the result.
  if (Utils::doPrint(Utils::OuterIteration)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "-- Newton Step " << niter << " -- \n";
    cout << "Residual Norm = " << Utils::sci(norm_k);
    cout << "  Step = " << Utils::sci(step);
    cout << "  Update Norm = " << Utils::sci(norm_newton);
    if (status > 0)
      cout << " (Converged!)";
    if (status < 0)
      cout << " (Failed!)";
    cout << "\n" << Utils::fill(72) << "\n" << endl;
  }
  
  if ((status != 0) && (Utils::doPrint(Utils::OuterIteration))) {
    cout << Utils::fill(72) << "\n";
    cout << "-- Final Status Test Results --\n";    
    testptr->print(cout);
    cout << Utils::fill(72) << "\n";
  }
}

// protected
void Newton::resetForcingTerm()
{
  // Reset the forcing term at the beginning on a nonlinear iteration,
  // based on the last iteration.

  if (niter == 0)
    return;

  if ((!iparamsptr->isParameter("Forcing Term Method")) ||
      (iparamsptr->isParameterEqual("Forcing Term Method", "None")))
    return;

  // Get forcing term parameters.
  const string method = iparamsptr->getParameter("Forcing Term Method", "");
  const double eta_min = iparamsptr->getParameter("Forcing Term Minimum Tolerance", 1.0e-6);
  const double eta_max = iparamsptr->getParameter("Forcing Term Maximum Tolerance", 0.01);

  // Get linear solver parameter list and current tolerance.
  Parameter::List& lsparams = iparamsptr->sublist("Linear Solver");
  const double eta_km1 = lsparams.getParameter("Tolerance", 0.0);

  // New forcing term.
  double eta_k;

  string indent = "       ";

  if (Utils::doPrint(Utils::Details)) {
    cout << indent << "CALCULATING FORCING TERM" << endl;
    cout << indent << "Method: " << method << endl;
  }

  if (method == "Constant") {    

    eta_k = eta_min;

  }        

  else if (method == "Type 1") {
    
    // Create a new vector to be the predicted RHS
    static Abstract::Vector* predrhs = oldsoln.getRHS().clone(CopyShape);

    // Compute predrhs = Jacobian * dir
    oldsoln.applyJacobian(dir, *predrhs);
    
    // Compute predrhs = RHSVector + step * predrhs (this is the predicted RHS)
    predrhs->update(1.0, oldsoln.getRHS(), step);

    // Return norm of predicted RHS
    const double normpredrhs = predrhs->norm();

    // Get other norms
    const double normrhs = solnptr->getNormRHS();
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
    
  else if (method == "Type 2") {  
    
    const double normrhs = solnptr->getNormRHS();
    const double normoldrhs = oldsoln.getNormRHS();
    const double alpha = iparamsptr->getParameter("Forcing Term Alpha", 1.5);
    const double gamma = iparamsptr->getParameter("Forcing Term Gamma", 0.9);
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

  else {

    cerr << "*** Warning: Invalid Forcing Term Method ***" << endl;
    return;
  }

  // Reset linear solver tolerance
  lsparams.setParameter("Tolerance", eta_k);

  if (Utils::doPrint(Utils::Details)) {
    cout << indent << "Forcing Term: " << eta_k << endl;
  }

}




