// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NOX_Solver_Newton.H"	// class definition

#include <iomanip>		// for setw
#include <cmath>		// for min, max
#include "NOX_Utils.H"		

using namespace NOX;
using namespace NOX::Solver;

Newton::Newton(Abstract::Group& xgrp, Status::Test& t, Parameter::List& p) :
  soln(xgrp),			// reference to xgrp
  oldsolnptr(xgrp.clone(DeepCopy)), // create via clone
  oldsoln(*oldsolnptr),		// reference to just-created pointer
  dirptr(xgrp.getX().clone(CopyShape)), // create via clone 
  dir(*dirptr),			// reference to just-created pointer
  test(t),			// reference to t
  iparams(p),			// reference to p
  oparams(),			// empty list
  linesearch(p.sublist("Line Search")), // initialize line search
  step(0.0),			// initialize to zero
  niter(0),			// initialize to zero
  status(Status::Unconverged)	// initialize convergence status
{
  // Set up utilities for printing, etc.
  Utils::setUtils(iparams);

  // Print out initialization information
  if (Utils::doPrint(1)) {

    cout << "\n" << Utils::fill(72) << "\n";
    cout << "\n-- Parameters Passed to Nonlinear Solver on Print Processor--\n\n";
    iparams.print(cout,5);
    cout << "\n-- Status Tests Passed to Nonlinear Solver on Print Processor--\n\n";
    test.print(cout, 5);
    cout <<"\n" << Utils::fill(72) << "\n";

  }

  // Compute RHS of initital guess
  soln.computeRHS();
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
  status = test(*this);
  return status;
}

Status::StatusType Newton::iterate()
{
  // Reset forcing term.
  resetForcingTerm();

  // Compute Jacobian at current solution.
  soln.computeJacobian();

  // Compute Newton direction for current solution.
  /* NOTE FROM TAMMY: Need to check the return status! */
  soln.computeNewton(iparams.sublist("Linear Solver"));

  // Set search direction.
  dir = oldsoln.getNewton();

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
  status = test(*this);
  printUpdate();

  // Iterate until converged or failed
  while (status == Status::Unconverged) {
    status = iterate();
    printUpdate();
  }

  return status;
}

Abstract::Group& Newton::getSolutionGroup() const
{
  return soln;
}

Abstract::Group& Newton::getPreviousSolutionGroup() const
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
  oparams.setParameter("2-Norm of Residual", soln.getNormRHS());
  return oparams;
}

// protected
void Newton::printUpdate() 
{
  double norm_k;
  double norm_newton;

  // All processors participate in the computation of these norms...
  if (Utils::doAllPrint(1)) {
    norm_k = soln.getNormRHS();
    norm_newton = (niter > 0) ? oldsoln.getNewton().norm() : 0;
  }

  // ...But only the print processors actually prints the result.
  if (Utils::doPrint(1)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "Newton Step " << niter;
    cout << " : Residual Norm = " << Utils::sci(norm_k);
    cout << "  Step = " << Utils::sci(step);
    cout << "  Update Norm = " << Utils::sci(norm_newton);
    cout << "\n" << Utils::fill(72) << "\n" << endl;
  }
  
  if ((status > 0) && (Utils::doPrint(1)))
    cout << "\n" << "Solution is CONVERGED!" << "\n" << endl;
  
  if ((status < 0) && (Utils::doPrint(1)))
    cout << "\n" << "Nonlinear solver failed." << "\n" << endl;
}

// protected
void Newton::resetForcingTerm()
{
  // Reset the forcing term at the beginning on a nonlinear iteration,
  // based on the last iteration.

  if (niter == 0)
    return;

  if ((!iparams.isParameter("Forcing Term Method")) ||
      (iparams.isParameterEqual("Forcing Term Method", "None")))
    return;

  // Get forcing term parameters.
  const string method = iparams.getParameter("Forcing Term Method", "");
  const double eta_min = iparams.getParameter("Forcing Term Minimum Tolerance", 1.0e-6);
  const double eta_max = iparams.getParameter("Forcing Term Maximum Tolerance", 0.01);

  // Get linear solver parameter list and current tolerance.
  Parameter::List& lsparams = iparams.sublist("Linear Solver");
  const double eta_km1 = lsparams.getParameter("Tolerance", 0.0);

  // New forcing term.
  double eta_k;

  string indent = "       ";

  if (Utils::doPrint(1)) {
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
    const double normrhs = soln.getNormRHS();
    const double normoldrhs = oldsoln.getNormRHS();

    // Compute forcing term
    eta_k = fabs(normrhs - normpredrhs) / normoldrhs;
     
    // Some output
    if (Utils::doPrint(1)) {
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
    
    const double normrhs = soln.getNormRHS();
    const double normoldrhs = oldsoln.getNormRHS();
    const double alpha = iparams.getParameter("Forcing Term Alpha", 1.5);
    const double gamma = iparams.getParameter("Forcing Term Gamma", 0.9);
    const double residual_ratio = normrhs / normoldrhs;
    
    eta_k = gamma * pow(residual_ratio, alpha);
     
    // Some output
    if (Utils::doPrint(1)) {
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

  if (Utils::doPrint(1)) {
    cout << indent << "Forcing Term: " << eta_k << endl;
  }

}




