// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_Newton.H"

NLS_Newton::NLS_Newton(NLS_Interface& i, NLS_ParameterList& p) :
  interface(i),
  soln(interface.getInitialGuess()),
  oldsoln(interface.createNewGroup()),
  niter(0)
{
  // Get initial guess for x and corresponding rhs
  soln.computeRHS();
}

int NLS_Newton::iterate()
{
  // compute Jacobian at current solution
  soln.computeJacobian();

  // compute Newton direction for current solution
  soln.computeNewtonDirection();

  // copy current solution to old solution
  oldsoln = soln;

  // update current solution
  soln.computeX(oldsoln, oldsoln.getNewton(), static_cast<double>(1.0));

  // compute RHS for new current solution
  soln.computeRHS();

  // compute norm of rhs
  normrhs = soln.getRHS().norm();

  // update iteration count
  niter ++;

  // completed successful iteration
  return true;
}

double NLS_Newton::getNormRHS() const
{
  return normrhs;
}


