// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_Newton.H"

NLS_Newton::NLS_Newton(NLS_Group& i, NLS_Group& s, NLS_ParameterList& p) :
  oldsoln(i),
  soln(s),
  params(p),
  niter(0)
{
  // Get initial guess for x and corresponding rhs
  soln.computeRHS();
}

NLS_Newton::~NLS_Newton() 
{
}

void NLS_Newton::resetParameters(NLS_ParameterList& p)
{
}

bool NLS_Newton::isConverged() 
{
  //Compute norms
  normrhs = soln.getRHS().norm();
  double normupdate = soln.getNewton().norm();

  // Output 
  if ((params.getParameter("MyPID",0)==0) && 
      (params.getParameter("Output Level",4) >= 2)) {
    cout << "**************************************************************"
	 << endl;
    cout << "Newton Step " << niter << " : Residual Norm = " << normrhs 
	 << "  Update Norm = " << normupdate << endl;
    cout << "**************************************************************"
	 << endl;
  }

  if ((normrhs < params.getParameter("Absolute Tolerance",1.0e-10))
      &&(normupdate < params.getParameter("Relative Tolerance",1.0e-6))) {
    if (params.getParameter("MyPID",0)==0) 
      cout << "Solution is CONVERGED!" << endl;
    return true;
  }
  return false;
}

      
int NLS_Newton::iterate()
{
  // compute Jacobian at current solution
  soln.computeJacobian();

  // compute Newton direction for current solution
  soln.computeNewton(params);

  // copy current group to the old group
  oldsoln.copy(soln);

  // update current solution
  soln.computeX(oldsoln, oldsoln.getNewton(), static_cast<double>(-1.0));

  // compute RHS for new current solution
  soln.computeRHS();

  // compute norm of rhs
  normrhs = soln.getRHS().norm();

  // update iteration count
  niter ++;

  // completed successful iteration
  return true;
}

int NLS_Newton::solve()
{
  cout << "Beginning nonlinear solve with Newtons method" << endl;
  // Check for convergence of initial guess
  isConverged();

  const NLS_Vector& RHS = soln.getX();

  // Get Parameters
  int maxit = params.getParameter("Max Nonlinear Iterations", 15);
  for (int i=0; i<maxit; i++) {
    iterate();
    if (isConverged()) break;
  }
}

NLS_Group& NLS_Newton::getSolutionGroup() const
{

}

bool NLS_Newton::getProfileInfo(string& name, NLS_Parameter& p) const
{

}

double NLS_Newton::getNormRHS() const
{
  return normrhs;
}


