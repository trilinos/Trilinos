// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_MethodManager.H"
#include "NLS_Newton.H"

NLS_Newton::NLS_Newton(NLS_Group& i, NLS_Group& s, NLS_ParameterList& p) :
  oldsoln(i),
  soln(s),
  params(p),
  niter(0)
{
  maxiter = p.getParameter("Max Nonlinear Iterations", 15);
  abstol = p.getParameter("Absolute Tolerance",1.0e-9);
  reltol = p.getParameter("Relative Tolerance",1.0e-4);
  soln.copy(oldsoln);
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
  if (NLS_MethodManager::util.isPrintProc() && 
      NLS_MethodManager::util.isOutput(1)) {
    cout << endl 
	 << "***********************************************************************"
	 << endl;
    cout << "Newton Step " << niter << " : Residual Norm = " << normrhs 
	 << "  Update Norm = " << normupdate << endl;
    cout << "***********************************************************************"
	 << endl << endl;;
  }

  if ((normrhs < abstol)
      &&(normupdate < reltol)) {
    if (NLS_MethodManager::util.isPrintProc() &&  
	NLS_MethodManager::util.isOutput(1)) 
      cout << endl << "Solution is CONVERGED!" << endl << endl;
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
  soln.computeX(oldsoln, soln.getNewton(), static_cast<double>(-1.0));

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
  if (NLS_MethodManager::util.isPrintProc() && 
      NLS_MethodManager::util.isOutput(2)) 
    cout << endl << "Beginning nonlinear solve with Newtons method!" << endl;

  // Check for convergence of initial guess
  isConverged();
  
  const NLS_Vector& RHS = soln.getX();

  // Get Parameters
  for (int i=0; i<maxiter; i++) {
    iterate();
    if (isConverged()) break;
  }
  return niter;
}

NLS_Group& NLS_Newton::getSolutionGroup() const
{
  return soln;
}

bool NLS_Newton::getProfileInfo(string& name, NLS_Parameter& p) const
{
  return false;
}

double NLS_Newton::getNormRHS() const
{
  return normrhs;
}


