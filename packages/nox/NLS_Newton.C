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

const string stars("***********************************************************************\n");

NLS_Newton::NLS_Newton(NLS_Group& initialguess, NLS_Group& workspace, NLS_ParameterList& p) :
  oldsoln(workspace),
  soln(initialguess),
  params(p),
  niter(0)
{
  NLS_Utilities::setUtilities(params);
  maxiter = params.getParameter("Max Nonlinear Iterations", 15);
  abstol = params.getParameter("Absolute Tolerance", 1.0e-9);
  reltol = params.getParameter("Relative Tolerance", 1.0e-4);
  soln.computeRHS();
}

NLS_Newton::~NLS_Newton() 
{
}

void NLS_Newton::resetParameters(NLS_ParameterList& p)
{
}

NLS_Method::STATUS NLS_Newton::isConverged() 
{
  // Compute norm of Newton step
  double normupdate = soln.getNewton().norm();

  // Output 
  if (NLS_Utilities::doPrint(1)) {

    cout << "\n" << stars;
    cout << "Newton Step " << niter 
	 << " : Residual Norm = " << soln.getNormRHS()
	 << "  Update Norm = " << normupdate;
    cout << "\n" << stars << endl;
  }

  if ((soln.getNormRHS() < abstol) && (normupdate < reltol)) {
    if (NLS_Utilities::doPrint(1)) 
      cout << "\n" << "Solution is CONVERGED!" << "\n" << endl;
    return NLS_Method::Converged;
  }

  // Check number of iterations
  if (niter >= maxiter) {
    if (NLS_Utilities::doPrint(1)) 
      cout << "\n" << "Max iterations exceeded in nonlinear solver." << "\n" << endl;
    return NLS_Method::MaxItersExceeded;
  }

  return NLS_Method::NotConverged;
}

      
NLS_Method::STATUS NLS_Newton::iterate()
{
  // compute the linear solver convergence criteria
  if ((niter > 0) && (params.isParameterEqual("Forcing Term Method", "None"))) {

    double tol;			// linear solver tolerance
    
    // Only compute the norm of the predicted RHS from the last iteration if necessary
    if (params.isParameterEqual("Forcing Term Method", "Type 1")) {
      double normPredRHS = oldsoln.computeNormPredictedRHS(oldsoln.getNewton());
      tol = forcingTerm(params, soln.getNormRHS(), oldsoln.getNormRHS(), normPredRHS);
    }
    else
      tol = forcingTerm(params, soln.getNormRHS(), oldsoln.getNormRHS());

    // Reset linear solver tolerance
    params.setParameter("Linear Solver Tolerance", tol);
  }

  // compute Jacobian at current solution
  soln.computeJacobian();

  // compute Newton direction for current solution
  soln.computeNewton(params);

  // copy current soln to the old soln
  oldsoln = soln;

  // compute new solution
  soln.computeX(oldsoln, oldsoln.getNewton(), 1.0);

  // compute RHS for new current solution
  soln.computeRHS();

  // update iteration count
  niter ++;

  // completed successful iteration
  return isConverged();
}

NLS_Method::STATUS NLS_Newton::solve()
{
  if (NLS_Utilities::doPrint(2)) 
    cout << "\n" << "Beginning nonlinear solve with Newtons method!" << endl;

  NLS_Method::STATUS status = isConverged();

  // Check for convergence of initial guess
  if (status != NotConverged)
    return status;

  // Iterate until converged or reach maxiter
  while (status != NotConverged) {
    status = iterate();
  }

  return status;
}

NLS_Group& NLS_Newton::getSolutionGroup() const
{
  return soln;
}

bool NLS_Newton::getProfileInfo(string& name, NLS_Parameter& p) const
{
  return false;
}



