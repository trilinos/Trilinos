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
  maxiter = params.getParameter("Max Iterations", 15);
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

NLS_Method::STATUS NLS_Newton::getStatus() 
{
  // Compute norm of Newton step 
  /* NOTE FROM TAMMY: This only works when we take full Newton
     steps. Need to change it if we do a linesearch. */
  double normupdate = soln.getNewton().norm();

  if (NLS_Utilities::doPrint(1)) {
    cout << "\n" << stars;
    cout << "Newton Step " << niter 
	 << " : Residual Norm = " << soln.getNormRHS()
	 << "  Update Norm = " << normupdate;
    cout << "\n" << stars << endl;
  }

  NLS_Method::STATUS status = NLS_Method::NotConverged;

  if (soln.getNormRHS() < abstol)
    status = NLS_Method::ConvergedAbsTol;

  if ((niter > 0) && (normupdate < reltol))
    status = NLS_Method::ConvergedRelTol;

  if (niter >= maxiter) {
    status = NLS_Method::MaxItersExceeded;
  }

  if ((status > 0) && (NLS_Utilities::doPrint(1)))
    cout << "\n" << "Solution is CONVERGED!" << "\n" << endl;

  if ((status < 0) && (NLS_Utilities::doPrint(1)))
    cout << "\n" << "Nonlinear solver failed." << "\n" << endl;

  return status;
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
  soln.computeNewton(params.sublist("Linear Solver Parameters"));

  // copy current soln to the old soln
  oldsoln = soln;

  // compute new solution
  soln.computeX(oldsoln, oldsoln.getNewton(), 1.0);

  // compute RHS for new current solution
  soln.computeRHS();

  // update iteration count
  niter ++;

  // completed successful iteration
  return getStatus();
}

NLS_Method::STATUS NLS_Newton::solve()
{
  if (NLS_Utilities::doPrint(2)) 
    cout << "\n" << "Beginning nonlinear solve with Newton's method!" << endl;

  NLS_Method::STATUS status = getStatus();

  // Check for convergence of initial guess
  if (status != NotConverged)
    return status;

  // Iterate until converged or reach maxiter
  while (status == NotConverged) {
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



