// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_Newton.H"		// class definition
#include "NLS_Utilities.H"	// for static doPrint function
#include "NLS_AdaptiveTolerance.H" // adaptive convergence tolerance
#include <iomanip>		// for setw

const string stars = NLS_Utilities::stars;

NLS_Newton::NLS_Newton(NLS_Group& initialguess, NLS_Group& workspace, NLS_ParameterList& p) :
  oldsoln(workspace),
  soln(initialguess),
  params(p),
  niter(0),
  linesearch(p.sublist("Line Search Parameters"))
{
  NLS_Utilities::setUtilities(params);

  if (NLS_Utilities::doPrint(4))
    cout << "Output Level " << NLS_Utilities::outputLevel << "." << endl; 

  if (NLS_Utilities::doPrint(5)) 
    cout << "NLS_Utilities: Processor " << NLS_Utilities::myPID 
	 << " is online." << endl;  

  //  maxiter = params.getParameter("Max Iterations", 15);
  //  abstol = params.getParameter("Absolute Tolerance", 1.0e-9);
  //  reltol = params.getParameter("Relative Tolerance", 1.0e-4);
  soln.computeRHS();
  statustest.setup(p.sublist("Convergence Tests"), soln.getNormRHS());
  step = 0;
  //NLS_ParameterList& tmp = params.sublist("Line Search Parameters");
  //defaultstep = tmp.getParameter("Default Step", 1.0);
  defaultstep = (params.sublist("Line Search Parameters")).getParameter("Default Step", 1.0);
}

NLS_Newton::~NLS_Newton() 
{
}

void NLS_Newton::resetParameters(NLS_ParameterList& p)
{
}

NLS_Method::STATUS NLS_Newton::getStatus() 
{
  double norm_k;
  double norm_km1;
  if (NLS_Utilities::doAllPrint(1)) {
    norm_k = soln.getNormRHS();
    norm_km1 = oldsoln.getNewton().norm();
  }
  if (NLS_Utilities::doPrint(1)) {
    cout.setf(ios::scientific);
    cout.precision(NLS_Utilities::precision);
    cout << "\n" << stars;
    cout << "Newton Step " << niter 
	 << " : Residual Norm = " << setw(NLS_Utilities::precision + 6) << norm_k
	 << "  Step = " << setw(NLS_Utilities::precision + 6) << step
	 << "  Update Norm = " << setw(NLS_Utilities::precision + 6) << norm_km1;
    cout << "\n" << stars << endl;
    cout.unsetf(ios::scientific);
  }
  
  NLS_Method::STATUS status = statustest(oldsoln, soln, oldsoln.getNewton(), step, niter);
  
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

    NLS_AdaptiveTolerance adapttol;  
    double tol;	 
    
    // Only compute the norm of the predicted RHS from the last iteration if necessary
    if (params.isParameterEqual("Forcing Term Method", "Type 1")) {
      double normPredRHS = oldsoln.computeNormPredictedRHS(oldsoln.getNewton());
      tol = adapttol(params, soln.getNormRHS(), oldsoln.getNormRHS(), normPredRHS);
    }
    else
      tol = adapttol(params, soln.getNormRHS(), oldsoln.getNormRHS());

    // Reset linear solver tolerance
    params.setParameter("Linear Solver Tolerance", tol);
  }

  // compute Jacobian at current solution
  soln.computeJacobian();

  // compute Newton direction for current solution
  soln.computeNewton(params.sublist("Linear Solver Parameters"));

  // copy current soln to the old soln
  oldsoln = soln;

  // Step default step
  step = defaultstep;

  // Do line search
  linesearch(soln, step, oldsoln, oldsoln.getNewton());

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



