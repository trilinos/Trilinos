// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_StatusTest.H"
#include "NLS_Utilities.H"

NLS_StatusTest::NLS_StatusTest()
{
}

NLS_StatusTest::~NLS_StatusTest()
{

}

void NLS_StatusTest::setup(const NLS_ParameterList& params, double normrhs)
{
  maxiter = params.getParameter("Max Iterations", 100);
  absresidtol = params.getParameter("Absolute Residual Tolerance", 1.0e-12);
  relresidtol = params.getParameter("Relative Residual Tolerance", 0.0);
  maxresidtol = params.getParameter("Max Norm Residual Tolerance", 0.0);
  ratioresidtol = params.getParameter("Minimum Convergence Rate", 0.0);
  absupdatetol = params.getParameter("Absolute Update Tolerance", 0.0);
  relupdate_epsilon_r = params.getParameter("Relative Update Epsilon r", 0.0);
  relupdate_epsilon_a = params.getParameter("Relative Update Epsilon a", 0.0);

  normrhszero = normrhs;

  if (NLS_Utilities::doPrint(1)) {
    
    cout.setf(ios::scientific);
    cout.precision(NLS_Utilities::precision);

    cout << "\n" << NLS_Utilities::stars;

    cout << "Nonlinear Solver Convergence (and Failure) Tests->\n";

    if (absresidtol > 0) 
      cout << "  2-norm of residual < " << absresidtol << "\n";

    if (maxresidtol > 0)
      cout << "  Max-norm of residual < " << maxresidtol << "\n";
  
    if (relresidtol > 0) 
      cout << "  Relative 2-norm of residual < " << relresidtol << "\n";
    
    if (absupdatetol > 0)
      cout << "  2-norm of update < " << absupdatetol << " \n";

    if ((relupdate_epsilon_r > 0) && (relupdate_epsilon_a > 0)) 
      cout << "  Relative update test with epsilon_a = " 
	   << relupdate_epsilon_a << " and epsilon_r = "
	   << relupdate_epsilon_r << "\n";
    
    if (ratioresidtol > 1)
      cout << "  Require residual convergence rate > " << ratioresidtol << "\n";

    if (maxiter > 0) 
      cout << "  Require max iterations < " << maxiter << "\n";

    cout << NLS_Utilities::stars << endl;

    cout.unsetf(ios::scientific);
  }

}

NLS_Method::STATUS NLS_StatusTest::operator()(const NLS_Group& oldsoln, const NLS_Group& soln, 
					    const NLS_Vector& update, double step, int niter)
{

  double normrhs = soln.getNormRHS();

  // 2-Norm of residual test
  if (absresidtol > 0) 
    if (normrhs < absresidtol)
      return NLS_Method::ConvAbsResidTol;

  // Max-Norm of residual test
  if (maxresidtol > 0)
    if (soln.getRHS().norm(NLS_Vector::INF) < maxresidtol)
      return NLS_Method::ConvMaxResidTol;
  
  // All other tests require at least one iteration
  if (niter == 0)
    return NLS_Method::NotConverged;

  // Relative residual norm (relative to initial guess)
  if (relresidtol > 0)
    if ((normrhs / normrhszero) < relresidtol)
      return NLS_Method::ConvRelResidTol;

  // 2-Norm of update
  /* NOTE FROM TAMMY: Note so sure about this test. Could lead to
     false convergence with the line search. */
  if (absupdatetol > 0)
    if ((step * update.norm(NLS_Vector::INF)) < absupdatetol)
      return NLS_Method::ConvAbsUpdateTol;

  // Compute special convergence test
  // (1/N) * SUM_i [ step * abs(update_i) / (epsilon_r * abs(soln_i) + epsilon_a) ] < 1
  /* NOTE FROM TAMMY: What is the citation for this test?? */
  if ((relupdate_epsilon_r > 0) && (relupdate_epsilon_a > 0)) {
    NLS_Vector* a = update.clone();
    NLS_Vector* b = update.clone();
    a->init(relupdate_epsilon_a);
    b->abs(soln.getX());
    a->update(relupdate_epsilon_r, *b, 1.0);
    b->reciprocal(*a);
    a->abs(update);
    a->scale(step);
    double result = a->dot(*b) / a->length();
    delete a,b;
    if (result < 1)
      return NLS_Method::ConvRelUpdateTol;
  }

  // Check for sufficient reduction
  if (ratioresidtol > 1)
    if ((oldsoln.getNormRHS() / normrhs) < ratioresidtol)
      return NLS_Method::TooSlowConv;

  // Maximum iterations test
  if (maxiter > 0) 
    if (niter >= maxiter)
      return NLS_Method::MaxItersExceeded;

  return NLS_Method::NotConverged;

}
  


