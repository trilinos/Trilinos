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

NOX::Direction::Newton::Newton(NOX::Parameter::List& p) :
  predRhs(NULL),
  stepDir(NULL)
  
{
  reset(p);
}

NOX::Direction::Newton::~Newton()
{
  delete predRhs;
  delete stepDir;
}

bool NOX::Direction::Newton::reset(NOX::Parameter::List& params)
{
  paramsPtr = &params;

  NOX::Parameter::List& p = params.sublist("Newton");

  doRescue = p.getParameter("Resuce Bad Newton Solve", true);
  if (!p.sublist("Linear Solver").isParameter("Tolerance"))
    p.sublist("Linear Solver").getParameter("Tolerance", 1.0e-10);
  return true;
}

bool NOX::Direction::Newton::compute(NOX::Abstract::Vector& dir, 
				     NOX::Abstract::Group& soln, 
				     const NOX::Solver::Generic& solver)
{
  NOX::Abstract::Group::ReturnType status;

  // Compute F at current solution.
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok) 
    NOX::Direction::Newton::throwError("compute", "Unable to compute F");

  // Reset the linear solver tolerance.
  resetForcingTerm(soln, solver.getPreviousSolutionGroup(), 
		   solver.getNumIterations(), solver.getParameterList());

  // Compute Jacobian at current solution.
  status = soln.computeJacobian();
  if (status != NOX::Abstract::Group::Ok) 
    NOX::Direction::Newton::throwError("compute", "Unable to compute Jacobian");
  
  // Compute the Newton direction
  status = soln.computeNewton(paramsPtr->sublist("Newton").sublist("Linear Solver"));
  
  // It didn't converge, but maybe we can recover. Otherwise, we throw an error.
  if (status == NOX::Abstract::Group::NotConverged)
  { 
    if (!NOX::Direction::Newton::rescueBadNewtonSolve(soln))
      return false;
  }
  else if (status != NOX::Abstract::Group::Ok) 
    NOX::Direction::Newton::throwError("compute", "Unable to solve Newton system");
    
  // Set search direction.
  dir = soln.getNewton();

  return true;
}


// protected
bool NOX::Direction::Newton::resetForcingTerm(const NOX::Abstract::Group& soln, 
			      const NOX::Abstract::Group& oldsoln, 
			      int niter,
			      const NOX::Parameter::List& solverParams)
{
  // Reset the forcing term at the beginning on a nonlinear iteration,
  // based on the last iteration.

  if ((!paramsPtr->sublist("Newton").isParameter("Forcing Term Method")) ||
      (paramsPtr->sublist("Newton").isParameterEqual("Forcing Term Method", "None")))
    return true;

  // Get forcing term parameters.
  const string method = paramsPtr->sublist("Newton").getParameter("Forcing Term Method", "Constant");
  const double eta_min = paramsPtr->sublist("Newton").getParameter("Forcing Term Minimum Tolerance", 1.0e-4);  const double eta_max = paramsPtr->sublist("Newton").getParameter("Forcing Term Maximum Tolerance", 0.9);

  // Get linear solver current tolerance.

  double eta_km1 = 0.0;
  if (((method == "Type 1") || (method == "Type 2")) 
	&& (solverParams.sublist("Line Search").isParameterDouble("Adjusted Tolerance"))) {
    
    // Tolerance may have been adjusted in a line search algorithm   
    eta_km1 = solverParams.sublist("Line Search").getParameter("Adjusted Tolerance", 0.0);
    
  }
  else {
    // Default to the old tolerance
    eta_km1 = paramsPtr->sublist("Newton").sublist("Linear Solver").getParameter("Tolerance", 0.0);
    
  }

  // New forcing term.
  double eta_k;

  const string indent = "       ";

  if (Utils::doPrint(Utils::Details)) {
    cout << indent << "CALCULATING FORCING TERM" << endl;
    cout << indent << "Method: " << method << endl;
  }

  if (method == "Constant") {    

    if (eta_km1 != 0.0) 
      eta_k = eta_km1;
    else
      eta_k = eta_min;

  }        

  else if (method == "Type 1") {
    
    if (niter == 0) {

      eta_k = paramsPtr->sublist("Newton").getParameter("Forcing Term Initial Tolerance", 0.01);

    }
    else {

      // Return norm of predicted F

      // do NOT use the following line!! This does NOT account for 
      // line search step length taken.
      //const double normpredf = oldsoln.getNormNewtonSolveResidual();
      
      // Create a new vector to be the predicted RHS
      if (predRhs == NULL) {
	predRhs = oldsoln.getF().clone(ShapeCopy);
      }
      if (stepDir == NULL) {
	stepDir = oldsoln.getF().clone(ShapeCopy);
      }
      
      // stepDir = X - oldX (i.e., the step times the direction)
      stepDir->update(1.0, soln.getX(), -1.0, oldsoln.getX(), 0);
      
      // Compute predRhs = Jacobian * step * dir
      oldsoln.applyJacobian(*stepDir, *predRhs);
      
      // Compute predRhs = RHSVector + predRhs (this is the predicted RHS)
      predRhs->update(1.0, oldsoln.getF(), 1.0);
      
      // Return norm of predicted RHS
      const double normpredf = predRhs->norm();
      
      if (normpredf < 0) {
	cerr << "NOX::Direction::Newton::resetForcingTerm " 
	     << "- getNormNewtonSolveResidual returned a negative value" << endl;
      }

      // Get other norms
      const double normf = soln.getNormF();
      const double normoldf = oldsoln.getNormF();
      
      // Compute forcing term
      eta_k = fabs(normf - normpredf) / normoldf;
      
      // Some output
      if (Utils::doPrint(Utils::Details)) {
	cout << indent << "Residual Norm k-1 =             " << normoldf << "\n";
	cout << indent << "Residual Norm Linear Model k =  " << normpredf << "\n";
	cout << indent << "Residual Norm k =               " << normf << "\n";
	cout << indent << "Calculated eta_k (pre-bounds) = " << eta_k << endl;
      }
      
      // Impose safeguard and constraints ...
      const double alpha = (1.0 + sqrt(5.0)) / 2.0;
      const double eta_km1_alpha = pow(eta_km1, alpha);
      eta_k = max(eta_k, eta_km1_alpha);
      eta_k = max(eta_k, eta_min);
      eta_k = min(eta_max, eta_k);
    }
  }
    
  else if (method == "Type 2") {  
    
    if (niter == 0) {

      eta_k = paramsPtr->sublist("Newton").getParameter("Forcing Term Initial Tolerance", 0.01);

    }
    else {

      const double normf = soln.getNormF();
      const double normoldf = oldsoln.getNormF();
      const double alpha = paramsPtr->sublist("Newton").getParameter("Forcing Term Alpha", 1.5);
      const double gamma = paramsPtr->sublist("Newton").getParameter("Forcing Term Gamma", 0.9);
      const double residual_ratio = normf / normoldf;
      
      eta_k = gamma * pow(residual_ratio, alpha);
      
      // Some output
      if (Utils::doPrint(Utils::Details)) {
	cout << indent << "Residual Norm k-1 =             " << normoldf << "\n";
	cout << indent << "Residual Norm k =               " << normf << "\n";
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
      cout << "NOX::Direction::Newton::resetForcingTerm "
	   << "- invalid forcing term method (" << method << ")" << endl;

    return false;
  }

  // Reset linear solver tolerance
  paramsPtr->sublist("Newton").sublist("Linear Solver").setParameter("Tolerance", eta_k);

  if (Utils::doPrint(Utils::Details)) 
    cout << indent << "Forcing Term: " << eta_k << endl;
  
  return true;
}


bool  NOX::Direction::Newton::rescueBadNewtonSolve(const NOX::Abstract::Group& grp) const
{
  //! Check if the "rescue" option has been selected
  if (!doRescue)
    return false;

  //! See if the group has compute the accuracy
  double accuracy;
  NOX::Abstract::Group::ReturnType status = grp.getNormLastLinearSolveResidual(accuracy);
    
  // If this functionality is not supported in the group, return false
  /* NOTE FROM TAMMY: We could later modify this to acutally caluclate
     the error itself if it's just a matter of the status being
     NotDefined. */
  if (status != NOX::Abstract::Group::Ok) 
    return false;

  // Check if there is any improvement in the relative residual
  double normF = grp.getNormF();

  // If we can't reduce the relative norm at all, we're not happy
  if (accuracy >= normF) 
    return false;

  // Otherwise, we just print a warning and keep going
  if (Utils::doPrint(Utils::Warning)) 
    cout << "WARNING: NOX::Direction::Newton::compute - Unable to achieve desired linear solve accuracy." << endl;

  return true;
}


void NOX::Direction::Newton::throwError(const string& functionName, const string& errorMsg)
{
    if (Utils::doPrint(Utils::Error))
      cerr << "NOX::Direction::Newton::" << functionName << " - " << errorMsg << endl;
    throw "NOX Error";
}
