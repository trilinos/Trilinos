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
#include "NOX_Parameter_UserNorm.H"

NOX::Direction::Newton::Newton(const NOX::Utils& u, NOX::Parameter::List& p) :
  utils(u),
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

  doRescue = p.getParameter("Rescue Bad Newton Solve", true);
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
  
  // It didn't converge, but maybe we can recover. 
  if ((status != NOX::Abstract::Group::Ok) &&
      (doRescue == false)) {
    NOX::Direction::Newton::throwError("compute", 
				       "Unable to solve Newton system");
  }
  else if ((status != NOX::Abstract::Group::Ok) &&
	   (doRescue == true)) {
    if (utils.isPrintProcessAndType(NOX::Utils::Warning))
      cout << "WARNING: NOX::Direction::Newton::compute() - Linear solve "
	   << "failed to achieve convergence - using the step anyway " 
	   << "since \"Rescue Bad Newton Solve\" is true " << endl;
  }

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

  if (utils.isPrintProcessAndType(Utils::Details)) {
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

      // do NOT use the following lines!! This does NOT account for 
      // line search step length taken.
      // const double normpredf = 0.0;
      // oldsoln.getNormLastLinearSolveResidual(normpredf);
      
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
      
      // Compute the norms
      double normpredf = 0.0;
      double normf = 0.0;
      double normoldf = 0.0;

      if (paramsPtr->sublist("Newton").isParameter("Forcing Term User Defined Norm")) {

	const NOX::Parameter::Arbitrary& arbitrary = paramsPtr->
	  sublist("Newton").getArbitraryParameter("Forcing Term User Defined Norm");
	const NOX::Parameter::UserNorm* userNorm = 0;
	userNorm = dynamic_cast<const NOX::Parameter::UserNorm*>(&arbitrary);
	
	if (userNorm != 0) {
	  if (utils.isPrintProcessAndType(Utils::Details)) {
	    cout << indent << "Forcing Term Norm: " << userNorm->getType()
		 << endl;
	  }
	  normpredf = userNorm->computeNorm(*predRhs);
	  normf = userNorm->computeNorm(soln.getF());
	  normoldf = userNorm->computeNorm(oldsoln.getF());
	}
	else {
	  if (utils.isPrintProcessAndType(Utils::Warning)) {
	    cout << "WARNING: NOX::Direction::Newton::resetForcingTerm() - "
		 << "\"Forcing Term User Defined Norm\" is not of type "
		 << "NOX::Parameter::UserNorm!\n" 
		 << "Defaulting to L-2 Norms!" << endl; 
	  }
	  normpredf = predRhs->norm();
	  normf = soln.getNormF();
	  normoldf = oldsoln.getNormF();
	}
      }
      else {
	if (utils.isPrintProcessAndType(Utils::Details)) {
	  cout << indent << "Forcing Term Norm: Using L-2 Norm."
	       << endl;
	}
	normpredf = predRhs->norm();
	normf = soln.getNormF();
	normoldf = oldsoln.getNormF();
      }      

      // Compute forcing term
      eta_k = fabs(normf - normpredf) / normoldf;
      
      // Some output
      if (utils.isPrintProcessAndType(Utils::Details)) {
	cout << indent << "Residual Norm k-1 =             " << normoldf << "\n";
	cout << indent << "Residual Norm Linear Model k =  " << normpredf << "\n";
	cout << indent << "Residual Norm k =               " << normf << "\n";
	cout << indent << "Calculated eta_k (pre-bounds) = " << eta_k << endl;
      }
      
      // Impose safeguard and constraints ...
      const double alpha = (1.0 + sqrt(5.0)) / 2.0;
      const double eta_km1_alpha = pow(eta_km1, alpha);
      if (eta_km1_alpha > 0.1) 
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

      double normf = 0.0;
      double normoldf = 0.0;

      if (paramsPtr->sublist("Newton").isParameter("Forcing Term User Defined Norm")) {

	const NOX::Parameter::Arbitrary& arbitrary = paramsPtr->
	  sublist("Newton").getArbitraryParameter("Forcing Term User Defined Norm");
	const NOX::Parameter::UserNorm* userNorm = 0;
	userNorm = dynamic_cast<const NOX::Parameter::UserNorm*>(&arbitrary);
	
	if (userNorm != 0) {
	  if (utils.isPrintProcessAndType(Utils::Details)) {
	    cout << indent << "Forcing Term Norm: " << userNorm->getType()
		 << endl;
	  }
	  normf = userNorm->computeNorm(soln.getF());
	  normoldf = userNorm->computeNorm(oldsoln.getF());
	}
	else {
	  if (utils.isPrintProcessAndType(Utils::Warning)) {
	    cout << "WARNING: NOX::Direction::Newton::resetForcingTerm() - "
		 << "\"Forcing Term User Defined Norm\" is not of type "
		 << "NOX::Parameter::UserNorm!\n" 
		 << "Defaulting to L-2 Norms!" << endl; 
	  }
	  normf = soln.getNormF();
	  normoldf = oldsoln.getNormF();
	}
      }
      else {
	if (utils.isPrintProcessAndType(Utils::Details)) {
	  cout << indent << "Forcing Term Norm: Using L-2 Norm."
	       << endl;
	}
	normf = soln.getNormF();
	normoldf = oldsoln.getNormF();
      }  

      const double alpha = paramsPtr->sublist("Newton").getParameter("Forcing Term Alpha", 1.5);
      const double gamma = paramsPtr->sublist("Newton").getParameter("Forcing Term Gamma", 0.9);
      const double residual_ratio = normf / normoldf;
      
      eta_k = gamma * pow(residual_ratio, alpha);
      
      // Some output
      if (utils.isPrintProcessAndType(Utils::Details)) {
	cout << indent << "Residual Norm k-1 =             " << normoldf << "\n";
	cout << indent << "Residual Norm k =               " << normf << "\n";
	cout << indent << "Calculated eta_k (pre-bounds) = " << eta_k << endl;
      }
      
      // Impose safeguard and constraints ... 
      const double eta_k_alpha = gamma * pow(eta_km1, alpha);
      if (eta_k_alpha > 0.1) 
	eta_k = max(eta_k, eta_k_alpha);
      eta_k = max(eta_k, eta_min);
      eta_k = min(eta_max, eta_k);
    }

  }

  else {

    if (utils.isPrintProcessAndType(Utils::Warning))
      cout << "NOX::Direction::Newton::resetForcingTerm "
	   << "- invalid forcing term method (" << method << ")" << endl;

    return false;
  }

  // Reset linear solver tolerance
  paramsPtr->sublist("Newton").sublist("Linear Solver").setParameter("Tolerance", eta_k);

  if (utils.isPrintProcessAndType(Utils::Details)) 
    cout << indent << "Forcing Term: " << eta_k << endl;
  
  return true;
}

void NOX::Direction::Newton::throwError(const string& functionName, const string& errorMsg)
{
    if (utils.isPrintProcessAndType(Utils::Error))
      cerr << "NOX::Direction::Newton::" << functionName << " - " << errorMsg << endl;
    throw "NOX Error";
}
