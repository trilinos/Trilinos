// $Id$ 
// $Source$ 

#ifdef WITH_PRERELEASE

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "NOX_Direction_Utils_InexactNewton.H" // class definition

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_LineSearchBased.H"
#include "NOX_Utils.H"
#include "NOX_Parameter_List.H"
#include "NOX_Parameter_UserNorm.H"

/* Some compilers (in particular the SGI and ASCI Red - TFLOP) 
 * fail to find the max and min function.  Therfore we redefine them 
 * here. 
 */ 
#ifdef max
#undef max
#endif
#define max(a,b) ((a)>(b)) ? (a) : (b);

#ifdef min
#undef min
#endif
#define min(a,b) ((a)<(b)) ? (a) : (b);

// **************************************************************************
// *** Constructor
// **************************************************************************
NOX::Direction::Utils::InexactNewton::
InexactNewton(const NOX::Utils& u, NOX::Parameter::List& directionSublist) :
  printing(0),
  paramsPtr(0),
  predRhs(0),
  stepDir(0),
  userNorm(0)
{
  reset(u, directionSublist);
}

// **************************************************************************
// *** Destructor
// **************************************************************************
NOX::Direction::Utils::InexactNewton::~InexactNewton()
{
  delete predRhs;
  delete stepDir;
}

// **************************************************************************
// *** reset
// **************************************************************************
bool NOX::Direction::Utils::InexactNewton::
reset(const NOX::Utils& u, NOX::Parameter::List& directionSublist)
{
  printing = &u;
  paramsPtr = &directionSublist;

  directionMethod = paramsPtr->getParameter("Method", "Newton");
  
  NOX::Parameter::List& p = paramsPtr->sublist(directionMethod);

  setTolerance = p.getParameter("Set Tolerance in Parameter List", true);

  method = p.getParameter("Forcing Term Method", "Constant");

  if (method == "Constant") {  
    forcingTermMethod = Constant;
    eta_k = p.sublist("Linear Solver").getParameter("Tolerance", 1.0e-4);
  }
  else {
    
    if (method == "Type 1") {
      forcingTermMethod = Type1;
    }
    else if (method == "Type 2") {
      forcingTermMethod = Type2;
    }
    else {
      throwError("reset", "\"Forcing Term Method\" is invalid!");
    }

    eta_min = p.getParameter("Forcing Term Minimum Tolerance", 1.0e-4);  
    eta_max = p.getParameter("Forcing Term Maximum Tolerance", 0.9);
    eta_initial = p.getParameter("Forcing Term Initial Tolerance", 0.01);
    alpha = p.getParameter("Forcing Term Alpha", 1.5);
    gamma = p.getParameter("Forcing Term Gamma", 0.9);
    eta_k = eta_min;
    
    userNorm = 0;
    if (p.isParameter("Forcing Term User Defined Norm")) {

      NOX::Parameter::Arbitrary& arbitrary = 
	const_cast<NOX::Parameter::Arbitrary&>
	(p.getArbitraryParameter("Forcing Term User Defined Norm"));

      userNorm = dynamic_cast<NOX::Parameter::UserNorm*>(&arbitrary);

      if (userNorm == 0) {
	if (printing->isPrintProcessAndType(NOX::Utils::Warning)) {
	  cout << "WARNING: NOX::InexactNewtonUtils::resetForcingTerm() - "
	       << "\"Forcing Term User Defined Norm\" is not of type "
	       << "NOX::Parameter::UserNorm!\n" 
	       << "Defaulting to L-2 Norms!" << endl; 
	}
      }

    }
  }

  return true;
}

// **************************************************************************
// *** computeForcingTerm
// **************************************************************************
double NOX::Direction::Utils::InexactNewton::
computeForcingTerm(const NOX::Abstract::Group& soln,
		   const NOX::Abstract::Group& oldsoln, 
		   int niter,
		   const NOX::Solver::Generic& solver,
		   double eta_last)
{
  const string indent = "       ";

  if (forcingTermMethod == Constant) {
    if (printing->isPrintProcessAndType(NOX::Utils::Details)) {
      cout << indent << "CALCULATING FORCING TERM" << endl;
      cout << indent << "Method: Constant" << endl;
      cout << indent << "Forcing Term: " << eta_k << endl;
    }
    if (setTolerance)
      paramsPtr->sublist(directionMethod).sublist("Linear Solver").
	setParameter("Tolerance", eta_k);

    return eta_k;
  }

  // Get linear solver current tolerance. 
  // NOTE: These values are changing at each nonlinear iteration and 
  // must either be updated from the parameter list each time a compute 
  // is called or supplied during the function call!
  double eta_km1 = 0.0;
  if (eta_last < 0.0)
    eta_km1 = paramsPtr->sublist(directionMethod).
      sublist("Linear Solver").getParameter("Tolerance", 0.0);
  else
    eta_km1 = eta_last;

  // Tolerance may have been adjusted in a line search algorithm so we 
  // have to account for this.
  const NOX::Solver::LineSearchBased* solverPtr = 0;
  solverPtr = dynamic_cast<const NOX::Solver::LineSearchBased*>(&solver);
  if (solverPtr != 0) {
    eta_km1 = 1.0 - solverPtr->getStepSize() * (1.0 - eta_km1);
  }

  if (printing->isPrintProcessAndType(NOX::Utils::Details)) {
    cout << indent << "CALCULATING FORCING TERM" << endl;
    cout << indent << "Method: " << method << endl;
  }


  if (forcingTermMethod == Type1) {
    
    if (niter == 0) {
      
      eta_k = eta_initial;

    }
    else {

      // Return norm of predicted F

      // do NOT use the following lines!! This does NOT account for 
      // line search step length taken.
      // const double normpredf = 0.0;
      // oldsoln.getNormLastLinearSolveResidual(normpredf);
      
      // Create a new vector to be the predicted RHS
      if (predRhs == 0) {
	predRhs = oldsoln.getF().clone(ShapeCopy);
      }
      if (stepDir == 0) {
	stepDir = oldsoln.getF().clone(ShapeCopy);
      }
      
      // stepDir = X - oldX (i.e., the step times the direction)
      stepDir->update(1.0, soln.getX(), -1.0, oldsoln.getX(), 0);
      
      // Compute predRhs = Jacobian * step * dir
      if (!(oldsoln.isJacobian())) {
	if (printing->isPrintProcessAndType(NOX::Utils::Details)) {
	  cout << "WARNING: NOX::InexactNewtonUtils::resetForcingTerm() - "
	       << "Jacobian is out of date! Recomputing Jacobian." << endl;
	}
	const_cast<NOX::Abstract::Group&>(oldsoln).computeJacobian();
      }
      oldsoln.applyJacobian(*stepDir, *predRhs);

      // Compute predRhs = RHSVector + predRhs (this is the predicted RHS)
      predRhs->update(1.0, oldsoln.getF(), 1.0);
      
      // Compute the norms
      double normpredf = 0.0;
      double normf = 0.0;
      double normoldf = 0.0;

      if (userNorm != 0) {
	if (printing->isPrintProcessAndType(NOX::Utils::Details)) {
	  cout << indent << "Forcing Term Norm: " << userNorm->getType()
	       << endl;
	}
	normpredf = userNorm->norm(*predRhs);
	normf = userNorm->norm(soln.getF());
	normoldf = userNorm->norm(oldsoln.getF());
      }
      else {
	if (printing->isPrintProcessAndType(NOX::Utils::Details)) {
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
      if (printing->isPrintProcessAndType(NOX::Utils::Details)) {
	cout << indent << "Residual Norm k-1 =             " 
	     << normoldf << "\n";
	cout << indent << "Residual Norm Linear Model k =  " 
	     << normpredf << "\n";
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
    
  else if (forcingTermMethod == Type2) {  
    
    if (niter == 0) {
      
      eta_k = eta_initial;
      
    }
    else {

      double normf = 0.0;
      double normoldf = 0.0;
      
      if (userNorm != 0) {
	if (printing->isPrintProcessAndType(NOX::Utils::Details)) {
	  cout << indent << "Forcing Term Norm: " << userNorm->getType()
	       << endl;
	}
	normf = userNorm->norm(soln.getF());
	normoldf = userNorm->norm(oldsoln.getF());
      }
      else {
	if (printing->isPrintProcessAndType(NOX::Utils::Details)) {
	  cout << indent << "Forcing Term Norm: Using L-2 Norm."
	       << endl;
	}
	normf = soln.getNormF();
	normoldf = oldsoln.getNormF();
      }  

      const double residual_ratio = normf / normoldf;
      
      eta_k = gamma * pow(residual_ratio, alpha);
      
      // Some output
      if (printing->isPrintProcessAndType(NOX::Utils::Details)) {
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
  
  // Set the new linear solver tolerance
  if (setTolerance) 
    paramsPtr->sublist(directionMethod).sublist("Linear Solver").
      setParameter("Tolerance", eta_k);

  if (printing->isPrintProcessAndType(NOX::Utils::Details)) 
    cout << indent << "Forcing Term: " << eta_k << endl;
  
  return eta_k;
}

// **************************************************************************
// *** throwError
// **************************************************************************
void NOX::Direction::Utils::InexactNewton::
throwError(const string& functionName, const string& errorMsg)
{
    if (printing->isPrintProcessAndType(NOX::Utils::Error))
      cerr << "NOX::InexactNewtonUtils::" << functionName << " - " 
	   << errorMsg << endl;
    throw "NOX Error";
}

// **************************************************************************
// **************************************************************************
// **************************************************************************

#endif
