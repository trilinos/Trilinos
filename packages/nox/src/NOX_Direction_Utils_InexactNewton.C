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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_Direction_Utils_InexactNewton.H" // class definition

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_LineSearchBased.H"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"
#include "Teuchos_ParameterList.hpp"

// **************************************************************************
// *** Constructor
// **************************************************************************
NOX::Direction::Utils::InexactNewton::
InexactNewton(const Teuchos::RCP<NOX::GlobalData>& gd, 
	      Teuchos::ParameterList& directionSublist) :
  paramsPtr(0)
{
  reset(gd, directionSublist);
}

// **************************************************************************
// *** Destructor
// **************************************************************************
NOX::Direction::Utils::InexactNewton::~InexactNewton()
{
}

// **************************************************************************
// *** reset
// **************************************************************************
bool NOX::Direction::Utils::InexactNewton::
reset(const Teuchos::RCP<NOX::GlobalData>& gd, 
      Teuchos::ParameterList& directionSublist)
{
  globalDataPtr = gd;
  printing = gd->getUtils();
  paramsPtr = &directionSublist;

  directionMethod = paramsPtr->get("Method", "Newton");
  
  Teuchos::ParameterList& p = paramsPtr->sublist(directionMethod);

  setTolerance = p.get("Set Tolerance in Parameter List", true);

  method = p.get("Forcing Term Method", "Constant");

  if (method == "Constant") {  
    forcingTermMethod = Constant;
    eta_k = p.sublist("Linear Solver").get("Tolerance", 1.0e-4);
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

    eta_min = p.get("Forcing Term Minimum Tolerance", 1.0e-4);  
    eta_max = p.get("Forcing Term Maximum Tolerance", 0.9);
    eta_initial = p.get("Forcing Term Initial Tolerance", 0.01);
    alpha = p.get("Forcing Term Alpha", 1.5);
    gamma = p.get("Forcing Term Gamma", 0.9);
    eta_k = eta_min;
    
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
  const std::string indent = "       ";

  if (forcingTermMethod == Constant) {
    if (printing->isPrintType(NOX::Utils::Details)) {
      printing->out() << indent << "CALCULATING FORCING TERM" << std::endl;
      printing->out() << indent << "Method: Constant" << std::endl;
      printing->out() << indent << "Forcing Term: " << eta_k << std::endl;
    }
    if (setTolerance)
      paramsPtr->sublist(directionMethod).sublist("Linear Solver").
	set("Tolerance", eta_k);

    return eta_k;
  }

  // Get linear solver current tolerance. 
  // NOTE: These values are changing at each nonlinear iteration and 
  // must either be updated from the parameter list each time a compute 
  // is called or supplied during the function call!
  double eta_km1 = 0.0;
  if (eta_last < 0.0)
    eta_km1 = paramsPtr->sublist(directionMethod).
      sublist("Linear Solver").get("Tolerance", 0.0);
  else
    eta_km1 = eta_last;

  // Tolerance may have been adjusted in a line search algorithm so we 
  // have to account for this.
  const NOX::Solver::LineSearchBased* solverPtr = 0;
  solverPtr = dynamic_cast<const NOX::Solver::LineSearchBased*>(&solver);
  if (solverPtr != 0) {
    eta_km1 = 1.0 - solverPtr->getStepSize() * (1.0 - eta_km1);
  }

  if (printing->isPrintType(NOX::Utils::Details)) {
    printing->out() << indent << "CALCULATING FORCING TERM" << std::endl;
    printing->out() << indent << "Method: " << method << std::endl;
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
      if (Teuchos::is_null(predRhs)) {
	predRhs = oldsoln.getF().clone(ShapeCopy);
      }
      if (Teuchos::is_null(stepDir)) {
	stepDir = oldsoln.getF().clone(ShapeCopy);
      }
      
      // stepDir = X - oldX (i.e., the step times the direction)
      stepDir->update(1.0, soln.getX(), -1.0, oldsoln.getX(), 0);
      
      // Compute predRhs = Jacobian * step * dir
      if (!(oldsoln.isJacobian())) {
	if (printing->isPrintType(NOX::Utils::Details)) {
	  printing->out() << "WARNING: NOX::InexactNewtonUtils::resetForcingTerm() - "
	       << "Jacobian is out of date! Recomputing Jacobian." << std::endl;
	}
	const_cast<NOX::Abstract::Group&>(oldsoln).computeJacobian();
      }
      oldsoln.applyJacobian(*stepDir, *predRhs);

      // Compute predRhs = RHSVector + predRhs (this is the predicted RHS)
      predRhs->update(1.0, oldsoln.getF(), 1.0);
      
      // Compute the norms
      double normpredf = predRhs->norm();
      double normf = soln.getNormF();
      double normoldf = oldsoln.getNormF();

      if (printing->isPrintType(NOX::Utils::Details)) {
	printing->out() << indent << "Forcing Term Norm: Using L-2 Norm."
			<< std::endl;
      }

      // Compute forcing term
      eta_k = fabs(normf - normpredf) / normoldf;
      
      // Some output
      if (printing->isPrintType(NOX::Utils::Details)) {
	printing->out() << indent << "Residual Norm k-1 =             " 
	     << normoldf << "\n";
	printing->out() << indent << "Residual Norm Linear Model k =  " 
	     << normpredf << "\n";
	printing->out() << indent << "Residual Norm k =               " << normf << "\n";
	printing->out() << indent << "Calculated eta_k (pre-bounds) = " << eta_k << std::endl;
      }
      
      // Impose safeguard and constraints ...
      const double tmp = (1.0 + sqrt(5.0)) / 2.0;
      const double eta_km1_alpha = pow(eta_km1, tmp);
      if (eta_km1_alpha > 0.1) 
	eta_k = NOX_MAX(eta_k, eta_km1_alpha);
      eta_k = NOX_MAX(eta_k, eta_min);
      eta_k = NOX_MIN(eta_max, eta_k);
    }
  }
    
  else if (forcingTermMethod == Type2) {  
    
    if (niter == 0) {
      
      eta_k = eta_initial;
      
    }
    else {

      double normf = soln.getNormF();
      double normoldf = oldsoln.getNormF();
      
      if (printing->isPrintType(NOX::Utils::Details)) {
	printing->out() << indent << "Forcing Term Norm: Using L-2 Norm."
			<< std::endl;
      }

      const double residual_ratio = normf / normoldf;
      
      eta_k = gamma * pow(residual_ratio, alpha);
      
      // Some output
      if (printing->isPrintType(NOX::Utils::Details)) {
	printing->out() << indent << "Residual Norm k-1 =             " 
			<< normoldf << "\n";
	printing->out() << indent << "Residual Norm k =               " 
			<< normf << "\n";
	printing->out() << indent << "Calculated eta_k (pre-bounds) = " 
			<< eta_k << std::endl;
      }
      
      // Impose safeguard and constraints ... 
      const double eta_k_alpha = gamma * pow(eta_km1, alpha);
      if (eta_k_alpha > 0.1) 
	eta_k = NOX_MAX(eta_k, eta_k_alpha);
      eta_k = NOX_MAX(eta_k, eta_min);
      eta_k = NOX_MIN(eta_max, eta_k);
    }
    
  }
  
  // Set the new linear solver tolerance
  if (setTolerance) 
    paramsPtr->sublist(directionMethod).sublist("Linear Solver").
      set("Tolerance", eta_k);

  if (printing->isPrintType(NOX::Utils::Details)) 
    printing->out() << indent << "Forcing Term: " << eta_k << std::endl;
  
  return eta_k;
}

// **************************************************************************
// *** throwError
// **************************************************************************
void NOX::Direction::Utils::InexactNewton::
throwError(const std::string& functionName, const std::string& errorMsg)
{
    if (printing->isPrintType(NOX::Utils::Error))
      printing->err() << "NOX::InexactNewtonUtils::" << functionName << " - " 
	   << errorMsg << std::endl;
    throw "NOX Error";
}

// **************************************************************************
// **************************************************************************
// **************************************************************************
