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

#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_ImplicitWeighting.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_LineSearchBased.H"
#include "NOX_Utils.H"

//#include "/home/rppawlo/Trilinos/packages/epetra/src/Epetra_Vector.h"
//#include "/home/rppawlo/nox/src-epetra/NOX_Epetra_Vector.H"

using namespace NOX::StatusTest;

NormWRMS::NormWRMS(double rtol_, double atol_, double BDFmult_, double tol_,
		   double alpha_, double beta_,
		   bool disable_implicit_weighting) :
  value(0.0),
  rtol(rtol_),
  atolIsScalar(true),
  atol(atol_),
  factor(BDFmult_),
  tolerance(tol_),
  alpha(alpha_),
  computedStepSize(1.0),
  beta(beta_),
  achievedTol(0.0),
  status(Unconverged),
  printCriteria2Info(false),
  printCriteria3Info(false),
  m_disable_implicit_weighting(disable_implicit_weighting)
{

}

NormWRMS::NormWRMS(double rtol_, 
		   const Teuchos::RCP<const NOX::Abstract::Vector>& atolVec_,
		   double BDFmult_, double tol_, double alpha_, double beta_,
		   bool disable_implicit_weighting) :
  value(0.0),
  rtol(rtol_),
  atolIsScalar(false),
  atol(0.0),
  atolVec(atolVec_),
  factor(BDFmult_),
  tolerance(tol_),
  alpha(alpha_),
  computedStepSize(1.0),
  beta(beta_),
  achievedTol(0.0),
  status(Unconverged),
  printCriteria2Info(false),
  printCriteria3Info(false),
  m_disable_implicit_weighting(disable_implicit_weighting)
{

}

NormWRMS::~NormWRMS()
{

}

StatusType NormWRMS::
checkStatus(const NOX::Solver::Generic& problem, 
	    NOX::StatusTest::CheckType checkType)
{
  if (checkType == NOX::StatusTest::None) {
    status = Unevaluated;
    value = 1.0e+12;
    return status;
  }

  status = Unconverged;

  const Abstract::Group& soln = problem.getSolutionGroup();
  const Abstract::Group& oldsoln = problem.getPreviousSolutionGroup();
  const Abstract::Vector& x = soln.getX();
  
  // On the first iteration, the old and current solution are the same so
  // we should return the test as unconverged until there is a valid 
  // old solution (i.e. the number of iterations is greater than zero).
  int niters = problem.getNumIterations();
  if (niters == 0) 
  {
    status = Unconverged;
    value = 1.0e+12;
    return status;
  } 

  // **** Begin check for convergence criteria #1 ****

  // Create the working vectors if this is the first time this
  // operator is called.
  if (Teuchos::is_null(u))
    u = x.clone(NOX::ShapeCopy);
  if (Teuchos::is_null(v))
    v = x.clone(NOX::ShapeCopy);
  
  // Create the weighting vector u = RTOL |x| + ATOL
  // |x| is evaluated at the old time step
  v->abs(oldsoln.getX());
  if (atolIsScalar) 
  {    
    u->init(1.0);
    u->update(rtol, *v, atol);
  }
  else 
  {
    u->update(rtol, *v, 1.0, *atolVec, 0.0);
  }

  // v = 1/u (elementwise)
  v->reciprocal(*u);

  // u = x - oldx (i.e., the update)
  u->update(1.0, x, -1.0, oldsoln.getX(), 0.0);

  // u = Cp * u @ v (where @ represents an elementwise multiply)
  u->scale(*v);

  // Turn off implicit scaling of norm if the vector supports it
  Teuchos::RCP<NOX::Abstract::ImplicitWeighting> iw_u;
  iw_u = Teuchos::rcp_dynamic_cast<NOX::Abstract::ImplicitWeighting>(u,false);
  bool saved_status = false;
  if (nonnull(iw_u) && m_disable_implicit_weighting) {
    saved_status = iw_u->getImplicitWeighting();
    iw_u->setImplicitWeighting(false);
  }

  // tmp = factor * sqrt (u * u / N)
  value = u->norm() * factor / sqrt(static_cast<double>(u->length()));

  // Set the implicit scaling back to original value
  if (nonnull(iw_u) && m_disable_implicit_weighting)
    iw_u->setImplicitWeighting(saved_status);

  StatusType status1 = Unconverged;
  if (value < tolerance)
    status1 = Converged;


  // **** Begin check for convergence criteria #2 ****
  StatusType status2 = Unconverged;
  
  // Determine if the Generic solver is a LineSearchBased solver
  // If it is not then return a "Converged" status
  const Solver::Generic* test = 0;
  test = dynamic_cast<const Solver::LineSearchBased*>(&problem);
  if (test == 0) 
  {
    status2 = Converged; 
  }
  else 
  {
    printCriteria2Info = true;
    computedStepSize = 
      (dynamic_cast<const Solver::LineSearchBased*>(&problem))->getStepSize();
    
    if (computedStepSize >= alpha)
      status2 = Converged;
  }

  // **** Begin check for convergence criteria #3 ****
  
  // First time through, make sure the output parameter list exists.
  // Since the list is const, a sublist call to a non-existent sublist 
  // throws an error.  Therefore we have to check the existence of each 
  // sublist before we call it.
  const Teuchos::ParameterList& p = problem.getList();
  if (niters == 1) {
    if (p.isSublist("Direction")) {
      if (p.sublist("Direction").isSublist("Newton")) {
	if (p.sublist("Direction").sublist("Newton").isSublist("Linear Solver")) {
	  if (p.sublist("Direction").sublist("Newton").sublist("Linear Solver").isSublist("Output")) {

	    const Teuchos::ParameterList& list = p.sublist("Direction").sublist("Newton").sublist("Linear Solver").sublist("Output");
	    
	    if (Teuchos::isParameterType<double>(list, "Achieved Tolerance")) {
	      
	      printCriteria3Info = true;
	      
	      
	    }
	  }
	}
      }
    }
  }
  
  StatusType status3 = Converged;
  if (printCriteria3Info) {
    achievedTol = const_cast<Teuchos::ParameterList&>(problem.getList()).
      sublist("Direction").sublist("Newton").sublist("Linear Solver").
      sublist("Output").get("Achieved Tolerance", -1.0);
    status3 = (achievedTol <= beta) ? Converged : Unconverged;
  }


  // Determine status of test
  if ((status1 == Converged) && 
      (status2 == Converged) &&
      (status3 == Converged))
    status = Converged;
  
  return status;
}

StatusType NormWRMS::getStatus() const
{
  return status;
}


std::ostream& NormWRMS::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "WRMS-Norm = " << Utils::sciformat(value, 3) << " < " << tolerance;
  if (printCriteria2Info) {
    stream << "\n";
    for (int j = 0; j < indent + 13; j ++)
      stream << ' ';
    stream << "(Min Step Size:  " << Utils::sciformat(computedStepSize, 3) << " >= " << alpha << ")";
  }
  if (printCriteria3Info) {
    stream << "\n";
    for (int j = 0; j < indent+ 13; j ++)
      stream << ' ';
    stream << "(Max Lin Solv Tol:  " << Utils::sciformat(achievedTol, 3) << " < " << beta << ")";
  }
  stream << std::endl;
  return stream;
}


double NormWRMS::getNormWRMS() const
{
  return value;
}   

double NormWRMS::getTolerance() const
{
  return tolerance;
}

double NormWRMS::getRTOL() const
{
  return rtol;
}

double NormWRMS::getATOL() const
{
  if (atolIsScalar)
    return atol;
  
  return (-1.0);
}

double NormWRMS::getBDFMultiplier() const
{
  return factor;
}
 
double NormWRMS::getAlpha() const
{
  return alpha;
}
 
double NormWRMS::getBeta() const
{
  return beta;
}

bool NormWRMS::getDisableImplicitWeighting() const
{
  return m_disable_implicit_weighting;
}
