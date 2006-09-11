// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_LineSearchBased.H"
#include "NOX_Utils.H"

//#include "/home/rppawlo/Trilinos/packages/epetra/src/Epetra_Vector.h"
//#include "/home/rppawlo/nox/src-epetra/NOX_Epetra_Vector.H"

using namespace NOX::StatusTest;

NormWRMS::NormWRMS(double rtol_, double atol_, double BDFmult_, double tol_,
		   double alpha_, double beta_) :
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
  printCriteria3Info(false)
{

}

NormWRMS::NormWRMS(double rtol_, 
		   const Teuchos::RefCountPtr<NOX::Abstract::Vector>& atolVec_,
		   double BDFmult_, double tol_, double alpha_, double beta_) :
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
  printCriteria3Info(false)
{

}

NormWRMS::~NormWRMS()
{

}

StatusType NormWRMS::
checkStatus(const Solver::Generic& problem, 
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

  // tmp = factor * sqrt (u * u / N)
  value = u->norm() * factor / sqrt(static_cast<double>(u->length()));

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


ostream& NormWRMS::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "WRMS-Norm = " << Utils::sciformat(value, 3) << " < " << tolerance;
  if (printCriteria2Info) {
    stream << "\n";
    for (int j = 0; j < indent + 13; j ++, 3)
      stream << ' ';
    stream << "(Min Step Size:  " << Utils::sciformat(computedStepSize, 3) << " >= " << alpha << ")";
  }
  if (printCriteria3Info) {
    stream << "\n";
    for (int j = 0; j < indent+ 13; j ++)
      stream << ' ';
    stream << "(Max Lin Solv Tol:  " << Utils::sciformat(achievedTol, 3) << " < " << beta << ")";
  }
  stream << endl;
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
