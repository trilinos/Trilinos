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

#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

//#include "/home/rppawlo/Trilinos/packages/epetra/src/Epetra_Vector.h"
//#include "/home/rppawlo/nox/src-epetra/NOX_Epetra_Vector.H"

using namespace NOX::StatusTest;

NormWRMS::NormWRMS(double a, double b, double c) :
  atolIsScalar(true),
  atoli(0),
  factor(c),
  u(0),
  v(0)
{
  rtol = a;
  atol = b;
  status = Unconverged;
}

NormWRMS::NormWRMS(double a, Abstract::Vector& b, double c) :
  atolIsScalar(false),
  atoli(0),
  factor(c),
  u(0),
  v(0)
{
  rtol = a;
  atoli = b.clone();
  status = Unconverged;
}

NormWRMS::~NormWRMS()
{
  delete atoli;
  delete u;
  delete v;
}

StatusType NormWRMS::checkStatus(const Solver::Generic& problem)
{
  status = Unconverged;

  const Abstract::Group& soln = problem.getSolutionGroup();
  const Abstract::Group& oldsoln = problem.getPreviousSolutionGroup();
  const Abstract::Vector& x = soln.getX();
  
  // On the first iteration, the old and current solution are the same so
  // we should return the test as unconverged until there is a valid 
  // old solution (i.e. the number of iterations is greater than zero).
  int niters = problem.getNumIterations();
  if (niters == 0) {
    status = Unconverged;
    value = -1.0;
    return status;
  } 

  // Create the working vectors if this is the first time this
  // operator is called.
  if (u == 0)
    u = x.clone(NOX::ShapeCopy);
  if (v == 0)
    v = x.clone(NOX::ShapeCopy);
  
  // Create the weighting vector u = RTOL |x| + ATOL
  // |x| is evaluated at the old time step
  v->abs(oldsoln.getX());
  if (atolIsScalar) {    
    u->init(1.0);
    u->update(rtol, *v, atol);
  }
  else {
    *u = *atoli;
    u->update(rtol, *v, 1.0);
  }

  // v = 1/u (elementwise)
  v->reciprocal(*u);

  // u = x - oldx (i.e., the update)
  u->update(1.0, soln.getX(), -1.0, oldsoln.getX(), 0.0);

  // u = Cp * u @ v (where @ represents an elementwise multiply)
  //u->multiply(factor, *u, *v, 0.0);
  u->scale(*v);
  u->scale(factor);

  // Compute the sum of u^2 then divide by vector length: tmp = u*u/N
  double tmp = u->dot(*u)/((double) u->length());

  // Finally, compute the WRMS norm value by taking the sqrt
  value = sqrt(tmp);

  if (value < 1.0)
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
    stream << "WRMS-Norm = " << Utils::sci(value) << " < 1.0";
  stream << endl;
  return stream;
}
