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

#include "NOX_Status_WRMS.H"
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

//#include "/home/rppawlo/Trilinos/packages/epetra/src/Epetra_Vector.h"
//#include "/home/rppawlo/nox/src-epetra/NOX_Epetra_Vector.H"

using namespace NOX::Status;

WRMS::WRMS(double a, double b, double c) :
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

WRMS::WRMS(double a, Abstract::Vector& b, double c) :
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

WRMS::~WRMS()
{
  delete atoli;
  delete u;
  delete v;
}

StatusType WRMS::operator()(const Solver::Generic& problem)
{
  status = Unconverged;

  const Abstract::Group& soln = problem.getSolutionGroup();
  const Abstract::Group& oldsoln = problem.getPreviousSolutionGroup();
  const Abstract::Vector& x = soln.getX();
  
  // Create the working vectors if this is the first time this
  // operator is called.
  if (u == 0)
    u = x.clone(NOX::CopyShape);
  if (v == 0)
    v = x.clone(NOX::CopyShape);
  
  // Create the weighting vector u = RTOL |x| + ATOL
  if (atolIsScalar) {    
    u->init(1.0);
    v->abs(oldsoln.getX());
    u->update(rtol, *v, atol);
  }
  else {
    *u = *atoli;
    v->abs(oldsoln.getX());
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

  /*
  cout << "Before Dynamic cast" << endl;
  NOX::Epetra::Vector* testVec = dynamic_cast<NOX::Epetra::Vector*>(u);
  Epetra_Vector* testVec2 = dynamic_cast<Epetra_Vector*>(&(testVec->getEpetraVector()));
  cout << "After elementwise multiply" << endl;
  testVec2->Print(cout);
  */

  // Compute the sum of u^2 then divide by vector length: tmp = u*u/N
  double tmp = u->dot(*u)/((double) u->length());

  // Finally, compute the ratio by taking the sqrt
  double ratio = sqrt(tmp);

  //cout << "WRMS = " << ratio << endl;
  
  if (ratio < 1.0)
    status = Converged;
  
  return status;
}

ostream& WRMS::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  if (atolIsScalar) {
    stream << "WRMS with RTOL = " << Utils::sci(rtol) 
	   << " and ATOL = " << Utils::sci(atol) << " < 1";
  }
  else {
    stream << "WRMS with RTOL = " << Utils::sci(rtol) 
	   << " and ATOL vector < 1";
  }
  stream << endl;
  return stream;
}
