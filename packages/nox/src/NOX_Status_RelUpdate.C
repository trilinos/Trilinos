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

#include "NOX_Status_RelUpdate.H"
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

using namespace NOX::Status;

RelUpdate::RelUpdate(double a, double b) :
  u(NULL),
  v(NULL)
{
  epsilon_r = a;
  epsilon_a = b;
  status = Unconverged;
}

RelUpdate::~RelUpdate()
{
  delete u;
  delete v;
}

StatusType RelUpdate::operator()(const Solver::Generic& problem)
{
  status = Unconverged;

  const Abstract::Group& soln = problem.getSolutionGroup();
  const Abstract::Group& oldsoln = problem.getPreviousSolutionGroup();
  const Abstract::Vector& x = soln.getX();
  
  // Create the working vectors if this is the first time this
  // operator is called.
  if (u == NULL)
    u = x.clone(NOX::CopyShape);
  if (v == NULL)
    v = x.clone(NOX::CopyShape);
  
  // u = 1 (elementwise)
  u->init(1.0);

  // v = |x| (elementwise)
  v->abs(x);

  // u = epsilon_r * v + epsilon_a * u
  u->update(epsilon_r, *v, epsilon_a);

  // v = 1/u (elementwise)
  v->reciprocal(*u);

  // u = x - oldx (i.e., the update)
  u->update(1.0, soln.getX(), -1.0, oldsoln.getX(), 0);
  
  // u = abs(u)
  u->abs(*u);

  // Finally, compute the ratio
  double ratio = u->dot(*v) / u->length();
  
  if (ratio < 1)
    status = Converged;
  
  return status;
}

ostream& RelUpdate::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Relative Update with e_r = " << Utils::sci(epsilon_r) 
	 << " and e_a = " << Utils::sci(epsilon_a) << " < 1";
  stream << endl;
  return stream;
}
