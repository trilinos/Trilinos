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

#include "NOX_Epetra_SharedJacobian.H" // class definition
// External include from Trilinos
#include "Epetra_CrsMatrix.h"

using namespace NOX::Epetra;

SharedJacobian::SharedJacobian(Epetra_RowMatrix& j)
{
  jacobian = &j;
}

SharedJacobian::SharedJacobian(Epetra_RowMatrix& j, Epetra_RowMatrix& p)
{
  jacobian = &j;
  prec = &p;
}

SharedJacobian::~SharedJacobian()
{
  // Do nothing
}

Epetra_RowMatrix& SharedJacobian::getJacobian(const Group* newowner)
{
  owner = newowner;
  return *jacobian;
}

const Epetra_RowMatrix& SharedJacobian::getJacobian() const
{
  return *jacobian;
}

bool SharedJacobian::isOwner(const Group* grp) const
{
  return (owner == grp);
}

Epetra_RowMatrix& SharedJacobian::getPrec(const Group* newowner)
{
  owner = newowner;
  return *prec;
}

const Epetra_RowMatrix& SharedJacobian::getPrec() const
{
  return *prec;
}

bool SharedJacobian::setJacobian(Epetra_RowMatrix& j)
{
  jacobian = & j;
}
