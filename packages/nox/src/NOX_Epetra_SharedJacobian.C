// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NOX_Epetra_SharedJacobian.H"

using namespace NOX::Epetra;

SharedJacobian::SharedJacobian(Epetra_RowMatrix& j)
{
  jacobian = &j;
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

bool SharedJacobian::isOwner(const Group* grp) const
{
  return (owner == grp);
}

const Epetra_RowMatrix& SharedJacobian::getJacobian() const
{
  return *jacobian;
}


