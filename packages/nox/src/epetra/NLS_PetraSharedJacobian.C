// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_PetraSharedJacobian.H"

NLS_PetraSharedJacobian::NLS_PetraSharedJacobian(Epetra_RowMatrix& j)
{
  jacobian = &j;
}

NLS_PetraSharedJacobian::~NLS_PetraSharedJacobian()
{
  // Do nothing
}

Epetra_RowMatrix& NLS_PetraSharedJacobian::getJacobian(const NLS_PetraGroup* newowner)
{
  owner = newowner;
  return *jacobian;
}


void NLS_PetraSharedJacobian::takeOwnership(const NLS_PetraGroup* newowner)
{
  owner = newowner;
}

bool NLS_PetraSharedJacobian::isOwner(const NLS_PetraGroup* grp) const
{
  return (owner == grp);
}

const Epetra_RowMatrix& NLS_PetraSharedJacobian::getJacobian() const
{
  return *jacobian;
}


