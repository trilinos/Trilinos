// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "FEApp_ConstantDirichletBC.hpp"

FEApp::ConstantDirichletBC::ConstantDirichletBC(unsigned int dof_GID, 
						double value) :
  dof(dof_GID),
  val(value)
{
}

FEApp::ConstantDirichletBC::~ConstantDirichletBC()
{
}

void
FEApp::ConstantDirichletBC::applyResidual(const Epetra_Vector& x, 
					  Epetra_Vector& f) const
{
  // Get map
  const Epetra_BlockMap& map = x.Map();

  if (map.MyGID(dof)) {
    int lid = map.LID(dof);

    // Set residual to f = x - a
    f[lid] = x[lid] - val;
  }
}

void
FEApp::ConstantDirichletBC::applyJacobian(const Epetra_Vector& x, 
					  Epetra_Vector& f,
					  Epetra_CrsMatrix& jac) const
{
  // Get map
  const Epetra_BlockMap& map = x.Map();

  if (map.MyGID(dof)) {
    int lid = map.LID(dof);

    // Set residual to f = x - a
    f[lid] = x[lid] - val;

    // Get view of row
    double *row_view;
    double done = 1.0;
    int num_entries;
    jac.ExtractGlobalRowView(dof, num_entries, row_view);

    // Set all entries to zero
    for (int k=0; k<num_entries; k++)
      row_view[k] = 0.0;

    // Set diagonal element to 1
    int idof = dof;  // remove const
    jac.ReplaceGlobalValues(idof, 1, &done, &idof);
  }
}
