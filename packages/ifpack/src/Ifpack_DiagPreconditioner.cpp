/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_DiagPreconditioner.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"

// ============================================================================ 
Ifpack_DiagPreconditioner::
Ifpack_DiagPreconditioner(const Epetra_Map& DomainMap,
                          const Epetra_Map& RangeMap,
                          const Epetra_Vector& diag) : 
  UseTranspose_(false),
  DomainMap_(DomainMap),
  RangeMap_(RangeMap),
  diag_(diag)
{ }

// ============================================================================ 
Ifpack_DiagPreconditioner::~Ifpack_DiagPreconditioner()
{ }


// ============================================================================ 
int Ifpack_DiagPreconditioner::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_RETURN(-1); // not defined
}

// ============================================================================ 
 int Ifpack_DiagPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-1); 

  for (int v = 0; v < X.NumVectors(); ++v)
    for (int i = 0; i < X.MyLength(); ++i)
      Y[v][i] = diag_[i] * X[v][i];
  ///Y.ReciprocalMultiply(1.0, diag_, X, 0.0);

  return(0);
}
