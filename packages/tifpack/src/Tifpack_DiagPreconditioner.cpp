/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_DiagPreconditioner.hpp"
#include "Tpetra_BlockMap.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Comm.hpp"

// ============================================================================ 
Tifpack_DiagPreconditioner::
Tifpack_DiagPreconditioner(const Tpetra_Map& DomainMap,
                          const Tpetra_Map& RangeMap,
                          const Tpetra_Vector& diag) : 
  UseTranspose_(false),
  DomainMap_(DomainMap),
  RangeMap_(RangeMap),
  diag_(diag)
{ }

// ============================================================================ 
Tifpack_DiagPreconditioner::~Tifpack_DiagPreconditioner()
{ }


// ============================================================================ 
int Tifpack_DiagPreconditioner::Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  TIFPACK_RETURN(-1); // not defined
}

// ============================================================================ 
 int Tifpack_DiagPreconditioner::ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  if (X.NumVectors() != Y.NumVectors())
    TIFPACK_CHK_ERR(-1); 

  for (int v = 0; v < X.NumVectors(); ++v)
    for (int i = 0; i < X.MyLength(); ++i)
      Y[v][i] = diag_[i] * X[v][i];
  ///Y.ReciprocalMultiply(1.0, diag_, X, 0.0);

  return(0);
}
