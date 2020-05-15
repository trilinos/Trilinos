/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
int Ifpack_DiagPreconditioner::Apply(const Epetra_MultiVector& /* X */, Epetra_MultiVector& /* Y */) const
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
