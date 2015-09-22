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

#include "Ifpack_OverlapFactor.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

//==============================================================================
Ifpack_OverlapFactor::Ifpack_OverlapFactor(const Ifpack_OverlapGraph * OverlapGraph) 
  : Factored_(false),
    Allocated_(false),
    ValuesInitialized_(false),
    OverlapGraph_(OverlapGraph),
    UserMatrix_(0)
{
}
//==============================================================================
Ifpack_OverlapFactor::Ifpack_OverlapFactor(const Epetra_RowMatrix * UserMatrix) 
  : Factored_(false),
    Allocated_(false),
    ValuesInitialized_(false),
    OverlapGraph_(0),
    UserMatrix_(UserMatrix)
{
}
//==============================================================================
Ifpack_OverlapFactor::Ifpack_OverlapFactor(const Ifpack_OverlapFactor & Source) 
  : Factored_(Source.Factored_),
    Allocated_(Source.Allocated_),
    ValuesInitialized_(Source.ValuesInitialized_),
    OverlapGraph_(Source.OverlapGraph_),
    UserMatrix_(Source.UserMatrix_)
{
}
//==============================================================================
int Ifpack_OverlapFactor::InitValues(const Epetra_RowMatrix * UserMatrix) {
  

  if (OverlapGraph_!=0) {

    Epetra_CrsMatrix * CrsMatrix = dynamic_cast<Epetra_CrsMatrix *>(UserMatrix);
    if (CrsMatrix!=0) 
  if (!Allocated()) EPETRA_CHK_ERR(-1); //Must be allocated
  if (ValuesInitialized()) EPETRA_CHK_ERR(1); // Values already init'ed, warn caller
  
  EPETRA_CHK_ERR(DerivedFactor()); // Call Derived class factorization
  SetValuesInitialized(false);
  SetFactored(true);
  return(0);
}
//==============================================================================
int Ifpack_OverlapFactor::Factor() {
  
  if (!ValuesInitialized()) EPETRA_CHK_ERR(-1); // Values must be initialized
  if (Factored()) EPETRA_CHK_ERR(1); // Return with a warning that factor already done
  
  EPETRA_CHK_ERR(DerivedFactor()); // Call Derived class factorization
  SetValuesInitialized(false);
  SetFactored(true);
  return(0);
}
