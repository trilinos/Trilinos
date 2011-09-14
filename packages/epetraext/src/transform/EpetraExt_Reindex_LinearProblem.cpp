//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#include <EpetraExt_Reindex_LinearProblem.h>

#include <EpetraExt_Reindex_CrsMatrix.h>
#include <EpetraExt_Reindex_MultiVector.h>

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Map.h>

namespace EpetraExt {

LinearProblem_Reindex::
~LinearProblem_Reindex()
{
  if( newObj_ ) delete newObj_;

  if( MatTrans_ ) delete MatTrans_;
  if( LHSTrans_ ) delete LHSTrans_;
  if( RHSTrans_ ) delete RHSTrans_;

  if( NewRowMapOwned_ ) delete NewRowMap_;
}

LinearProblem_Reindex::NewTypeRef
LinearProblem_Reindex::
operator()( OriginalTypeRef orig )
{
  Epetra_CrsMatrix * OldMatrix = dynamic_cast<Epetra_CrsMatrix*>( orig.GetMatrix() );
  Epetra_MultiVector * OldRHS = orig.GetRHS();
  Epetra_MultiVector * OldLHS = orig.GetLHS();
  const Epetra_BlockMap & OldRowMap = OldMatrix->Map();

  //If new map not passed in create one with lex ordering
  if( !NewRowMap_ )
  {
    int NumMyElements = OldRowMap.NumMyElements();
    int NumGlobalElements = OldRowMap.NumGlobalElements();

    NewRowMap_ = new Epetra_Map( NumGlobalElements, NumMyElements, 0, OldRowMap.Comm() );
    NewRowMapOwned_ = true;
  }

  MatTrans_ = new CrsMatrix_Reindex( *NewRowMap_ );
  LHSTrans_ = new MultiVector_Reindex( *NewRowMap_ );
  RHSTrans_ = new MultiVector_Reindex( *NewRowMap_ );

  Epetra_CrsMatrix * NewMatrix = &((*MatTrans_)( *OldMatrix ));
  Epetra_MultiVector * NewLHS = &((*LHSTrans_)( *OldLHS ));
  Epetra_MultiVector * NewRHS = &((*RHSTrans_)( *OldRHS ));

  newObj_ = new Epetra_LinearProblem( NewMatrix, NewLHS, NewRHS );

  return *newObj_;
}

} //namespace EpetraExt

