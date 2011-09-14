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

#include <EpetraExt_LPTrans_From_GraphTrans.h>

#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

namespace EpetraExt {

LinearProblem_GraphTrans::
~LinearProblem_GraphTrans()
{
  if( MatExporter_ ) delete MatExporter_;
  if( VecExporter_ ) delete VecExporter_;
  if( Importer_ ) delete Importer_;

  if( NewProblem_ ) delete NewProblem_;
  if( NewRHS_ ) delete NewRHS_;
  if( NewLHS_ ) delete NewLHS_;
  if( NewMatrix_ ) delete NewMatrix_;
}

LinearProblem_GraphTrans::NewTypeRef
LinearProblem_GraphTrans::
operator()( OriginalTypeRef orig )
{
  OldProblem_ = &orig;
  OldMatrix_ = dynamic_cast<Epetra_CrsMatrix*>( orig.GetMatrix() );
  OldGraph_ = const_cast<Epetra_CrsGraph*>(&OldMatrix_->Graph());
  OldRHS_ = orig.GetRHS();
  OldLHS_ = orig.GetLHS();
  OldRowMap_ = const_cast<Epetra_Map*>(&OldMatrix_->RowMap());

  int ierr = 0;


  if( !OldMatrix_ ) ierr = -2;
  if( !OldRHS_ )    ierr = -3;
  if( !OldLHS_ )    ierr = -4;

  Epetra_CrsGraph & NewGraph = graphTrans_( *OldGraph_ );
  NewMatrix_ = new Epetra_CrsMatrix( Copy, NewGraph );

  Epetra_BlockMap & NewRowMap = const_cast<Epetra_BlockMap&>(NewGraph.RowMap());

  NewRHS_ = new Epetra_MultiVector( NewRowMap, 1 );
  NewLHS_ = new Epetra_MultiVector( NewRowMap, 1 );

  MatExporter_ = new Epetra_Export( *OldRowMap_, NewRowMap );
  VecExporter_ = new Epetra_Export( *OldRowMap_, NewRowMap );
  Importer_ = new Epetra_Import( *OldRowMap_, NewRowMap );

  NewProblem_ = new Epetra_LinearProblem( NewMatrix_, NewLHS_, NewRHS_ );

  return *NewProblem_;
}

bool
LinearProblem_GraphTrans::
fwd()
{
  NewLHS_->Export( *OldLHS_, *VecExporter_, Insert );
  NewRHS_->Export( *OldRHS_, *VecExporter_, Insert );
  NewMatrix_->Export( *OldMatrix_, *MatExporter_, Insert );

  return true;
}

bool
LinearProblem_GraphTrans::
rvs()
{
  OldLHS_->Import( *NewLHS_, *Importer_, Insert );
//  OldRHS_->Import( *NewRHS_, *Importer_, Insert );
//  OldMatrix_->Import( *NewMatrix_, *Importer_, Insert );

  return true;
}

} //namespace EpetraExt

