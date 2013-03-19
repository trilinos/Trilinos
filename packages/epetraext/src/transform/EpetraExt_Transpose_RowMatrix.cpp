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

#include <EpetraExt_Transpose_RowMatrix.h>

#include <Epetra_Export.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>

#include <Epetra_Comm.h>    //DEBUG

#include <vector>

namespace EpetraExt {

RowMatrix_Transpose::
~RowMatrix_Transpose()
{
  if( TransposeMatrix_ ) delete TransposeMatrix_;
  if( TransposeExporter_ ) delete TransposeExporter_;

  if( !OrigMatrixIsCrsMatrix_ )
  {
    delete [] Indices_;
    delete [] Values_;
  }

  for( int i = 0; i < NumMyCols_; ++i )
    if( TransNumNz_[i] )
    {
      delete [] TransIndices_[i];
      delete [] TransValues_[i];
    }

  delete [] TransNumNz_;
  delete [] TransIndices_;
  delete [] TransValues_;
  delete [] TransMyGlobalEquations_;
}

//=========================================================================
RowMatrix_Transpose::NewTypeRef
RowMatrix_Transpose::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  int i, j, err;

  if( !TransposeRowMap_ )
  {
    if( IgnoreNonLocalCols_ )
      TransposeRowMap_ = (Epetra_Map *) &(orig.OperatorRangeMap()); // Should be replaced with refcount =
    else
      TransposeRowMap_ = (Epetra_Map *) &(orig.OperatorDomainMap()); // Should be replaced with refcount =
  }

  // This routine will work for any RowMatrix object, but will attempt cast the matrix to a CrsMatrix if
  // possible (because we can then use a View of the matrix and graph, which is much cheaper).

  // First get the local indices to count how many nonzeros will be in the 
  // transpose graph on each processor
  Epetra_CrsMatrix * OrigCrsMatrix = dynamic_cast<Epetra_CrsMatrix*>(&orig);

  OrigMatrixIsCrsMatrix_ = (OrigCrsMatrix!=0); // If this pointer is non-zero, the cast to CrsMatrix worked

  NumMyRows_ = orig.NumMyRows();
  NumMyCols_ = orig.NumMyCols();
  TransNumNz_ = new int[NumMyCols_];
  TransIndices_ = new int*[NumMyCols_];
  TransValues_ = new double*[NumMyCols_];
  TransMyGlobalEquations_ = new int[NumMyCols_];

  int NumIndices;

  if (OrigMatrixIsCrsMatrix_)
  {
    const Epetra_CrsGraph & OrigGraph = OrigCrsMatrix->Graph(); // Get matrix graph

    for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0;
    for (i=0; i<NumMyRows_; i++)
    {
      err = OrigGraph.ExtractMyRowView(i, NumIndices, Indices_); // Get view of ith row
      if (err != 0) throw OrigGraph.ReportError("ExtractMyRowView failed",err);
      for (j=0; j<NumIndices; j++) ++TransNumNz_[Indices_[j]];
    }
  }
  else // Original is not a CrsMatrix
  {
    MaxNumEntries_ = 0;
    int NumEntries;
    for (i=0; i<NumMyRows_; i++)
    {
      orig.NumMyRowEntries(i, NumEntries);
      MaxNumEntries_ = EPETRA_MAX(MaxNumEntries_, NumEntries);
    }
    Indices_ = new int[MaxNumEntries_];
    Values_ = new double[MaxNumEntries_];

    for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0;
    for (i=0; i<NumMyRows_; i++)
    {
      err = orig.ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_); 
      if (err != 0) {
        std::cerr << "ExtractMyRowCopy failed."<<std::endl;
        throw err;
      }
      for (j=0; j<NumIndices; j++) ++TransNumNz_[Indices_[j]];
    }
  }

  // Most of remaining code is common to both cases
  for(i=0; i<NumMyCols_; i++)
  {
    NumIndices = TransNumNz_[i];
    if (NumIndices>0)
    {
      TransIndices_[i] = new int[NumIndices];
      TransValues_[i] = new double[NumIndices];
    }
  }

  // Now copy values and global indices into newly create transpose storage

  for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0; // Reset transpose NumNz counter
  for (i=0; i<NumMyRows_; i++)
  {
    if (OrigMatrixIsCrsMatrix_)
      err = OrigCrsMatrix->ExtractMyRowView(i, NumIndices, Values_, Indices_);
    else
      err = orig.ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_);
    if (err != 0) {
      std::cerr << "ExtractMyRowCopy failed."<<std::endl;
      throw err;
    }

    int ii = orig.RowMatrixRowMap().GID(i);
    for (j=0; j<NumIndices; j++)
    {
      int TransRow = Indices_[j];
      int loc = TransNumNz_[TransRow];
      TransIndices_[TransRow][loc] = ii;
      TransValues_[TransRow][loc] = Values_[j];
      ++TransNumNz_[TransRow]; // increment counter into current transpose row
    }
  }

  //  Build Transpose matrix with some rows being shared across processors.
  //  We will use a view here since the matrix will not be used for anything else

  const Epetra_Map & TransMap = orig.RowMatrixColMap();

  Epetra_CrsMatrix * TempTransA1 = new Epetra_CrsMatrix(Copy, TransMap, TransNumNz_);
  TransMap.MyGlobalElements(TransMyGlobalEquations_);
  
  for (i=0; i<NumMyCols_; i++) {
    err = TempTransA1->InsertGlobalValues(TransMyGlobalEquations_[i], 
                     TransNumNz_[i], TransValues_[i], TransIndices_[i]);
    if (err < 0) throw TempTransA1->ReportError("InsertGlobalValues failed.",err);
  }

  err = TempTransA1->FillComplete(orig.OperatorRangeMap(),*TransposeRowMap_); 
  if (err != 0) {
    throw TempTransA1->ReportError("FillComplete failed.",err);
  }

  if(!TempTransA1->Exporter()) {
    // Short Circuit: There is no need to make another matrix since TransposeRowMap_== TransMap
    newObj_ = TransposeMatrix_ = TempTransA1;
    return *newObj_;    
  }
  else
    TransposeExporter_ = new Epetra_Export(*TempTransA1->Exporter());

  // Now that transpose matrix with shared rows is entered, create a new matrix that will
  // get the transpose with uniquely owned rows (using the same row distribution as A).
  TransposeMatrix_ = new Epetra_CrsMatrix(*TempTransA1,*TransposeExporter_,TransposeRowMap_);

  newObj_ = TransposeMatrix_;
  return *newObj_;
}

//=========================================================================
bool EpetraExt::RowMatrix_Transpose::fwd()
{
  int i, j, NumIndices, err;

  Epetra_CrsMatrix * OrigCrsMatrix = dynamic_cast<Epetra_CrsMatrix*>(origObj_);

  // Now copy values and global indices into newly create transpose storage

  for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0; // Reset transpose NumNz counter
  for (i=0; i<NumMyRows_; i++)
  {
    if(OrigMatrixIsCrsMatrix_)
      err = OrigCrsMatrix->ExtractMyRowView(i, NumIndices, Values_, Indices_);
    else
      err = origObj_->ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_);
    if (err != 0) {
      std::cerr << "ExtractMyRowCopy/View failed."<<std::endl;
      throw err;
    }

    int ii = origObj_->RowMatrixRowMap().GID(i);
    for (j=0; j<NumIndices; j++)
    {
      int TransRow = Indices_[j];
      int loc = TransNumNz_[TransRow];
      TransIndices_[TransRow][loc] = ii;
      TransValues_[TransRow][loc] = Values_[j];
      ++TransNumNz_[TransRow]; // increment counter into current transpose row
    }
  }

  //  Build Transpose matrix with some rows being shared across processors.
  //  We will use a view here since the matrix will not be used for anything else
  const Epetra_Map & TransMap = origObj_->RowMatrixColMap();

  Epetra_CrsMatrix TempTransA1(View, TransMap, TransNumNz_);
  TransMap.MyGlobalElements(TransMyGlobalEquations_);
  
  for (i=0; i<NumMyCols_; i++)
      EPETRA_CHK_ERR(TempTransA1.InsertGlobalValues(TransMyGlobalEquations_[i], 
                     TransNumNz_[i], TransValues_[i], TransIndices_[i]));
 
  // Note: The following call to FillComplete is currently necessary because
  //     some global constants that are needed by the Export () are computed in this routine
  EPETRA_CHK_ERR(TempTransA1.FillComplete(false));

  // Now that transpose matrix with shared rows is entered, update values of target transpose matrix
  TransposeMatrix_->PutScalar(0.0);  // Zero out all values of the matrix

  EPETRA_CHK_ERR(TransposeMatrix_->Export(TempTransA1, *TransposeExporter_, Add));

  return(0);
}

//=========================================================================
bool EpetraExt::RowMatrix_Transpose::rvs()
{
  EPETRA_CHK_ERR(-1); //Not Implemented Yet
  return false;
}

} // namespace EpetraExt
