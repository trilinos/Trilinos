// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
// @HEADER

#include <EpetraExt_Transpose_RowMatrix.h>

#include <Epetra_Export.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>

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

  int i, j;

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

  int NumIndices;

  if (OrigMatrixIsCrsMatrix_)
  {
    const Epetra_CrsGraph & OrigGraph = OrigCrsMatrix->Graph(); // Get matrix graph

    for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0;
    for (i=0; i<NumMyRows_; i++)
    {
      assert(OrigGraph.ExtractMyRowView(i, NumIndices, Indices_)==0); // Get view of ith row
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
      assert(orig.ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_)==0); 
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
      assert(OrigCrsMatrix->ExtractMyRowView(i, NumIndices, Values_, Indices_)==0);
    else
      assert(orig.ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_)==0);

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

  Epetra_CrsMatrix TempTransA1(View, TransMap, TransNumNz_);
  TransMyGlobalEquations_ = new int[NumMyCols_];
  TransMap.MyGlobalElements(TransMyGlobalEquations_);
  
  for (i=0; i<NumMyCols_; i++)
      assert(TempTransA1.InsertGlobalValues(TransMyGlobalEquations_[i], 
                     TransNumNz_[i], TransValues_[i], TransIndices_[i])>=0);
 
  // Note: The following call to TransformToLocal is currently necessary because
  //      some global constants that are needed by the Export () are computed in this routine
  assert(TempTransA1.TransformToLocal()==0);

  // Now that transpose matrix with shared rows is entered, create a new matrix that will
  // get the transpose with uniquely owned rows (using the same row distribution as A).
  if( IgnoreNonLocalCols_ )
    TransposeMatrix_ = new Epetra_CrsMatrix(Copy, *TransposeRowMap_, *TransposeRowMap_, 0);
  else
    TransposeMatrix_ = new Epetra_CrsMatrix(Copy, *TransposeRowMap_,0);

  // Create an Export object that will move TempTransA around
  TransposeExporter_ = new Epetra_Export(TransMap, *TransposeRowMap_);

  assert(TransposeMatrix_->Export(TempTransA1, *TransposeExporter_, Add)==0);
  
  assert(TransposeMatrix_->TransformToLocal()==0);

  if (MakeDataContiguous_)
    assert(TransposeMatrix_->MakeDataContiguous()==0);

  newObj_ = TransposeMatrix_;

  return *newObj_;
}

//=========================================================================
bool EpetraExt::RowMatrix_Transpose::fwd()
{
  int i, j, NumIndices;

  Epetra_CrsMatrix * OrigCrsMatrix = dynamic_cast<Epetra_CrsMatrix*>(origObj_);

  // Now copy values and global indices into newly create transpose storage

  for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0; // Reset transpose NumNz counter
  for (i=0; i<NumMyRows_; i++)
  {
    if(OrigMatrixIsCrsMatrix_)
      assert(OrigCrsMatrix->ExtractMyRowView(i, NumIndices, Values_, Indices_)==0);
    else
      assert(origObj_->ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_)==0);

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
  TransMyGlobalEquations_ = new int[NumMyCols_];
  TransMap.MyGlobalElements(TransMyGlobalEquations_);
  
  for (i=0; i<NumMyCols_; i++)
      EPETRA_CHK_ERR(TempTransA1.InsertGlobalValues(TransMyGlobalEquations_[i], 
                     TransNumNz_[i], TransValues_[i], TransIndices_[i]));
 
  // Note: The following call to TransformToLocal is currently necessary because
  //     some global constants that are needed by the Export () are computed in this routine
  EPETRA_CHK_ERR(TempTransA1.TransformToLocal());

  // Now that transpose matrix with shared rows is entered, update values of target transpose matrix
  TransposeMatrix_->PutScalar(0.0);  // Zero out all values of the matrix

  EPETRA_CHK_ERR(TransposeMatrix_->Export(TempTransA1, *TransposeExporter_, Add));

  return(0);
}

//=========================================================================
bool EpetraExt::RowMatrix_Transpose::rvs()
{
  EPETRA_CHK_ERR(-1); //Not Implemented Yet
}

} // namespace EpetraExt
