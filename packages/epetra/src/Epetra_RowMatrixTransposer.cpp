
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
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
// ************************************************************************
//@HEADER

#include "Epetra_RowMatrixTransposer.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"
#include "Epetra_Export.h"
//=============================================================================
Epetra_RowMatrixTransposer::Epetra_RowMatrixTransposer(Epetra_RowMatrix * OrigMatrix)
  : OrigMatrix_(OrigMatrix),
    TransposeMatrix_(0),
    TransposeExporter_(0),
    TransposeRowMap_(0),
    TransposeCreated_(false),
    MakeDataContiguous_(false),
    NumMyRows_(0),
    NumMyCols_(0),
    MaxNumEntries_(0),
    Indices_(NULL),
    Values_(NULL),
    TransNumNz_(NULL),
    TransIndices_(NULL),
    TransValues_(NULL),
    TransMyGlobalEquations_(NULL),
    OrigMatrixIsCrsMatrix_(false)
{
}
//=============================================================================
Epetra_RowMatrixTransposer::Epetra_RowMatrixTransposer(const Epetra_RowMatrixTransposer& Source)
  :OrigMatrix_(Source.OrigMatrix_),
   TransposeMatrix_(0),
   TransposeExporter_(0),
   TransposeRowMap_(0),
   TransposeCreated_(Source.TransposeCreated_),
   MakeDataContiguous_(Source.MakeDataContiguous_),
   NumMyRows_(0),
   NumMyCols_(0),
   MaxNumEntries_(0),
   Indices_(NULL),
   Values_(NULL),
   TransNumNz_(NULL),
   TransIndices_(NULL),
   TransValues_(NULL),
   TransMyGlobalEquations_(NULL),
   OrigMatrixIsCrsMatrix_(false)
{
  TransposeMatrix_ = new Epetra_CrsMatrix(*Source.TransposeMatrix_);
  if (MakeDataContiguous_) TransposeMatrix_->MakeDataContiguous();
  TransposeExporter_ = new Epetra_Export(*Source.TransposeExporter_);
}
//=========================================================================
Epetra_RowMatrixTransposer::~Epetra_RowMatrixTransposer(){

  DeleteData();

}

//=========================================================================
void Epetra_RowMatrixTransposer::DeleteData (){

  int i;

  if (TransposeMatrix_!=0) {delete TransposeMatrix_; TransposeMatrix_=0;}
  if (TransposeExporter_!=0) {delete TransposeExporter_; TransposeExporter_=0;}

  // Delete any intermediate storage

  if (!OrigMatrixIsCrsMatrix_) {
    delete [] Indices_;
    delete [] Values_;
  }
	
	
  for(i=0; i<NumMyCols_; i++) {
    int NumIndices = TransNumNz_[i];
    if (NumIndices>0) {
      delete [] TransIndices_[i];
      delete [] TransValues_[i];
    }
  }
  delete [] TransNumNz_;
  delete [] TransIndices_;
  delete [] TransValues_;
  delete [] TransMyGlobalEquations_;
}

//=========================================================================
int Epetra_RowMatrixTransposer::CreateTranspose (const bool MakeDataContiguous, 
						 Epetra_CrsMatrix *& TransposeMatrix, 
						 Epetra_Map * TransposeRowMap) {

  int i, j;

  if (TransposeCreated_) DeleteData(); // Get rid of existing data first

  if (TransposeRowMap==0)
    TransposeRowMap_ = (Epetra_Map *) &(OrigMatrix_->OperatorDomainMap()); // Should be replaced with refcount =
  else
    TransposeRowMap_ = TransposeRowMap; 

  // This routine will work for any RowMatrix object, but will attempt cast the matrix to a CrsMatrix if
  // possible (because we can then use a View of the matrix and graph, which is much cheaper).

  // First get the local indices to count how many nonzeros will be in the 
  // transpose graph on each processor


  Epetra_CrsMatrix * OrigCrsMatrix = dynamic_cast<Epetra_CrsMatrix *>(OrigMatrix_);

  OrigMatrixIsCrsMatrix_ = (OrigCrsMatrix!=0); // If this pointer is non-zero, the cast to CrsMatrix worked

  NumMyRows_ = OrigMatrix_->NumMyRows();
  NumMyCols_ = OrigMatrix_->NumMyCols();
  NumMyRows_ = OrigMatrix_->NumMyRows();
  TransNumNz_ = new int[NumMyCols_];
  TransIndices_ = new int*[NumMyCols_];
  TransValues_ = new double*[NumMyCols_];


  int NumIndices;

  if (OrigMatrixIsCrsMatrix_) {


    const Epetra_CrsGraph & OrigGraph = OrigCrsMatrix->Graph(); // Get matrix graph

    for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0;
    for (i=0; i<NumMyRows_; i++) {
      EPETRA_CHK_ERR(OrigGraph.ExtractMyRowView(i, NumIndices, Indices_)); // Get view of ith row
      for (j=0; j<NumIndices; j++) ++TransNumNz_[Indices_[j]];
    }
  }
  else { // OrigMatrix is not a CrsMatrix

    MaxNumEntries_ = 0;
    int NumEntries;
    for (i=0; i<NumMyRows_; i++) {
      OrigMatrix_->NumMyRowEntries(i, NumEntries);
      MaxNumEntries_ = EPETRA_MAX(MaxNumEntries_, NumEntries);
    }
    Indices_ = new int[MaxNumEntries_];
    Values_ = new double[MaxNumEntries_];

    for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0;
    for (i=0; i<NumMyRows_; i++) {
      // Get ith row
      EPETRA_CHK_ERR(OrigMatrix_->ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_)); 
      for (j=0; j<NumIndices; j++) ++TransNumNz_[Indices_[j]];
    }
  }


  // Most of remaining code is common to both cases

  for(i=0; i<NumMyCols_; i++) {
    NumIndices = TransNumNz_[i];
    if (NumIndices>0) {
      TransIndices_[i] = new int[NumIndices];
      TransValues_[i] = new double[NumIndices];
    }
  }

  // Now copy values and global indices into newly create transpose storage

  for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0; // Reset transpose NumNz counter
  for (i=0; i<NumMyRows_; i++) {
    if (OrigMatrixIsCrsMatrix_) {
      EPETRA_CHK_ERR(OrigCrsMatrix->ExtractMyRowView(i, NumIndices, Values_, Indices_));
    }
    else {
      EPETRA_CHK_ERR(OrigMatrix_->ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_));
    }

    int ii = OrigMatrix_->RowMatrixRowMap().GID(i);
    for (j=0; j<NumIndices; j++) {
      int TransRow = Indices_[j];
      int loc = TransNumNz_[TransRow];
      TransIndices_[TransRow][loc] = ii;
      TransValues_[TransRow][loc] = Values_[j];
      ++TransNumNz_[TransRow]; // increment counter into current transpose row
    }
  }

  //  Build Transpose matrix with some rows being shared across processors.
  //  We will use a view here since the matrix will not be used for anything else

  const Epetra_Map & TransMap = OrigMatrix_->RowMatrixColMap();

  Epetra_CrsMatrix TempTransA1(View, TransMap, TransNumNz_);
  TransMyGlobalEquations_ = new int[NumMyCols_];
  TransMap.MyGlobalElements(TransMyGlobalEquations_);
  
  /* Add  rows one-at-a-time */

  for (i=0; i<NumMyCols_; i++)
    {
      EPETRA_CHK_ERR(TempTransA1.InsertGlobalValues(TransMyGlobalEquations_[i], 
						    TransNumNz_[i], TransValues_[i], TransIndices_[i]));
    }
 
  // Note: The following call to TransformToLocal is currently necessary because
  //       some global constants that are needed by the Export () are computed in this routine
  EPETRA_CHK_ERR(TempTransA1.TransformToLocal());

  // Now that transpose matrix with shared rows is entered, create a new matrix that will
  // get the transpose with uniquely owned rows (using the same row distribution as A).

  TransposeMatrix_ = new Epetra_CrsMatrix(Copy, *TransposeRowMap_,0);

  // Create an Export object that will move TempTransA around

  TransposeExporter_ = new Epetra_Export(TransMap, *TransposeRowMap_);

  EPETRA_CHK_ERR(TransposeMatrix_->Export(TempTransA1, *TransposeExporter_, Add));
  
  EPETRA_CHK_ERR(TransposeMatrix_->TransformToLocal());

  if (MakeDataContiguous) {
    EPETRA_CHK_ERR(TransposeMatrix_->MakeDataContiguous());
  }

  TransposeMatrix = TransposeMatrix_;
  TransposeCreated_ = true;

  return(0);
}
//=========================================================================
int Epetra_RowMatrixTransposer::UpdateTransposeValues(Epetra_RowMatrix * MatrixWithNewValues){

  int i, j, NumIndices;

  if (!TransposeCreated_) EPETRA_CHK_ERR(-1); // Transpose must be already created

  // Sanity check of incoming matrix.  Perform some tests to see if it is compatible with original input matrix
  if (OrigMatrix_!=MatrixWithNewValues) { // Check if pointer of new matrix is same as previous input matrix
    OrigMatrix_ = MatrixWithNewValues; // Reset this pointer if not, then check for other attributes
    if (NumMyRows_ != OrigMatrix_->NumMyRows() ||
	NumMyCols_ != OrigMatrix_->NumMyCols() ||
	NumMyRows_ != OrigMatrix_->NumMyRows()) {
      EPETRA_CHK_ERR(-2); // New matrix not compatible with previous
    }
  }

  Epetra_CrsMatrix * OrigCrsMatrix = dynamic_cast<Epetra_CrsMatrix *>(MatrixWithNewValues);

	
  OrigMatrixIsCrsMatrix_ = (OrigCrsMatrix!=0); // If this pointer is non-zero, the cast to CrsMatrix worked


  // Now copy values and global indices into newly create transpose storage

  for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0; // Reset transpose NumNz counter
  for (i=0; i<NumMyRows_; i++) {
    if (OrigMatrixIsCrsMatrix_) {
      EPETRA_CHK_ERR(OrigCrsMatrix->ExtractMyRowView(i, NumIndices, Values_, Indices_));
    }
    else {
      EPETRA_CHK_ERR(OrigMatrix_->ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_));
    }

    int ii = OrigMatrix_->RowMatrixRowMap().GID(i);
    for (j=0; j<NumIndices; j++) {
      int TransRow = Indices_[j];
      int loc = TransNumNz_[TransRow];
      TransIndices_[TransRow][loc] = ii;
      TransValues_[TransRow][loc] = Values_[j];
      ++TransNumNz_[TransRow]; // increment counter into current transpose row
    }
  }

  //  Build Transpose matrix with some rows being shared across processors.
  //  We will use a view here since the matrix will not be used for anything else

  const Epetra_Map & TransMap = OrigMatrix_->RowMatrixColMap();

  Epetra_CrsMatrix TempTransA1(View, TransMap, TransNumNz_);
  TransMyGlobalEquations_ = new int[NumMyCols_];
  TransMap.MyGlobalElements(TransMyGlobalEquations_);
  
  /* Add  rows one-at-a-time */

  for (i=0; i<NumMyCols_; i++)
    {
      EPETRA_CHK_ERR(TempTransA1.InsertGlobalValues(TransMyGlobalEquations_[i], 
						    TransNumNz_[i], TransValues_[i], TransIndices_[i]));
    }
 
  // Note: The following call to TransformToLocal is currently necessary because
  //       some global constants that are needed by the Export () are computed in this routine
  EPETRA_CHK_ERR(TempTransA1.TransformToLocal());

  // Now that transpose matrix with shared rows is entered, update values of target transpose matrix

  TransposeMatrix_->PutScalar(0.0);  // Zero out all values of the matrix

  EPETRA_CHK_ERR(TransposeMatrix_->Export(TempTransA1, *TransposeExporter_, Add));

  return(0);
}
