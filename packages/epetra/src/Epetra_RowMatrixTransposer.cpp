
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#include "Epetra_RowMatrixTransposer.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_CrsGraph.h"
//=============================================================================
Epetra_RowMatrixTransposer::Epetra_RowMatrixTransposer(Epetra_RowMatrix * OrigMatrix)
  : OrigMatrix_(OrigMatrix),
    TransposeMatrix_(0),
    MakeDataContiguous_(false)
{
}
//=============================================================================
Epetra_RowMatrixTransposer::Epetra_RowMatrixTransposer(const Epetra_RowMatrixTransposer& Source)
  :OrigMatrix_(Source.OrigMatrix_),
    TransposeMatrix_(0),
    MakeDataContiguous_(Source.MakeDataContiguous_) 
{
	TransposeMatrix_ = new Epetra_CrsMatrix(*Source.TransposeMatrix_);
	if (MakeDataContiguous_) TransposeMatrix_->MakeDataContiguous();
}
//=========================================================================
Epetra_RowMatrixTransposer::~Epetra_RowMatrixTransposer(){

	DeleteData();

}

//=========================================================================
void Epetra_RowMatrixTransposer::Allocate (Epetra_CrsMatrix * CrsMatrix){

	NumMyRows_ = OrigMatrix_->NumMyRows();
	NumMyCols_ = OrigMatrix_->NumMyCols();
	if (CrsMatrix!=0) {
		MaxNumEntries_ = CrsMatrix_->MaxNumEntries();
	}
	else {
		MaxNumEntries_ = 0;
		for (i=0; i<NumMyRows_; i++) EPETRA_MAX(MaxNumEntries_, NumMyRowEntries(i));
		Indices_ = new int[MaxNumEntries_];
		Values_ = new double[MaxNumEntries_];
	}
		TransNumNz = new int[NumMyCols];
}

//=========================================================================
void Epetra_RowMatrixTransposer::DeleteData (){

	if (TransposeMatrix_!=0) {delete TransposeMatrix_; TransposeMatrix_=0;}
}

//=========================================================================
int Epetra_RowMatrixTransposer::CreateTranspose (const bool MakeDataContiguous, Epetra_CrsMatrix *& TransposeMatrix
																								 Epetra_Map * TransposeRowMap){

	int i, j;

  if (TransposeCreated_) DeleteData(); // Get rid of existing data first

	// This routine will work for any RowMatrix object, but will attempt cast the matrix to a CrsMatrix if
	// possible (because we can then use a View of the matrix and graph, which is much cheaper).

  // First get the local indices to count how many nonzeros will be in the 
  // transpose graph on each processor


	Epetra_CrsMatrix * CrsMatrix = dynamic_cast<Epetra_CrsMatrix *>(OrigMatrix_);

	bool IsCrsMatrix = (CrsMatrix!=0); // If this pointer is non-zero, the cast to CrsMatrix worked

	// Allocate internal storage
	Allocate(CrsMatrix);

	if (IsCrsMatrix) {


		const Epetra_CrsGraph & AG = A->Graph(); // Get matrix graph

		for (i=0;i<NumMyCols; i++) TransNumNz[i] = 0;
		for (i=0; i<NumMyEquations; i++) {
			EPETRA_CHK_ERR(AG.ExtractMyRowView(i, NumIndices, Indices)); // Get view of ith row
			for (j=0; j<NumIndices; j++) ++TransNumNz[Indices[j]];
		}
	}
	else { // OrigMatrix is not a CrsMatrix



		for (i=0;i<NumMyCols; i++) TransNumNz[i] = 0;
		for (i=0; i<NumMyEquations; i++) {
			// Get view of ith row
			EPETRA_CHK_ERR(OrigMatrix_->ExtractMyRowCopy(i, MaxNumEntries, NumIndices, Values, Indices)); 
			for (j=0; j<NumIndices; j++) ++TransNumNz[Indices[j]];
		}
	}

	// Most of remaining code is common to both cases

  int ** TransIndices = new int*[NumMyCols];
  double ** TransValues = new double*[NumMyCols];

  for(i=0; i<NumMyCols; i++) {
    NumIndices = TransNumNz[i];
    if (NumIndices>0) {
      TransIndices[i] = new int[NumIndices];
      TransValues[i] = new double[NumIndices];
    }
  }

  // Now copy values and global indices into newly create transpose storage

  for (i=0;i<NumMyCols; i++) TransNumNz[i] = 0; // Reset transpose NumNz counter
  for (i=0; i<NumMyEquations; i++) {
		if (IsCrsMatrix) {
			EPETRA_CHK_ERR(CrsMatrix->ExtractMyRowView(i, NumIndices, Values, Indices));
		}
		else {
			EPETRA_CHK_ERR(OrigMatrix_->ExtractMyRowCopy(i, MaxNumEntries, NumIndices, Values, Indices));

    int ii = OrigMatrix_->RowMatrixRowMap().GID(i);
    for (j=0; j<NumIndices; j++) {
      int TransRow = Indices[j];
      int loc = TransNumNz[TransRow];
      TransIndices[TransRow][loc] = ii;
      TransValues[TransRow][loc] = Values[j];
      ++TransNumNz[TransRow]; // increment counter into current transpose row
    }
  }

  //  Build Transpose matrix with some rows being shared across processors.
  //  We will use a view here since the matrix will not be used for anything else

  const Epetra_Map & TransMap = OrigMatrix_->RowMatrixColMap();

  Epetra_CrsMatrix TempTransA1(View, TransMap, TransNumNz);
  int * TransMyGlobalEquations = new int[NumMyCols];
  TransMap.MyGlobalElements(TransMyGlobalEquations);
  
  /* Add  rows one-at-a-time */

  for (i=0; i<NumMyCols; i++)
    {
     EPETRA_CHK_ERR(TempTransA1.InsertGlobalValues(TransMyGlobalEquations[i], 
				      TransNumNz[i], TransValues[i], TransIndices[i]));
    }
 
  // Note: The following call to TransformToLocal is currently necessary because
  //       some global constants that are needed by the Export () are computed in this routine
  EPETRA_CHK_ERR(TempTransA1.TransformToLocal());

  // Now that transpose matrix with shared rows is entered, create a new matrix that will
  // get the transpose with uniquely owned rows (using the same row distribution as A).

   TransposeMatrix_ = new Epetra_CrsMatrix(Copy, TransposeRowMap,0);

  // Create an Export object that will move TempTransA around

  Epetra_Export Export(TransMap, TransposeRowMap);

  EPETRA_CHK_ERR(Transposematrix_->Export(TempTransA1, Export, Add));
  
  EPETRA_CHK_ERR(Transposematrix_->TransformToLocal());

	if (MakeDataContiguous) {
		EPETRA_CHK_ERR(Transposematrix_->MakeDataContiguous());
	}

	// Delete any temp storage

	if (!IsCrsMatrix) {
		delete [] Indices;
		delete [] Values;
	}
	
  for(i=0; i<NumMyCols; i++) {
    NumIndices = TransNumNz[i];
    if (NumIndices>0) {
      delete [] TransIndices[i];
      delete [] TransValues[i];
    }
  }
	delete [] TransNumNz;
	delete [] TransIndices;
	delete [] TransValues;
	delete [] TransMyGlobalEquations
	TransposeCreated_ = true;

  return(0);
}
