//@HEADER
// ***********************************************************************
// 
//                     New_Package Example Package
//                 Copyright (2004) Sandia Corporation
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

#include "EpetraExt_MatlabEngine.h"
#include "EpetraExt_PutMultiVector.h"
#include "EpetraExt_PutRowMatrix.h"

#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_IntVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_DataAccess.h"
#include "Epetra_Import.h"

#include "Epetra_CrsMatrix.h"

using namespace EpetraExt;
namespace EpetraExt {

//=============================================================================
EpetraExt_MatlabEngine::EpetraExt_MatlabEngine (const Epetra_Comm& Comm):Comm_(Comm) {

    // MATLAB engOpen, to construct the MATLAB engine

    if (Comm_.MyPID() == 0) {
		// Gentlemen, start your engines ...
		Engine_ = engOpen(NULL);
    }
} 
 
//=============================================================================
EpetraExt_MatlabEngine::~EpetraExt_MatlabEngine (void) {

    // MATLAB engClose, to destruct the MATLAB engine

    if (Comm_.MyPID() == 0) {
		engClose(Engine_);
    }
}

//=============================================================================
int EpetraExt_MatlabEngine::EvalString (char* command, char* outputBuffer, int outputBufferSize) {
    // send a string command to the MATLAB engine
    if (Comm_.MyPID() == 0) {
		if (outputBuffer != NULL) {
		  if (engOutputBuffer(Engine_, outputBuffer, outputBufferSize))
			return(-4);
		}
		if (engEvalString(Engine_, command))
		  return(-3);
    }
	
	return(0);
}

//=============================================================================
int EpetraExt_MatlabEngine::PutMultiVector(const Epetra_MultiVector& A, const char * variableName) {
    mxArray* matlabA = 0;
    double* pr = 0;
    if (Comm_.MyPID() == 0) {
      matlabA = mxCreateDoubleMatrix(A.GlobalLength(), A.NumVectors(), mxREAL);
      pr = mxGetPr(matlabA);
    }

	if (Matlab::CopyMultiVector(&pr, A)) {
	  mxDestroyArray(matlabA);
	  return(-2);
	}

	if (Comm_.MyPID() == 0)
	  if (PutIntoMatlab(variableName, matlabA)) {
		mxDestroyArray(matlabA);
		return(-1);
	  }
	
	mxDestroyArray(matlabA);
	return(0);
}

//=============================================================================
int EpetraExt_MatlabEngine::PutRowMatrix(const Epetra_RowMatrix& A, const char* variableName, bool transA) {
    mxArray* matlabA = 0;
    if (Comm_.MyPID() == 0)
	  // since matlab uses column major for matrices, switch row and column numbers
	  matlabA = mxCreateSparse(A.RowMatrixColMap().MaxAllGID() - A.RowMatrixColMap().MinAllGID()+1, A.RowMatrixRowMap().MaxAllGID() - A.RowMatrixRowMap().MinAllGID() + 1, A.NumGlobalNonzeros(), mxREAL);
	
	//cout << "numrows: " << A.RowMatrixColMap().NumGlobalElements() << " numCols: " << A.NumGlobalRows() << "numNonZeros: " << A.NumGlobalNonzeros() << "\n";

	cout << "calling CopyRowMatrix\n";
	if (Matlab::CopyRowMatrix(matlabA, A)) {
	  mxDestroyArray(matlabA);
	  return(-2);
	}

	cout << "done doing CopyRowMatrix\n";
	if (Comm_.MyPID() == 0) {

	  /*cout << "printing matlabA pointers\n";
		double* matlabAvaluesPtr = mxGetPr(matlabA);
		int* matlabAcolumnIndicesPtr = mxGetJc(matlabA);
		int* matlabArowIndicesPtr = mxGetIr(matlabA);
		for(int i=0; i < A.NumGlobalNonzeros(); i++) {
		  cout << "*matlabAvaluesPtr: " << *matlabAvaluesPtr++ << " *matlabAcolumnIndicesPtr: " << *matlabAcolumnIndicesPtr++ << " *matlabArowIndicesPtr" << *matlabArowIndicesPtr++ << "\n";
		}
		cout << "done printing matlabA pointers\n";
	  */
	  if (PutIntoMatlab(variableName, matlabA)) {
		mxDestroyArray(matlabA);
		return(-1);
	  }

	  if (!transA) {
		char* buff = new char[128];;
		sprintf(buff, "%s = %s'", variableName, variableName);
		if (EvalString(buff)) {
		  mxDestroyArray(matlabA);
		  return(-3);
		}
	  }
	}

	cout << "done with everything in PutRowMatrix, going to destroy matlabA\n" << "matlabA=" << matlabA << "\n";
	mxDestroyArray(matlabA);
	cout << "done destroying matlabA\n";
	return(0);
}

//=============================================================================
int EpetraExt_MatlabEngine::PutCrsGraph(const Epetra_CrsGraph& A, const char* variableName, bool transA) {
	return(-1);
}

int EpetraExt_MatlabEngine::PutSerialDenseMatrix(const Epetra_SerialDenseMatrix& A, const char* variableName, int proc) {

  if (proc == 0) {
	if (Comm_.MyPID() == 0) {
	  mxArray* matlabA = 0;

	  int numRows = A.M();
	  int numCols = A.N();
	  
	  matlabA = mxCreateDoubleMatrix(numRows, numCols, mxREAL);

	  int row;
	  int col;
	  double* targetPtr = 0;
	  double* sourcePtr = 0;
	  double* source = (double*)A.A();
	  double* target = (double*)mxGetPr(matlabA);
	  int source_LDA = A.LDA();
	  int target_LDA = A.LDA();
	  for (col = 0; col < numCols; col++) {
		targetPtr = target + (col * target_LDA);
		sourcePtr = source + (col * source_LDA);
		for (row = 0; row < numRows; row++) {
			*targetPtr++ = *sourcePtr++;
		}
	  }

	  if (PutIntoMatlab(variableName, matlabA)) {
		mxDestroyArray(matlabA);
		return(-1);
	  }

	  mxDestroyArray(matlabA);
	}
  }
  else {
	Epetra_MultiVector* tempMultiVector = 0;
	if (proc == Comm_.MyPID()) {
	  int* numVectors = new int[1];
	  int* temp = new int[1];
	  temp[0] = A.N();
	  Comm_.MaxAll (temp, numVectors, 1);
	  const Epetra_BlockMap tempBlockMap (-1, A.LDA(), 1, 0, Comm_);
	  tempMultiVector = new Epetra_MultiVector (View, tempBlockMap, A.A(), A.LDA(), A.N());
	}
	else {
	  int* numVectors = new int[1];
	  int* temp = new int[1];
	  temp[0] = 0;
	  Comm_.MaxAll (temp, numVectors, 1);
	  const Epetra_BlockMap tempBlockMap (-1, 0, 1, 0, Comm_);
	  tempMultiVector = new Epetra_MultiVector (tempBlockMap, numVectors[0], false);
	}

	return(PutMultiVector(*tempMultiVector, variableName));	
  }
}

//=============================================================================
int EpetraExt_MatlabEngine::PutIntSerialDenseMatrix(const Epetra_IntSerialDenseMatrix& A, const char* variableName, int proc) {
  mxArray* matlabA = 0;
  
    if (proc == 0) {
	if (Comm_.MyPID() == 0) {
	  int numRows = A.M();
	  int numCols = A.N();
	  
	  matlabA = mxCreateDoubleMatrix(numRows, numCols, mxREAL);

	  int row;
	  int col;
	  double* targetPtr = 0;
	  int* sourcePtr = 0;
	  int* source = (int*)A.A();
	  double* target = (double*)mxGetPr(matlabA);
	  int source_LDA = A.LDA();
	  int target_LDA = A.LDA();
	  for (col = 0; col < numCols; col++) {
		targetPtr = target + (col * target_LDA);
		sourcePtr = source + (col * source_LDA);
		for (row = 0; row < numRows; row++) {
			*targetPtr++ = *sourcePtr++;
		}
	  }

	  if (PutIntoMatlab(variableName, matlabA)) {
		mxDestroyArray(matlabA);
		return(-1);
	  }
	}
  }
  else {
  int* numVectors = new int[2];
  if (proc == Comm_.MyPID()) {
	int* temp = new int[2];
	temp[0] = A.N();
	temp[1] = A.M();
	Comm_.MaxAll (temp, numVectors, 2);
  }
  else {
	int* temp = new int[2];
	temp[0] = 0;
	temp[1] = 0;
	Comm_.MaxAll (temp, numVectors, 2);
  }

  int numRows = numVectors[1];
  int numCols = numVectors[0];
  double* targetPtr = 0;
  int* sourcePtr = 0;
  int row;
  double* target = 0;
  const Epetra_BlockMap* tgBlockMap = 0;
  if (Comm_.MyPID() == 0) {	
	matlabA = mxCreateDoubleMatrix(numRows, numCols, mxREAL);
	target = (double*)mxGetPr(matlabA);
	tgBlockMap = new Epetra_BlockMap(-1, numRows, 1, 0, Comm_);
  }
  else {
	tgBlockMap = new Epetra_BlockMap(-1, 0, 1, 0, Comm_);
  }

  const Epetra_BlockMap* srcBlockMap = 0;
  Epetra_IntVector* srcIntVector = 0;
  Epetra_IntVector tgIntVector (*tgBlockMap, false);
  if (proc == Comm_.MyPID()) {
	srcBlockMap = new Epetra_BlockMap(-1, A.LDA(), 1, 0, Comm_);
  }
  else {
	srcBlockMap = new Epetra_BlockMap(-1, 0, 1, 0, Comm_);
  }

  Epetra_Import importer (*tgBlockMap, *srcBlockMap);
  
  for(int i=0; i < numVectors[0]; i++) {
	if (proc == Comm_.MyPID()) {
	  srcIntVector = new Epetra_IntVector(View, *srcBlockMap, (int*)A[i]);
	}
	else {
	  srcIntVector = new Epetra_IntVector(*srcBlockMap, false);
	}
	// need to add some error checking for this!
	tgIntVector.Import(*srcIntVector, importer, Insert);
	if (Comm_.MyPID() == 0) {
	  targetPtr = target + (i * numRows);
	  sourcePtr = (int*)tgIntVector.Values();
	  for (row = 0; row < numRows; row++) {
		*targetPtr++ = *sourcePtr++;
	  }
	}
  }

  if (PutIntoMatlab(variableName, matlabA)) {
	mxDestroyArray(matlabA);
	return(-1);
  }
}
  mxDestroyArray(matlabA);
  return(0);
}

//=============================================================================
int EpetraExt_MatlabEngine::PutBlockMap(const Epetra_BlockMap& blockMap, const char* variableName) {
	return(-1);
}

int EpetraExt_MatlabEngine::PutIntoMatlab(const char* variableName, mxArray* matlabA) {
  if (Comm_.MyPID() != 0)
    return(0);
#ifdef USE_ENGPUTARRAY
    // for matlab versions < 6.5 (release 13)
    mxSetName(matlabA, variableName);
    return engPutArray(Engine_, matlabA);
#else
    // for matlab versions >= 6.5 (release 13)
    return engPutVariable(Engine_, variableName, matlabA);
#endif
}

} // namespace EpetraExt
