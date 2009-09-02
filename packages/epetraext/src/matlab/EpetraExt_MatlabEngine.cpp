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
#include "EpetraExt_PutBlockMap.h"

#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_IntVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_DataAccess.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CombineMode.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

using namespace EpetraExt;
namespace EpetraExt {

//=============================================================================
EpetraExt_MatlabEngine::EpetraExt_MatlabEngine (const Epetra_Comm& Comm):Comm_(Comm) {
  if (Comm_.MyPID() == 0) {
		// construct the MATLAB engine
		Engine_ = engOpen(NULL);
  }
} 
 
//=============================================================================
EpetraExt_MatlabEngine::~EpetraExt_MatlabEngine (void) {
  if (Comm_.MyPID() == 0) {
    // to destruct the MATLAB engine
		engClose(Engine_);
  }
}

//=============================================================================
int EpetraExt_MatlabEngine::EvalString (char* command, char* outputBuffer, int outputBufferSize) {
  // send a string command to the MATLAB engine
  if (Comm_.MyPID() == 0) {
    if (outputBuffer != NULL) {
      if (engOutputBuffer(Engine_, outputBuffer, outputBufferSize)) {
        return(-4);
      }
    }
    if (engEvalString(Engine_, command)) {
      return(-3);
    }
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
  if (Comm_.MyPID() == 0) {
	  // since matlab uses column major for matrices, switch row and column numbers
	  matlabA = mxCreateSparse(A.OperatorDomainMap().MaxAllGID() - A.OperatorDomainMap().MinAllGID()+1, A.OperatorRangeMap().MaxAllGID() - A.OperatorRangeMap().MinAllGID() + 1, A.NumGlobalNonzeros(), mxREAL);
	}
	//cout << "numrows: " << A.RowMatrixColMap().NumGlobalElements() << " numCols: " << A.NumGlobalRows() << "numNonZeros: " << A.NumGlobalNonzeros() << "\n";

	//cout << "calling CopyRowMatrix\n";
	if (Matlab::CopyRowMatrix(matlabA, A)) {
	  mxDestroyArray(matlabA);
	  return(-2);
	}

	//cout << "done doing CopyRowMatrix\n";
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

	//cout << "done with everything in PutRowMatrix, going to destroy matlabA\n" << "matlabA=" << matlabA << "\n";
	mxDestroyArray(matlabA);
	//cout << "done destroying matlabA\n";
	return(0);
}

//=============================================================================
int EpetraExt_MatlabEngine::PutCrsGraph(const Epetra_CrsGraph& A, const char* variableName, bool transA) {
	return(-9999);
}

//=============================================================================
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

	  return(0);
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
int EpetraExt_MatlabEngine::PutBlockMap(const Epetra_BlockMap& blockMap, const char* variableName, bool transA) {
  mxArray* matlabA = 0;
  if (Comm_.MyPID() == 0) {		 
    int M = blockMap.NumGlobalElements();
    int N = 1;
  
    if (blockMap.MaxElementSize()>1) N = 2; // Non-trivial block map, store element sizes in second column
  
    matlabA = mxCreateSparse(N, M, M*N, mxREAL);
    int* matlabAcolumnIndicesPtr = mxGetJc(matlabA);
    matlabAcolumnIndicesPtr += M;
    *matlabAcolumnIndicesPtr = M*N; // set max cap.
  }
  
  if (Matlab::CopyBlockMap(matlabA, blockMap)) {
    mxDestroyArray(matlabA);
    return(-2);
  }

  if (Comm_.MyPID() == 0) {
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

  mxDestroyArray(matlabA);
  return(0);
}

//=============================================================================
int EpetraExt_MatlabEngine::PutIntoMatlab(const char* variableName, mxArray* matlabA) {
  if (Comm_.MyPID() != 0) {
    return(0);
  }
#ifdef USE_ENGPUTARRAY
  // for matlab versions < 6.5 (release 13)
  mxSetName(matlabA, variableName);
  return(engPutArray(Engine_, matlabA));
#else
  // for matlab versions >= 6.5 (release 13)
  return(engPutVariable(Engine_, variableName, matlabA));
#endif
}

//=============================================================================
int EpetraExt_MatlabEngine::GetMultiVector(const char* variableName, Epetra_MultiVector& A) {
  mxArray* matlabA = 0;
  int ierr = 0;
  if (ierr = GetmxArray(variableName, &matlabA)) {
    mxDestroyArray(matlabA);
    return(ierr);
  }

  bool isSparse = false;
  int numRows = 0;
  int numCols = 0;
  int numNonZeros = 0;

  if (ierr = GetmxArrayDimensions(matlabA, isSparse, numRows, numCols, numNonZeros)) {
    mxDestroyArray(matlabA);
    return(ierr);
  }

  if ((Comm_.MyPID() == 0) && isSparse) {
    mxDestroyArray(matlabA);
    return(-7);
  }
  
  // tempMap is nontrivial only on PE0
  Epetra_Map tempMap (-1, numRows, 0, Comm_);
  
  int numVectors = 0;
	Comm_.MaxAll (&numCols, &numVectors, 1);
  double* tempValues = 0;
  if (Comm_.MyPID() == 0) {
    tempValues = mxGetPr(matlabA);
  }
  Epetra_MultiVector tempMV (View, tempMap, tempValues, numRows, numVectors);
  Epetra_Export exporter (tempMap, A.Map());
  A.Export(tempMV, exporter, Insert);
  
  mxDestroyArray(matlabA);
  return(0);
}

//=============================================================================
int EpetraExt_MatlabEngine::GetSerialDenseMatrix(const char* variableName, Epetra_SerialDenseMatrix& A, int proc) {
  int ierr = 0;
  int numMyGIDs = 0;
  if (Comm_.MyPID() == proc) {
    numMyGIDs = A.M();
  }
  Epetra_Map tempMap (-1, numMyGIDs, 0, Comm_);
  Epetra_MultiVector tempMV (View, tempMap, A.A(), numMyGIDs, A.N());
  
  if (ierr = GetMultiVector(variableName, tempMV)) {
    return(ierr);
  }
  
  if (Comm_.MyPID() == proc) {
    double* aValues = A.A();
    double* mvValues = tempMV.Values();
    for(int i=0; i < A.M() * A.N(); i++) {
      *aValues++ = *mvValues++;
    }
  }
  
  return(0);
}

//=============================================================================
int EpetraExt_MatlabEngine::GetIntSerialDenseMatrix(const char* variableName, Epetra_IntSerialDenseMatrix& A, int proc) {
  int ierr = 0;
  int numMyGIDs = 0;
  double* values = 0;
  if (Comm_.MyPID() == proc) {
    numMyGIDs = A.M();
    int* aValues = A.A();
    values = new double[A.M() * A.N()];
    for (int i=0; i < A.M() * A.N(); i++) {
      *values++ = *aValues++;
    }
  }
  Epetra_Map tempMap (-1, numMyGIDs, 0, Comm_);
  Epetra_MultiVector tempMV (View, tempMap, values, numMyGIDs, A.N());
  
  if (ierr = GetMultiVector(variableName, tempMV)) {
    return(ierr);
  }
  
  if (Comm_.MyPID() == proc) {
    int* aValues = A.A();
    double* mvValues = tempMV.Values();
    for(int i=0; i < A.M() * A.N(); i++) {
      *aValues++ = *mvValues++;
    }
  }
  
  return(0);
}

//=============================================================================
int EpetraExt_MatlabEngine::GetCrsMatrix(const char* inVariableName, Epetra_CrsMatrix& A, bool getTrans) {
  const char* variableName;
  
  if (!getTrans) {
    int inVariableNameLength = strlen(inVariableName);
    char* buff = new char[inVariableNameLength*2 + 11];
    sprintf(buff, "TRANS_%s = %s'", inVariableName, inVariableName);
    if (EvalString(buff)) {
      return(-3);
    }
    char* tempStr = new char[inVariableNameLength + 7];
    sprintf(tempStr, "TRANS_%s", inVariableName);
    variableName = tempStr;
  }
  else {
    variableName = inVariableName;
  }
  
  mxArray* matlabA = 0;
  int ierr = 0;
  if (ierr = GetmxArray(variableName, &matlabA)) {
    mxDestroyArray(matlabA);
    return(ierr);
  }
  
  if (!getTrans) {
    char* buff = new char[strlen(variableName) + 7];
    sprintf(buff, "clear %s", variableName);
    if (EvalString(buff)) {
      return(-3);
    }
  }
  
  bool isSparse = false;
  int numRows = 0;
  int numCols = 0;
  int numNonZeros = 0;

  // note that numCols and numRows are transposed on purpose here
  // we will treat the column major mxArray as a row major mxArray
  if (ierr = GetmxArrayDimensions(matlabA, isSparse, numCols, numRows, numNonZeros)) {
    mxDestroyArray(matlabA);
    return(ierr);
  }

  if ((Comm_.MyPID() == 0) && !isSparse) {
    mxDestroyArray(matlabA);
    return(-8);
  }
  
  // tempMap is nontrivial only on PE0
  Epetra_Map tempMap (-1, numRows, 0, Comm_);
  
  int numVectors = 0;
  Comm_.MaxAll (&numCols, &numVectors, 1);
  Epetra_CrsMatrix tempCRSM (View, tempMap, numVectors);
  if (Comm_.MyPID() == 0) {
    int* rowIndex = mxGetJc(matlabA);
    int* colIndex = mxGetIr(matlabA);
    double* values = mxGetPr(matlabA);
    int numCols = 0;
    for(int row = 0; row < numRows; row++) {
      numCols = *(rowIndex+1) - *rowIndex;
      if (numCols > 0) {
        int* indices = new int[numCols];
        int* indicesPtr = indices;
        for(int col=0; col < numCols; col++) {
          *indicesPtr++ = *colIndex++;
        }
        tempCRSM.InsertGlobalValues(row, numCols, values, indices);
      }
      values += numCols;
      rowIndex++;
    }
  }
  
  tempCRSM.FillComplete();
  Epetra_Export exporter (tempMap, A.RowMatrixRowMap());
  A.Export(tempCRSM, exporter, Insert);
  
  mxDestroyArray(matlabA);
  return(0);
}

//=============================================================================
int EpetraExt_MatlabEngine::GetmxArrayDimensions(mxArray* matlabA, bool& isSparse, int& numRows, int& numCols, int& numNonZeros) {
  if (Comm_.MyPID() == 0) {
    if (mxGetNumberOfDimensions(matlabA) > 2) {
      return(-6);
    }

    const int* dimensions = mxGetDimensions(matlabA);
    numRows = dimensions[0];
    numCols = dimensions[1];
    isSparse = mxIsSparse(matlabA);

    if (isSparse) {
      numNonZeros = mxGetNzmax(matlabA);
    }
    else {
      numNonZeros = numRows * numCols;
    }
  }
  
  return(0);
}

//=============================================================================
int EpetraExt_MatlabEngine::GetmxArray(const char* variableName, mxArray** matlabA) {
  if (Comm_.MyPID() == 0) {
#ifdef USE_ENGPUTARRAY
  // for matlab versions < 6.5 (release 13)
  *matlabA = engGetArray(Engine_, variableName);
#else
  // for matlab versions >= 6.5 (release 13)
  *matlabA = engGetVariable(Engine_, variableName);
#endif
    if (matlabA == NULL) {
      return(-5);
    }
  }
  
  return(0);
}

} // namespace EpetraExt
