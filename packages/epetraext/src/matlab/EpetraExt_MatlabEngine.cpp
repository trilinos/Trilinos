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
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"

using namespace EpetraExt;
namespace EpetraExt {

//=============================================================================
MatlabEngine::MatlabEngine (const Epetra_Comm& Comm):Comm_(Comm) {

    // MATLAB engOpen, to construct the MATLAB engine

    MyPID_ = Comm_.MyPID();
    if (MyPID_ == 0) {
		// Gentlemen, start your engines ...
		Engine_ = engOpen(NULL);
    }
} 

//=============================================================================
MatlabEngine::~MatlabEngine (void) {

    // MATLAB engClose, to destruct the MATLAB engine

    if (MyPID_ == 0) {
		int result = engClose(Engine_);
		if (result != 0) {
	    	// need to handle a nonzero result somehow
		}
    }
}

//=============================================================================
int MatlabEngine::EvalString (char* command) {
	
	// call this.EvalString without an ouputBuffer defined
	return EvalString(command, NULL, -1);
}

int MatlabEngine::EvalString (char* command, char* outputBuffer, int outputBufferSize) {
    // send a string command to the MATLAB engine
    if (MyPID_ == 0) {
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
int MatlabEngine::PutMultiVector(const Epetra_MultiVector & multiVector, const char * variableName) {
    mxArray* matlabA = 0;
    if (comm.MyPID() == 0)
      matlabA = mxCreateDoubleMatrix(A.GlobalLength(), A.NumVectors(), mxREAL);
	
	if (CopyMultiVector(matlabA, multiVector)) {
	  mxDestroyArray(matlabA);
	  return(-2);
	}

	if (comm.MyPID() == 0)
	  if (engPutVariable(Engine_, variableName, matlabA)) {
		mxDestroyArray(matlabA)
		return(-1);
	  }
	
	mxDestroyArray(matlabA);
	return(0);
}

//=============================================================================
int MatlabEngine::PutRowMatrix(const Epetra_RowMatrix& A, const char* variableName, bool transA) {
    mxArray* matlabA = 0;
    if (comm.MyPID() == 0)
	  // since matlab uses column major for matrices, switch row and column numbers
	  mablabA = createSparse(A.NumGlobalCols(), A.NumGlobalRows(), A.NumGlobalNonzeros(), mxREAL);

	if (CopyRowMatrix(matlabA, A)) {
	  mxDestroyArray(matlabA);
	  return(-2);
	}

	if (comm.MyPID() == 0) {
	  if (engPutVariable(Engine_, variableName, matlabA)) {
		mxDestroyArray(matlabA);
		return(-1);
	  }

	  if (!tranA)
		if (EvalString(sprintf("%s = %s'", variableName, variableName))) {
		  mxDestroyArray(matlabA);
		  return(-3);
		}
	}

	mxDestroyArray(matlabA);
	return(0);
}

//=============================================================================
int MatlabEngine::PutCrsGraph(const Epetra_CrsGraph& A, const char* variableName, bool transA) {
	return(-1);
}

int MatlabEngine::PutSerialDenseMatrix(const Epetra_SerialDenseMatrix& A, const char* variableName) {
  int ierr = 0;

	if (comm.MyPID() == 0) {
	  int numRows = A.M();
	  int numCols = A.N();
	  mxArray* matlabA = 0;
	  // need to add global opp to send remote sdMatrix to proc 0
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

	  if (engPutVariable(Engine_, variableName, matlabA))
		ierr = -1;
	}
	
	mxDestroyArray(matlabA);
	return(ierr);
}

//=============================================================================
int MatlabEngine::PutIntSerialDenseMatrix(const Epetra_IntSerialDenseMatrix& A, const char* variableName) {
	int ierr = 0;
	if (comm.MyPID() == 0) {
	  int numRows = A.M();
	  int numCols = A.N();
	  mxArray* matlabA = 0;
	  // need to add global opp to send remote sdMatrix to proc 0
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
	
	  if (engPutVariable(Engine_, variableName, matlabA))
		ierr = -1;
	}

	mxDestroyArray(matlabA);
	return(ierr);
}

//=============================================================================
int MatlabEngine::PutBlockMap(const Epetra_BlockMap& blockMap, const char* variableName) {
	return(-1);
}


} // namespace EpetraExt
