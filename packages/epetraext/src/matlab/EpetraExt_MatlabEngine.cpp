
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
	int err = 0;
    // send a string command to the MATLAB engine
    if (MyPID_ == 0) {
		if (outputBuffer != NULL) {
			err = engOutputBuffer(Engine_, outputBuffer, outputBufferSize);
			if (err != 0) {
				return(-1);
			}
		}
		err = engEvalString(Engine_, command);
		return err;
    }
	
	return(0);
}

//=============================================================================
int MatlabEngine::PutMultiVector(const Epetra_MultiVector & multiVector, const char * variableName) {
	mxArray * matlabMatrix = 0;
	matlabMatrix = mxCreateDoubleMatrix(multiVector.MyLength(), multiVector.NumVectors(), mxREAL);
	multiVector.ExtractCopy((double *)mxGetPr(matlabMatrix), multiVector.MyLength());
	engPutVariable(Engine_, variableName, matlabMatrix);
	
	mxDestroyArray(matlabMatrix);
	return(0);
}

//=============================================================================
int MatlabEngine::PutRowMatrix(const Epetra_RowMatrix& rowMatrix, const char* variableName) {
	return(-1);
}

//=============================================================================
int MatlabEngine::PutCrsGraph(const Epetra_CrsGraph& crsGraph, const char* variableName) {
	return(-1);
}

int MatlabEngine::PutSerialDenseMatrix(const Epetra_SerialDenseMatrix& sdMatrix, const char* variableName) {
	int numRows = sdMatrix.M();
	int numCols = sdMatrix.N();
	mxArray* matlabMatrix = 0;
	// need to add global opp to send remote sdMatrix to proc 0
	matlabMatrix = mxCreateDoubleMatrix(numRows, numCols, mxREAL);

	int row;
	int col;
	double* targetPtr = 0;
	double* sourcePtr = 0;
	double* source = (double*)sdMatrix.A();
    double* target = (double*)mxGetPr(matlabMatrix);
	int source_LDA = sdMatrix.LDA();
	int target_LDA = sdMatrix.LDA();
	for (col = 0; col < numCols; col++) {
		targetPtr = target + (col * target_LDA);
		sourcePtr = source + (col * source_LDA);
		for (row = 0; row < numRows; row++) {
			*targetPtr++ = *sourcePtr++;
		}
	}
	
	int err = engPutVariable(Engine_, variableName, matlabMatrix);

	mxDestroyArray(matlabMatrix);
	return(err);
}

//=============================================================================
int MatlabEngine::PutIntSerialDenseMatrix(const Epetra_IntSerialDenseMatrix& isdMatrix, const char* variableName) {
	int numRows = isdMatrix.M();
	int numCols = isdMatrix.N();
	mxArray* matlabMatrix = 0;
	// need to add global opp to send remote sdMatrix to proc 0
	matlabMatrix = mxCreateDoubleMatrix(numRows, numCols, mxREAL);

	int row;
	int col;
	double* targetPtr = 0;
	int* sourcePtr = 0;
	int* source = (int*)isdMatrix.A();
    double* target = (double*)mxGetPr(matlabMatrix);
	int source_LDA = isdMatrix.LDA();
	int target_LDA = isdMatrix.LDA();
	for (col = 0; col < numCols; col++) {
		targetPtr = target + (col * target_LDA);
		sourcePtr = source + (col * source_LDA);
		for (row = 0; row < numRows; row++) {
			*targetPtr++ = *sourcePtr++;
		}
	}
	
	int err = engPutVariable(Engine_, variableName, matlabMatrix);

	mxDestroyArray(matlabMatrix);
	return(err);
}

//=============================================================================
int MatlabEngine::PutBlockMap(const Epetra_BlockMap& blockMap, const char* variableName) {
	return(-1);
}


} // namespace EpetraExt
