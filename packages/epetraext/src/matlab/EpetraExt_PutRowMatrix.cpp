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

#include "EpetraExt_PutRowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

using namespace Matlab;
namespace Matlab {

int CopyRowMatrix(mxArray* matlabA, const Epetra_RowMatrix& A) {
  int valueCount = 0;
  //int* valueCount = &temp;

  Epetra_Map map = A.RowMatrixRowMap();
  const Epetra_Comm & comm = map.Comm();
  int numProc = comm.NumProc();

  if (numProc==1) 
    DoCopyRowMatrix(matlabA, valueCount, A);
  else {
    int numRows = map.NumMyElements();
    
    //cout << "creating allGidsMap\n";
    Epetra_Map allGidsMap(-1, numRows, 0,comm);
    //cout << "done creating allGidsMap\n";
    
    Epetra_IntVector allGids(allGidsMap);
    for (int i=0; i<numRows; i++) allGids[i] = map.GID(i);
    
    // Now construct a RowMatrix on PE 0 by strip-mining the rows of the input matrix A.
    int numChunks = numProc;
    int stripSize = allGids.GlobalLength()/numChunks;
    int remainder = allGids.GlobalLength()%numChunks;
    int curStart = 0;
    int curStripSize = 0;
    Epetra_IntSerialDenseVector importGidList;
    int numImportGids = 0;
    if (comm.MyPID()==0) 
      importGidList.Size(stripSize+1); // Set size of vector to max needed
    for (int i=0; i<numChunks; i++) {
      if (comm.MyPID()==0) { // Only PE 0 does this part
	curStripSize = stripSize;
	if (i<remainder) curStripSize++; // handle leftovers
	for (int j=0; j<curStripSize; j++) importGidList[j] = j + curStart;
	curStart += curStripSize;
      }
      // The following import map will be non-trivial only on PE 0.
      //cout << "creating importGidMap\n";
      Epetra_Map importGidMap(-1, curStripSize, importGidList.Values(), 0, comm);
      //cout << "done creating importGidMap\n";
      Epetra_Import gidImporter(importGidMap, allGidsMap);
      Epetra_IntVector importGids(importGidMap);
      if (importGids.Import(allGids, gidImporter, Insert)) return(-1); 

      // importGids now has a list of GIDs for the current strip of matrix rows.
      // Use these values to build another importer that will get rows of the matrix.

      // The following import map will be non-trivial only on PE 0.
      //cout << "creating importMap\n";
      //cout << "A.RowMatrixRowMap().MinAllGID: " << A.RowMatrixRowMap().MinAllGID() << "\n";
      Epetra_Map importMap(-1, importGids.MyLength(), importGids.Values(), A.RowMatrixRowMap().MinAllGID(), comm);
      //cout << "done creating importMap\n";
      Epetra_Import importer(importMap, map);
      Epetra_CrsMatrix importA(Copy, importMap, 0);
      if (importA.Import(A, importer, Insert)) return(-1); 
      if (importA.FillComplete()) return(-1);

      // Finally we are ready to write this strip of the matrix to ostream
      if (DoCopyRowMatrix(matlabA, valueCount, importA)) return(-1);
    }
  }

  if (A.RowMatrixRowMap().Comm().MyPID() == 0) {
	// set max cap
	int* matlabAcolumnIndicesPtr = mxGetJc(matlabA);
	matlabAcolumnIndicesPtr[A.NumGlobalRows()] = valueCount;
  }

  return(0);
}

int DoCopyRowMatrix(mxArray* matlabA, int& valueCount, const Epetra_RowMatrix& A) {
  //cout << "doing DoCopyRowMatrix\n";
  int ierr = 0;
  int numRows = A.NumGlobalRows();
  //cout << "numRows: " << numRows << "\n";
  Epetra_Map rowMap = A.RowMatrixRowMap();
  Epetra_Map colMap = A.RowMatrixColMap();
  int minAllGID = rowMap.MinAllGID();

  const Epetra_Comm & comm = rowMap.Comm();
  //cout << "did global setup\n";
  if (comm.MyPID()!=0) {
    if (A.NumMyRows()!=0) ierr = -1;
    if (A.NumMyCols()!=0) ierr = -1;
  }
  else {
	// declare and get initial values of all matlabA pointers
	double* matlabAvaluesPtr = mxGetPr(matlabA);
	int* matlabAcolumnIndicesPtr = mxGetJc(matlabA);
	int* matlabArowIndicesPtr = mxGetIr(matlabA);

	// set all matlabA pointers to the proper offset
	matlabAvaluesPtr += valueCount;
	matlabArowIndicesPtr += valueCount;

    if (numRows!=A.NumMyRows()) ierr = -1;
    Epetra_SerialDenseVector values(A.MaxNumEntries());
    Epetra_IntSerialDenseVector indices(A.MaxNumEntries());
    //cout << "did proc0 setup\n";
    for (int i=0; i<numRows; i++) {
	  //cout << "extracting a row\n";
	  int I = rowMap.GID(i);
      int numEntries = 0;
      if (A.ExtractMyRowCopy(i, values.Length(), numEntries, 
	  		     values.Values(), indices.Values())) return(-1);
	  matlabAcolumnIndicesPtr[I - minAllGID] = valueCount;  // set the starting index of column I
	  double* serialValuesPtr = values.Values();
      for (int j=0; j<numEntries; j++) {
		int J = colMap.GID(indices[j]);
		*matlabAvaluesPtr = *serialValuesPtr++;
		*matlabArowIndicesPtr = J;
		// increment matlabA pointers
		matlabAvaluesPtr++;
		matlabArowIndicesPtr++;
		valueCount++;
      }
    }
    //cout << "proc0 row extraction for this chunck is done\n";
  }

/*
  if (comm.MyPID() == 0) {
  cout << "printing matlabA pointers\n";
	double* matlabAvaluesPtr = mxGetPr(matlabA);
	int* matlabAcolumnIndicesPtr = mxGetJc(matlabA);
	int* matlabArowIndicesPtr = mxGetIr(matlabA);
  for(int i=0; i < numRows; i++) {
	for(int j=0; j < A.MaxNumEntries(); j++) {
	  cout << "*matlabAvaluesPtr: " << *matlabAvaluesPtr++ << " *matlabAcolumnIndicesPtr: " << *matlabAcolumnIndicesPtr++ << " *matlabArowIndicesPtr" << *matlabArowIndicesPtr++ << "\n";
	}
  }
  
  cout << "done printing matlabA pointers\n";
  }
  */
  
  int ierrGlobal;
  comm.MinAll(&ierr, &ierrGlobal, 1); // If any processor has -1, all return -1
  return(ierrGlobal);
}

} // namespace Matlab
