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

#include "EpetraExt_PutMultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"

using namespace Matlab;
namespace Matlab {
int CopyMultiVector(double** matlabApr, const Epetra_MultiVector& A) {

  Epetra_BlockMap bmap = A.Map();
  const Epetra_Comm & comm = bmap.Comm();
  int numProc = comm.NumProc();

  if (numProc==1)
    DoCopyMultiVector(matlabApr, A);
  else {

    // In the case of more than one column in the multivector, and writing to MatrixMarket
    // format, we call this function recursively, passing each vector of the multivector
    // individually so that we can get all of it written to file before going on to the next 
    // multivector
    if (A.NumVectors() > 1) {
      for (int i=0; i < A.NumVectors(); i++)
	if (CopyMultiVector(matlabApr, *(A(i)))) return(-1);
      return(0);
    }

    Epetra_Map map(-1, bmap.NumMyPoints(), 0, comm);
    // Create a veiw of this multivector using a map (instead of block map)
    Epetra_MultiVector A1(View, map, A.Pointers(), A.NumVectors());
    int numRows = map.NumMyElements();
    
    Epetra_Map allGidsMap(-1, numRows, 0,comm);
    
    Epetra_IntVector allGids(allGidsMap);
    for (int i=0; i<numRows; i++) allGids[i] = map.GID(i);
    
    // Now construct a MultiVector on PE 0 by strip-mining the rows of the input matrix A.
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
      Epetra_Map importGidMap(-1, curStripSize, importGidList.Values(), 0, comm);
      Epetra_Import gidImporter(importGidMap, allGidsMap);
      Epetra_IntVector importGids(importGidMap);
      if (importGids.Import(allGids, gidImporter, Insert)) return(-1); 

      // importGids now has a list of GIDs for the current strip of matrix rows.
      // Use these values to build another importer that will get rows of the matrix.

      // The following import map will be non-trivial only on PE 0.
      Epetra_Map importMap(-1, importGids.MyLength(), importGids.Values(), 0, comm);
      Epetra_Import importer(importMap, map);
      Epetra_MultiVector importA(importMap, A1.NumVectors());
      if (importA.Import(A1, importer, Insert)) return(-1); 

      // Finally we are ready to write this strip of the matrix to ostream
      if (DoCopyMultiVector(matlabApr, importA)) return(-1);
    }
  }
  return(0);
}

int DoCopyMultiVector(double** matlabApr, const Epetra_MultiVector& A) {

  int ierr = 0;
  int length = A.GlobalLength();
  int numVectors = A.NumVectors();
  const Epetra_Comm & comm = A.Map().Comm();
  if (comm.MyPID()!=0) {
    if (A.MyLength()!=0) ierr = -1;
  }
  else {
    if (length!=A.MyLength()) ierr = -1;
    double* matlabAvalues = *matlabApr;
    double* Aptr = A.Values();
    memcpy((void *)matlabAvalues, (void *)Aptr, sizeof(*Aptr) * length * numVectors);
    *matlabApr += length;   
  }
  int ierrGlobal;
  comm.MinAll(&ierr, &ierrGlobal, 1); // If any processor has -1, all return -1
  return(ierrGlobal);
}
} // namespace Matlab
