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

#include "EpetraExt_PutBlockMap.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"

using namespace Matlab;
namespace Matlab {

int CopyBlockMap(mxArray* matlabA, const Epetra_BlockMap& map) {

  int valueCount = 0;
  const Epetra_Comm & comm = map.Comm();
  int numProc = comm.NumProc();
  bool doSizes = !map.ConstantElementSize();

  if (numProc==1) {
    int * myElements = map.MyGlobalElements();
    int * elementSizeList = 0;
    if (doSizes) elementSizeList = map.ElementSizeList();
    return(DoCopyBlockMap(matlabA, valueCount, map.NumGlobalElements(), myElements, elementSizeList, doSizes));
  }

  int numRows = map.NumMyElements();
  
  Epetra_Map allGidsMap(-1, numRows, 0,comm);
  
  Epetra_IntVector allGids(allGidsMap);
  for (int i=0; i<numRows; i++) allGids[i] = map.GID(i);
  
  Epetra_IntVector allSizes(allGidsMap);
  for (int i=0; i<numRows; i++) allSizes[i] = map.ElementSize(i);
  
  // Now construct a Map on PE 0 by strip-mining the rows of the input matrix map.
  int numChunks = numProc;
  int stripSize = allGids.GlobalLength()/numChunks;
  int remainder = allGids.GlobalLength()%numChunks;
  int curStart = 0;
  int curStripSize = 0;
  Epetra_IntSerialDenseVector importGidList;
  Epetra_IntSerialDenseVector importSizeList;
  int numImportGids = 0;
  if (comm.MyPID()==0) {
    importGidList.Size(stripSize+1); // Set size of vector to max needed
    if (doSizes) importSizeList.Size(stripSize+1); // Set size of vector to max needed
  }
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
    Epetra_IntVector importSizes(importGidMap);
    if (doSizes) if (importSizes.Import(allSizes, gidImporter, Insert)) return(-1); 
    
    // importGids (and importSizes, if non-trivial block map)
    // now have a list of GIDs (and sizes, respectively) for the current strip of map.
    
    int * myElements = importGids.Values();
    int * elementSizeList = 0;
    if (doSizes) elementSizeList = importSizes.Values();
    // Finally we are ready to write this strip of the map to file
    if (comm.MyPID()==0) {
      DoCopyBlockMap(matlabA, valueCount, importGids.MyLength(), myElements, elementSizeList, doSizes);
    }
  }
  return(0);
}

int DoCopyBlockMap(mxArray* matlabA, int& valueCount, int length, const int * v1, const int * v2, bool doSizes) {
  // declare and get initial values of all matlabA pointers
  double* matlabAvaluesPtr = mxGetPr(matlabA);
  int* matlabAcolumnIndicesPtr = mxGetJc(matlabA);
  int* matlabArowIndicesPtr = mxGetIr(matlabA);

  // set all matlabA pointers to the proper offset
  matlabAvaluesPtr += valueCount;
  matlabArowIndicesPtr += valueCount;
  int numGidsDone = valueCount;
  if (doSizes) {
    numGidsDone /= 2;
  }
  
  matlabAcolumnIndicesPtr += numGidsDone;
  
  for (int i=0; i<length; i++) {
    *matlabAcolumnIndicesPtr++ = valueCount;
    *matlabArowIndicesPtr++ = 0;
    *matlabAvaluesPtr++ = v1[i]; // row GID of block map

    if (doSizes) {
      *matlabAvaluesPtr++ = v2[i]; // size of block
      valueCount += 2;
      *matlabArowIndicesPtr++ = 1;
    }
    else {
      valueCount++;
    }
  }
  
  return(0);
}
} // namespace Matlab
