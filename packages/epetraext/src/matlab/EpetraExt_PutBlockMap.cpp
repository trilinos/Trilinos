//@HEADER
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

  const Epetra_Comm & comm = map.Comm();
  int numProc = comm.NumProc();
  bool doSizes = !map.ConstantElementSize();

  if (numProc==1) {
    int * myElements = map.MyGlobalElements();
    int * elementSizeList = 0;
    if (doSizes) elementSizeList = map.ElementSizeList();
    return(DoCopyBlockMap(matlabA, map.NumGlobalElements(), myElements, elementSizeList, doSizes));
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
    DoCopyBlockMap(matlabA, importGids.MyLength(), myElements, elementSizeList, doSizes);
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
  matlabAcolumnIndicesPtr += valueCount;
  
  for (int i=0; i<length; i++) {
    *matlabAcolumnIndicesPtr++ = valueCount;
    *matlabArowIndicesPtr++ = 0;
    *matlabAvaluesPtr++ = v1[i]; # row GID of block map
    if (doSizes) {
      *matlabAvaluesPtr++ = v2[i]; # size of block
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
