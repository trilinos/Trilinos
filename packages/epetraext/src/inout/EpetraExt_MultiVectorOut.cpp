//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_mmio.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"

using namespace EpetraExt;
namespace EpetraExt {

int MultiVectorToMatlabFile( const char *filename, const Epetra_MultiVector & A) {

  // Simple wrapper to make it clear what can be used to write to Matlab format
  return(MultiVectorToMatrixMarketFile(filename, A, 0, 0, false));
}

int MultiVectorToMatrixMarketFile( const char *filename, const Epetra_MultiVector & A, 
				 const char * matrixName,
				 const char *matrixDescription, 
				 bool writeHeader) {
  int M = A.NumVectors();
  int N = A.GlobalLength();

  FILE * handle = 0;

  if (A.Map().Comm().MyPID()==0) { // Only PE 0 does this section

    handle = fopen(filename,"w");
    if (!handle) return(-1);
    MM_typecode matcode;
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_array(&matcode);
    mm_set_real(&matcode);

    if (writeHeader==true) { // Only write header if requested (true by default)
    
      if (mm_write_banner(handle, matcode)) return(-1);
      
      if (matrixName!=0) fprintf(handle, "%% \n%% %s\n", matrixName);
      if (matrixDescription!=0) fprintf(handle, "%% %s\n%% \n", matrixDescription);
      
      if (mm_write_mtx_array_size(handle, M, N)) return(-1);
    }
  }
    
  if (MultiVectorToHandle(handle, A)) return(-1); // Everybody calls this routine

  if (A.Map().Comm().MyPID()==0) // Only PE 0 opened a file
    if (fclose(handle)) return(-1);
  return(0);
}

int MultiVectorToHandle(FILE * handle, const Epetra_MultiVector & A) {

  Epetra_BlockMap bmap = A.Map();
  const Epetra_Comm & comm = bmap.Comm();
  int numProc = comm.NumProc();

  if (numProc==1)
    writeMultiVector(handle, A);
  else {

    Epetra_Map map(-1, bmap.NumMyPoints(), 0, comm);
    Epetra_MultiVector A1(View, map, A.Pointers(), A.NumVectors()); // Create a veiw of this multivector using a map (instead of block map)
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
      if (writeMultiVector(handle, importA)) return(-1);
    }
  }
  return(0);
}
int writeMultiVector(FILE * handle, const Epetra_MultiVector & A) {

  // The multivector will be written to file in transpose form, so that
 //  we can write all contributions from a processor at once
  int ierr = 0;
  int length = A.GlobalLength();
  int numVectors = A.NumVectors();
  const Epetra_Comm & comm = A.Map().Comm();
  if (comm.MyPID()!=0) {
    if (A.MyLength()!=0) ierr = -1;
  }
  else {
    if (length!=A.MyLength()) ierr = -1;
    for (int i=0; i<length; i++) {
      for (int j=0; j<numVectors; j++) {
	double val = A[j][i];
	fprintf(handle, "%22.16e\n", val);
      }
    }
  }
  int ierrGlobal;
  comm.MinAll(&ierr, &ierrGlobal, 1); // If any processor has -1, all return -1
  return(ierrGlobal);
}
} // namespace EpetraExt
