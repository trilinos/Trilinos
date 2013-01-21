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
#include "EpetraExt_BlockMapIn.h"
#include "Epetra_Comm.h"
#include "Epetra_Util.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"
#include "EpetraExt_mmio.h"

using namespace EpetraExt;
namespace EpetraExt {

int MatrixMarketFileToMap( const char *filename, const Epetra_Comm & comm, Epetra_Map * & map) {

  Epetra_BlockMap * bmap;
  if (MatrixMarketFileToBlockMap(filename, comm, bmap)) return(-1);
  map = dynamic_cast<Epetra_Map *>(bmap);
  return(0);
}

int MatrixMarketFileToBlockMap( const char *filename, const Epetra_Comm & comm, Epetra_BlockMap * & map) {

  const int lineLength = 1025;
  char line[lineLength];
  char token[lineLength];
  int M, N, numProc, MaxElementSize, MinElementSize, NumMyElements, IndexBase, NumGlobalElements, firstGid;

  FILE * handle = 0;

  bool inHeader = true;

  handle = fopen(filename,"r");
  if (handle == 0)
    EPETRA_CHK_ERR(-1); // file not found

  while (inHeader) {
    if(fgets(line, lineLength, handle)==0) return(-1);
    if(sscanf(line, "%s", token)==0) return(-1);
    if (!strcmp(token, "%NumProc:")) inHeader = false;
  }

  if(fgets(line, lineLength, handle)==0) return(-1); // numProc value
  if(sscanf(line, "%s %d", token, &numProc)==0) return(-1);

  if(fgets(line, lineLength, handle)==0) return(-1); // MaxElementSize header line
  if(fgets(line, lineLength, handle)==0) return(-1); // MaxElementSize value
  if(sscanf(line, "%s %d", token, &MaxElementSize)==0) return(-1);

  if(fgets(line, lineLength, handle)==0) return(-1); // MinElementSize header line
  if(fgets(line, lineLength, handle)==0) return(-1); // MinElementSize value
  if(sscanf(line, "%s %d", token, &MinElementSize)==0) return(-1);

  if(fgets(line, lineLength, handle)==0) return(-1); // IndexBase header line
  if(fgets(line, lineLength, handle)==0) return(-1); // IndexBase value
  if(sscanf(line, "%s %d", token, &IndexBase)==0) return(-1);

  if(fgets(line, lineLength, handle)==0) return(-1); // NumGlobalElements header line
  if(fgets(line, lineLength, handle)==0) return(-1); // NumGlobalElements value
  if(sscanf(line, "%s %d", token, &NumGlobalElements)==0) return(-1);

  int ierr = 0;
  (void) ierr; // mfh 13 Jan 2013: Forestall compiler warning for unused var.
  if (comm.NumProc()==numProc) {
    if(fgets(line, lineLength, handle)==0) return(-1); // NumMyElements header line
    firstGid = 0;
    for (int i=0; i<comm.MyPID(); i++) {
      if(fgets(line, lineLength, handle)==0) return(-1); // ith NumMyElements value
      if(sscanf(line, "%s %d", token, &NumMyElements)==0) return(-1);
      firstGid += NumMyElements;
    }
 
    if(fgets(line, lineLength, handle)==0) return(-1); // This PE's NumMyElements value
    if(sscanf(line, "%s %d", token, &NumMyElements)==0) return(-1);

    for (int i=comm.MyPID()+1; i<numProc; i++) {
      if(fgets(line, lineLength, handle)==0) return(-1); // ith NumMyElements value (dump these)
    }
  }
  else {
    ierr = 1; // Warning error, different number of processors.

    if(fgets(line, lineLength, handle)==0) return(-1); // NumMyElements header line
    for (int i=0; i<numProc; i++) {
      if(fgets(line, lineLength, handle)==0) return(-1); // ith NumMyElements value (dump these)
    }

    NumMyElements = NumGlobalElements/comm.NumProc();
    firstGid = comm.MyPID()*NumMyElements;
    int remainder = NumGlobalElements%comm.NumProc();
    if (comm.MyPID()<remainder) NumMyElements++;
    int extra = remainder;
    if (comm.MyPID()<remainder) extra = comm.MyPID();
    firstGid += extra;
  }
  if(fgets(line, lineLength, handle)==0) return(-1); // Number of rows, columns
  if(sscanf(line, "%d %d", &M, &N)==0) return(-1);

  bool doSizes = (N>1);
  Epetra_IntSerialDenseVector v1(NumMyElements);
  Epetra_IntSerialDenseVector v2(NumMyElements);
  for (int i=0; i<firstGid; i++) {
    if(fgets(line, lineLength, handle)==0) return(-1); // dump these
  }

  if (doSizes) {
    for (int i=0; i<NumMyElements; i++) {
      if(fgets(line, lineLength, handle)==0) return(-1);
      if(sscanf(line, "%d %d", &v1[i], &v2[i])==0) return(-1); // load v1, v2
    }
  }
  else {
    for (int i=0; i<NumMyElements; i++) {
      if(fgets(line, lineLength, handle)==0) return(-1);
      if(sscanf(line, "%d", &v1[i])==0) return(-1); // load v1
      v2[i] = MinElementSize; // Fill with constant size
    }
  }
  if (fclose(handle)) return(-1);

  comm.Barrier();

  if (MinElementSize==1 && MaxElementSize==1)
    map = new Epetra_Map(-1, NumMyElements, v1.Values(), IndexBase, comm);
  else
    map = new Epetra_BlockMap(-1, NumMyElements, v1.Values(), v2.Values(), IndexBase, comm);
  return(0);
}

int MatrixMarketFileToRowMap(const char* filename,
                             const Epetra_Comm& comm,
                             Epetra_BlockMap*& rowmap)
{
  FILE* infile = fopen(filename, "r");
  MM_typecode matcode;

  int err = mm_read_banner(infile, &matcode);
  if (err != 0) return(err);

  if (!mm_is_matrix(matcode) || !mm_is_coordinate(matcode) ||
      !mm_is_real(matcode)   || !mm_is_general(matcode)) {
    return(-1);
  }

  int numrows, numcols;
  err = mm_read_mtx_array_size(infile, &numrows, &numcols);
  if (err != 0) return(err);

  fclose(infile);

  rowmap = new Epetra_BlockMap(numrows, 1, 0, comm);
  return(0);
}

int MatrixMarketFileToBlockMaps(const char* filename,
                                const Epetra_Comm& comm,
                                Epetra_BlockMap*& rowmap,
                                Epetra_BlockMap*& colmap,
                                Epetra_BlockMap*& rangemap,
                                Epetra_BlockMap*& domainmap)
{
  FILE* infile = fopen(filename, "r");
  if (infile == NULL) {
    return(-1);
  }

  MM_typecode matcode;

  int err = mm_read_banner(infile, &matcode);
  if (err != 0) return(err);

  if (!mm_is_matrix(matcode) || !mm_is_coordinate(matcode) ||
      !mm_is_real(matcode)   || !mm_is_general(matcode)) {
    return(-1);
  }

  int numrows, numcols, nnz;
  err = mm_read_mtx_crd_size(infile, &numrows, &numcols, &nnz);
  if (err != 0) return(err);

  //for this case, we'll assume that the row-map is the same as
  //the range-map.
  //create row-map and range-map with linear distributions.

  rowmap = new Epetra_BlockMap(numrows, 1, 0, comm);
  rangemap = new Epetra_BlockMap(numrows, 1, 0, comm);

  int I, J;
  double val, imag;

  int num_map_cols = 0, insertPoint, foundOffset;
  int allocLen = numcols;
  int* map_cols = new int[allocLen];

  //read through all matrix data and construct a list of the column-
  //indices that occur in rows that are local to this processor.
 
  for(int i=0; i<nnz; ++i) {
    err = mm_read_mtx_crd_entry(infile, &I, &J, &val,
                                &imag, matcode);

    if (err == 0) {
      --I;
      --J;
      if (rowmap->MyGID(I)) {
        foundOffset = Epetra_Util_binary_search(J, map_cols, num_map_cols,
                                                insertPoint);
        if (foundOffset < 0) {
          Epetra_Util_insert(J, insertPoint, map_cols,
                             num_map_cols, allocLen);
        }
      }
    } 
  }

  //create colmap with the list of columns associated with rows that are
  //local to this processor.
  colmap = new Epetra_Map(-1, num_map_cols, map_cols, 0, comm);

  //create domainmap which has a linear distribution
  domainmap = new Epetra_BlockMap(numcols, 1, 0, comm);

  delete [] map_cols;

  return(0);
}

} // namespace EpetraExt

