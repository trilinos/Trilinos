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
#include "EpetraExt_CrsMatrixIn.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_Util.h"

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// Functions to read a MatrixMarket file and load it into a matrix.
// Adapted from 
//    Trilinos/packages/epetraext/src/inout/EpetraExt_CrsMatrixIn.cpp
// Modified by Jon Berry and Karen Devine to make matrix reallocations 
// more efficient; rather than insert each non-zero separately, we 
// collect rows of non-zeros for insertion.
// Modified by Karen Devine and Steve Plimpton to prevent all processors 
// from reading the same file at the same time (ouch!); chunks of the file 
// are read and broadcast by processor zero; each processor then extracts 
// its nonzeros from the chunk, sorts them by row, and inserts them into 
// the matrix.
// The variable "chunk_read" can be changed to modify the size of the chunks
// read from the file.  Larger values of chunk_read lead to faster execution
// and greater memory use.
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

using namespace EpetraExt;
namespace EpetraExt {

static void sort_three(int* list, int *parlista, double *parlistb, 
                       int start, int end);
//////////////////////////////////////////////////////////////////////////////
int MatrixMarketFileToCrsMatrix(const char *filename, 
    const Epetra_Comm & comm, Epetra_CrsMatrix * & A)
{
  EPETRA_CHK_ERR(MatrixMarketFileToCrsMatrixHandle(filename, comm, A));
  return(0);
}

//////////////////////////////////////////////////////////////////////////////
int MatrixMarketFileToCrsMatrix(const char *filename, 
    const Epetra_Comm & comm, Epetra_CrsMatrix * & A, const bool transpose)
{
  EPETRA_CHK_ERR(MatrixMarketFileToCrsMatrixHandle(filename, comm, A,
                                                   0, 0, 0, 0, transpose));
  return(0);
}
//////////////////////////////////////////////////////////////////////////////
  
int MatrixMarketFileToCrsMatrix(const char *filename, 
    const Epetra_Comm & comm, Epetra_CrsMatrix * & A, const bool transpose,
    const bool verbose) 
{
  EPETRA_CHK_ERR(MatrixMarketFileToCrsMatrixHandle(filename, comm, A,
                                                   0, 0, 0, 0,
                                                   transpose, verbose));
  return(0);
}
  
//////////////////////////////////////////////////////////////////////////////
int MatrixMarketFileToCrsMatrix(const char *filename, 
    const Epetra_Map & rowMap, const Epetra_Map& rangeMap, 
    const Epetra_Map& domainMap, Epetra_CrsMatrix * & A, const bool transpose,
    const bool verbose) 
{
  EPETRA_CHK_ERR(MatrixMarketFileToCrsMatrixHandle(filename, 
                                                   rowMap.Comm(), A,
                                                   &rowMap, 0, 
                                                   &rangeMap, &domainMap,
                                                   transpose, verbose));
  return(0);
}

//////////////////////////////////////////////////////////////////////////////
int MatrixMarketFileToCrsMatrix(const char *filename, 
    const Epetra_Map & rowMap, Epetra_CrsMatrix * & A, const bool transpose,
    const bool verbose) 
{
  EPETRA_CHK_ERR(MatrixMarketFileToCrsMatrixHandle(filename, 
                                                   rowMap.Comm(), A,
                                                   &rowMap, 0, 0, 0, 
                                                   transpose, verbose));
  return(0);
}

//////////////////////////////////////////////////////////////////////////////
int MatrixMarketFileToCrsMatrix(const char *filename, 
    const Epetra_Map & rowMap, const Epetra_Map & colMap, 
    Epetra_CrsMatrix * & A, const bool transpose,
    const bool verbose) 
{
  EPETRA_CHK_ERR(MatrixMarketFileToCrsMatrixHandle(filename, 
                                                   rowMap.Comm(), A,
                                                   &rowMap, &colMap, 0, 0,
                                                   transpose, verbose));
  return(0);
}

//////////////////////////////////////////////////////////////////////////////
int MatrixMarketFileToCrsMatrix(const char *filename, 
    const Epetra_Map & rowMap, const Epetra_Map & colMap,
    const Epetra_Map& rangeMap, const Epetra_Map& domainMap, 
    Epetra_CrsMatrix * & A, const bool transpose,
    const bool verbose) 
{
  EPETRA_CHK_ERR(MatrixMarketFileToCrsMatrixHandle(filename, 
                                                   rowMap.Comm(), A, 
                                                   &rowMap, &colMap, 
                                                   &rangeMap, &domainMap,
                                                   transpose, verbose));
  return(0);
}

//////////////////////////////////////////////////////////////////////////////
int MatrixMarketFileToCrsMatrixHandle(const char *filename,
  const Epetra_Comm & comm,
  Epetra_CrsMatrix * & A,
  const Epetra_Map * rowMap,
  const Epetra_Map * colMap,
  const Epetra_Map * rangeMap,
  const Epetra_Map * domainMap,
  const bool transpose,
  const bool verbose
)
{
  const int chunk_read = 500000;  //  Modify this variable to change the
                                  //  size of the chunks read from the file.
  const int headerlineLength = 257;
  const int lineLength = 81;
  const int tokenLength = 35;
  char line[headerlineLength];
  char token1[tokenLength];
  char token2[tokenLength];
  char token3[tokenLength];
  char token4[tokenLength];
  char token5[tokenLength];
  int M, N, NZ;      // Matrix dimensions
  int i;
  int me = comm.MyPID();

  Epetra_Time timer(comm);

  // Make sure domain and range maps are either both defined or undefined 
  if ((domainMap!=0 && rangeMap==0) || (domainMap==0 && rangeMap!=0)) {
    EPETRA_CHK_ERR(-3);
  }

  // check maps to see if domain and range are 1-to-1

  if (domainMap!=0) {
    if (!domainMap->UniqueGIDs()) {EPETRA_CHK_ERR(-2);}
    if (!rangeMap->UniqueGIDs()) {EPETRA_CHK_ERR(-2);}
  }
  else {
    // If domain and range are not specified, 
    // rowMap becomes both and must be 1-to-1 if specified
    if (rowMap!=0) {
      if (!rowMap->UniqueGIDs()) {EPETRA_CHK_ERR(-2);}
    }
  }
      
  FILE * handle = 0;
  if (me == 0) {
    if (verbose) cout << "Reading MatrixMarket file " << filename << endl;
    handle = fopen(filename,"r");  // Open file
    if (handle == 0)
      EPETRA_CHK_ERR(-1); // file not found

    // Check first line, which should be 
    // %%MatrixMarket matrix coordinate real general
    if(fgets(line, headerlineLength, handle)==0) {
      if (handle!=0) fclose(handle); 
      EPETRA_CHK_ERR(-1);
    }
    if(sscanf(line, "%s %s %s %s %s", token1,token2,token3,token4,token5)==0) {
      if (handle!=0) fclose(handle); 
      EPETRA_CHK_ERR(-1);
    }
    if (strcmp(token1, "%%MatrixMarket") ||
        strcmp(token2, "matrix") ||
        strcmp(token3, "coordinate") ||
        strcmp(token4, "real") ||
        strcmp(token5, "general")) {
      if (handle!=0) fclose(handle); 
      EPETRA_CHK_ERR(-1);
    }

    // Next, strip off header lines (which start with "%")
    do {
      if(fgets(line, headerlineLength, handle)==0) {
        if (handle!=0) fclose(handle); 
        EPETRA_CHK_ERR(-1);
      }
    } while (line[0] == '%');

    // Next get problem dimensions: M, N, NZ
    if(sscanf(line, "%d %d %d", &M, &N, &NZ)==0) {
      if (handle!=0) fclose(handle); 
      EPETRA_CHK_ERR(-1);
    }
  }
  comm.Broadcast(&M, 1, 0);
  comm.Broadcast(&N, 1, 0);
  comm.Broadcast(&NZ, 1, 0);

  // Now create matrix using user maps if provided.


  // Now read in chunks of triplets and broadcast them to other processors.
  char *buffer = new char[chunk_read*lineLength];
  int nchunk; 
  int nmillion = 0;
  int nread = 0;
  int rlen;

  // Storage for this processor's nonzeros.
  const int localblock = 100000;
  int localsize = NZ / comm.NumProc() + localblock;
  int *iv = (int *) malloc(localsize * sizeof(int));
  int *jv = (int *) malloc(localsize * sizeof(int));
  double *vv = (double *) malloc(localsize * sizeof(double));
  int lnz = 0;   //  Number of non-zeros on this processor.

  if (!iv || !jv || !vv) 
    EPETRA_CHK_ERR(-1);

  Epetra_Map *rowMap1;
  bool allocatedHere=false;
  if (rowMap != 0) 
    rowMap1 = const_cast<Epetra_Map *>(rowMap);
  else {
    rowMap1 = new Epetra_Map(M, 0, comm);
    allocatedHere = true;
  }
  int ioffset = rowMap1->IndexBase()-1;
  int joffset = (colMap != 0 ? colMap->IndexBase()-1 : ioffset);

  int rowmajor = 1;  // Assume non-zeros are listed in row-major order, 
  int prevrow = -1;  // but test to detect otherwise.  If non-zeros
                     // are row-major, we can avoid the sort.

  while (nread < NZ) {
    if (NZ-nread > chunk_read) nchunk = chunk_read;
    else nchunk = NZ - nread;

    if (me == 0) {
      char *eof;
      rlen = 0;
      for (int i = 0; i < nchunk; i++) {
        eof = fgets(&buffer[rlen],lineLength,handle);
        if (eof == NULL) {
          fprintf(stderr, "%s", "Unexpected end of matrix file.");
          EPETRA_CHK_ERR(-1);
        }
        rlen += strlen(&buffer[rlen]);
      }
      buffer[rlen++]='\n';
    }
    comm.Broadcast(&rlen, 1, 0);
    comm.Broadcast(buffer, rlen, 0);

    buffer[rlen++] = '\0';
    nread += nchunk;

    // Process the received data, saving non-zeros belonging on this proc.
    char *lineptr = buffer;
    for (rlen = 0; rlen < nchunk; rlen++) {
      char *next = strchr(lineptr,'\n');
      int I = atoi(strtok(lineptr," \t\n")) + ioffset;
      int J = atoi(strtok(NULL," \t\n")) + joffset;
      double V = atof(strtok(NULL," \t\n"));
      lineptr = next + 1;
      if (transpose) {
        // swap I and J indices.
        int tmp = I;
        I = J;
        J = tmp;
      }
      if (rowMap1->MyGID(I)) {
        //  This processor keeps this non-zero.
        if (lnz >= localsize) {  
          // Need more memory to store nonzeros.
          localsize += localblock;
          iv = (int *) realloc(iv, localsize * sizeof(int));
          jv = (int *) realloc(jv, localsize * sizeof(int));
          vv = (double *) realloc(vv, localsize * sizeof(double));
        }
        iv[lnz] = I;
        jv[lnz] = J;
        vv[lnz] = V;
        lnz++;
        if (I < prevrow) rowmajor = 0;
        prevrow = I;
      }
    }

    // Status check
    if (nread / 1000000 > nmillion) {
      nmillion++;
      if (verbose && me == 0) cout << nmillion << "M ";
    }
  }

  delete [] buffer;

  if (!rowmajor) {
    // Sort into row-major order (by iv) so can insert entire rows at once.
    // Reorder jv and vv to parallel iv.
    if (verbose && me == 0) cout << endl << "   Sorting local nonzeros" << endl;
    sort_three(iv, jv, vv, 0, lnz-1);
  }

  //  Compute number of entries per local row for use in constructor;
  //  saves reallocs in FillComplete.

  //  Now construct the matrix.
  //  Count number of entries in each row so can use StaticProfile=2.
  if (verbose && me == 0) cout << endl << "   Constructing the matrix" << endl;
  int numRows = rowMap1->NumMyElements();
  int *numNonzerosPerRow = new int[numRows];
  for (i = 0; i < numRows; i++) numNonzerosPerRow[i] = 0;
  for (i = 0; i < lnz; i++) 
    numNonzerosPerRow[rowMap1->LID(iv[i])]++;

  if (rowMap!=0 && colMap !=0) 
    A = new Epetra_CrsMatrix(Copy, *rowMap, *colMap, numNonzerosPerRow);
  else if (rowMap!=0) {
    // Construct with StaticProfile=true since we know numNonzerosPerRow.
    // Less memory will be needed in FillComplete.
    A = new Epetra_CrsMatrix(Copy, *rowMap, numNonzerosPerRow, true);
  }
  else {
    // Construct with StaticProfile=true since we know numNonzerosPerRow.
    // Less memory will be needed in FillComplete.
    A = new Epetra_CrsMatrix(Copy, *rowMap1, numNonzerosPerRow, true);
  }
  A->SetTracebackMode(1);

  // Rows are inserted in ascending global number, and the insertion uses numNonzerosPerRow.
  // Therefore numNonzerosPerRow must be permuted, as it was created in ascending local order.
  int *gidList = new int[numRows];
  for (i = 0; i < numRows; i++) gidList[i] = rowMap1->GID(i);
  Epetra_Util::Sort(true,numRows,gidList,0,0,1,&numNonzerosPerRow);
  delete [] gidList;

  //  Insert the global values into the matrix row-by-row.
  if (verbose && me == 0) cout << "   Inserting global values" << endl;
  for (int sum = 0, i = 0; i < numRows; i++) {
    if (numNonzerosPerRow[i]) {
      int ierr = A->InsertGlobalValues(iv[sum], numNonzerosPerRow[i], 
                                       &vv[sum], &jv[sum]);
      if (ierr<0) EPETRA_CHK_ERR(ierr);
      sum += numNonzerosPerRow[i];
    }
  }

  delete [] numNonzerosPerRow;
  free(iv);
  free(jv);
  free(vv);
    
  if (verbose && me == 0) cout << "   Completing matrix fill" << endl;
  if (rangeMap != 0 && domainMap != 0) {
    EPETRA_CHK_ERR(A->FillComplete(*domainMap, *rangeMap));
  }
  else if (M!=N) {
    Epetra_Map newDomainMap(N,rowMap1->IndexBase(), comm);
    EPETRA_CHK_ERR(A->FillComplete(newDomainMap, *rowMap1));
  }
  else {
    EPETRA_CHK_ERR(A->FillComplete());
  }

  if (allocatedHere) delete rowMap1;
  
  if (handle!=0) fclose(handle);
  double dt = timer.ElapsedTime();
  if (verbose && me == 0) cout << "File Read time (secs):  " << dt << endl;
  return(0);
}

///////////////////////////////////////////////////////////////////////////
// Sorting values in array list in increasing order. Criteria is int. 
// Also rearrange values in arrays parlista and partlistb
// to match the new order of list. 

static void quickpart_list_inc_int (
  int *list, int *parlista, double *parlistb,
  int start, int end, int *equal, int *larger)
{
int i, key, itmp;
double dtmp;

  key = list ? list[(end+start)/2] : 1;

  *equal = *larger = start;
  for (i = start; i <= end; i++)
    if (list[i] < key) {
      itmp                = parlista[i];
      parlista[i]         = parlista[*larger];
      parlista[(*larger)] = parlista[*equal];
      parlista[(*equal)]  = itmp;
      dtmp                = parlistb[i];
      parlistb[i]         = parlistb[*larger];
      parlistb[(*larger)] = parlistb[*equal];
      parlistb[(*equal)]  = dtmp;
      itmp                = list[i];
      list[i]             = list[*larger];
      list[(*larger)++]   = list[*equal];
      list[(*equal)++]    = itmp;
    }
    else if (list[i] == key) {
      itmp                = parlista[i];
      parlista[i]         = parlista[*larger];
      parlista[(*larger)] = itmp;
      dtmp                = parlistb[i];
      parlistb[i]         = parlistb[*larger];
      parlistb[(*larger)] = dtmp;
      list[i]             = list[*larger];
      list[(*larger)++]   = key;
    }
}

///////////////////////////////////////////////////////////////////////////
static void sort_three(int* list, int *parlista, double *parlistb, 
                       int start, int end)
{
int  equal, larger;

  if (start < end) {
    quickpart_list_inc_int(list, parlista, parlistb, start, end, 
                           &equal, &larger);
    sort_three(list, parlista, parlistb, start,  equal-1);
    sort_three(list, parlista, parlistb, larger, end);
  }
}

///////////////////////////////////////////////////////////////////////////
int MatlabFileToCrsMatrix(const char *filename,
				const Epetra_Comm & comm,
				Epetra_CrsMatrix * & A)
{
  const int lineLength = 1025;
  char line[lineLength];
  int I, J;
  double V;

  FILE * handle = 0;

  handle = fopen(filename,"r");  // Open file
  if (handle == 0)
    EPETRA_CHK_ERR(-1); // file not found

  int numGlobalRows = 0;
  int numGlobalCols = 0;
  while(fgets(line, lineLength, handle)!=0) {
    if(sscanf(line, "%d %d %lg\n", &I, &J, &V)==0) {if (handle!=0) fclose(handle); EPETRA_CHK_ERR(-1);}
    if (I>numGlobalRows) numGlobalRows = I;
    if (J>numGlobalCols) numGlobalCols = J;
  }

  if (handle!=0) fclose(handle);
  Epetra_Map rangeMap(numGlobalRows, 0, comm);
  Epetra_Map domainMap(numGlobalCols, 0, comm);
  A = new Epetra_CrsMatrix(Copy, rangeMap, 0);

  // Now read in each triplet and store to the local portion of the matrix if the row is owned.
  const Epetra_Map & rowMap1 = A->RowMap();
  
  handle = 0;

  handle = fopen(filename,"r");  // Open file
  if (handle == 0)
    EPETRA_CHK_ERR(-1); // file not found

  while (fgets(line, lineLength, handle)!=0) {
    if(sscanf(line, "%d %d %lg\n", &I, &J, &V)==0) {if (handle!=0) fclose(handle); EPETRA_CHK_ERR(-1);}
    I--; J--; // Convert to Zero based
    if (rowMap1.MyGID(I)) {
      int ierr = A->InsertGlobalValues(I, 1, &V, &J);
      if (ierr<0) EPETRA_CHK_ERR(ierr);
    }
  }
    
  EPETRA_CHK_ERR(A->FillComplete(domainMap, rangeMap));

  if (handle!=0) fclose(handle);
  return(0);
}

int HypreFileToCrsMatrix(const char *filename, const Epetra_Comm &comm, Epetra_CrsMatrix *&Matrix){
  int MyPID = comm.MyPID();
  // This double will be in the format we want for the extension besides the leading zero
  double filePID = (double)MyPID/(double)100000;
  std::ostringstream stream;
  // Using setprecision() puts it in the string
  stream << std::setiosflags(std::ios::fixed) << std::setprecision(5) << filePID;
  // Then just ignore the first character
  std::string fileName(filename);
  fileName += stream.str().substr(1,7);
  // Open the file
  std::ifstream file(fileName.c_str());
  string line;
  if(file.is_open()){
    std::getline(file, line);
    int ilower, iupper;
    std::istringstream istream(line);
    // The first line of the file has the beginning and ending rows
    istream >> ilower;
    istream >> iupper;
    // Using those we can create a row map
    Epetra_Map RowMap(-1, iupper-ilower+1, 0, comm);
    Matrix = new Epetra_CrsMatrix(Copy, RowMap, 0);
    int currRow = -1;
    int counter = 0;
    std::vector<int> indices;
    std::vector<double> values;
    while(!file.eof()){
      std::getline(file, line);
      std::istringstream lineStr(line);
      int row, col;
      double val;
      lineStr >> row;
      lineStr >> col;
      lineStr >> val;
      if(currRow == -1) currRow = row; // First line
      if(row == currRow){
        // add to the vector
        counter = counter + 1;
        indices.push_back(col);
        values.push_back(val);
      } else {
        Matrix->InsertGlobalValues(currRow, counter, &values[0], &indices[0]);
        indices.clear();
        values.clear();
        counter = 0;
        currRow = row;
        // make a new vector
        indices.push_back(col);
        values.push_back(val);
        counter = counter + 1;
      }
    }
    Matrix->InsertGlobalValues(currRow, counter, &values[0], &indices[0]);
    Matrix->Comm().Barrier();
    Matrix->FillComplete();
    file.close();
    return 0;
  } else {
    cout << "\nERROR:\nCouldn't open " << fileName << ".\n";
    return -1;
  }
}

} // namespace EpetraExt

