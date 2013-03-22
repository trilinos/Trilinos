#ifndef __BUILD_SIMPLE_MATRICES_HPP
#define __BUILD_SIMPLE_MATRICES_HPP

#include <stdio.h>


#ifdef EPETRA_MPI
#include <mpi.h>
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
 
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

#ifdef KDDKDD_DEBUG
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#endif

#include "build_maps.hpp"

////////////////////////////////////////////////////////////////////////////

// Build the following Laplacian matrix with nGlobalRows rows.
//
//    | 1 -1          |
//    |-1  2 -1       |
//    |   -1  2 -1    |
//    |        ...    |
//    |           -1 1|
//
// Distribute the matrix according to the provided row and column maps.

template <typename itype>
void build_simple_matrix(
  Epetra_Comm &comm,         // Communicator to use
  Epetra_CrsMatrix *&A,      // OUTPUT:  Matrix returned
  itype nGlobalRows,         // Number of global matrix rows and columns
  bool testEpetra64,         // if true, add 2*INT_MAX to each global ID
                             // to exercise Epetra64
  bool verbose               // if true, print out matrix information
)
{
  Epetra_Map *rowMap = NULL;        // Row map for the created matrix
  Epetra_Map *colMap = NULL;        // Col map for the created matrix
  Epetra_Map *vectorMap = NULL;     // Range/Domain map for the created matrix

  long long offsetEpetra64;

  build_maps(nGlobalRows, testEpetra64, comm, 
             &vectorMap, &rowMap, &colMap, offsetEpetra64, verbose);

  // Create an integer vector nnzPerRow that is used to build the Epetra Matrix.
  // nnzPerRow[i] is the number of entries for the ith global equation
  int nMyRows = rowMap->NumMyElements();
  std::vector<int> nnzPerRow(nMyRows+1, 0);

  // Also create lists of the nonzeros to be assigned to processors.
  // To save programming time and complexity, these vectors are allocated 
  // bigger than they may actually be needed.
  std::vector<itype> iv(3*nMyRows+1);
  std::vector<itype> jv(3*nMyRows+1);
  std::vector<double> vv(3*nMyRows+1);

  itype nMyNonzeros = 0;
  for (itype i = 0, myrowcnt = 0; i < nGlobalRows; i++) {
    if (rowMap->MyGID(i+offsetEpetra64)) { 
      // This processor owns part of this row; see whether it owns the nonzeros
      if (i > 0 && (!colMap || colMap->MyGID(i-1+offsetEpetra64))) {
        iv[nMyNonzeros] = i + offsetEpetra64;
        jv[nMyNonzeros] = i-1 + offsetEpetra64;
        vv[nMyNonzeros] = -1;
        nMyNonzeros++;
        nnzPerRow[myrowcnt]++;
      }
      if (!colMap || colMap->MyGID(i+offsetEpetra64)) {
        iv[nMyNonzeros] = i + offsetEpetra64;
        jv[nMyNonzeros] = i + offsetEpetra64;
        vv[nMyNonzeros] = ((i == 0 || i == nGlobalRows-1) ? 1. : 2.);
        nMyNonzeros++;
        nnzPerRow[myrowcnt]++;
      }
      if (i < nGlobalRows - 1 && (!colMap ||  colMap->MyGID(i+1+offsetEpetra64))) {
        iv[nMyNonzeros] = i + offsetEpetra64;
        jv[nMyNonzeros] = i+1 + offsetEpetra64;
        vv[nMyNonzeros] = -1;
        nMyNonzeros++;
        nnzPerRow[myrowcnt]++;
      }
      myrowcnt++;
    }
  }

  // Create an Epetra_Matrix
  A = new Epetra_CrsMatrix(Copy, *rowMap, &nnzPerRow[0], true);

  int info;

  for (int sum = 0, i=0; i < nMyRows; i++) {
    if (nnzPerRow[i]) {
      info = A->InsertGlobalValues(iv[sum],nnzPerRow[i],&vv[sum],&jv[sum]);
      assert(info==0);
      sum += nnzPerRow[i];
    }
  }

  // Finish up
  if (vectorMap)
    info = A->FillComplete(*vectorMap, *vectorMap);
  else
    info = A->FillComplete();

  assert(info==0);

}

#endif
