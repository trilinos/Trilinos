// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// @HEADER
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
#include "Teuchos_Assert.hpp"

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
  A = new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *rowMap, &nnzPerRow[0], true);

  int info;

  for (int sum = 0, i=0; i < nMyRows; i++) {
    if (nnzPerRow[i]) {
      info = A->InsertGlobalValues(iv[sum],nnzPerRow[i],&vv[sum],&jv[sum]);
      TEUCHOS_ASSERT( info==0);
      sum += nnzPerRow[i];
    }
  }

  // Finish up
  if (vectorMap)
    info = A->FillComplete(*vectorMap, *vectorMap);
  else
    info = A->FillComplete();

  TEUCHOS_ASSERT( info==0);

}

#endif
