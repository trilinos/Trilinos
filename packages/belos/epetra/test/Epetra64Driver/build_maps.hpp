// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//  Build maps for 1D or 2D matrix distribution
//  Karen Devine, SNL
//

#ifndef __BUILD_MAPS_HPP
#define __BUILD_MAPS_HPP

#include <cstdio>
#include <cstdlib>
#include <iostream>
#ifndef __cplusplus
#define __cplusplus
#endif

#include "Epetra_Comm.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"

#define TWODPROW(p, npRows, npCols) ((p) % (npRows))
#define TWODPCOL(p, npRows, npCols) ((p) / (npRows))
#define TWODPRANK(i, j, npRows, npCols) ((j) * (npRows) + ((j)+(i)) % (npRows))
#define MIN(i,j) ((i) < (j) ? (i) : (j))

///////////////////////////////////////////////////////////////////////////////
template <typename itype>
void random_distribution_1D(
  itype nrows,          // Number of global matrix rows
  Epetra_Comm &comm,    // Epetra communicator to be used in maps
  Epetra_Map **rowMap,  // OUTPUT: pointer to row map to be created
  long long offsetEpetra64
)
{
  // Randomly assign matrix rows to processor's row Map.

  int me = comm.MyPID();
  int np = comm.NumProc();

  std::vector<itype> myGlobalElements(1.2 * (nrows / np) + 1);
  int nMyRows = 0;
  srandom(1);
  double denom = (double) RAND_MAX + 1.;
  for (itype i = 0; i < nrows; i++) {
    int p = (int) ((double) np * (double) random() / denom);
    if (p == me) {
      if (nMyRows >= myGlobalElements.size())
        myGlobalElements.resize(1.5*myGlobalElements.size());
      myGlobalElements[nMyRows] = i + offsetEpetra64;
      nMyRows++;
    }
  }
  *rowMap = new Epetra_Map(nrows, nMyRows, &myGlobalElements[0], 0, comm);
}

///////////////////////////////////////////////////////////////////////////////
template <typename itype>
void random_distribution_2D(
  int mypRow,      // processor's row (in 2D)
  int mypCol,      // processor's col (in 2D)
  int npRows,      // number of processor rows
  int npCols,      // number of processor cols
  itype nrows,     // Number of global matrix rows
  Epetra_Comm &comm,    // Epetra communicator to be used in maps
  Epetra_Map **vectorMap,  // OUTPUT: Map to be used for the vector
  Epetra_Map **rowMap,     // OUTPUT: Map to be used for the matrix rows
  Epetra_Map **colMap,     // OUTPUT: Map to be used for the matrix cols
  long long offsetEpetra64
)
{
  // Randomly assign matrix rows to processor's vector Map.
  // Build appropriate GlobalElements lists for row and column maps at the
  // same time.

  int me = comm.MyPID();
  int np = comm.NumProc();

  int nMyEntries = 0;
  int nMyRows = 0;
  int nMyCols = 0;
  std::vector<itype> myGlobalElements(1.2 * nrows / np + 1);
  std::vector<itype> myGlobalRowElements(1.2 * nrows / npRows + 1);
  std::vector<itype> myGlobalColElements(1.2 * nrows / npCols + 1);

  srandom(1);
  double denom = (double) RAND_MAX + 1.;

  for (itype i = 0; i < nrows; i++) {
    // Compute rank to receive the vector entry i
    int p = (int) ((double) np * (double) random() / denom);

    if (p == me) {
      // Add entry i to my vector map
      if (nMyEntries >= myGlobalElements.size())
        myGlobalElements.resize(1.5*myGlobalElements.size());
      myGlobalElements[nMyEntries] = i + offsetEpetra64;
      nMyEntries++;
    }

    if (mypRow == TWODPROW(p, npRows, npCols)) {
      // Add entry i to my row map
      if (nMyRows >= myGlobalRowElements.size())
        myGlobalRowElements.resize(1.5*myGlobalRowElements.size());
      myGlobalRowElements[nMyRows] = i + offsetEpetra64;
      nMyRows++;
    }

    if (mypCol == TWODPCOL(p, npRows, npCols)) {
      // Add entry i to my col map
      if (nMyCols >= myGlobalColElements.size())
        myGlobalColElements.resize(1.5*myGlobalColElements.size());
      myGlobalColElements[nMyCols] = i + offsetEpetra64;
      nMyCols++;
    }
  }

  *vectorMap = new Epetra_Map(-1, nMyEntries, &myGlobalElements[0], 0, comm);
  *rowMap = new Epetra_Map(-1, nMyRows, &myGlobalRowElements[0], 0, comm);
  *colMap = new Epetra_Map(-1, nMyCols, &myGlobalColElements[0], 0, comm);
}

///////////////////////////////////////////////////////////////////////////////
template <typename itype>
void linear_distribution_2D(
  int mypRow,      // processor's row (in 2D)
  int mypCol,      // processor's col (in 2D)
  int npRows,      // number of processor rows
  int npCols,      // number of processor cols
  itype nrows,     // Number of global matrix rows
  Epetra_Comm &comm,    // Epetra communicator to be used in maps
  Epetra_Map **vectorMap,  // OUTPUT: Map to be used for the vector
  Epetra_Map **rowMap,     // OUTPUT: Map to be used for the matrix rows
  Epetra_Map **colMap,     // OUTPUT: Map to be used for the matrix cols
  long long offsetEpetra64
)
{
  // The vector will be distributed linearly:
  //  [0 1 2 3 4 5 6 7 8]
  //
  // If nrows is not divisible by np, extra rows will be distributed among
  // the processors.  (See below.)

  int me = comm.MyPID();
  int np = comm.NumProc();

  // Create vector map first

  std::vector<itype> entries(np+1, nrows / np); // Initial # entries per proc
  int nExtraEntries = nrows % np;

  // Distribute the extra entries evenly among processors.
  // To evenly distribute them extra entries among processor rows and
  // columns, we distribute them along diagonals of the matrix distribution.
  // For our example, assume we have seven extra values (the max possible
  // with np=8).  Then we give one extra entry each to ranks
  // [0, 3, 4, 7, 1, 2, 5].  For fewer extra entries, we follow the same
  // order of assignment, and just stop early.

  for (int cnt = 0, i = 0; (cnt < nExtraEntries) && (i < npRows); i++) {
    for (int j = 0; (cnt < nExtraEntries) && (j < npCols); cnt++, j++) {
      int rankForExtra = TWODPRANK(i, j, npRows, npCols);
      entries[rankForExtra+1]++;  // Store in rankForExtra+1 to simply later
                                   // prefix sum.
    }
  }

  // Perform prefix sum of entries.
  entries[0] = 0;
  for (int i = 1; i <= np; i++)
    entries[i] = entries[i-1] + entries[i];
  // Now entries contains the first vector entry for each rank.

  // Create the global elements for the vector.
  int nMyGlobalElements = entries[me+1]-entries[me];
  std::vector<itype> myGlobalElements(nMyGlobalElements+1);

  for (int i = 0; i < nMyGlobalElements; i++)
    myGlobalElements[i] = entries[me] + i + offsetEpetra64;

  *vectorMap = new Epetra_Map(-1, nMyGlobalElements, &myGlobalElements[0],
                               0, comm);

  // Column map:  Easy; consecutive entries for all ranks in column.
  int firstRank = mypCol * npRows;  // First rank in my column
  nMyGlobalElements = 0;
  for (int i = firstRank; i < firstRank + npRows; i++)
    nMyGlobalElements += entries[i+1] - entries[i];

  itype myFirstCol = entries[firstRank];
  myGlobalElements.resize(nMyGlobalElements+1);
  for (int i = 0; i < nMyGlobalElements; i++)
    myGlobalElements[i] = myFirstCol + i + offsetEpetra64;

  *colMap = new Epetra_Map(-1, nMyGlobalElements, &myGlobalElements[0],
                            0, comm);


  // Row map:  trickier since corresponding vector entries are not
  //           consecutive
  firstRank = mypRow;  // First rank in my row
  nMyGlobalElements = 0;
  for (int i = 0; i < npCols; i++) {
    int rank = firstRank + i * npRows;
    nMyGlobalElements += entries[rank+1] - entries[rank];
  }
  myGlobalElements.resize(nMyGlobalElements+1);
  for (int cnt = 0, i = 0; i < npCols; i++) {
    int rank = firstRank + i * npRows;
    for (itype j = entries[rank]; j < entries[rank+1]; j++)
      myGlobalElements[cnt++] = j + offsetEpetra64;
  }

  *rowMap = new Epetra_Map(-1, nMyGlobalElements, &myGlobalElements[0],
                            0, comm);
}

///////////////////////////////////////////////////////////////////////////////
template <typename itype>
void build_maps(
  itype nrows,      // Number of global matrix rows
  bool testEpetra64,// Flag indicating whether to adjust global row/column
                    // indices to exercise Epetra64 capability.
  Epetra_Comm &comm,       // Epetra communicator to be used in maps
  Epetra_Map **vectorMap,  // OUTPUT: Map to be used for the vector
  Epetra_Map **rowMap,     // OUTPUT: Map to be used for the matrix rows
  Epetra_Map **colMap,     // OUTPUT: Map to be used for the matrix cols
  long long &offsetEpetra64, // OUTPUT for testing Epetra64: add offsetEpetra64
                             // to all row/column indices.
  bool verbose             // print out generated maps
)
{
  // Function to build the maps for 1D or 2D matrix distribution.
  // Output for 1D includes rowMap and NULL colMap and vectorMap.
  // Output for 2D includes rowMap, colMap and vectorMap.

  int me = comm.MyPID();
  int np = comm.NumProc();

  *rowMap = NULL;
  *colMap = NULL;
  *vectorMap = NULL;

//  offsetEpetra64 = (testEpetra64 ? (long long) INT_MAX - (long long) 5 : 0);
  offsetEpetra64 = (testEpetra64 ? (long long) 2 * INT_MAX : 0);

  // Generate 1D row-based decomposition.

  if ((me == 0) && verbose)
    std::cout << std::endl
         << "1D Distribution: " << std::endl
         << "    np     = " << np << std::endl;

  // Linear map similar to Trilinos default.
  itype nMyRows = nrows / np + (nrows % np > me);
  itype myFirstRow = me * (nrows / np) + MIN(nrows % np, me);
  itype *myGlobalRows = new itype[nMyRows];
  for (itype i = 0; i < nMyRows; i++)
    myGlobalRows[i] = i + myFirstRow + offsetEpetra64;
  *rowMap = new Epetra_Map(nrows, nMyRows, &myGlobalRows[0], 0, comm);
  delete [] myGlobalRows;
}

#endif
