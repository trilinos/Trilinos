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
#include "EpetraExt_CrsMatrixIn.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"

using namespace EpetraExt;
namespace EpetraExt {

int MatrixMarketFileToCrsMatrix( const char *filename, const Epetra_Map & rowMap, Epetra_CrsMatrix * & A) {

  A = new Epetra_CrsMatrix(Copy, rowMap, 0);
  return(MatrixMarketFileToCrsMatrixHandle(filename, A));
}

int MatrixMarketFileToCrsMatrix( const char *filename, const Epetra_Map & rowMap, const Epetra_Map & colMap, Epetra_CrsMatrix * & A) {

  A = new Epetra_CrsMatrix(Copy, rowMap, colMap, 0);
  return(MatrixMarketFileToCrsMatrixHandle(filename, A));
}

int MatrixMarketFileToCrsMatrixHandle( const char *filename, Epetra_CrsMatrix * A) {

  const int lineLength = 1025;
  const int tokenLength = 35;
  char line[lineLength];
  char token[tokenLength];
  char token1[tokenLength];
  char token2[tokenLength];
  char token3[tokenLength];
  char token4[tokenLength];
  char token5[tokenLength];
  int M, N, NZ;

  FILE * handle = 0;

  handle = fopen(filename,"r");  // Open file

  // Check first line, which should be "%%MatrixMarket matrix coordinate real general" (without quotes)
  if(fgets(line, lineLength, handle)==0) return(-1);
  if(sscanf(line, "%s", token1, token2, token3, token4, token5 )==0) return(-1);
  if (strcmp(token1, "%%MatrixMarket") ||
      strcmp(token2, "matrix") ||
      strcmp(token3, "coordinate") ||
      strcmp(token4, "real") ||
      strcmp(token5, "general")) return(-1);


  bool inHeader = true;

  // Next, strip off header lines (which start with "%")
  while (inHeader) {
    if(fgets(line, lineLength, handle)==0) return(-1);
    if(sscanf(line, "%c", token)==0) return(-1);
    if (!strcmp(token, "%")) inHeader = false;
  }

  // Next get problem dimensions: M, N, NZ
  if(fgets(line, lineLength, handle)==0) return(-1); // numProc value
  if(sscanf(line, "%d %d %d", &M, &N, &NZ)==0) return(-1);

  // Now read in each triplet and store to the local portion of the matrix if the row is owned.
  int I, J;
  double V;
  const Epetra_Map & map = A->RowMap();
  for (int i=0; i<NZ; i++) {
    if(fgets(line, lineLength, handle)==0) return(-1); // MaxElementSize header line
    if(sscanf(line, "%d %d %f", &I, &J, &V)==0) return(-1);
    if (map.MyGID(I))
      A->InsertGlobalValues(I, 1, &V, &J);
  }

  A->FillComplete();

  return(0);
}
} // namespace EpetraExt
