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
#include "EpetraExt_CrsMatrixIn.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"
#include "Epetra_Util.h"

using namespace EpetraExt;
namespace EpetraExt {

int MatrixMarketFileToCrsMatrix( const char *filename, const Epetra_Comm & comm, Epetra_CrsMatrix * & A) {
    
  EPETRA_CHK_ERR(MatrixMarketFileToCrsMatrixHandle(filename, comm, A));
  return(0);
}
  
int MatrixMarketFileToCrsMatrix(const char *filename, const Epetra_Map & rowMap,
				  const Epetra_Map& rangeMap, const Epetra_Map& domainMap, Epetra_CrsMatrix * & A) {
  EPETRA_CHK_ERR(MatrixMarketFileToCrsMatrixHandle(filename, rowMap.Comm(), A, &rowMap, 0, &rangeMap, &domainMap));
  return(0);
}

int MatrixMarketFileToCrsMatrix( const char *filename, const Epetra_Map & rowMap, Epetra_CrsMatrix * & A) {

  EPETRA_CHK_ERR(MatrixMarketFileToCrsMatrixHandle(filename, rowMap.Comm(), A, &rowMap));
  return(0);
}

int MatrixMarketFileToCrsMatrix( const char *filename, const Epetra_Map & rowMap, const Epetra_Map & colMap, Epetra_CrsMatrix * & A) {

  EPETRA_CHK_ERR(MatrixMarketFileToCrsMatrixHandle(filename, rowMap.Comm(), A, &rowMap, &colMap));
  return(0);
}

int MatrixMarketFileToCrsMatrix(const char *filename, const Epetra_Map & rowMap, const Epetra_Map & colMap,
                                const Epetra_Map& rangeMap, const Epetra_Map& domainMap, Epetra_CrsMatrix * & A) {
  EPETRA_CHK_ERR(MatrixMarketFileToCrsMatrixHandle(filename, rowMap.Comm(), A, &rowMap, &colMap, &rangeMap, &domainMap));
  return(0);
}

int MatrixMarketFileToCrsMatrixHandle(const char *filename,
				      const Epetra_Comm & comm,
                                      Epetra_CrsMatrix * & A,
				      const Epetra_Map * rowMap,
				      const Epetra_Map * colMap,
                                      const Epetra_Map * rangeMap,
                                      const Epetra_Map * domainMap)
{

  const int lineLength = 1025;
  const int tokenLength = 35;
  char line[lineLength];
  char token1[tokenLength];
  char token2[tokenLength];
  char token3[tokenLength];
  char token4[tokenLength];
  char token5[tokenLength];
  int M, N, NZ;

  // Make sure domain and range maps are either both defined or undefined 
  if ((domainMap!=0 && rangeMap==0) || (domainMap==0 && rangeMap!=0)) {EPETRA_CHK_ERR(-3);}

  // check maps to see if domain and range are 1-to-1

  if (domainMap!=0) {
    if (!domainMap->UniqueGIDs()) {EPETRA_CHK_ERR(-2);}
    if (!rangeMap->UniqueGIDs()) {EPETRA_CHK_ERR(-2);}
  }
  else {
    // If domain and range are not specified, rowMap becomes both and must be 1-to-1 if specified
    if (rowMap!=0) {
      if (!rowMap->UniqueGIDs()) {EPETRA_CHK_ERR(-2);}
    }
  }
      

  FILE * handle = 0;

  handle = fopen(filename,"r");  // Open file
  if (handle == 0)
    EPETRA_CHK_ERR(-1); // file not found

  // Check first line, which should be "%%MatrixMarket matrix coordinate real general" (without quotes)
  if(fgets(line, lineLength, handle)==0) {if (handle!=0) fclose(handle); EPETRA_CHK_ERR(-1);}
  if(sscanf(line, "%s %s %s %s %s", token1, token2, token3, token4, token5 )==0) {if (handle!=0) fclose(handle); EPETRA_CHK_ERR(-1);}
  if (strcmp(token1, "%%MatrixMarket") ||
      strcmp(token2, "matrix") ||
      strcmp(token3, "coordinate") ||
      strcmp(token4, "real") ||
      strcmp(token5, "general")) {if (handle!=0) fclose(handle); EPETRA_CHK_ERR(-1);}

  // Next, strip off header lines (which start with "%")
  do {
    if(fgets(line, lineLength, handle)==0) {if (handle!=0) fclose(handle); EPETRA_CHK_ERR(-1);}
  } while (line[0] == '%');

  // Next get problem dimensions: M, N, NZ
  if(sscanf(line, "%d %d %d", &M, &N, &NZ)==0) {if (handle!=0) fclose(handle); EPETRA_CHK_ERR(-1);}

  // Now create matrix using user maps if provided.

  if (rowMap!=0 && colMap !=0)
    A = new Epetra_CrsMatrix(Copy, *rowMap, *colMap, 0);
  else if (rowMap!=0)
    A = new Epetra_CrsMatrix(Copy, *rowMap, 0);
  else {
    Epetra_Map newRowMap(M,0, comm);
    A = new Epetra_CrsMatrix(Copy, newRowMap, 0);
  }
  // Now read in each triplet and store to the local portion of the matrix if the row is owned.
  int I, J;
  double V;
  const Epetra_Map & rowMap1 = A->RowMap();
  const Epetra_Map & colMap1 = A->ColMap();
  int ioffset = rowMap1.IndexBase()-1;
  int joffset = colMap1.IndexBase()-1;
  
  for (int i=0; i<NZ; i++) {
    if(fgets(line, lineLength, handle)==0) {if (handle!=0) fclose(handle); EPETRA_CHK_ERR(-1);};
    if(sscanf(line, "%d %d %lg\n", &I, &J, &V)==0) {if (handle!=0) fclose(handle); EPETRA_CHK_ERR(-1);}
    I+=ioffset; J+=joffset; // Convert to Zero based
    if (rowMap1.MyGID(I)) {
      int ierr = A->InsertGlobalValues(I, 1, &V, &J);
      if (ierr<0) EPETRA_CHK_ERR(ierr);
    }
  }
    
  if (rangeMap != 0 && domainMap != 0) {
    A->FillComplete(*domainMap, *rangeMap);
  }
  else {
    A->FillComplete();
  }
  
  if (handle!=0) fclose(handle);
  return(0);
}

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
} // namespace EpetraExt

