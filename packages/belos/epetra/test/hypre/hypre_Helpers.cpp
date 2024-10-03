// @HEADER
// *****************************************************************************
//                 Didasko: Tutorial Package
//
// Copyright 2005 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "hypre_Helpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include <iostream>
#include <fstream>
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#include "Ifpack_Hypre.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <sstream>
#include <map>

using Teuchos::RCP;
using Teuchos::rcp;
EpetraExt_HypreIJMatrix* newHypreMatrix(const int N)
{
  HYPRE_IJMatrix Matrix;
  int ierr = 0;
  int i;
  int numprocs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int ilower = (N/numprocs)*rank;
  int iupper = (N/numprocs)*(rank+1);
  if(rank == numprocs-1){ /*Last processor */
    iupper = N;
  }

  // Create
  ierr += HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper-1, ilower, iupper-1, &Matrix);
  ierr += HYPRE_IJMatrixSetObjectType(Matrix,HYPRE_PARCSR);
  // Initialize
  ierr += HYPRE_IJMatrixInitialize(Matrix);

  if(Matrix == NULL){
    printf("Error! Matrix is NULL!\n");
    std::cin >> ierr;
  }

  srand(time(NULL));
  int rows[1];
  for(i = ilower; i < iupper; i++) {
    int ncols = (int)(1+( (double)rand()/(double)RAND_MAX ) * 0.5*(N-1));
    TEUCHOS_TEST_FOR_EXCEPTION(ncols <= 0, std::logic_error, "ncols is negative");
    Teuchos::Array<int> cols; cols.resize(ncols);
    Teuchos::Array<double> values; values.resize(ncols);
    for(int j = 0; j < ncols; j++){
      int index = 0;
      if(i-(ncols/2) >= 0 && i+(ncols/2) < N){
        index = j + i - (ncols/2);
      } else if (i-(ncols/2) < 0){
        index = j + i;
      } else{
        index = j + i - (ncols-1);
        }

      cols[j] = index;
      double currVal =  ( (double)rand()/(double)RAND_MAX ) * 100;
      values[j] = currVal;
    }
    rows[0] = i;
    ierr += HYPRE_IJMatrixAddToValues(Matrix, 1, &ncols, rows, &cols[0], &values[0]);
  }
  // Assemble
  ierr += HYPRE_IJMatrixAssemble(Matrix);
  EpetraExt_HypreIJMatrix* RetMat = new EpetraExt_HypreIJMatrix(Matrix);
  //Don't HYPRE_IJMatrixDestroy(Matrix);
  return RetMat;
}

Epetra_CrsMatrix* newCrsMatrix(int N){

  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  // pointer to the matrix to be created
  Epetra_CrsMatrix*      Matrix;
  //Epetra_CrsMatrix* Matrix;
  // container for parameters
  Teuchos::ParameterList GaleriList;
  int nx = N * Comm.NumProc();
  int ny = N * Comm.NumProc();
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);

  // amk Nov 24, 2015: This map described a non-contiguous distribution before.
  // hypre didn't like that at all, so I changed it
  Epetra_Map Map(nx*ny,0,Comm);

  Matrix   = Galeri::CreateCrsMatrix("Laplace2D", &Map, GaleriList);
  return Matrix;
}

Epetra_CrsMatrix* GetCrsMatrix(EpetraExt_HypreIJMatrix *Matrix)
{
  int N = Matrix->NumGlobalRows();
  Epetra_CrsMatrix* TestMat = new Epetra_CrsMatrix(Copy, Matrix->RowMatrixRowMap(), Matrix->RowMatrixColMap(), N, false);

  for(int i = 0; i < Matrix->NumMyRows(); i++){
    int entries;
    Matrix->NumMyRowEntries(i,entries);
    Teuchos::Array<double> Values; Values.resize(entries);
    Teuchos::Array<int> Indices; Indices.resize(entries);
    int NumEntries;
    Matrix->ExtractMyRowCopy(i,entries,NumEntries,&Values[0], &Indices[0]);
    for(int j = 0; j < NumEntries; j++){
      double currVal[1];
      currVal[0] = Values[j];
      int currInd[1];
      currInd[0] = Matrix->RowMatrixColMap().GID(Indices[j]);
      TestMat->InsertGlobalValues(Matrix->RowMatrixRowMap().GID(i), 1, currVal, currInd);
    }
  }
  TestMat->FillComplete();
  return TestMat;
}

bool EquivalentVectors(const Epetra_MultiVector &Y1, const Epetra_MultiVector &Y2, const double tol){

  bool retVal = true;

  int num_vectors = Y1.NumVectors();
  if(Y2.NumVectors() != num_vectors){
    printf("Multivectors do not have same number of vectors.\n");
    return false;
  }

  for(int j = 0; j < num_vectors; j++){
    if(Y1.MyLength() != Y2.MyLength()){
      printf("Vectors are not same size on local processor.\n");
      return false;
    }
    for(int i = 0; i < Y1.GlobalLength(); i++){
      int Y1_local = Y1.Map().LID(i);
      int Y2_local = Y2.Map().LID(i);
      if(Y1_local < 0 || Y2_local < 0){
        continue;
      }
      if(fabs((*Y1(j))[Y1_local] - (*Y2(j))[Y2_local]) > tol){
        printf("Vector number[%d] ", j);
        printf("Val1[%d] = %f != Val2[%d] = %f\n", i, (*Y1(j))[Y1_local], i, (*Y2(j))[Y2_local]);
        retVal = false;
      }
    }
  }
  Teuchos::Array<int> vals; vals.resize(Y1.Comm().NumProc());
  int my_vals[1]; my_vals[0] = (int)retVal;
  Y1.Comm().GatherAll(my_vals, &vals[0], 1);
  for(int i = 0; i < Y1.Comm().NumProc(); i++){
    if(vals[i] == false){
      retVal = false;
    }
  }
  if(retVal == false){
    printf("[%d]Failed vector equivalency test.\n", Y1.Comm().MyPID());
  }
  return retVal;
}



bool EquivalentMatrices(const Epetra_RowMatrix &HypreMatrix, const Epetra_RowMatrix &CrsMatrix, const double tol){
  bool retVal = true;
  int MyPID = HypreMatrix.Comm().MyPID();
  if(HypreMatrix.NumMyRows() != CrsMatrix.NumMyRows()){
    printf("Different number of local rows.");
    return false;
  }
  for(int j = 0; j < HypreMatrix.NumGlobalRows(); j++){
    int hyp_j = HypreMatrix.RowMatrixRowMap().LID(j);
    int crs_j = CrsMatrix.RowMatrixRowMap().LID(j);
    if(hyp_j < 0 || crs_j < 0){
      continue;
    }

    int NumEntries = HypreMatrix.NumMyCols();
    int entries1;
    int entries2;

    Teuchos::Array<double> Y1_vals; Y1_vals.resize(NumEntries);
    Teuchos::Array<double> Y2_vals; Y2_vals.resize(NumEntries);
    Teuchos::Array<int> indices1; indices1.resize(NumEntries);
    Teuchos::Array<int> indices2; indices2.resize(NumEntries);

    HypreMatrix.ExtractMyRowCopy(hyp_j,NumEntries, entries1, &Y1_vals[0], &indices1[0]);
    CrsMatrix.ExtractMyRowCopy(crs_j,NumEntries, entries2, &Y2_vals[0], &indices2[0]);
    if(entries1 != entries2){
      printf("Row[%d], Differing number of entries %d != %d\n", j, entries1, entries2);
    }
    for(int i = 0; i < entries1; i++){
      indices1[i] = HypreMatrix.RowMatrixColMap().GID(indices1[i]);
      indices2[i] = CrsMatrix.RowMatrixColMap().GID(indices2[i]);
    }
    for(int i = 1; i < entries1; ++i){
      int value = indices1[i];
      double my_val = Y1_vals[i];
      int jj = i-1;
      while(jj >= 0 && indices1[jj] > value){
        indices1[jj+1] = indices1[jj];
        Y1_vals[jj+1] = Y1_vals[jj];
        jj = jj-1;
      }
      indices1[jj+1] = value;
      Y1_vals[jj+1] = my_val;
    }
    for(int i = 1; i < entries2; ++i){
      int value = indices2[i];
      double my_val = Y2_vals[i];
      int jj = i-1;
      while(jj >= 0 && indices2[jj] > value){
        indices2[jj+1] = indices2[jj];
        Y2_vals[jj+1] = Y2_vals[jj];
        jj = jj-1;
      }
      indices2[jj+1] = value;
      Y2_vals[jj+1] = my_val;
    }

    for(int i = 0; i < entries1; i++){
      if(indices1[i] != indices2[i]){
        printf("indices[%d], %d != %d\n", i, indices1[i], indices2[i]);
        retVal = false;
      }
      if(fabs(Y1_vals[i] - Y2_vals[i]) > tol){
        printf("Failed at %d\n", i);
        retVal = false;
      }
    }
  }
  Teuchos::Array<int> vals; vals.resize(HypreMatrix.Comm().NumProc());
  int my_vals[1]; my_vals[0] = (int)retVal;
  HypreMatrix.Comm().GatherAll(my_vals, &vals[0], 1);
  for(int i = 0; i < HypreMatrix.Comm().NumProc(); i++){
    if(vals[i] == false){
      retVal = false;
    }
  }
  if(retVal == false){
    printf("[%d]Failed matrix equivalency test.\n", MyPID);
  }
  return retVal;
}
