//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2009) Sandia Corporation
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

#include "hypre_Helpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include <iostream>
#include <fstream>

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <map>

EpetraExt_HypreIJMatrix::EpetraExt_HypreIJMatrix* newHypreMatrix(const int N, const int type)
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
  
  if(type == 0){
    // Set values
    int rows[1];
    int cols[1];
    double values[1];
    int ncols = 1;
    for(i = ilower; i < iupper; i++) {
      rows[0] = i;
      cols[0] = i;
      values[0] = 1.0;
      ierr += HYPRE_IJMatrixSetValues(Matrix, 1, &ncols, rows, cols, values);
    }
  } else if(type == 1){
    srand(time(NULL));
    // Set values
    int rows[1];
    int cols[1];
    double values[1];
    int ncols = 1;
    for(i = ilower; i < iupper; i++) {
      rows[0] = i;
      cols[0] = i;
      values[0] =  ( (double)rand()/(double)RAND_MAX ) * 100;
      ierr += HYPRE_IJMatrixSetValues(Matrix, 1, &ncols, rows, cols, values);
    }
    
  } else if(type == 2){
    // Set values
    int rows[1];
    Teuchos::Array<int> cols; cols.resize(N);
    Teuchos::Array<double> values; values.resize(N);
    int ncols = N;
    for(i = ilower; i < iupper; i++) {
      for(int j = 0; j < N; j++){
        cols[j] = j;
        values[j] = j;
      }
      rows[0] = i;
      ierr += HYPRE_IJMatrixSetValues(Matrix, 1, &ncols, rows, &cols[0], &values[0]);
    }
  } else if(type == 3){
    srand(time(NULL));
    int rows[1];
    Teuchos::Array<int> cols; cols.resize(N);
    Teuchos::Array<double> values; values.resize(N);
    int ncols = N;
    for(i = ilower; i < iupper; i++) {
      for(int j = 0; j < N; j++){
        cols[j] = j;
        double currVal =  ( (double)rand()/(double)RAND_MAX ) * 100;
        values[j] = currVal;
      }
      rows[0] = i;
      ierr += HYPRE_IJMatrixSetValues(Matrix, 1, &ncols, rows, &cols[0], &values[0]);
    }
  } else {
    srand(time(NULL));
    int rows[1];
    for(i = ilower; i < iupper; i++) {
      int ncols = (int)(1+( (double)rand()/(double)RAND_MAX ) * 0.5*(N-1));
      TEST_FOR_EXCEPTION(ncols <= 0, std::logic_error, "ncols is negative");
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
      ierr += HYPRE_IJMatrixSetValues(Matrix, 1, &ncols, rows, &cols[0], &values[0]);
    }
  }
  // Assemble
  ierr += HYPRE_IJMatrixAssemble(Matrix);
  EpetraExt_HypreIJMatrix* RetMat = new EpetraExt_HypreIJMatrix(Matrix);
  //Don't HYPRE_IJMatrixDestroy(Matrix);
  return RetMat;
}

Epetra_CrsMatrix::Epetra_CrsMatrix* newCrsMatrix(EpetraExt_HypreIJMatrix &Matrix)
{
  int N = Matrix.NumGlobalRows();
  Epetra_CrsMatrix* TestMat = new Epetra_CrsMatrix(Copy, Matrix.RowMatrixRowMap(), Matrix.RowMatrixColMap(), N, false);
  
  for(int i = 0; i < Matrix.NumMyRows(); i++){
    int entries;
    Matrix.NumMyRowEntries(i,entries);
    Teuchos::Array<double> Values; Values.resize(entries);
    Teuchos::Array<int> Indices; Indices.resize(entries);
    int NumEntries;
    Matrix.ExtractMyRowCopy(i,entries,NumEntries,&Values[0], &Indices[0]);
    for(int j = 0; j < NumEntries; j++){
      double currVal[1];
      currVal[0] = Values[j];
      int currInd[1];
      currInd[0] = Matrix.RowMatrixColMap().GID(Indices[j]);
      TestMat->InsertGlobalValues(Matrix.RowMatrixRowMap().GID(i), 1, currVal, currInd);
    }
  }
  TestMat->FillComplete();
  return TestMat;
}

bool EquivalentVectors(Epetra_MultiVector &Y1, Epetra_MultiVector &Y2, const double tol){
  
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
    Teuchos::Array<double> Y1_vals; Y1_vals.resize(Y1.MyLength());
    Teuchos::Array<double> Y2_vals; Y2_vals.resize(Y2.MyLength());
    (*Y1(j)).ExtractCopy(&Y1_vals[0]);
    (*Y2(j)).ExtractCopy(&Y2_vals[0]);
    
    for(int i = 0; i < Y1.MyLength(); i++){
      if(fabs(Y1_vals[i] - Y2_vals[i]) > tol){
        printf("Vector number[%d] ", j);
        printf("Val1[%d] = %f != Val2[%d] = %f\n", i, Y1_vals[i], i, Y2_vals[i]);
        retVal = false;
      }
    }
  }
  if(retVal == false){
    printf("Failed equivalent vectors.\n");
  }
  return retVal;
}



bool EquivalentMatrices(Epetra_RowMatrix &HypreMatrix, Epetra_RowMatrix &CrsMatrix, const double tol){
  bool retVal = true;
  for(int j = 0; j < HypreMatrix.NumMyRows(); j++){
    
    int NumEntries;
    int entries1;
    int entries2;
    HypreMatrix.NumMyRowEntries(j,NumEntries);

    Teuchos::Array<double> Y1_vals; Y1_vals.resize(NumEntries);
    Teuchos::Array<double> Y2_vals; Y2_vals.resize(NumEntries);
    Teuchos::Array<int> indices1; indices1.resize(NumEntries);
    Teuchos::Array<int> indices2; indices2.resize(NumEntries);
     
    HypreMatrix.ExtractMyRowCopy(j,NumEntries, entries1, &Y1_vals[0], &indices1[0]);
    CrsMatrix.ExtractMyRowCopy(j,NumEntries, entries2, &Y2_vals[0], &indices2[0]);

    std::map<int,double> Y1map;
    std::map<int,double> Y2map;
    for (int i=0; i < NumEntries ; ++i) {
      Y1map[indices1[i]] = Y1_vals[i]; 
      Y2map[indices2[i]] = Y2_vals[i]; 
    }
    retVal = retVal && (Y1map == Y2map);
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
    printf("[%d]Failed matrix equivalency test.\n", HypreMatrix.Comm().MyPID());
    if(HypreMatrix.Comm().MyPID() == 0){
      //std::ofstream outfile("Matrices.txt");
      
      //HypreMatrix.Print(outfile);
      //CrsMatrix.Print(outfile);
      //outfile.close();
    }
  }
  return retVal;
}
