
//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "ml_config.h"

// This example reads a matrix in AIJ format, and creates
// the corresponding Epetra_CrsMatrix. Then matrix is then
// distributed over all the available processes.
//
// The matrix is stored in an ASCII file (specified by the first
// argument of the command line), which contains the following lines:
// ----(file begins line below)---
// <NumRows>
// <NumElements>
// <Offset>
// i j A_{i,j}
// ...
// ---(file ends line above)---
// In the file, NumRows is the number of rows, <NumCols> the number of
// columns, and <NumElements> the number of nonzero elements of the matrix.
// <Offset> is the offset of the first row (1 for MATLAB and FORTRAN matrices,
// 0 for C/C++ matrices). Elements can be stored in any order.
//
// An example of use can be as:
// $ mpirun -np 2 ./ml_read_matrix_AIJ.exe <matrix-name>
// (a small example of matrix.aij is contained in the
// ml/examples/ExampleMatrices // subdirectory).
//
// NOTE: this code is not efficient for serial runs,
// as it creates the distributed matrix anyway (which in that
// case is simply a copy of the serial matrix).

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "AztecOO.h"

// includes required by ML
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra.h"
#include <fstream>

using namespace Teuchos;
using namespace Trilinos_Util;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int NumRows;
  int Offset;
  int NumElements;
  std::ifstream data_file;

  if (Comm.MyPID() == 0) {

    // proc 0 reads the number of rows, columns, nonzero elements and offset
    char *FileName = argv[1];
    assert (FileName != 0);
    string Title;

    data_file.open(FileName);

    if (!data_file.good()) {
      std::cerr << "Error opening file `" << FileName << "'" << endl;
      exit(EXIT_FAILURE);
    }

    data_file >> NumRows;
    data_file >> NumElements;
    data_file >> Offset;

    cout << "Number of rows             = " << NumRows << endl;
    cout << "Number of nonzero elements = " << NumElements << endl;
    cout << "Offset                     = " << Offset << endl;
  }
  else
    NumRows = 0;
  
  // creates a map with all elements on proc 0
  Epetra_Map* SerialMap = new Epetra_Map(-1,NumRows,0,Comm);
  Epetra_CrsMatrix* SerialMatrix;
  SerialMatrix = new Epetra_CrsMatrix(Copy,*SerialMap,0);

  if (Comm.MyPID() == 0) {

    // now proc 0 read the actual matrix, element by element
    for (int i = 0 ; i < NumElements ; ++i) {
      int row;
      int col;
      double val;
      data_file >> row;
      data_file >> col;
      data_file >> val;
      row -= Offset;
      col -= Offset;
      SerialMatrix->InsertGlobalValues(row,1,&val,&col);
    }

  }
  SerialMatrix->FillComplete();

  data_file.close();

  // need to create the distributed map, this
  // is for simplicity linear
  Comm.Broadcast(&NumRows,1,0);
  Epetra_Map DistributedMap(NumRows, 0, Comm);

  Epetra_CrsMatrix DistributedMatrix(Copy, DistributedMap,0);

  // creates the import 
  Epetra_Import Importer(DistributedMap,*SerialMap);

  ML_CHK_ERR(DistributedMatrix.Import(*SerialMatrix,
				      Importer, Insert));
  
  ML_CHK_ERR(DistributedMatrix.FillComplete());
  
  // can delete serial objects, no longer needed
  delete SerialMap;
  delete SerialMatrix;

  // =========================== begin of ML part ===========================

  // create a parameter list for ML options
  ParameterList MLList;

  ML_Epetra::SetDefaults("SA",MLList);
 
  MLList.set("smoother: type", "symmetric Gauss-Seidel");

  vector<double> x_coord(NumRows);
  vector<double> y_coord(NumRows);
  
  char *CoordFileName = argv[2];
  assert (CoordFileName != 0);

  data_file.open(CoordFileName);

  if (!data_file.good()) {
	  std::cerr << "Error opening file `" << CoordFileName << "'" << endl;
	  exit(EXIT_FAILURE);
  }

  for (int i = 0 ; i < NumRows ; ++i) {
	  int itmp;
	  double dtmp;
	  data_file >> itmp;
	  assert (itmp == i);
	  data_file >> x_coord[i];
	  data_file >> y_coord[i];
	  data_file >> dtmp;
  }

  data_file.close();

  for (int i = 0 ; i < NumRows ; ++i) 
    cout << i << ": " << x_coord[i] << ", " << y_coord[i] << endl;

  vector<double> NullSpace;
  
  char *NullSpaceFileName = argv[3];
  assert (NullSpaceFileName != 0);

  data_file.open(NullSpaceFileName);

  if (!data_file.good()) {
	  std::cerr << "Error opening file `" << NullSpaceFileName 
                    << "'" << endl;
	  exit(EXIT_FAILURE);
  }

  int size;
  int NullSpaceSize;
  data_file >> size;
  data_file >> NullSpaceSize;
  NullSpace.resize(NullSpaceSize * NumRows);

  for (int i = 0 ; i < NullSpaceSize * NumRows ; ++i) {
	  data_file >> NullSpace[i];
  }

  data_file.close();

  for (int i = 0 ; i < NumRows ; ++i) 
    for (int j = 0 ; j < NullSpaceSize ; ++j) 
	    cout << NullSpace[i + j * NumRows] << endl;

  // set the read null space
  MLList.set("null space: type", "pre-computed");
  MLList.set("null space: dimension", NullSpaceSize);
  MLList.set("null space: vectors", &NullSpace[0]);

  // number of relaxation sweeps
  MLList.set("adaptive: max sweeps", 10);
  // number of additional null space vectors to compute
  MLList.set("adaptive: num vectors",2);

  MLList.set("adaptive: visualize", true);
  MLList.set("viz: x-coordinates", &x_coord[0]);
  MLList.set("viz: y-coordinates", &y_coord[0]);
#if 1
  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(DistributedMatrix, MLList, false);
  MLPrec->ComputeAdaptivePreconditioner(NullSpaceSize,&NullSpace[0]);
#else
  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(DistributedMatrix, MLList);
#endif

  // =========================== end of ML part =============================

  Epetra_Vector LHS(DistributedMap);       // solution vector
  Epetra_Vector LHSexact(DistributedMap);  // exact solution, check later
  Epetra_Vector RHS(DistributedMap);       // right-hand side
  LHS.PutScalar(0.0);                      // zero starting solution
  LHSexact.Random();                       // random exact solution
  DistributedMatrix.Multiply(false,LHSexact,RHS);
  
  Epetra_LinearProblem Problem(&DistributedMatrix,&LHS,&RHS);
  AztecOO solver(Problem);

  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  solver.SetPrecOperator(MLPrec);

  // solve with 500 iterations and 1e-12 as tolerance on the
  // relative residual  
  solver.Iterate(500, 1e-12);

  // delete the preconditioner. Do it BEFORE MPI_Finalize
  delete MLPrec;

  // check the error
  LHSexact.Update(1.0,LHS,-1.0);
  double Norm;
  LHSexact.Norm2(&Norm);

  if (Norm > 1e-3) {
    cerr << "TEST FAILED" << endl;
    exit(EXIT_FAILURE);
  }
  
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
  
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with --enable-epetra --enable-teuchos");
  puts("--enable-aztecoo --enable-triutils");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */
