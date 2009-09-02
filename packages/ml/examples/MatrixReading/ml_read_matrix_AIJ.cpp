
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_config.h"

// This example reads a SERIAL matrix in AIJ format, and creates
// the corresponding SERIAL Epetra_CrsMatrix. Then matrix is then
// distributed over all the available processes, to obtain a
// DISTRIBUTED matrix that can be solved with all the available
// processes. 
//
// The matrix must be stored in an ASCII file (specified by the first
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
// A small example of matrix.aij is contained in the
// ml/examples/ExampleMatrices // subdirectory. For example, run the
// code as follows:
// $ mpirun -np 2 ./ml_read_matrix_AIJ.exe ~/Trilinos/packages/ml/example/ExampleMatrices/simple_matrix.aij
// This matrix is very small and it is included only to explain how
// to write the .aij file.
//

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_AZTECOO)

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

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int NumRows;
  int NumElements;
  int Offset;
  std::ifstream data_file;

  if (Comm.MyPID() == 0) 
  {
    // proc 0 reads the number of rows, columns, nonzero elements
    // The matrix is supposed to be square (otherwise ML doesn't work)
    
    char *FileName = argv[1];
    if (FileName == 0)
    {
      cerr << "Usage: <executable name> <matrix name>" << endl;
#ifdef HAVE_MPI
      MPI_Finalize();
#endif
      exit(EXIT_SUCCESS);
    }

    data_file.open(FileName);

    if (!data_file.good()) {
      std::cerr << "Error opening file `" << FileName << "'" << endl;
      exit(EXIT_FAILURE);
    }

    data_file >> NumRows;
    data_file >> NumElements;
    data_file >> Offset;

    cout << "Matrix name                = " << FileName << endl;
    cout << "Number of rows             = " << NumRows << endl;
    cout << "Number of nonzero elements = " << NumElements << endl;
    cout << "Offset                     = " << Offset << endl;
  }
  else
    NumRows = 0;
  
  // creates a map with all elements on proc 0
  Epetra_Map* SerialMap = new Epetra_Map(-1,NumRows,0,Comm);
  Epetra_CrsMatrix* SerialMatrix = new Epetra_CrsMatrix(Copy,*SerialMap,0);

  if (Comm.MyPID() == 0) 
  {
    // now proc 0 read the actual matrix, element by element
    for (int i = 0 ; i < NumElements ; ++i) 
    {
      int row;
      int col;
      double val;
      data_file >> row;
      data_file >> col;

      if (row < Offset || col < Offset || row >= NumRows + Offset || col >= NumRows + Offset)
      {
        cout << "Something wrong at element " << i << endl;
        cout << "row = " << row;
        cout << ", col = " << col << ", while NumRows = " << NumRows << endl;
        exit(EXIT_FAILURE);
      }
             
      data_file >> val;
      row -= Offset;
      col -= Offset;
      int ierr;
      ierr = SerialMatrix->InsertGlobalValues(row,1,&val,&col);
      if (ierr < 0)
      {
        cout << "Error at element " << i << endl;
        ML_CHK_ERR(ierr);
      }
#if 0
      // If only one half of a symmetric matrix is stored, then this
      // part of the code will insert the other half as well.
      if (row != col)
      {
        ierr = SerialMatrix->InsertGlobalValues(col,1,&val,&row);
        if (ierr < 0)
          ML_CHK_ERR(ierr);
      }
#endif
    }
#if 0
    // The matrix can still be modified here, for example this is to
    // add or insert a diagonal value
    for (int i = 0 ; i < NumRows ; ++i) 
    {
      double value = 3.99476e+16;
      if (SerialMatrix->SumIntoGlobalValues(i, 1, &value, &i))
        SerialMatrix->InsertGlobalValues(i, 1, &value, &i);
    }
#endif
  }

  SerialMatrix->FillComplete();

  Epetra_Map* DistributedMap = 0;
  Epetra_CrsMatrix* DistributedMatrix = 0;

  // Distributes the matrix but only if necessary
  
  if (Comm.NumProc() > 1)
  {
    // need to create the distributed map, this
    // is for simplicity linear
    Comm.Broadcast(&NumRows,1,0);
    DistributedMap = new Epetra_Map(NumRows, 0, Comm);

    DistributedMatrix = new Epetra_CrsMatrix(Copy, *DistributedMap,0);

    // creates the import 
    Epetra_Import Importer(*DistributedMap,*SerialMap);

    ML_CHK_ERR(DistributedMatrix->Import(*SerialMatrix, Importer, Insert));

    ML_CHK_ERR(DistributedMatrix->FillComplete());

    // can delete serial objects, no longer needed
    delete SerialMap;
    delete SerialMatrix;
  }
  else
  {
    DistributedMap = SerialMap;
    DistributedMatrix = SerialMatrix;
  }
  DistributedMatrix->OptimizeStorage();

  // =========================== begin of ML part ===========================

  // create a parameter list for ML options
  ParameterList MLList;

  ML_Epetra::SetDefaults("SA", MLList);          // Use NSSA for highly nonsymmetric
  MLList.set("smoother: type", "Chebyshev");
  MLList.set("smoother: sweeps", 3);
  MLList.set("eigen-analysis: type", "cg");     // use power-method and 15
  MLList.set("eigen-analysis: iterations", 10); // iterations for nonsymmetric
                                                // systems.
  MLList.set("aggregation: threshold", 0.0);

  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*DistributedMatrix, MLList);

  // =========================== end of ML part =============================

  Epetra_Vector LHS(*DistributedMap);       // solution vector
  Epetra_Vector LHSexact(*DistributedMap);  // exact solution, check later
  Epetra_Vector RHS(*DistributedMap);       // right-hand side
  LHS.PutScalar(0.0);                       // zero starting solution
  LHSexact.Random();                        // random exact solution
  DistributedMatrix->Multiply(false,LHSexact,RHS);
  
  Epetra_LinearProblem Problem(DistributedMatrix,&LHS,&RHS);
  AztecOO solver(Problem);

  solver.SetAztecOption(AZ_solver, AZ_cg);  // Change (e.g. AZ_gmres) for nonsymmetric
  solver.SetAztecOption(AZ_output, 32);
  solver.SetPrecOperator(MLPrec);

  // solve with 500 iterations and 1e-12 as tolerance on the
  // relative residual  
  solver.Iterate(500, 1e-8);

  // delete the preconditioner. Do it BEFORE calling MPI_Finalize
  delete MLPrec;

  // check the error
  LHSexact.Update(1.0,LHS,-1.0);
  double Norm;
  LHSexact.Norm2(&Norm);

  delete DistributedMatrix;
  delete DistributedMap;

#ifdef HAVE_MPI
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

  puts("Please configure ML with:");
  puts("--enable-epetra");
  puts("--enable-teuchos");
  puts("--enable-aztecoo");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(EXIT_SUCCESS);
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) */
