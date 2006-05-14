#include "Amesos_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Amesos.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
using namespace Teuchos;
using namespace Galeri;

// ==================== //
// M A I N  D R I V E R //
// ==================== //
//
// This example will:
// 1.- create a linear system, stored as an
//     Epetra_LinearProblem. The matrix corresponds
//     to a 5pt Laplacian (2D on Cartesian grid).
//     The user can change the global size of the problem 
//     by modifying variable NumGlobalRows.
// 2.- The linear system matrix, solution and rhs
//     are distributed among the available processors,
//     using a linear distribution. This is for 
//     simplicity only! Amesos can support any Epetra_Map.
// 3.- Once the linear problem is created, we
//     create an Amesos Factory object.
// 4.- Using the Factory, we create the required Amesos_BaseSolver
//     solver. Any supported (and compiled) Amesos
//     solver can be used. If the selected solver
//     is not available (that is, if Amesos has *not*
//     been configured with support for this solver),
//     the factory returns 0. Usually, Amesos_Klu
//     is always available.
// 5.- At this point we can factorize the matrix,
//     and solve the linear system. Only three methods
//     should be used for any Amesos_BaseSolver object:
//     1.- NumericFactorization();
//     2.- SymbolicFactorization();
//     3.- Solve();
// 6.- We note that the header files of Amesos-supported
//     libraries are *not* required in this file. They are
//     actually needed to compile the Amesos library only.
//
// NOTE: this example can be run with any number of processors.
//
// Author: Marzio Sala, SNL 9214
// Last modified: Oct-05.

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int NumGlobalRows = 10000; // must be a square for the
                             // matrix generator.
  int NumVectors = 1;        // number of rhs's. Amesos
                             // supports single or
			     // multiple RHS.

  // Initializes an Gallery object.
  // NOTE: this example uses the Trilinos package Galeri
  // to define in an easy way the linear system matrix.
  // The user can easily change the matrix type; consult the 
  // Galeri documentation for mode details.
  //
  // Here the problem has size nx x ny, and the 2D Cartesian
  // grid is divided into mx x my subdomains.
  ParameterList GaleriList;
  GaleriList.set("nx", 100);
  GaleriList.set("ny", 100 * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* Matrix = CreateCrsMatrix("Laplace2D", Map, GaleriList);

  // Creates vectors for right-hand side and solution, and the
  // linear problem container.

  Epetra_Vector LHS(*Map); LHS.PutScalar(0.0); // zero solution
  Epetra_Vector RHS(*Map); RHS.Random();       // random rhs
  Epetra_LinearProblem Problem(Matrix, &LHS, &RHS);

  // ===================================================== //
  // B E G I N N I N G   O F  T H E   AM E S O S   P A R T //
  // ===================================================== //

  // Initializes the Amesos solver. This is the base class for
  // Amesos. It is a pure virtual class (hence objects of this
  // class cannot be allocated, and can exist only as pointers 
  // or references).
  //
  Amesos_BaseSolver* Solver;

  // Initializes the Factory. Factory is a function class (a
  // class that contains methods only, no data). Factory
  // will be used to create Amesos_BaseSolver derived objects.
  //
  Amesos Factory;

  Solver = Factory.Create("Klu", Problem);

  // Parameters for all Amesos solvers are set through
  // a call to SetParameters(List). List is a Teuchos
  // parameter list (Amesos requires Teuchos to compile).
  // In most cases, users can proceed without calling
  // SetParameters(). Please refer to the Amesos guide
  // for more details.
  // NOTE: you can skip this call; then the solver will
  // use default parameters.
  //
  // Parameters in the list are set using 
  // List.set("parameter-name", ParameterValue);
  // In this example, we specify that we want more output.
  //
  Teuchos::ParameterList List;
  List.set("PrintTiming", true);
  List.set("PrintStatus", true);
  
  Solver->SetParameters(List);
  
  // Now we are ready to solve. Generally, users will
  // call SymbolicFactorization(), then NumericFactorization(),
  // and finally Solve(). Note that:
  // - the numerical values of the linear system matrix
  //   are *not* required before NumericFactorization();
  // - solution and rhs are *not* required before calling
  //   Solve().
  if (Comm.MyPID() == 0)
    cout << "Starting symbolic factorization..." << endl;
  Solver->SymbolicFactorization();
  
  // you can change the matrix values here
  if (Comm.MyPID() == 0)
    cout << "Starting numeric factorization..." << endl;
  Solver->NumericFactorization();
  
  // you can change LHS and RHS here
  if (Comm.MyPID() == 0)
    cout << "Starting solution phase..." << endl;
  Solver->Solve();

  // =========================================== //
  // E N D   O F   T H E   A M E S O S   P A R T //
  // =========================================== //

  // delete Solver. MPI calls can occur.
  delete Solver;
    
  // delete the objects created by Galeri
  delete Matrix;
  delete Map;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);

} // end of main()
