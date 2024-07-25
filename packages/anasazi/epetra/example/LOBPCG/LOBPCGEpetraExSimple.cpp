// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/// \example LOBPCGEpetraExSimple.cpp
/// \brief Use "Simple LOBPCG" with Epetra test problem (computed here).
///
/// This example computes the eigenvalues of largest magnitude of an
/// eigenvalue problem $A x = \lambda x$, using Anasazi's "simple"
/// implementation of the LOBPCG method, with Epetra linear algebra.
/// It constructs the test problem within the example itself.

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Assert.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int
main (int argc, char *argv[])
{
  using namespace Anasazi;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init (&argc, &argv);
#endif // HAVE_MPI

  // Create an Epetra communicator
#ifdef HAVE_MPI
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif // HAVE_MPI

  // Create an Anasazi output manager
  BasicOutputManager<double> printer;
  printer.stream(Errors) << Anasazi_Version() << std::endl << std::endl;

  // Get the sorting std::string from the command line
  std::string which ("SM");
  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  if (cmdp.parse (argc, argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize ();
#endif // HAVE_MPI
    return -1;
  }

  // Dimension of the matrix
  //
  // Discretization points in any one direction.
  const int nx = 10;
  // Size of matrix nx*nx
  const int NumGlobalElements = nx*nx;

  // Construct a Map that puts approximately the same number of
  // equations on each process.
  Epetra_Map Map (NumGlobalElements, 0, Comm);

  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map.NumMyElements ();

  std::vector<int> MyGlobalElements (NumMyElements);
  Map.MyGlobalElements (&MyGlobalElements[0]);

  // Create an integer vector NumNz that is used to build the Petra
  // matrix.  NumNz[i] is the number of OFF-DIAGONAL terms for the ith
  // global equation on this process.
  std::vector<int> NumNz(NumMyElements);

  /* We are building a matrix of block structure:

      | T -I          |
      |-I  T -I       |
      |   -I  T       |
      |        ...  -I|
      |           -I T|

   where each block is dimension nx by nx and the matrix is on the order of
   nx*nx.  The block T is a tridiagonal matrix.
  */
  for (int i=0; i<NumMyElements; ++i) {
    if (MyGlobalElements[i] == 0 || MyGlobalElements[i] == NumGlobalElements-1 ||
        MyGlobalElements[i] == nx-1 || MyGlobalElements[i] == nx*(nx-1) ) {
      NumNz[i] = 3;
    }
    else if (MyGlobalElements[i] < nx || MyGlobalElements[i] > nx*(nx-1) ||
             MyGlobalElements[i]%nx == 0 || (MyGlobalElements[i]+1)%nx == 0) {
      NumNz[i] = 4;
    }
    else {
      NumNz[i] = 5;
    }
  }

  // Create an Epetra_Matrix
  RCP<Epetra_CrsMatrix> A = rcp (new Epetra_CrsMatrix (Epetra_DataAccess::Copy, Map, &NumNz[0]));

  // Compute coefficients for discrete convection-diffution operator
  const double one = 1.0;
  std::vector<double> Values(4);
  std::vector<int> Indices(4);
  double rho = 0.0;
  double h = one /(nx+1);
  double h2 = h*h;
  double c = 5.0e-01*rho/ h;
  Values[0] = -one/h2 - c; Values[1] = -one/h2 + c; Values[2] = -one/h2; Values[3]= -one/h2;
  double diag = 4.0 / h2;
  int NumEntries;

  for (int i=0; i<NumMyElements; ++i) {
    if (MyGlobalElements[i]==0) {
      Indices[0] = 1;
      Indices[1] = nx;
      NumEntries = 2;
      int info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if (MyGlobalElements[i] == nx*(nx-1)) {
      Indices[0] = nx*(nx-1)+1;
      Indices[1] = nx*(nx-2);
      NumEntries = 2;
      int info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if (MyGlobalElements[i] == nx-1) {
      Indices[0] = nx-2;
      NumEntries = 1;
      int info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
      Indices[0] = 2*nx-1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1) {
      Indices[0] = NumGlobalElements-2;
      NumEntries = 1;
      int info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
      Indices[0] = nx*(nx-1)-1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if (MyGlobalElements[i] < nx) {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]+nx;
      NumEntries = 3;
      int info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if (MyGlobalElements[i] > nx*(nx-1)) {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]-nx;
      NumEntries = 3;
      int info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if (MyGlobalElements[i]%nx == 0) {
      Indices[0] = MyGlobalElements[i]+1;
      Indices[1] = MyGlobalElements[i]-nx;
      Indices[2] = MyGlobalElements[i]+nx;
      NumEntries = 3;
      int info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if ((MyGlobalElements[i]+1)%nx == 0) {
      Indices[0] = MyGlobalElements[i]-nx;
      Indices[1] = MyGlobalElements[i]+nx;
      NumEntries = 2;
      int info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
      Indices[0] = MyGlobalElements[i]-1;
      NumEntries = 1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]-nx;
      Indices[3] = MyGlobalElements[i]+nx;
      NumEntries = 4;
      int info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    // Put in the diagonal entry
    int info = A->InsertGlobalValues(MyGlobalElements[i], 1, &diag, &MyGlobalElements[i]);
    TEUCHOS_ASSERT( info==0 );
  }

  // Finish up
  int info = A->FillComplete ();
  TEUCHOS_ASSERT( info==0 );
  A->SetTracebackMode (1); // Shutdown Epetra Warning tracebacks

  // Create a identity matrix for the temporary mass matrix
  RCP<Epetra_CrsMatrix> M = rcp (new Epetra_CrsMatrix (Epetra_DataAccess::Copy, Map, 1));
  for (int i=0; i<NumMyElements; i++) {
    Values[0] = one;
    Indices[0] = i;
    NumEntries = 1;
    info = M->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
    TEUCHOS_ASSERT( info==0 );
  }
  // Finish up
  info = M->FillComplete ();
  TEUCHOS_ASSERT( info==0 );
  M->SetTracebackMode (1); // Shutdown Epetra Warning tracebacks

  //************************************
  // Call the LOBPCG solver manager
  //***********************************
  //
  // Variables used for the LOBPCG Method
  const int nev       = 10;
  const int blockSize = 5;
  const int maxIters  = 500;
  const double tol    = 1.0e-8;
  const int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;
 
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef MultiVecTraits<double, Epetra_MultiVector> MVT;

  // Create an Epetra_MultiVector for an initial vector to start the
  // solver.  Note: This needs to have the same number of columns as
  // the blocksize.
  RCP<Epetra_MultiVector> ivec = rcp (new Epetra_MultiVector (Map, blockSize));
  ivec->Random (); // fill the initial vector with random values

  // Create the eigenproblem.
  RCP<BasicEigenproblem<double, MV, OP> > MyProblem =
    rcp (new BasicEigenproblem<double, MV, OP> (A, ivec));

  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->setHermitian (true);

  // Set the number of eigenvalues requested
  MyProblem->setNEV (nev);

  // Inform the eigenproblem that you are finishing passing it information
  const bool success = MyProblem->setProblem ();
  if (! success) {
    printer.print (Errors, "Anasazi::BasicEigenproblem::setProblem() reported an error.\n");
#ifdef HAVE_MPI
    MPI_Finalize ();
#endif
    return -1;
  }

  // Create parameter list to pass into the solver manager
  //
  Teuchos::ParameterList MyPL;
  MyPL.set ("Which", which);
  MyPL.set ("Block Size", blockSize);
  MyPL.set ("Maximum Iterations", maxIters);
  MyPL.set ("Convergence Tolerance", tol);
  MyPL.set ("Verbosity", verbosity);

  // Create the solver manager
  SimpleLOBPCGSolMgr<double, MV, OP> MySolverMan (MyProblem, MyPL);

  // Solve the problem
  ReturnType returnCode = MySolverMan.solve ();

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Eigensolution<double, MV> sol = MyProblem->getSolution ();
  std::vector<Value<double> > evals = sol.Evals;
  RCP<MV> evecs = sol.Evecs;

  // Compute residuals
  std::vector<double> normR (sol.numVecs);
  if (sol.numVecs > 0) {
    Teuchos::SerialDenseMatrix<int,double> T (sol.numVecs, sol.numVecs);
    Epetra_MultiVector tempAevec (Map, sol.numVecs );
    T.putScalar (0.0);
    for (int i = 0; i < sol.numVecs; ++i) {
      T(i,i) = evals[i].realpart;
    }
    A->Apply (*evecs, tempAevec);
    MVT::MvTimesMatAddMv (-1.0, *evecs, T, 1.0, tempAevec);
    MVT::MvNorm (tempAevec, normR);
  }

  // Print the results
  std::ostringstream os;
  os.setf (std::ios_base::right, std::ios_base::adjustfield);
  os << "Solver manager returned "
     << (returnCode == Converged ? "converged." : "unconverged.") << endl;
  os << endl;
  os << "------------------------------------------------------" << endl;
  os << std::setw(16) << "Eigenvalue"
     << std::setw(18) << "Direct Residual"
     << endl;
  os << "------------------------------------------------------" << endl;
  for (int i = 0; i < sol.numVecs; ++i) {
    os << std::setw(16) << evals[i].realpart
       << std::setw(18) << normR[i] / evals[i].realpart
       << endl;
  }
  os << "------------------------------------------------------" << endl;
  printer.print (Errors, os.str ());

#ifdef HAVE_MPI
  MPI_Finalize ();
#endif // HAVE_MPI
  return 0;
}
