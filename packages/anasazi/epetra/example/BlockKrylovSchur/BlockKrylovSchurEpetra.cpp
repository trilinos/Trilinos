// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
// This example computes the eigenvalues of largest magnitude of an
// eigenvalue problem $A x = \lambda x$, using Anasazi's
// implementation of the Block Krylov-Schur method.

// Include header for block Davidson eigensolver
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
// Include header to define eigenproblem Ax = \lambda*x
#include "AnasaziBasicEigenproblem.hpp"
// Include header to provide Anasazi with Epetra adapters.  If you
// plan to use Tpetra objects instead of Epetra objects, include
// AnasaziTpetraAdapter.hpp instead; do analogously if you plan to use
// Thyra objects instead of Epetra objects.
#include "AnasaziEpetraAdapter.hpp"
// Include header for Epetra sparse matrix and multivector.
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
// The Trilinos package Galeri has many example problems.
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
// Include selected communicator class required by Epetra objects
#ifdef EPETRA_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif // EPETRA_MPI

// ****************************************************************************
// BEGIN MAIN ROUTINE
// ****************************************************************************

int
main (int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cerr;
  using std::cout;
  using std::endl;
  // Anasazi solvers have the following template parameters:
  //
  //   - Scalar: The type of dot product results.
  //   - MV: The type of (multi)vectors.
  //   - OP: The type of operators (functions from multivector to
  //     multivector).  A matrix (like Epetra_CrsMatrix) is an example
  //     of an operator; an Ifpack preconditioner is another example.
  //
  // Here, Scalar is double, MV is Epetra_MultiVector, and OP is
  // Epetra_Operator.
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;

#ifdef EPETRA_MPI
  MPI_Init (&argc, &argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif // EPETRA_MPI

  const int MyPID = Comm.MyPID ();

  //
  // Set up the test problem.
  //
  // We use Trilinos' Galeri package to construct a test problem.
  // Here, we use a discretization of the 2-D Laplacian operator.
  // The global mesh size is nx * nx.
  //
  Teuchos::ParameterList GaleriList;
  const int nx = 30;
  GaleriList.set ("n", nx * nx);
  GaleriList.set ("nx", nx);
  GaleriList.set ("ny", nx);
  RCP<Epetra_Map> Map = rcp (Galeri::CreateMap ("Linear", Comm, GaleriList));
  RCP<Epetra_RowMatrix> A =
    rcp (Galeri::CreateCrsMatrix ("Laplace2D", &*Map, GaleriList));

  // Set eigensolver parameters.
  const double tol = 1.0e-8; // convergence tolerance
  const int nev = 4; // number of eigenvalues for which to solve
  const int blockSize = 5; // block size (number of eigenvectors processed at once)
  const int numBlocks = 8; // restart length
  const int maxRestarts = 100; // maximum number of restart cycles

  // Create a set of initial vectors to start the eigensolver.
  // This needs to have the same number of columns as the block size.
  RCP<MV> ivec = rcp (new MV (*Map, blockSize));
  ivec->Random ();

  // Create the eigenproblem.  This object holds all the stuff about
  // your problem that Anasazi will see.  In this case, it knows about
  // the matrix A and the inital vectors.
  RCP<Anasazi::BasicEigenproblem<double, MV, OP> > problem =
    rcp (new Anasazi::BasicEigenproblem<double, MV, OP> (A, ivec));

  // Tell the eigenproblem that the operator A is symmetric.
  problem->setHermitian (true);

  // Set the number of eigenvalues requested
  problem->setNEV (nev);

  // Tell the eigenproblem that you are finishing passing it information.
  const bool boolret = problem->setProblem();
  if (boolret != true) {
    if (MyPID == 0) {
      cerr << "Anasazi::BasicEigenproblem::setProblem() returned an error." << endl;
    }
#ifdef EPETRA_MPI
    MPI_Finalize ();
#endif // EPETRA_MPI
    return -1;
  }

  // Create a ParameterList, to pass parameters into the Block
  // Davidson eigensolver.
  Teuchos::ParameterList anasaziPL;
  anasaziPL.set ("Which", "LM");
  anasaziPL.set ("Block Size", blockSize);
  anasaziPL.set ("Num Blocks", numBlocks);
  anasaziPL.set ("Maximum Restarts", maxRestarts);
  anasaziPL.set ("Convergence Tolerance", tol);
  anasaziPL.set ("Verbosity", Anasazi::Errors + Anasazi::Warnings +
                 Anasazi::TimingDetails + Anasazi::FinalSummary);

  // Create the Block Davidson eigensolver.
  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> anasaziSolver (problem, anasaziPL);

  // Solve the eigenvalue problem.
  //
  // Note that creating the eigensolver is separate from solving it.
  // After creating the eigensolver, you may call solve() multiple
  // times with different parameters or initial vectors.  This lets
  // you reuse intermediate state, like allocated basis vectors.
  Anasazi::ReturnType returnCode = anasaziSolver.solve ();
  if (returnCode != Anasazi::Converged && MyPID == 0) {
    cout << "Anasazi eigensolver did not converge." << endl;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem.
  Anasazi::Eigensolution<double,MV> sol = problem->getSolution ();
  // Anasazi returns eigenvalues as Anasazi::Value, so that if
  // Anasazi's Scalar type is real-valued (as it is in this case), but
  // some eigenvalues are complex, you can still access the
  // eigenvalues correctly.  In this case, there are no complex
  // eigenvalues, since the matrix pencil is symmetric.
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  RCP<MV> evecs = sol.Evecs;

  // Compute residuals.
  std::vector<double> normR (sol.numVecs);
  if (sol.numVecs > 0) {
    Teuchos::SerialDenseMatrix<int,double> T (sol.numVecs, sol.numVecs);
    MV tempAevec (*Map, sol.numVecs);
    T.putScalar (0.0);
    for (int i=0; i<sol.numVecs; ++i) {
      T(i,i) = evals[i].realpart;
    }
    A->Apply (*evecs, tempAevec);
    MVT::MvTimesMatAddMv (-1.0, *evecs, T, 1.0, tempAevec);
    MVT::MvNorm (tempAevec, normR);
  }

  // Print the results on MPI process 0.
  if (MyPID == 0) {
    cout << "Solver manager returned "
         << (returnCode == Anasazi::Converged ? "converged." : "unconverged.")
         << endl << endl
         << "------------------------------------------------------" << endl
         << std::setw(16) << "Eigenvalue"
         << std::setw(18) << "Direct Residual"
         << endl
         << "------------------------------------------------------" << endl;
    for (int i=0; i<sol.numVecs; ++i) {
      cout << std::setw(16) << evals[i].realpart
           << std::setw(18) << normR[i] / evals[i].realpart
           << endl;
    }
    cout << "------------------------------------------------------" << endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize () ;
#endif // EPETRA_MPI

  return 0;
}
