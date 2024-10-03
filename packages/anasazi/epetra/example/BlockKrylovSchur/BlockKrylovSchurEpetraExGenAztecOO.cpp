// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
// This example computes the eigenvalues of smallest magnitude of a
// generalized eigenvalue problem $K x = \lambda M x$, using Anasazi's
// implementation of the block Krylov-Schur method.  It implements
// inverse iteration using an AztecOO iterative linear solver with an
// Ifpack preconditioner.
//
// Anasazi computes the smallest-magnitude eigenvalues using a
// shift-and-invert strategy.  For simplicity, the code below uses a
// shift of zero.  It illustrates the general pattern for using
// Anasazi for this problem:
//
//   1. Construct an "operator" A such that $Az = K^{-1} M z$.
//   2. Use Anasazi to solve $Az = \sigma z$, which is a spectral
//      transformation of the original problem $K x = \lambda M x$.
//   3. The eigenvalues $\lambda$ of the original problem are the
//      inverses of the eigenvalues $\sigma$ of the transformed
//      problem.
//
// The "operator A such that $A z = K^{-1} M z$" is a subclass of
// Epetra_Operator.  The Apply method of that operator takes the
// vector b, and computes $x = K^{-1} M b$.  It does so first by
// applying the matrix M, and then by solving the linear system $K x =
// M b$ for x.  Trilinos implements many different ways to solve
// linear systems.  The example below uses an iterative linear solver
// provided by the AztecOO package, with an Ifpack preconditioner, to
// solve linear systems.

// Include header for block Krylov-Schur solver
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
// Include header to define basic eigenproblem Ax = \lambda*Bx
#include "AnasaziBasicEigenproblem.hpp"
// Include header to provide Anasazi with Epetra adapters.  If you
// plan to use Tpetra objects instead of Epetra objects, include
// AnasaziTpetraAdapter.hpp instead; do analogously if you plan to use
// Thyra objects instead of Epetra objects.
#include "AnasaziEpetraAdapter.hpp"
// Include header for Epetra sparse matrix, Map (representation of
// parallel distributions), and linear problem.  AztecOO uses the
// latter to encapsulate linear systems to solve.
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
// Include header for AztecOO iterative linear solver, and
// AztecOO_Operator.  The latter wraps an AztecOO solver in an
// Epetra_Operator.
#include "AztecOO.h"
#include "AztecOO_Operator.h"
// Include header for Ifpack preconditioner factory
#include "Ifpack.h"
// Include header for Teuchos serial dense matrix, which we use to
// compute the eigenvectors.
#include "Teuchos_SerialDenseMatrix.hpp"
// Include header for the problem definition
#include "ModeLaplace2DQ2.h"

// Include selected communicator class and map required by Epetra objects
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
  typedef Anasazi::MultiVecTraits<double, MV> MVT;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init (&argc, &argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif // EPETRA_MPI

  const int MyPID = Comm.MyPID ();

  //
  // Set up the test problem
  //

  // Dimensionality of the spatial domain to discretize
  const int space_dim = 2;

  // Size of each of the dimensions of the (discrete) domain
  std::vector<double> brick_dim (space_dim);
  brick_dim[0] = 1.0;
  brick_dim[1] = 1.0;

  // Number of elements in each of the dimensions of the domain
  std::vector<int> elements (space_dim);
  elements[0] = 10;
  elements[1] = 10;

  // Create the test problem
  RCP<ModalProblem> testCase =
    rcp (new ModeLaplace2DQ2 (Comm, brick_dim[0], elements[0],
                              brick_dim[1], elements[1]));

  // Get the stiffness and mass matrices.
  //
  // rcp (T*, false) returns a nonowning (doesn't deallocate) RCP.
  RCP<Epetra_CrsMatrix> K =
    rcp (const_cast<Epetra_CrsMatrix* > (testCase->getStiffness ()), false);
  RCP<Epetra_CrsMatrix> M =
    rcp (const_cast<Epetra_CrsMatrix*> (testCase->getMass ()), false);

  //
  // Create linear solver for linear systems with K
  //
  // Anasazi uses shift and invert, with a "shift" of zero, to find
  // the eigenvalues of least magnitude.  In this example, we
  // implement the "invert" part of shift and invert by using an
  // AztecOO iterative linear solver with an Ifpack preconditioner.
  //

  //
  // Construct Ifpack preconditioner
  //

  // An Ifpack "factory" knows how to create Ifpack preconditioners.
  Ifpack Factory;

  // Use the Factory to create the preconditioner.  Factory.Create()
  // returns a raw Ifpack2_Preconditioner pointer.  The caller (that's
  // us!) is responsible for deallocation, so we use an owning RCP to
  // deallocate it automatically.
  //
  // The Factory needs a string name of the preconditioner, and the
  // overlap level.  (Almost) all Ifpack preconditioners are local to
  // each MPI process.  Ifpack automatically does an additive Schwarz
  // decomposition across processes.  The default overlap level for
  // additive Schwarz is zero, but you can use a different overlap
  // level here if you want.  The overlap level must be nonnegative,
  // and only matters if running with more than one MPI process.
  std::string PrecType = "ICT"; // incomplete Cholesky
  const int OverlapLevel = 0;
  // Create the preconditioner.
  RCP<Ifpack_Preconditioner> Prec =
    rcp (Factory.Create (PrecType, &*K, OverlapLevel));
  if (Prec.is_null ()) {
    throw std::logic_error ("Ifpack's factory returned a NULL preconditioner!");
  }

  //
  // Set Ifpack preconditioner parameters.
  //
  // Set parameters after creating the preconditioner.  Please refer
  // to Ifpack's documentation for a list of valid parameters.
  //
  Teuchos::ParameterList ifpackList;
  ifpackList.set ("fact: drop tolerance", 1e-4);
  ifpackList.set ("fact: ict level-of-fill", 0.0);
  // The "combine mode" describes how to combine contributions from
  // other MPI processes.  It only matters if the overlap level is
  // nonzero.  See Epetra_CombineMode.h for documentation of of the
  // various combine modes.
  ifpackList.set ("schwarz: combine mode", "Add");

  // Set the parameters.
  IFPACK_CHK_ERR(Prec->SetParameters(ifpackList));

  // Initialize the preconditioner.  This only looks at the structure
  // of the matrix, not the values.  Nevertheless, the matrix must
  // generally be fill complete at this point.  Compare Initialize()
  // to a symbolic factorization, and Compute() (see below) to a
  // numeric factorization.
  IFPACK_CHK_ERR(Prec->Initialize());

  // Compute the preconditioner.  This looks at both the structure and
  // the values of the matrix.  Compare Initialize() (see above) to a
  // symbolic factorization, and Compute() (see below) to a numeric
  // factorization.
  IFPACK_CHK_ERR(Prec->Compute());

  //
  // Set up AztecOO iterative solver for solving linear systems with K.
  //

  // Create Epetra linear problem class for solving linear systems
  // with K.  This implements the inverse operator for shift and
  // invert.
  Epetra_LinearProblem precProblem;
  // Tell the linear problem about the matrix K.  Epetra_LinearProblem
  // doesn't know about RCP, so we have to give it a raw pointer.
  precProblem.SetOperator (K.getRawPtr ());

  // Create AztecOO iterative solver for solving linear systems with K.
  AztecOO precSolver (precProblem);
  // Tell the solver to use the Ifpack preconditioner we created above.
  precSolver.SetPrecOperator (Prec.get ());
  // Set AztecOO solver options.
  precSolver.SetAztecOption (AZ_output, AZ_none); // Don't print output
  precSolver.SetAztecOption (AZ_solver, AZ_cg); // Use CG

  // Use the above AztecOO solver to create the AztecOO_Operator.
  // This is the place where we tell the AztecOO solver the maximum
  // number of iterations (here, we use the matrix dimension; in
  // practice, you'll want a smaller number) and the convergence
  // tolerance (here, 1e-12).
  RCP<AztecOO_Operator> precOperator =
    rcp (new AztecOO_Operator (&precSolver, K->NumGlobalRows (), 1e-12));

  // Create an Operator that computes y = K^{-1} M x.
  //
  // This Operator object is the operator we give to Anasazi.  Thus,
  // Anasazi just sees an operator that computes y = K^{-1} M x.  The
  // matrix K got absorbed into precOperator via precProblem (the
  // Epetra_LinearProblem object).  Later, when we set up the Anasazi
  // eigensolver, we will need to tell it about M, so that it can
  // orthogonalize basis vectors with respect to the inner product
  // defined by M (since it is symmetric positive definite).
  RCP<Anasazi::EpetraGenOp> Aop = rcp (new Anasazi::EpetraGenOp (precOperator, M));

  //
  // Set parameters for the block Krylov-Schur eigensolver
  //

  double tol = 1.0e-8;
  int nev = 10;
  int blockSize = 3;
  int numBlocks = 3 * nev / blockSize;
  int maxRestarts = 10;

  // We're looking for the largest-magnitude eigenvalues of the
  // _inverse_ operator, thus, the smallest-magnitude eigenvalues of
  // the original operator.
  std::string which = "LM";
  int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;

  // Create ParameterList to pass into eigensolver
  Teuchos::ParameterList MyPL;
  MyPL.set ("Verbosity", verbosity);
  MyPL.set ("Which", which);
  MyPL.set ("Block Size", blockSize);
  MyPL.set ("Num Blocks", numBlocks);
  MyPL.set ("Maximum Restarts", maxRestarts);
  MyPL.set ("Convergence Tolerance", tol);

  // Create an initial set of vectors to start the eigensolver.  Note:
  // This needs to have the same number of columns as the block size.
  RCP<MV> ivec = rcp (new MV (K->Map (), blockSize));
  MVT::MvRandom (*ivec);

  // This object holds all the stuff about your problem that Anasazi
  // will see.
  //
  // Anasazi only needs M so that it can orthogonalize basis vectors
  // with respect to the M inner product.  Wouldn't it be nice if
  // Anasazi didn't require M in two different places?  Alas, this is
  // not currently the case.
  RCP<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem =
    rcp (new Anasazi::BasicEigenproblem<double,MV,OP> (Aop, M, ivec));

  // Tell the eigenproblem that the matrix pencil (K,M) is symmetric.
  MyProblem->setHermitian (true);

  // Set the number of eigenvalues requested
  MyProblem->setNEV (nev);

  // Tell the eigenproblem that you are finished passing it information.
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (MyPID == 0) {
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
    }
#ifdef EPETRA_MPI
    MPI_Finalize ();
#endif // EPETRA_MPI
    return -1;
  }

  // Create the Block Krylov-Schur eigensolver.
  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMgr (MyProblem, MyPL);

  // Solve the eigenvalue problem.
  //
  // Note that creating the eigensolver is separate from solving it.
  // After creating the eigensolver, you may call solve() multiple
  // times with different parameters or initial vectors.  This lets
  // you reuse intermediate state, like allocated basis vectors.
  Anasazi::ReturnType returnCode = MySolverMgr.solve ();
  if (returnCode != Anasazi::Converged && MyPID == 0) {
    cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem.
  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution ();
  // Anasazi returns eigenvalues as Anasazi::Value, so that if
  // Anasazi's Scalar type is real-valued (as it is in this case), but
  // some eigenvalues are complex, you can still access the
  // eigenvalues correctly.  In this case, there are no complex
  // eigenvalues, since the matrix pencil is symmetric.
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  if (numev > 0) {
    // Reconstruct the eigenvalues.  The ones that Anasazi gave back
    // are the inverses of the original eigenvalues.  Reconstruct the
    // eigenvectors too.
    Teuchos::SerialDenseMatrix<int,double> dmatr(numev,numev);
    MV tempvec (K->Map (), MVT::GetNumberVecs (*evecs));
    K->Apply (*evecs, tempvec);
    MVT::MvTransMv (1.0, tempvec, *evecs, dmatr);

    if (MyPID == 0) {
      double compeval = 0.0;
      cout.setf (std::ios_base::right, std::ios_base::adjustfield);
      cout << "Actual Eigenvalues (obtained by Rayleigh quotient) : " << endl;
      cout << "------------------------------------------------------" << endl;
      cout << std::setw(16) << "Real Part"
           << std::setw(16) << "Rayleigh Error" << endl;
      cout << "------------------------------------------------------" << endl;
      for (int i = 0; i < numev; ++i) {
        compeval = dmatr(i,i);
        cout << std::setw(16) << compeval
             << std::setw(16)
             << std::fabs (compeval - 1.0/evals[i].realpart)
             << endl;
      }
      cout << "------------------------------------------------------" << endl;
    }

  }

#ifdef EPETRA_MPI
  MPI_Finalize ();
#endif // EPETRA_MPI

  return 0;
}
