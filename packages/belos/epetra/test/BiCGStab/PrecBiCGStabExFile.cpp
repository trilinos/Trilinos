// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver reads a problem from a file, which can be in Harwell-Boeing (*.hb),
// Matrix Market (*.mtx), or triplet format (*.triU, *.triS).  The right-hand side
// from the problem, if it exists, will be used instead of multiple random
// right-hand-sides.  The initial guesses are all set to zero.  An ILU preconditioner
// is constructed using the Ifpack factory.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBiCGStabSolMgr.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"
#include "Epetra_Map.h"
#ifdef EPETRA_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Ifpack.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main (int argc, char *argv[])
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cout;
  using std::endl;
  typedef double                            ST;
  typedef Teuchos::ScalarTraits<ST>        SCT;
  typedef SCT::magnitudeType                MT;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  typedef Belos::MultiVecTraits<ST,MV>     MVT;
  typedef Belos::OperatorTraits<ST,MV,OP>  OPT;

  int MyPID = 0;
#ifdef EPETRA_MPI
  MPI_Init (&argc, &argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
  MyPID = Comm.MyPID ();
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = false;
  bool success = true;

  // This "try" relates to TEUCHOS_STANDARD_CATCH_STATEMENTS near the
  // bottom of main().  That macro has the corresponding "catch".
  try {
    bool proc_verbose = false;
    bool leftprec = false;     // left preconditioning or right.
    int frequency = -1;        // frequency of status test output.
    int numrhs = 1;            // number of right-hand sides to solve for
    int maxiters = -1;         // maximum number of iterations allowed per linear system
    std::string filename ("cage4.hb");
    MT tol = 1.0e-5;           // relative residual tolerance

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption ("verbose", "quiet", &verbose, "Whether to print messages "
                    "and results.");
    cmdp.setOption ("left-prec", "right-prec", &leftprec, "Left or right "
                    "preconditioning.");
    cmdp.setOption ("frequency", &frequency, "Frequency of solver output "
                    "(1 means every iteration; -1 means never).");
    cmdp.setOption ("filename", &filename, "Test matrix filename.  "
                    "Allowed file extensions: *.hb, *.mtx, *.triU, *.triS");
    cmdp.setOption ("tol", &tol, "Relative residual tolerance for solver.");
    cmdp.setOption ("num-rhs", &numrhs, "Number of right-hand sides to solve.");
    cmdp.setOption ("max-iters", &maxiters, "Maximum number of iterations per "
                    "linear system (-1 = adapted to problem/block size).");
    if (cmdp.parse (argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (! verbose) {
      frequency = -1;  // reset frequency if test is not verbose
    }

    //
    // *************Get the problem*********************
    //
    RCP<Epetra_Map> Map;
    RCP<Epetra_CrsMatrix> A;
    RCP<Epetra_MultiVector> B, X;
    RCP<Epetra_Vector> vecB, vecX;
    EpetraExt::readEpetraLinearSystem (filename, Comm, &A, &Map, &vecX, &vecB);
    A->OptimizeStorage ();
    proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */

    // Check to see if the number of right-hand sides is the same as requested.
    if (numrhs > 1) {
      X = rcp (new Epetra_MultiVector (*Map, numrhs));
      B = rcp (new Epetra_MultiVector (*Map, numrhs));
      X->Random ();
      OPT::Apply (*A, *X, *B);
      X->PutScalar (0.0);
    }
    else {
      X = Teuchos::rcp_implicit_cast<Epetra_MultiVector> (vecX);
      B = Teuchos::rcp_implicit_cast<Epetra_MultiVector> (vecB);
    }

    // FIXME (mfh 05 Aug 2015) Doesn't this negate what we just did
    // above with numrhs > 1 ?
    B->PutScalar (1.0);

    //
    // Construct preconditioner
    //

    // allocates an IFPACK factory. No data is associated
    // to this object (only method Create()).
    Ifpack Factory;

    // create the preconditioner. For valid PrecType values,
    // please check the documentation
    std::string PrecType = "point relaxation stand-alone"; // incomplete Cholesky
    int OverlapLevel = 0; // must be >= 0. If Comm.NumProc() == 1,
    // it is ignored.

    RCP<Ifpack_Preconditioner> Prec;
    try {
      Prec = rcp (Factory.Create (PrecType, &*A, OverlapLevel));
    } catch (std::exception& e) {
      std::ostringstream os;
      os << "Ifpack preconditioner construction threw an exception: "
         << e.what ();
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, os.str ());
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (Prec.is_null (), std::runtime_error, "Failed to create Ifpack "
       "preconditioner!");

    try {
      ParameterList ifpackList;
      ifpackList.set ("relaxation: type", "Jacobi");
      IFPACK_CHK_ERR( Prec->SetParameters (ifpackList) );
      IFPACK_CHK_ERR( Prec->Initialize () );
      IFPACK_CHK_ERR( Prec->Compute () );
    } catch (std::exception& e) {
      std::ostringstream os;
      os << "Ifpack preconditioner setup threw an exception: " << e.what ();
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, os.str ());
    }

    // Create the Belos preconditioned operator from the Ifpack
    // preconditioner.  This is necessary because Belos expects an
    // operator to apply the preconditioner with Apply(), NOT with
    // ApplyInverse().
    RCP<Belos::EpetraPrecOp> belosPrec = rcp (new Belos::EpetraPrecOp (Prec));

    //
    // Create parameter list for the Belos solver
    //
    const int NumGlobalElements = B->GlobalLength ();
    if (maxiters == -1) {
      maxiters = NumGlobalElements - 1; // maximum number of iterations to run
    }
    RCP<ParameterList> belosList (new ParameterList ("Belos"));
    belosList->set ("Maximum Iterations", maxiters);
    belosList->set ("Convergence Tolerance", tol);
    if (verbose) {
      belosList->set ("Verbosity", Belos::Errors + Belos::Warnings +
                     Belos::TimingDetails + Belos::StatusTestDetails);
      if (frequency > 0) {
        belosList->set ("Output Frequency", frequency);
      }
    }
    else {
      belosList->set ("Verbosity", Belos::Errors + Belos::Warnings +
                      Belos::FinalSummary);
    }
    //
    // Construct a preconditioned linear problem
    //
    RCP<Belos::LinearProblem<double,MV,OP> > problem
      = rcp (new Belos::LinearProblem<double,MV,OP> (A, X, B));
    if (leftprec) {
      problem->setLeftPrec (belosPrec);
    }
    else {
      problem->setRightPrec (belosPrec);
    }
    bool set = problem->setProblem ();
    if (! set) {
      if (proc_verbose) {
        cout << endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << endl;
      }
      return -1;
    }

    // Create a Belos solver.
    RCP<Belos::SolverManager<double,MV,OP> > solver
      = rcp (new Belos::BiCGStabSolMgr<double,MV,OP> (problem, belosList));

    if (proc_verbose) {
      cout << endl << endl;
      cout << "Dimension of matrix: " << NumGlobalElements << endl;
      cout << "Number of right-hand sides: " << numrhs << endl;
      cout << "Max number of CG iterations: " << maxiters << endl;
      cout << "Relative residual tolerance: " << tol << endl;
      cout << endl;
    }

    // Ask Belos to solve the linear system.
    Belos::ReturnType ret = solver->solve();

    //
    // Compute actual residuals.
    //
    bool badRes = false;
    std::vector<double> actual_resids (numrhs);
    std::vector<double> rhs_norm (numrhs);
    Epetra_MultiVector resid (*Map, numrhs);
    OPT::Apply (*A, *X, resid);
    MVT::MvAddMv (-1.0, resid, 1.0, *B, resid);
    MVT::MvNorm (resid, actual_resids);
    MVT::MvNorm (*B, rhs_norm);
    if (proc_verbose) {
      cout << "---------- Actual Residuals (normalized) ----------" << endl << endl;
      for (int i = 0; i < numrhs; ++i) {
        double actRes = actual_resids[i] / rhs_norm[i];
        cout <<"Problem " << i << " : \t" << actRes << endl;
        if (actRes > tol) {
          badRes = true;
        }
      }
    }

    if (ret != Belos::Converged || badRes) {
      success = false;
      if (proc_verbose) {
        cout << endl << "ERROR:  Belos did not converge!" << endl
             << "End Result: TEST FAILED" << endl;
      }
    } else {
      success = true;
      if (proc_verbose) {
        cout << endl << "SUCCESS:  Belos converged!" << endl
             << "End Result: TEST PASSED" << endl;
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
