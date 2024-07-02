// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Purpose
// The example tests the successive right-hand sides capabilities of ML
// and Belos on a heat flow u_t = u_xx problem.
//
// A sequence of linear systems with the same coefficient matrix and
// different right-hand sides is solved.  A seed space is generated dynamically,
// and a deflated linear system is solved.  After each solves, the first
// few Krylov vectors are saved, and used to reduce the number of iterations
// for later solves.
// The optimal numbers of vectors to deflate and save are not known.
// Presently, the maximum number of vectors to deflate (seed space dimension)
// and to save are user paraemters.
// The seed space dimension is less than or equal to total number of vectors
// saved. The difference between the seed space dimension and the total number
// of vectors, is the number of vectors used to update the seed space after each
// solve. I guess that a seed space whose dimension is a small fraction of the
// total space will be best.
//
// maxSave=1 and maxDeflate=0 uses no recycling (not tested ).
//
// TODO: Instrument with timers, so that we can tell what is going on besides
//       by counting the numbers of iterations.
//
//
// Adapted from PCPGEpetraExFile.cpp by David M. Day (with original comments)

// All preconditioning has been commented out

// Tpetra
#include <TpetraExt_MatrixMatrix.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix_fwd.hpp>
#include <Tpetra_Map_fwd.hpp>
#include <Tpetra_Vector_fwd.hpp>

// MueLu
// #include <MueLu_TpetraOperator.hpp>
// #include <MueLu_CreateTpetraPreconditioner.hpp>

// Teuchos
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

// Belos
#include "BelosLinearProblem.hpp"
#include "BelosPCPGSolMgr.hpp"
#include "BelosTpetraAdapter.hpp"

template <typename ScalarType>
int run(int argc, char *argv[]) {
  //
  // Laplace's equation, homogeneous Dirichlet boundary conditions, [0,1]^2
  // regular mesh, Q1 finite elements
  //

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;

  using ST = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO = typename Tpetra::Vector<>::local_ordinal_type;
  using GO = typename Tpetra::Vector<>::global_ordinal_type;
  using NT = typename Tpetra::Vector<>::node_type;

  using V = typename Tpetra::Vector<ST, LO, GO, NT>;
  using MV = typename Tpetra::MultiVector<ST, LO, GO, NT>;
  using OP = typename Tpetra::Operator<ST, LO, GO, NT>;
  using MVT = typename Belos::MultiVecTraits<ST, MV>;
  using OPT = typename Belos::OperatorTraits<ST, MV, OP>;
  using MAP = typename Tpetra::Map<LO, GO, NT>;
  using MAT = typename Tpetra::CrsMatrix<ST, LO, GO, NT>;
  using SCT = typename Teuchos::ScalarTraits<ST>;
  using MT = typename SCT::magnitudeType;

  using LinearProblem = typename Belos::LinearProblem<ST, MV, OP>;
  using PCPGSolMgr = ::Belos::PCPGSolMgr<ST, MV, OP>;

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

  bool verbose = false;
  bool success = true;

  try {
    bool procVerbose = false;
    int frequency = -1;  // frequency of status test output.
    int blocksize = 1;   // blocksize, PCPGIter
    int numRhs = 1;      // number of right-hand sides to solve for
    int maxIters = 30;   // maximum number of iterations allowed per linear system

    int maxDeflate = 4;  // maximum number of vectors deflated from the linear system;
    // There is no overhead cost assoc with changing maxDeflate between solves
    int maxSave = 8;  // maximum number of vectors saved from current and previous .");
    // If maxSave changes between solves, then re-initialize (setSize).

    // Hypothesis: seed vectors are conjugate.
    // Initial versions allowed users to supply a seed space et cetera, but no
    // longer.

    // The documentation it suitable for certain tasks, like defining a modules
    // grammar,
    std::string ortho("ICGS");  // The Belos documentation obscures the fact that
    // IMGS is Iterated Modified Gram Schmidt,
    // ICGS is Iterated Classical Gram Schmidt, and
    // DKGS is another Iterated Classical Gram Schmidt.
    // Mathematical issues, such as the difference between ICGS and DKGS, are
    // not documented at all. UH tells me that Anasazi::SVQBOrthoManager is
    // available;  I need it for Belos
    MT tol = sqrt(std::numeric_limits<ST>::epsilon());  // relative residual tolerance

    // How do command line parsers work?
    Teuchos::CommandLineProcessor cmdp(false, true);
    cmdp.setOption("verbose", "quiet", &verbose, "Print messages and results");
    cmdp.setOption("frequency", &frequency, "Solvers frequency for printing residuals (#iters)");
    cmdp.setOption("tol", &tol, "Relative residual tolerance used by PCPG solver");
    cmdp.setOption("num-rhs", &numRhs, "Number of right-hand sides to be solved for");
    cmdp.setOption("max-iters", &maxIters,
                   "Maximum number of iterations per linear system (-1 = "
                   "adapted to problem/block size)");
    cmdp.setOption("num-deflate", &maxDeflate, "Number of vectors deflated from the linear system");
    cmdp.setOption("num-save", &maxSave, "Number of vectors saved from old Krylov subspaces");
    cmdp.setOption("ortho-type", &ortho, "Orthogonalization type, either DGKS, ICGS or IMGS");
    if (cmdp.parse(argc, argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose

    //
    // *************Form the problem****************
    //
    int numTimeStep = 4;
    GO numElePerDirection = 14 * comm->getSize();  // 5 -> 20
    size_t numNodes = (numElePerDirection - 1) * (numElePerDirection - 1);

    // By the way, either matrix has (3*numElePerDirection - 2)^2 nonzeros.
    RCP<MAP> map = rcp(new MAP(numNodes, 0, comm));
    RCP<MAT> stiff = rcp(new MAT(map, numNodes));
    RCP<MAT> mass = rcp(new MAT(map, numNodes));
    RCP<V> vecLHS = rcp(new V(map));
    RCP<V> vecRHS = rcp(new V(map));
    RCP<MV> LHS, RHS;

    ST ko = 8.0 / 3.0, k1 = -1.0 / 3.0;

    ST h = 1.0 / static_cast<ST>(numElePerDirection);  // x=(iX,iY)h
    ST mo = h * h * 4.0 / 9.0, m1 = h * h / 9.0, m2 = h * h / 36.0;

    ST pi = 4.0 * atan(1.0), valueLHS;
    GO node, iX, iY;

    for (LO lid = map->getMinLocalIndex(); lid <= map->getMaxLocalIndex(); lid++) {
      node = map->getGlobalElement(lid);
      iX = node % (numElePerDirection - 1);
      iY = (node - iX) / (numElePerDirection - 1);
      GO pos = node;
      stiff->insertGlobalValues(node, tuple(pos), tuple(ko));
      mass->insertGlobalValues(node, tuple(pos),
                               tuple(mo));  // init guess violates hom Dir bc
      valueLHS = sin(pi * h * ((ST)iX + 1)) * cos(2.0 * pi * h * ((ST)iY + 1));
      vecLHS->replaceGlobalValue(1, valueLHS);

      if (iY > 0) {
        pos = iX + (iY - 1) * (numElePerDirection - 1);
        stiff->insertGlobalValues(node, tuple(pos), tuple(k1));  // North
        mass->insertGlobalValues(node, tuple(pos), tuple(m1));
      }
      if (iY < numElePerDirection - 2) {
        pos = iX + (iY + 1) * (numElePerDirection - 1);
        stiff->insertGlobalValues(node, tuple(pos), tuple(k1));  // South
        mass->insertGlobalValues(node, tuple(pos), tuple(m1));
      }

      if (iX > 0) {
        pos = iX - 1 + iY * (numElePerDirection - 1);
        stiff->insertGlobalValues(node, tuple(pos), tuple(k1));  // West
        mass->insertGlobalValues(node, tuple(pos), tuple(m1));
        if (iY > 0) {
          pos = iX - 1 + (iY - 1) * (numElePerDirection - 1);
          stiff->insertGlobalValues(node, tuple(pos), tuple(k1));  // North West
          mass->insertGlobalValues(node, tuple(pos), tuple(m2));
        }
        if (iY < numElePerDirection - 2) {
          pos = iX - 1 + (iY + 1) * (numElePerDirection - 1);
          stiff->insertGlobalValues(node, tuple(pos), tuple(k1));  // South West
          mass->insertGlobalValues(node, tuple(pos), tuple(m2));
        }
      }

      if (iX < numElePerDirection - 2) {
        pos = iX + 1 + iY * (numElePerDirection - 1);
        stiff->insertGlobalValues(node, tuple(pos), tuple(k1));  // East
        mass->insertGlobalValues(node, tuple(pos), tuple(m1));
        if (iY > 0) {
          pos = iX + 1 + (iY - 1) * (numElePerDirection - 1);
          stiff->insertGlobalValues(node, tuple(pos), tuple(k1));  // North East
          mass->insertGlobalValues(node, tuple(pos), tuple(m2));
        }
        if (iY < numElePerDirection - 2) {
          pos = iX + 1 + (iY + 1) * (numElePerDirection - 1);
          stiff->insertGlobalValues(node, tuple(pos), tuple(k1));  // South East
          mass->insertGlobalValues(node, tuple(pos), tuple(m2));
        }
      }
    }

    stiff->fillComplete();
    mass->fillComplete();

    const ST one = SCT::one();
    ST hdt = .00005;  // half time step

    // A = Mass+Stiff*dt/2
    RCP<MAT> A = Tpetra::MatrixMatrix::add(one, false, *mass, hdt, false, *stiff);

    // B = Mass-Stiff*dt/2
    hdt = -hdt;
    RCP<MAT> B = Tpetra::MatrixMatrix::add(one, false, *mass, hdt, false, *stiff);

    B->apply(*vecLHS, *vecRHS);

    procVerbose = verbose && (comm->getRank() == 0);  // Only print on the zero processor

    LHS = Teuchos::rcp_implicit_cast<MV>(vecLHS);
    RHS = Teuchos::rcp_implicit_cast<MV>(vecRHS);

    //
    // **********Construct preconditioner***********
    //

    // Teuchos::ParameterList MLList; // Set MLList for Smoothed Aggregation
    // ML_Tpetra::SetDefaults("SA",MLList); // reset parameters ML User's Guide
    // MLList.set("smoother: type","CHEBYSHEV"); // Chebyshev smoother  ... aztec??
    // MLList.set("smoother: sweeps",3);
    // MLList.set("smoother: pre or post", "both"); // both pre- and post-smoothing

    // #ifdef HAVE_MUELU_AMESOS2
    //     MueLuList.set("coarse: type", "KLU2"); // solve with serial direct
    //     solver KLU
    // #else
    //     MLList.set("coarse: type", "Jacobi");     // not recommended
    //     puts("Warning: Iterative coarse grid solve");
    // #endif

    //
    // RCP<OP> prec = MueLu::CreateTpetraPreconditioner(A, MLList );
    // assert(prec != Teuchos::null);

    // Create the Belos preconditioned operator from the preconditioner.
    // NOTE:  This is necessary because Belos expects an operator to apply the
    //        preconditioner with Apply() NOT ApplyInverse().
    // RCP<Belos::EpetraPrecOp> belosPrec = rcp( new Belos::EpetraPrecOp( Prec )
    // );

    // Create parameter list for the PCPG solver manager
    const int numGlobalElements = RHS->getGlobalLength();
    if (maxIters == -1)
      maxIters = numGlobalElements / blocksize - 1;  // maximum number of iterations to run

    RCP<ParameterList> belosList = rcp(new ParameterList());
    belosList->set("Block Size",
                   blocksize);  // Blocksize to be used by iterative solver
    belosList->set("Maximum Iterations",
                   maxIters);  // Maximum number of iterations allowed
    belosList->set("Convergence Tolerance",
                   tol);  // Relative convergence tolerance requested
    belosList->set("Num Deflated Blocks",
                   maxDeflate);  // Number of vectors in seed space
    belosList->set("Num Saved Blocks",
                   maxSave);                     // Number of vectors saved from old spaces
    belosList->set("Orthogonalization", ortho);  // Orthogonalization type

    if (numRhs > 1) {
      belosList->set("Show Maximum Residual Norm Only",
                     true);  // although numRhs = 1.
    }
    if (verbose) {
      belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails +
                                      Belos::FinalSummary + Belos::StatusTestDetails);
      if (frequency > 0)
        belosList->set("Output Frequency", frequency);
    } else
      belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::FinalSummary);

    // Construct a preconditioned linear problem
    RCP<LinearProblem> problem = rcp(new LinearProblem(A, LHS, RHS));

    // problem->setLeftPrec( prec );

    bool set = problem->setProblem();
    if (set == false) {
      if (procVerbose)
        std::cout << std::endl
                  << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    // Create an iterative solver manager.
    RCP<PCPGSolMgr> solver = rcp(new PCPGSolMgr(problem, belosList));

    //
    // *******************************************************************
    // ************************* Iterate PCPG ****************************
    // *******************************************************************
    //
    if (procVerbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << numGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numRhs << std::endl;
      std::cout << "Block size used by solver: " << blocksize << std::endl;
      std::cout << "Maximum number of iterations allowed: " << maxIters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }
    bool badRes;
    for (int timeStep = 0; timeStep < numTimeStep; timeStep++) {
      if (timeStep) {
        // old epetra B->multiply(false, *vecLHS, *vecRHS); // rhs_new :=
        // B*lhs_old,
        B->apply(*LHS, *RHS);  // rhs_new := B*lhs_old,
        set = problem->setProblem(LHS, RHS);
        if (set == false) {
          if (procVerbose)
            std::cout << std::endl
                      << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
          return -1;
        }
      }  // if timeStep
      std::vector<ST> rhs_norm(numRhs);
      MVT::MvNorm(*RHS, rhs_norm);
      std::cout << "                  RHS norm is ... " << rhs_norm[0] << std::endl;

      // Perform solve
      Belos::ReturnType ret = solver->solve();

      // Compute actual residuals.
      badRes = false;
      std::vector<ST> actual_resids(numRhs);
      MV resid(map, numRhs);
      OPT::Apply(*A, *LHS, resid);
      MVT::MvAddMv(-1.0, resid, 1.0, *RHS, resid);
      MVT::MvNorm(resid, actual_resids);
      MVT::MvNorm(*RHS, rhs_norm);
      std::cout << "                    RHS norm is ... " << rhs_norm[0] << std::endl;

      if (procVerbose) {
        std::cout << "---------- Actual Residuals (normalized) ----------" << std::endl
                  << std::endl;
        for (int i = 0; i < numRhs; i++) {
          double actRes = actual_resids[i] / rhs_norm[i];
          std::cout << "Problem " << i << " : \t" << actRes << std::endl;
          if (actRes > tol)
            badRes = true;
        }
      }
      if (ret != Belos::Converged || badRes) {
        success = false;
        break;
      }
    }  // for timeStep

    if (procVerbose) {
      if (success)
        std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
      else
        std::cout << std::endl << "ERROR:  Belos did not converge!" << std::endl;
    }
  }  // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

int main(int argc, char *argv[]) { 
  return run<double>(argc, argv);
  // return run<float>(argc, argv); }
}

// end PCPGTpetraExFile.cpp
