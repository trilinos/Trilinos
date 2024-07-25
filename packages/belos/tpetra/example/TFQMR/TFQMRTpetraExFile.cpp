// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver reads a problem from a Harwell-Boeing (HB) file.
// Multiple right-hand-sides are created randomly.
// The initial guesses are all set to zero.
//
// Adapted from TFQMREpetraExFile.cpp (with original comments)

// All preconditioning has been commented out

// Ifpack2
// #include <Ifpack2_Factory.hpp>
// #include <Ifpack2_Preconditioner.hpp>

// Teuchos
#include <Teuchos_Assert.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Tpetra
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MatrixIO.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTFQMRSolMgr.hpp"
#include "BelosTpetraAdapter.hpp"

template <typename ScalarType>
int run(int argc, char *argv[]) {
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  using ST = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO = typename Tpetra::Vector<>::local_ordinal_type;
  using GO = typename Tpetra::Vector<>::global_ordinal_type;
  using NT = typename Tpetra::Vector<>::node_type;

  using MV = typename Tpetra::MultiVector<ST, LO, GO, NT>;
  using OP = typename Tpetra::Operator<ST, LO, GO, NT>;
  using MAP = typename Tpetra::Map<LO, GO, NT>;
  using MAT = typename Tpetra::CrsMatrix<ST, LO, GO, NT>;

  using MVT = typename Belos::MultiVecTraits<ST, MV>;
  using OPT = typename Belos::OperatorTraits<ST, MV, OP>;

  using MT = typename Teuchos::ScalarTraits<ST>::magnitudeType;
  using STM = typename Teuchos::ScalarTraits<MT>;

  using LinearProblem = typename Belos::LinearProblem<ST, MV, OP>;
  // using Preconditioner = typename Ifpack2::Preconditioner<ST, LO, GO, NT>;
  using TFQMRSolMgr = ::Belos::TFQMRSolMgr<ST, MV, OP>;

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  bool verbose = false;
  bool success = true;

  try {
    bool procVerbose = true;
    bool leftPrec = true;  // use left preconditioning to solve these linear systems
    int frequency = -1;    // how often residuals are printed by solver
    int numRhs = 1;
    int maxiters = -1;  // maximum iterations allowed
    std::string filename("osrirr1.hb");
    MT tol = STM::squareroot (STM::eps()); // 1.e-5  // relative residual tolerance

    if (procVerbose)
      std::cout << "Init cli processor" << std::endl;

    Teuchos::CommandLineProcessor cmdp(false, true);
    cmdp.setOption("verbose", "quiet", &verbose, "Print messages and results.");
    cmdp.setOption("left-prec", "right-prec", &leftPrec, "Left preconditioning or right.");
    cmdp.setOption("frequency", &frequency, "Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename", &filename, "Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("tol", &tol, "Relative residual tolerance used by GMRES solver.");
    cmdp.setOption("num-rhs", &numRhs, "Number of right-hand sides to be solved for.");
    cmdp.setOption("maxiters", &maxiters,
                   "Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    if (cmdp.parse(argc, argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose

    // Get the problem
    if (procVerbose)
      std::cout << "Read Matrix from HB file" << std::endl;
    RCP<MAT> A;
    Tpetra::Utils::readHBMatrix(filename, comm, A);
    if (procVerbose)
      std::cout << "Read Matrix from HB file OK" << std::endl;
    RCP<const MAP> map = A->getRowMap();
    if (procVerbose)
      std::cout << "Read Map OK" << std::endl;

    procVerbose = verbose && (comm->getRank() == 0); // Only print on zero processor

    // Create the preconditioner
    // if (procVerbose)
    //   std::cout << "Create preconditioner" << std::endl;
    // std::string precType = "RILUK";
    // RCP<Preconditioner> prec = Ifpack2::Factory::create<MAT>(precType, A);
    // assert(prec != Teuchos::null);

    // // Specify parameters for the preconditioner
    // int lFill = 2;     // if (argc > 2) lFill = atoi(argv[2]);
    // int overlap = 2;   // if (argc > 3) overlap = atoi(argv[3]);
    // ST absThresh = 0.0;  // if (argc > 4) AThresh = atof(argv[4]);
    // ST relThresh = 1.0;  // if (argc >5) RThresh = atof(argv[5]);

    // if (procVerbose) {
    //   std::cout << std::endl << std::endl;
    //   std::cout << "Constructing ILU preconditioner" << std::endl;
    //   std::cout << "Using Level of fill = " << lFill << std::endl;
    //   std::cout << "Using Level Overlap = " << overlap << std::endl;
    //   std::cout << "Using Absolute Threshold Value of " << absThresh << std::endl;
    //   std::cout << "Using Relative Threshold Value of " << relThresh << std::endl;
    // }

    // if (procVerbose)
    //   std::cout << "Init preconditioner" << std::endl;

    // ParameterList precParams;
    // precParams.set("fact: iluk level-of-fill", lFill);
    // precParams.set("fact: iluk level-of-overlap", overlap);
    // precParams.set("fact: absolute threshold", absThresh);
    // precParams.set("fact: relative threshold", relThresh);
    // prec->setParameters(precParams);

    // // Initialize and build the preconditioner
    // prec->initialize();
    // prec->compute();

    // Specify parameters for the TFQMR solver manager
    const int numGlobalElements = map->getGlobalNumElements();
    if (maxiters == -1)
      maxiters = numGlobalElements - 1;  // maximum number of iterations to run
    RCP<ParameterList> belosList = rcp(new ParameterList());
    belosList->set("Maximum Iterations", maxiters);  // Maximum number of iterations allowed
    belosList->set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
    if (leftPrec)
      belosList->set("Explicit Residual Test", true);  // Need to check for the explicit residual before returning
    if (verbose) {
      belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);
      if (frequency > 0)
        belosList->set("Output Frequency", frequency);
    } else
      belosList->set("Verbosity", Belos::Errors + Belos::Warnings);

    // Construct solution std::vector and random right-hand-sides
    RCP<MV> X = rcp(new MV(map, numRhs));
    X->putScalar(0.0);
    RCP<MV> B = rcp(new MV(map, numRhs));
    B->randomize();
    RCP<LinearProblem> problem = rcp(new LinearProblem(A, X, B));
    // if (leftPrec)
    //   problem->setLeftPrec(prec);
    // else
    //   problem->setRightPrec(prec);

    bool set = problem->setProblem();
    if (set == false) {
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    // Create the TFQMR solver manager
    RCP<TFQMRSolMgr> solver = rcp(new TFQMRSolMgr(problem, belosList));

    // Print out information about problem
    if (procVerbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << numGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numRhs << std::endl;
      std::cout << "Max number of TFQMR iterations: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }

    // Perform solve
    if (procVerbose) {
      std::cout << "Solving..." << std::endl;
    }
    Belos::ReturnType ret = solver->solve();
    if (procVerbose) {
      std::cout << "Solve end" << std::endl;
    }

    // Compute actual residuals.
    bool badRes = false;
    std::vector<ST> actualResids(numRhs);
    std::vector<ST> rhsNorm(numRhs);
    MV resids(map, numRhs);
    OPT::Apply(*A, *X, resids);
    MVT::MvAddMv(-1.0, resids, 1.0, *B, resids);
    MVT::MvNorm(resids, actualResids);
    MVT::MvNorm(*B, rhsNorm);
    if (procVerbose) {
      std::cout << "---------- Actual Residuals (normalized) ----------" << std::endl << std::endl;
      for (int i = 0; i < numRhs; i++) {
        ST actRes = actualResids[i] / rhsNorm[i];
        std::cout << "Problem " << i << " : \t" << actRes << std::endl;
        if (actRes > tol)
          badRes = true;
      }
    }

    if (ret != Belos::Converged || badRes) {
      success = false;
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Belos did not converge!" << std::endl;
    } else {
      success = true;
      if (procVerbose)
        std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

int main(int argc, char *argv[]) {
  return run<double>(argc, argv);
  // return run<float>(argc, argv);
}
// end TFQMRpetraExFile.cpp
