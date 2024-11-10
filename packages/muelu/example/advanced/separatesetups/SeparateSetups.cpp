// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Note: use --help to list available options.

/*
   This example demonstrates:
   1) How to create the hierarchy in a separate step from the smoothers.
   2) How grid transfers can be generated using a different matrix than the smoothers.
   3) How grid transfers P & R can be generated using matrix A, and coarse grid operators
   can be generated via R*B*P.
   */

#include <iostream>

// Teuchos
#include "Teuchos_StandardCatchMacros.hpp"

// MueLu
#include "MueLu.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_RAPFactory.hpp"

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;
  using Teuchos::rcp;

  //
  // MPI initialization
  //

  Teuchos::oblackholestream blackhole;

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    //
    // Process command line arguments
    //
    Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 81);  // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);                  // manage parameters of xpetra

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
      default:;
    }

    if (comm->getRank() == 0) std::cout << xpetraParameters << matrixParameters;

    //
    // Setup test case (Ax = b)
    //

    // Distribution
    RCP<const Map> map = MapFactory::Build(lib, matrixParameters.GetNumGlobalElements(), 0, comm);

    // Matrix
    RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
    RCP<Matrix> A = Pr->BuildMatrix();

    // User defined nullspace
    RCP<MultiVector> nullSpace = VectorFactory::Build(map, 1);
    nullSpace->putScalar((SC)1.0);

    // Define B
    RCP<Vector> X = VectorFactory::Build(map, 1);
    RCP<Vector> B = VectorFactory::Build(map, 1);
    X->setSeed(846930886);
    X->randomize();
    A->apply(*X, *B, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // X = 0
    X->putScalar((SC)0.0);

    //
    // Create a multigrid configuration
    //

    // Transfer operators
    RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
    RCP<SaPFactory> SaPFact               = rcp(new SaPFactory());
    RCP<TransPFactory> RFact              = rcp(new TransPFactory());

    FactoryManager M;
    M.SetFactory("Ptent", TentativePFact);
    M.SetFactory("P", SaPFact);
    M.SetFactory("R", RFact);

    M.SetFactory("Smoother", Teuchos::null);      // skips smoother setup
    M.SetFactory("CoarseSolver", Teuchos::null);  // skips coarsest solve setup

    //
    // Multigrid setup phase
    //

    int startLevel = 0;
    int maxLevels  = 10;

    std::cout << "=============== Setup transfers only ====================" << std::endl;

    Hierarchy H;
    H.SetDefaultVerbLevel(MueLu::Medium);

    RCP<Level> finestLevel = H.GetLevel();
    finestLevel->Set("A", A);
    finestLevel->Set("Nullspace", nullSpace);

    // Indicate which Hierarchy operators we want to keep
    H.Keep("P", SaPFact.get());             // SaPFact is the generating factory for P.
    H.Keep("R", RFact.get());               // RFact is the generating factory for R.
    H.Keep("Ptent", TentativePFact.get());  // SaPFact is the generating factory for P.

    H.Setup(M, startLevel, maxLevels);

    std::cout << "=============== Setup smoothers only ====================" << std::endl;

    // Create a new A.
    RCP<Matrix> newA = Pr->BuildMatrix();
    finestLevel->Set("A", newA);

    // Create Gauss-Seidel smoother.
    std::string ifpackType = "RELAXATION";
    Teuchos::ParameterList ifpackList;
    ifpackList.set("relaxation: sweeps", (LO)3);
    ifpackList.set("relaxation: damping factor", (SC)1.0);
    RCP<SmootherPrototype> smootherPrototype = rcp(new TrilinosSmoother(ifpackType, ifpackList));

    M.SetFactory("Smoother", rcp(new SmootherFactory(smootherPrototype)));

    // Create coarsest solver.
    RCP<SmootherPrototype> coarseSolverPrototype = rcp(new DirectSolver());
    RCP<SmootherFactory> coarseSolverFact        = rcp(new SmootherFactory(coarseSolverPrototype, Teuchos::null));
    M.SetFactory("CoarseSolver", coarseSolverFact);

    // Note that we pass the number of levels back in.
    H.Setup(M, startLevel, H.GetNumLevels());

    std::cout << "=============== Solve ====================" << std::endl;

    //
    // Solve Ax = B
    //

    LO nIts = 9;
    H.Iterate(*B, *X, nIts);

    //
    // Print relative residual norm
    //

    typename Teuchos::ScalarTraits<SC>::magnitudeType residualNorms = Utilities::ResidualNorm(*A, *X, *B)[0];
    if (comm->getRank() == 0) {
      std::ios::fmtflags f(std::cout.flags());
      std::cout << "||Residual|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorms << std::endl;
      std::cout.flags(f);
    }

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
