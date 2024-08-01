// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This example demonstrates how to clone a Hierarchy from one type of compute
// node to another.

#include <iostream>

#include <Xpetra_MultiVectorFactory.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_ParameterListInterpreter.hpp>  // TODO: move into MueLu.hpp
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Ifpack2Smoother.hpp"
#include <MueLu_Utilities.hpp>

#include <MueLu_MutuallyExclusiveTime.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>  // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>   // => This header defines Belos::MueLuOp
#endif

int main(int argc, char* argv[]) {
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;  // reference count pointers
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();

  // Instead of checking each time for rank, create a rank 0 stream
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& fancyout  = *fancy;
  fancyout.setOutputToRootOnly(0);

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  Teuchos::CommandLineProcessor clp(false);

  GO nx = 100, ny = 100, nz = 100;
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace2D");  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                                       // manage parameters of Xpetra

  std::string xmlFileName = "scalingTest.xml";
  clp.setOption("xml", &xmlFileName, "read parameters from a file. Otherwise, this example uses by default 'scalingTest.xml'");
  int amgAsPrecond = 1;
  clp.setOption("precond", &amgAsPrecond, "apply multigrid as preconditioner");
  int amgAsSolver = 0;
  clp.setOption("fixPoint", &amgAsSolver, "apply multigrid as solver");
  bool printTimings = true;
  clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
  int writeMatricesOPT = -2;
  clp.setOption("write", &writeMatricesOPT, "write matrices to file (-1 means all; i>=0 means level i)");
  double tol = 1e-12;
  clp.setOption("tol", &tol, "solver convergence tolerance");
  std::string krylovMethod = "cg";
  clp.setOption("krylov", &krylovMethod, "outer Krylov method");
  std::string optSmooType = "cheby";
  clp.setOption("smooType", &optSmooType, "smoother type ('l1-sgs', 'sgs 'or 'cheby')");
  int optSweeps = 2;
  clp.setOption("sweeps", &optSweeps, "sweeps to be used in SGS (or Chebyshev degree)");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  fancyout << "========================================================\n"
           << xpetraParameters << matrixParameters;

  // =========================================================================
  // Problem construction
  // =========================================================================
  RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time"))), tm;

  comm->barrier();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));

  RCP<const Map> map;
  RCP<MultiVector> coordinates;

  // Retrieve matrix parameters (they may have been changed on the command line), and pass them to Galeri.
  // Galeri will attempt to create a square-as-possible distribution of subdomains di, e.g.,
  //                                 d1  d2  d3
  //                                 d4  d5  d6
  //                                 d7  d8  d9
  //                                 d10 d11 d12
  // A perfect distribution is only possible when the #processors is a perfect square.
  // This *will* result in "strip" distribution if the #processors is a prime number or if the factors are very different in
  // size. For example, np=14 will give a 7-by-2 distribution.
  // If you don't want Galeri to do this, specify mx or my on the galeriList.
  Teuchos::ParameterList pl = matrixParameters.GetParameterList();
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", pl.get("nx", nx));
  galeriList.set("ny", pl.get("ny", ny));
  galeriList.set("nz", pl.get("nz", nz));
  // galeriList.set("mx", comm->getSize());
  // galeriList.set("my", 1);

  // Create map and coordinates
  // In the future, we hope to be able to first create a Galeri problem, and then request map and coordinates from it
  // At the moment, however, things are fragile as we hope that the Problem uses same map and coordinates inside
  if (matrixParameters.GetMatrixType() == "Laplace1D") {
    map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("1D", map, matrixParameters.GetParameterList());
  } else if (matrixParameters.GetMatrixType() == "Laplace2D" || matrixParameters.GetMatrixType() == "Star2D" || matrixParameters.GetMatrixType() == "Elasticity2D") {
    map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("2D", map, matrixParameters.GetParameterList());
  } else if (matrixParameters.GetMatrixType() == "Laplace3D" || matrixParameters.GetMatrixType() == "Elasticity3D") {
    map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("3D", map, matrixParameters.GetParameterList());
  }
  // Expand map to do multiple DOF per node for block problems
  if (matrixParameters.GetMatrixType() == "Elasticity2D")
    map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 2);
  if (matrixParameters.GetMatrixType() == "Elasticity3D")
    map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 3);

  if (comm->getRank() == 0) {
    GO mx = galeriList.get("mx", -1);
    GO my = galeriList.get("my", -1);
    GO mz = galeriList.get("mz", -1);
    fancyout << "Processor subdomains in x direction: " << mx << std::endl
             << "Processor subdomains in y direction: " << my << std::endl
             << "Processor subdomains in z direction: " << mz << std::endl
             << "========================================================" << std::endl;
  }

  Teuchos::ParameterList matrixParams = matrixParameters.GetParameterList();
  matrixParams.set("mx", galeriList.get("mx", -1));
  matrixParams.set("my", galeriList.get("my", -1));
  matrixParams.set("mz", galeriList.get("mz", -1));
  if (matrixParameters.GetMatrixType() == "Elasticity2D" || matrixParameters.GetMatrixType() == "Elasticity3D") {
    // Our default test case for elasticity: all boundaries of a square/cube have Neumann b.c. except left which has Dirichlet
    matrixParams.set("right boundary", "Neumann");
    matrixParams.set("bottom boundary", "Neumann");
    matrixParams.set("top boundary", "Neumann");
    matrixParams.set("front boundary", "Neumann");
    matrixParams.set("back boundary", "Neumann");
  }

  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixParameters.GetMatrixType(), map, matrixParams);
  RCP<Matrix> A = Pr->BuildMatrix();

  RCP<MultiVector> nullspace = MultiVectorFactory::Build(map, 1);
  if (matrixParameters.GetMatrixType() == "Elasticity2D" ||
      matrixParameters.GetMatrixType() == "Elasticity3D") {
    nullspace = Pr->BuildNullspace();
    A->SetFixedBlockSize((matrixParameters.GetMatrixType() == "Elasticity2D") ? 2 : 3);

  } else {
    nullspace->putScalar(one);
  }

  fancyout << "Galeri complete.\n========================================================" << std::endl;

  // =========================================================================
  // Preconditioner construction
  // =========================================================================

  // Multigrid Hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy(A));
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  FactoryManager M;

  // Smoothers
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO)optSweeps);
  ifpackList.set("relaxation: damping factor", (SC)1.0);
  if (optSmooType == "sgs") {
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
  } else if (optSmooType == "l1-sgs") {
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
    ifpackList.set("relaxation: use l1", true);
  } else if (optSmooType == "cheby") {
    ifpackType = "CHEBYSHEV";
    ifpackList.set("chebyshev: degree", (LO)optSweeps);

    if (matrixParameters.GetMatrixType() == "Laplace1D") {
      ifpackList.set("chebyshev: ratio eigenvalue", (SC)3);
    } else if (matrixParameters.GetMatrixType() == "Laplace2D") {
      ifpackList.set("chebyshev: ratio eigenvalue", (SC)7);
    } else if (matrixParameters.GetMatrixType() == "Laplace3D") {
      ifpackList.set("chebyshev: ratio eigenvalue", (SC)20);
    }
  }
  RCP<SmootherPrototype> smootherPrototype = rcp(new Ifpack2Smoother(ifpackType, ifpackList));
  M.SetFactory("Smoother", rcp(new SmootherFactory(smootherPrototype)));

  // create coarsest smoother
  RCP<SmootherPrototype> coarsestSmooProto;
  Teuchos::ParameterList coarsestSmooList;
  coarsestSmooProto                     = rcp(new Ifpack2Smoother("RILUK", coarsestSmooList));
  RCP<SmootherFactory> coarsestSmooFact = rcp(new SmootherFactory(coarsestSmooProto, Teuchos::null));
  M.SetFactory("CoarseSolver", coarsestSmooFact);

  int startLevel   = 0;
  int optMaxLevels = 10;
  H->Setup(M, startLevel, optMaxLevels);

  // Print out the hierarchy stats. We should not need this line, but for some reason the
  // print out in the hierarchy construction does not work.
  H->print(fancyout);

  // =========================================================================
  // System solution (Ax = b)
  // =========================================================================

  RCP<MultiVector> X = MultiVectorFactory::Build(map, 1);
  RCP<MultiVector> B = MultiVectorFactory::Build(map, 1);

  {
    // we set seed for reproducibility
    X->setSeed(846930886);
    X->randomize();
    A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

    Teuchos::Array<ST::magnitudeType> norms(1);
    B->norm2(norms);
    B->scale(1.0 / norms[0]);
    X->putScalar(zero);
  }

  if (amgAsSolver) {
    H->IsPreconditioner(false);
    H->Iterate(*B, *X, 25);

  } else if (amgAsPrecond) {
#ifdef HAVE_MUELU_BELOS
    // Operator and Multivector type that will be used with Belos
    typedef MultiVector MV;
    typedef Belos::OperatorT<MV> OP;
    H->IsPreconditioner(true);

    // Define Operator and Preconditioner
    Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A));  // Turns a Xpetra::Matrix object into a Belos operator
    Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO>(H));   // Turns a MueLu::Hierarchy object into a Belos operator

    // Construct a Belos LinearProblem object
    RCP<Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
    belosProblem->setLeftPrec(belosPrec);

    bool set = belosProblem->setProblem();
    if (set == false) {
      fancyout << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }

    // Belos parameter list
    int maxIts = 2000;
    Teuchos::ParameterList belosList;
    belosList.set("Maximum Iterations", maxIts);  // Maximum number of iterations allowed
    belosList.set("Convergence Tolerance", tol);  // Relative convergence tolerance requested
    belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList.set("Output Frequency", 1);
    belosList.set("Output Style", Belos::Brief);

    // Create an iterative solver manager
    RCP<Belos::SolverManager<SC, MV, OP> > solver;
    if (krylovMethod == "cg") {
      solver = rcp(new Belos::BlockCGSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
    } else if (krylovMethod == "gmres") {
      solver = rcp(new Belos::BlockGmresSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Invalid Krylov method.  Options are \"cg\" or \" gmres\".");
    }

    // Perform solve
    Belos::ReturnType ret = Belos::Unconverged;
    try {
      ret = solver->solve();

      // Get the number of iterations for this solve.
      fancyout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

    } catch (...) {
      fancyout << std::endl
               << "ERROR:  Belos threw an error! " << std::endl;
    }

    // Check convergence
    if (ret != Belos::Converged)
      fancyout << std::endl
               << "ERROR:  Belos did not converge! " << std::endl;
    else
      fancyout << std::endl
               << "SUCCESS:  Belos converged!" << std::endl;

    // Clone the preconditioner to ThrustGPU node type
    typedef KokkosClassic::ThrustGPUNode NO2;
    typedef MueLu::Hierarchy<SC, LO, GO, NO2> Hierarchy2;
    typedef Xpetra::MultiVector<SC, LO, GO, NO2> MV2;
    typedef Belos::OperatorT<MV2> OP2;

    ParameterList plClone;
    plClone.set<LocalOrdinal>("Verbose", 1);
    RCP<NO2> node           = rcp(new NO2(plClone));
    RCP<Hierarchy2> clonedH = H->clone<NO2>(node);

    // Clone A, X, B to new node type
    RCP<Xpetra::Matrix<SC, LO, GO, NO2> > clonedA = Xpetra::clone(*A, node);
    RCP<MV2> clonedX                              = Xpetra::clone(*X, node);
    clonedX->putScalar(zero);
    RCP<MV2> clonedB = Xpetra::clone(*B, node);
    clonedH->IsPreconditioner(true);

    // Define Operator and Preconditioner
    Teuchos::RCP<OP2> belosOp2   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO2>(clonedA));  // Turns a Xpetra::Matrix object into a Belos operator
    Teuchos::RCP<OP2> belosPrec2 = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO2>(clonedH));   // Turns a MueLu::Hierarchy object into a Belos operator

    // Construct a Belos LinearProblem object
    RCP<Belos::LinearProblem<SC, MV2, OP2> > belosProblem2 = rcp(new Belos::LinearProblem<SC, MV2, OP2>(belosOp2, clonedX, clonedB));
    belosProblem2->setLeftPrec(belosPrec2);

    bool set2 = belosProblem2->setProblem();
    if (set2 == false) {
      fancyout << "\nERROR:  Belos::LinearProblem for new node failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }
    // Create an iterative solver manager
    RCP<Belos::SolverManager<SC, MV2, OP2> > solver2;
    if (krylovMethod == "cg") {
      solver2 = rcp(new Belos::BlockCGSolMgr<SC, MV2, OP2>(belosProblem2, rcp(&belosList, false)));
    } else if (krylovMethod == "gmres") {
      solver2 = rcp(new Belos::BlockGmresSolMgr<SC, MV2, OP2>(belosProblem2, rcp(&belosList, false)));
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Invalid Krylov method.  Options are \"cg\" or \" gmres\".");
    }
    // Perform solve
    Belos::ReturnType ret2 = Belos::Unconverged;
    try {
      ret2 = solver2->solve();
      // Get the number of iterations for this solve.
      fancyout << "Number of iterations performed for this solve: " << solver2->getNumIters() << std::endl;
    } catch (...) {
      fancyout << std::endl
               << "ERROR:  Belos threw an error! " << std::endl;
    }
    // Check convergence
    if (ret2 != Belos::Converged)
      fancyout << std::endl
               << "ERROR:  Belos did not converge! " << std::endl;
    else
      fancyout << std::endl
               << "SUCCESS:  Belos converged!" << std::endl;

    // Determine if example passed
    RCP<KokkosClassic::DefaultNode::DefaultNodeType> defaultNode =
        rcp(new Tpetra::KokkosClassic::DefaultNode::DefaultNodeType(pl));
    RCP<MV> clonedXcpu = Xpetra::clone(*clonedX, defaultNode);
    clonedXcpu->update(1.0, *X, -1.0);
    Scalar norm;
    clonedXcpu->norm2(Teuchos::arrayView(&norm, 1));
    std::cout << "\nNorm of serial node soln - ThrustGPU node soln = "
              << norm << std::endl;

    bool passed = false;
    if (norm <= Scalar(1e-10))
      passed = true;
    if (passed)
      std::cout << "Example Passed!" << std::endl;
    else
      std::cout << "Example Failed!" << std::endl;
  }
#endif  // ifdef HAVE_MUELU_BELOS

  return 0;
}  // main
