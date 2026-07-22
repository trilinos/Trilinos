// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

#include <Amesos2.hpp>
#include <Amesos2_config.h>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_IO.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosStatusTestCombo.hpp>
#include <BelosXpetraStatusTestGenResSubNorm.hpp>

#include <Galeri_XpetraMaps.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

#include <MueLu.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_ConfigDefs.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RepartitionHeuristicFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"

#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

#include <Xpetra_Map.hpp>

int main(int argc, char* argv[]) {
  using Scalar        = Tpetra::MultiVector<>::scalar_type;
  using LocalOrdinal  = Tpetra::MultiVector<>::local_ordinal_type;
  using GlobalOrdinal = Tpetra::MultiVector<>::global_ordinal_type;
  using Node          = Tpetra::MultiVector<>::node_type;
#include "MueLu_UseShortNames.hpp"
  bool success = false;
  try {
    Tpetra::ScopeGuard tpetraScope(&argc, &argv);
    {
      using Teuchos::ParameterList;
      using Teuchos::RCP;
      using Teuchos::rcp;

      RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

      // ================================
      // Problem construction
      // ================================
      Teuchos::ParameterList galeriList;
      galeriList.set("nx", 200);
      galeriList.set("ny", 400);
      galeriList.set("mx", comm->getSize());
      galeriList.set("my", 1);

      // Create map
      RCP<const Map> xpetra_map = Galeri::Xpetra::CreateMap<LO, GO, NO>(Xpetra::UseTpetra, "Cartesian2D", comm, galeriList);
      // Create coordinates
      RCP<MultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("2D", xpetra_map, galeriList);

      //! [2DLaplacianOperator begin]
      // Create the matrix
      RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> galeriProblem =
          Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", xpetra_map, galeriList);
      RCP<Matrix> A = galeriProblem->BuildMatrix();
      //! [2DLaplacianOperator end]

      //! [DefineNearNullSpace begin]
      // define near null space
      RCP<MultiVector> nullspace = MultiVectorFactory::Build(xpetra_map, 1);
      nullspace->putScalar((SC)1.0);
      //! [DefineNearNullSpace end]

      Teuchos::Array<typename Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);

      //! [CreateNewHierarchy begin]
      // create new hierarchy
      RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H;
      //! [CreateNewHierarchy end]

      //! [InstantiateNewHierarchyObject begin]
      // Instantiate new Hierarchy object
      H = rcp(new Hierarchy());
      H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
      H->SetMaxCoarseSize((GO)50);
      //! [InstantiateNewHierarchyObject end]

      //! [CreateFineLevelObject begin]
      // create a fine level object
      RCP<Level> Finest = H->GetLevel();
      Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
      Finest->Set("A", A);
      Finest->Set("Nullspace", nullspace);
      Finest->Set("Coordinates", coordinates);
      //! [CreateFineLevelObject end]

      //! [DefineFactoryManager begin]     // define a factory manager
      MueLu::FactoryManager M;
      M.SetKokkosRefactor(false);
      //! [DefineFactoryManager end]

      RCP<UncoupledAggregationFactory> AggregationFact = rcp(new UncoupledAggregationFactory());
      AggregationFact->SetMinNodesPerAggregate(2);
      AggregationFact->SetMaxNeighAlreadySelected(0);
      AggregationFact->SetOrdering("natural");
      // AggregationFact->SetPhase3AggCreation(0.5);
      M.SetFactory("Aggregates", AggregationFact);

      //! [DeclareSomeFactories begin]
      // declare some factories (potentially overwrite default factories)
      RCP<SaPFactory> PFact = rcp(new SaPFactory());
      PFact->SetParameter("sa: damping factor", Teuchos::ParameterEntry(4. / 3));

      RCP<Factory> RFact = rcp(new TransPFactory());

      RCP<RAPFactory> AcFact = rcp(new RAPFactory());
      AcFact->setVerbLevel(Teuchos::VERB_HIGH);

      H->SetImplicitTranspose(true);
      Teuchos::ParameterList Aclist = *(AcFact->GetValidParameterList());
      Aclist.set("transpose: use implicit", true);
      AcFact->SetParameterList(Aclist);
      //! [DeclareSomeFactories end]

      //! [ConfigureFactoryManager begin]
      // configure factory manager (no repartitioning)
      M.SetFactory("P", PFact);
      M.SetFactory("R", RFact);
      M.SetFactory("A", AcFact);
      //! [ConfigureFactoryManager end]

      //! [DefineSmootherObject begin]
      // define smoother object
      std::string ifpackType;
      Teuchos::ParameterList ifpackList;
      ifpackList.set("relaxation: sweeps", (LO)2);
      ifpackList.set("relaxation: damping factor", (SC)1.0);
      ifpackType = "RELAXATION";
      ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
      //! [DefineSmootherObject end]

      //! [CreateSmootherFactory begin]
      // create smoother factory
      RCP<SmootherPrototype> smootherPrototype = rcp(new TrilinosSmoother(ifpackType, ifpackList));
      M.SetFactory("Smoother", rcp(new SmootherFactory(smootherPrototype)));
      //! [CreateSmootherFactory end]

      //! [SetupMultigridHierarchy begin]
      // setup multigrid hierarchy
      int startLevel = 0;
      H->Setup(M, startLevel, 10);
      //! [SetupMultigridHierarchy end]

      //! [DefineXB begin]
      // Define X, B
      RCP<MultiVector> X = MultiVectorFactory::Build(xpetra_map, 1);
      RCP<MultiVector> B = MultiVectorFactory::Build(xpetra_map, 1);

      X->setSeed(846930886);
      X->randomize();
      A->apply(*X, *B, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);
      B->norm2(norms);
      B->scale(1.0 / norms[0]);
      //! [DefineXB end]

      //! [OperatorAndMultivectorTypeBelos begin]
      // Operator and Multivector type that will be used with Belos
      typedef Belos::OperatorT<MultiVector> OP;
      H->IsPreconditioner(true);

      // Define Operator and Preconditioner
      Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A));  // Turns a Xpetra::Operator object into a Belos operator
      Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO>(H));   // Turns a MueLu::Hierarchy object into a Belos operator

      // Construct a Belos LinearProblem object
      RCP<Belos::LinearProblem<SC, MultiVector, OP>> belosProblem = rcp(new Belos::LinearProblem<SC, MultiVector, OP>(belosOp, X, B));
      belosProblem->setLeftPrec(belosPrec);

      bool set = belosProblem->setProblem();
      if (set == false) {
        if (comm->getRank() == 0)
          std::cout << std::endl
                    << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
        return EXIT_FAILURE;
      }
      //! [OperatorAndMultivectorTypeBelos end]

      //! [BelosParameterList begin]
      // Belos parameter list
      Teuchos::ParameterList belosList;
      belosList.set("Maximum Iterations", 100);      // Maximum number of iterations allowed
      belosList.set("Convergence Tolerance", 1e-7);  // Relative convergence tolerance requested
      belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
      belosList.set("Output Frequency", 1);
      belosList.set("Output Style", Belos::Brief);

      // Create an iterative solver manager
      RCP<Belos::SolverManager<SC, MultiVector, OP>> solver = rcp(new Belos::BlockCGSolMgr<SC, MultiVector, OP>(belosProblem, rcp(&belosList, false)));
      //! [BelosParameterList end]

      // Perform solve
      Belos::ReturnType ret = Belos::Unconverged;
      try {
        //! [SolveLinearSystem begin]
        // solve linear system
        ret = solver->solve();
        //! [SolveLinearSystem end]

        // Get the number of iterations for this solve.
        if (comm->getRank() == 0)
          std::cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;
      } catch (...) {
        if (comm->getRank() == 0)
          std::cout << std::endl
                    << "ERROR:  Belos threw an error! " << std::endl;
      }

      //! [CheckConvergence begin]
      // Check convergence
      if (ret != Belos::Converged) {
        if (comm->getRank() == 0) std::cout << std::endl
                                            << "ERROR:  Belos did not converge! " << std::endl;
      } else {
        if (comm->getRank() == 0) std::cout << std::endl
                                            << "SUCCESS:  Belos converged!" << std::endl;
      }
      //! [CheckConvergence end]
    }  // end of Tpetra::ScopeGuard
    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}  // main
