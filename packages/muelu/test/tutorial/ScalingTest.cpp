// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <unistd.h>
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultComm.hpp>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_IO.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RepartitionHeuristicFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_Zoltan2Interface.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"

// Belos
#ifdef HAVE_MUELU_BELOS
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosXpetraAdapter.hpp" // this header defines Belos::XpetraOp()
#include "BelosMueLuAdapter.hpp"  // this header defines Belos::MueLuOp()
#endif

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib& lib, int argc, char *argv[])
{
#include "MueLu_UseShortNames.hpp"

  using Teuchos::RCP; using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  //using Galeri::Xpetra::CreateCartesianCoordinates;

  Teuchos::oblackholestream blackhole;

  // USER GUIDE // define communicator
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  // USER GUIDE // create fancy output stream
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);
  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  #ifndef HAVE_XPETRA_INT_LONG_LONG
  *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
  #endif

  //
  // SET TEST PARAMETERS
  //

  // Default is Laplace1D with nx = 8748.
  // It's a nice size for 1D and perfect aggregation. (6561 = 3^8)
  //Nice size for 1D and perfect aggregation on small numbers of processors. (8748 = 4*3^7)
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748); // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                   // manage parameters of xpetra

  // Custom command line parameters
  // - Debug
  int optDump    = 0;                     clp.setOption("dump",           &optDump,               "write matrix to file");
  int optTimings = 0;                     clp.setOption("timings",        &optTimings,            "print timings to screen");

  // - Levels
  LO  optMaxLevels     = 2;               clp.setOption("maxLevels",      &optMaxLevels,          "maximum number of levels allowed");
  int optMaxCoarseSize = 50;              clp.setOption("maxCoarseSize",  &optMaxCoarseSize,      "maximum #dofs in coarse operator"); //FIXME clp doesn't like long long int

  // - Smoothed-Aggregation
  double optSaDamping = 4. / 3;           clp.setOption("saDamping",      &optSaDamping,          "prolongator damping factor");

  // - Aggregation
  std::string optAggOrdering = "natural"; clp.setOption("aggOrdering",    &optAggOrdering,        "aggregation ordering strategy (natural, random, graph)");
  int optMinPerAgg = 2;                   clp.setOption("minPerAgg",      &optMinPerAgg,          "minimum #DOFs per aggregate");
  int optMaxNbrSel = 0;                   clp.setOption("maxNbrSel",      &optMaxNbrSel,          "maximum # of nbrs allowed to be in other aggregates");

  // - R
  int optExplicitR = 1;                   clp.setOption("explicitR",      &optExplicitR,          "restriction will be explicitly stored as transpose of prolongator");

  // - Smoothers
  std::string optSmooType = "sgs";        clp.setOption("smooType",       &optSmooType,           "smoother type ('l1-sgs', 'sgs 'or 'cheby')");
  int optSweeps = 2;                      clp.setOption("sweeps",         &optSweeps,             "sweeps to be used in SGS (or Chebyshev degree)");

  // - Repartitioning
#if defined(HAVE_MPI) && defined(HAVE_MUELU_ZOLTAN2)
  int optRepartition = 1;                 clp.setOption("repartition",    &optRepartition,        "enable repartitioning (0=no repartitioning, 1=Zoltan2 RCB");
  LO optMinRowsPerProc = 2000;            clp.setOption("minRowsPerProc", &optMinRowsPerProc,     "min #rows allowable per proc before repartitioning occurs");
  double optNnzImbalance = 1.2;           clp.setOption("nnzImbalance",   &optNnzImbalance,       "max allowable nonzero imbalance before repartitioning occurs");
#else
  int optRepartition = 0;
#endif // HAVE_MPI && HAVE_MUELU_ZOLTAN2

  // - Solve
  int    optFixPoint = 1;                 clp.setOption("fixPoint",       &optFixPoint,           "apply multigrid as solver");
  int    optPrecond  = 1;                 clp.setOption("precond",        &optPrecond,            "apply multigrid as preconditioner");
  LO     optIts      = 10;                clp.setOption("its",            &optIts,                "number of multigrid cycles");
  double optTol      = 1e-7;              clp.setOption("tol",            &optTol,                "stopping tolerance for Krylov method");

  switch (clp.parse(argc, argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  RCP<TimeMonitor> globalTimeMonitor = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time")));

  matrixParameters.check();
  xpetraParameters.check();
  // TODO: check custom parameters
  std::transform(optSmooType.begin(), optSmooType.end(), optSmooType.begin(), ::tolower);

  if (comm->getRank() == 0) {
    std::cout << xpetraParameters << matrixParameters;
    // TODO: print custom parameters // Or use paramList::print()!
  }

  //
  // CREATE INITIAL MATRIX                                                          */
  //
  RCP<const Map> map;
  RCP<Matrix> A;

  RCP<MultiVector> coordinates;
  {
    TimeMonitor tm(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build"));

    map = MapFactory::Build(lib, matrixParameters.GetNumGlobalElements(), 0, comm);
    Teuchos::RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Matrix vs. CrsMatrixWrap
    A = Pr->BuildMatrix();

    if (matrixParameters.GetMatrixType() == "Laplace1D") {
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("1D", map, matrixParameters.GetParameterList());
    }
    else if (matrixParameters.GetMatrixType() == "Laplace2D") {
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("2D", map, matrixParameters.GetParameterList());
    }
    else if (matrixParameters.GetMatrixType() == "Laplace3D") {
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("3D", map, matrixParameters.GetParameterList());
    }
  }

  //
  //
  //

  // dump matrix to file
  if (optDump) {
    std::string fileName = "Amat.mm";
    Xpetra::IO<SC,LO,GO,NO>::Write(fileName, *A);
  }

  //! [DefineNearNullSpace begin]
  // define near null space
  RCP<MultiVector> nullspace = MultiVectorFactory::Build(map, 1);
  nullspace->putScalar( (SC) 1.0);
 //! [DefineNearNullSpace end]
  Teuchos::Array<typename Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);

  nullspace->norm1(norms);
  if (comm->getRank() == 0)
    std::cout << "||NS|| = " << norms[0] << std::endl;

  //! [CreateNewHierarchy begin]
  // create new hierarchy
  RCP<MueLu::Hierarchy<SC, LO, GO, NO> > H;
  //! [CreateNewHierarchy end]

  //
  //
  // SETUP
  //
  //

  {
    TimeMonitor tm(*TimeMonitor::getNewTimer("ScalingTest: 2 - MueLu Setup"));

    //
    // Hierarchy
    //

    //! [InstantiateNewHierarchyObject begin]
    // instantiate new Hierarchy object
    H = rcp(new Hierarchy());
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    H->SetMaxCoarseSize((GO) optMaxCoarseSize);
    //! [InstantiateNewHierarchyObject end]

    //
    // Finest level
    //

    //! [CreateFineLevelObject begin]
    // create a fine level object
    RCP<Level> Finest = H->GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A",           A);
    Finest->Set("Nullspace",   nullspace);
    Finest->Set("Coordinates", coordinates); //! [CreateFineLevelObject end] //FIXME: XCoordinates, YCoordinates, ..
       

    //
    // FactoryManager
    //

    //! [DefineFactoryManager begin]     // define a factory manager
    FactoryManager M;
    M.SetKokkosRefactor(false);
    //! [DefineFactoryManager end]

    //
    //
    // Aggregation
    //

    {
      RCP<UncoupledAggregationFactory> AggregationFact = rcp(new UncoupledAggregationFactory());
      *out << "========================= Aggregate option summary =========================" << std::endl;
      *out << "min DOFs per aggregate :                " << optMinPerAgg << std::endl;
      *out << "min # of root nbrs already aggregated : " << optMaxNbrSel << std::endl;
      AggregationFact->SetMinNodesPerAggregate(optMinPerAgg);  //TODO should increase if run anything othpermRFacter than 1D
      AggregationFact->SetMaxNeighAlreadySelected(optMaxNbrSel);
      std::transform(optAggOrdering.begin(), optAggOrdering.end(), optAggOrdering.begin(), ::tolower);
      if (optAggOrdering == "natural" || optAggOrdering == "random" || optAggOrdering == "graph") {
        *out << "aggregate ordering :                    " << optAggOrdering << std::endl;
        AggregationFact->SetOrdering(optAggOrdering);
      } else {
        std::string msg = "main: bad aggregation option """ + optAggOrdering + """.";
        throw(MueLu::Exceptions::RuntimeError(msg));
      }
      //AggregationFact->SetPhase3AggCreation(0.5);
      M.SetFactory("Aggregates", AggregationFact);

    *out << "=============================================================================" << std::endl;
    }

    //
    // Transfer
    //

    {
      //
      // Non rebalanced factories
      //

      //! [DeclareSomeFactories begin]       
      // declare some factories (potentially overwrite default factories)
      RCP<SaPFactory> PFact = rcp(new SaPFactory());
      PFact->SetParameter("sa: damping factor", Teuchos::ParameterEntry(optSaDamping));

      RCP<Factory>    RFact = rcp(new TransPFactory());

      RCP<RAPFactory> AcFact = rcp(new RAPFactory());
      AcFact->setVerbLevel(Teuchos::VERB_HIGH);
      //! [DeclareSomeFactories end]

      if (!optExplicitR) {
        H->SetImplicitTranspose(true);
        Teuchos::ParameterList Aclist = *(AcFact->GetValidParameterList());
        Aclist.set("transpose: use implicit", true);
        AcFact->SetParameterList(Aclist);
        if (comm->getRank() == 0) std::cout << "\n\n* ***** USING IMPLICIT RESTRICTION OPERATOR ***** *\n" << std::endl;
      }

      //
      // Repartitioning (if needed)
      //

      if (optRepartition == 0) {
        // No repartitioning

        //! [ConfigureFactoryManager begin]         
        // configure factory manager
        M.SetFactory("P", PFact);
        M.SetFactory("R", RFact);
        M.SetFactory("A", AcFact);
        //! [ConfigureFactoryManager end] 

      } else {
#if defined(HAVE_MPI) && defined(HAVE_MUELU_ZOLTAN2)
        // Repartitioning

        // The Factory Manager will be configured to return the rebalanced versions of P, R, A by default.
        // Everytime we want to use the non-rebalanced versions, we need to explicitly define the generating factory.
        RFact->SetFactory("P", PFact);
        //
        AcFact->SetFactory("P", PFact);
        AcFact->SetFactory("R", RFact);

        // Transfer coordinates
        RCP<CoordinatesTransferFactory> TransferCoordinatesFact = rcp(new CoordinatesTransferFactory());
        //AcFact->AddTransferFactory(TransferCoordinatesFact); // FIXME REMOVE

        // Repartitioning heuristic (decides whether to rebalance based on params)
        RCP<RepartitionHeuristicFactory> RepHeuFact = Teuchos::rcp(new RepartitionHeuristicFactory());
        RepHeuFact->SetFactory("A", AcFact);
        RepHeuFact->SetParameter("repartition: start level", Teuchos::ParameterEntry(0));
        RepHeuFact->SetParameter("repartition: min rows per proc", Teuchos::ParameterEntry(optMinRowsPerProc));
        RepHeuFact->SetParameter("repartition: max imbalance", Teuchos::ParameterEntry(optNnzImbalance));

        // Repartitioning (creates repartition vector)
        RCP<Zoltan2Interface> ZoltanFact = rcp(new Zoltan2Interface());
        ZoltanFact->SetFactory("A", AcFact);
        ZoltanFact->SetFactory("number of partitions", RepHeuFact);
        ZoltanFact->SetFactory("Coordinates", TransferCoordinatesFact);

        // Repartitioning (creates "Importer" from "Partition")
        RCP<Factory> RepartitionFact = rcp(new RepartitionFactory());
        RepartitionFact->SetFactory("A", AcFact);
        RepartitionFact->SetFactory("number of partitions", RepHeuFact);
        RepartitionFact->SetFactory("Partition", ZoltanFact);

        // Reordering of the transfer operators
        RCP<Factory> RebalancedPFact = rcp(new RebalanceTransferFactory());
        RebalancedPFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Interpolation")));
        RebalancedPFact->SetFactory("P", PFact);
        RebalancedPFact->SetFactory("Coordinates", TransferCoordinatesFact);
        RebalancedPFact->SetFactory("Nullspace", M.GetFactory("Ptent")); // TODO

        RCP<Factory> RebalancedRFact = rcp(new RebalanceTransferFactory());
        RebalancedRFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Restriction")));
        RebalancedRFact->SetFactory("R", RFact);

        // Compute Ac from rebalanced P and R
        RCP<Factory> RebalancedAFact = rcp(new RebalanceAcFactory());
        RebalancedAFact->SetFactory("A", AcFact);

        // Configure FactoryManager
        M.SetFactory("A", RebalancedAFact);
        M.SetFactory("P", RebalancedPFact);
        M.SetFactory("R", RebalancedRFact);
        M.SetFactory("Nullspace",   RebalancedPFact);
        M.SetFactory("Coordinates", RebalancedPFact);
        M.SetFactory("Importer",    RepartitionFact);

#else
        TEUCHOS_TEST_FOR_EXCEPT(true);
#endif
      } // optRepartition

    } // Transfer

    //
    // Smoothers
    //

    {
      //! [DefineSmootherObject begin]
      // define smoother object
      std::string ifpackType;
      Teuchos::ParameterList ifpackList;
      ifpackList.set("relaxation: sweeps", (LO) optSweeps);
      ifpackList.set("relaxation: damping factor", (SC) 1.0);
      if (optSmooType == "sgs") {
        ifpackType = "RELAXATION";
        ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
      }
      //! [DefineSmootherObject end]
      else if (optSmooType == "l1-sgs") {
        ifpackType = "RELAXATION";
        ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
        ifpackList.set("relaxation: use l1", true);
      } else if (optSmooType == "cheby") {
        ifpackType = "CHEBYSHEV";
        ifpackList.set("chebyshev: degree", (LO) optSweeps);

        if (matrixParameters.GetMatrixType() == "Laplace1D") {
          ifpackList.set("chebyshev: ratio eigenvalue", (SC) 3);
        }
        else if (matrixParameters.GetMatrixType() == "Laplace2D") {
          ifpackList.set("chebyshev: ratio eigenvalue", (SC) 7);
        }
        else if (matrixParameters.GetMatrixType() == "Laplace3D") {
          ifpackList.set("chebyshev: ratio eigenvalue", (SC) 20);
        }
        // ifpackList.set("chebyshev: max eigenvalue", (double) -1.0);
        // ifpackList.set("chebyshev: min eigenvalue", (double) 1.0);
      }

      //! [CreateSmootherFactory begin]
      // create smoother factory
      RCP<SmootherPrototype> smootherPrototype = rcp(new TrilinosSmoother(ifpackType, ifpackList));
      M.SetFactory("Smoother", rcp(new SmootherFactory(smootherPrototype)));
      //! [CreateSmootherFactory end]
    }

    //
    // Setup preconditioner
    //

    //! [SetupMultigridHierarchy begin]
    // setup multigrid hierarchy
    int startLevel = 0;
    H->Setup(M, startLevel, optMaxLevels);
    //! [SetupMultigridHierarchy end]

  } // end of Setup TimeMonitor

  /*{ // some debug output
    // print out content of levels
    std::cout << "FINAL CONTENT of multigrid levels" << std::endl;
    for(LO l = 0; l < H->GetNumLevels(); l++) {
      RCP<Level> coarseLevel = H->GetLevel(l);
      coarseLevel->print(*out);
    }
    std::cout << "END FINAL CONTENT of multigrid levels" << std::endl;
  } // end debug output*/

  //
  //
  // SOLVE
  //
  //

  //! [DefineXB begin] 
  // Define X, B
  RCP<MultiVector> X = MultiVectorFactory::Build(map, 1);
  RCP<MultiVector> B = MultiVectorFactory::Build(map, 1);

  X->setSeed(846930886);
  X->randomize();
  A->apply(*X, *B, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);
  B->norm2(norms);
  B->scale(1.0/norms[0]);
  //! [DefineXB end] 

  //
  // Use AMG directly as an iterative method
  //

  if (optFixPoint) {

    X->putScalar( (SC) 0.0);

    TimeMonitor tm(*TimeMonitor::getNewTimer("ScalingTest: 3 - Fixed Point Solve"));

    H->IsPreconditioner(false);
    H->Iterate(*B, *X, optIts);

  } // optFixedPt

  //
  // Use AMG as a preconditioner in Belos
  //

#ifdef HAVE_MUELU_BELOS

  if (optPrecond) {

    RCP<TimeMonitor> tm;
    tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 5 - Belos Solve")));

    //! [OperatorAndMultivectorTypeBelos begin] 
    // Operator and Multivector type that will be used with Belos
    typedef MultiVector          MV;
    typedef Belos::OperatorT<MV> OP;
    H->IsPreconditioner(true);

    // Define Operator and Preconditioner
    Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A)); // Turns a Xpetra::Operator object into a Belos operator
    Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO>(H));  // Turns a MueLu::Hierarchy object into a Belos operator

    // Construct a Belos LinearProblem object
    RCP< Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
    belosProblem->setLeftPrec(belosPrec);

    bool set = belosProblem->setProblem();
    if (set == false) {
      if (comm->getRank() == 0)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }
    //! [OperatorAndMultivectorTypeBelos end] 

    //! [BelosParameterList begin]
    // Belos parameter list
    int maxIts = 100;
    Teuchos::ParameterList belosList;
    belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
    belosList.set("Convergence Tolerance", optTol);    // Relative convergence tolerance requested
    //belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);
    belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList.set("Output Frequency", 1);
    belosList.set("Output Style", Belos::Brief);

    // Create an iterative solver manager
    RCP< Belos::SolverManager<SC, MV, OP> > solver = rcp(new Belos::BlockCGSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
    //! [BelosParameterList end]

    // Perform solve
    Belos::ReturnType ret = Belos::Unconverged;
    try {
      {
        TimeMonitor tm2(*TimeMonitor::getNewTimer("ScalingTest: 5bis - Belos Internal Solve"));
        //! [SolveLinearSystem begin] 
        // solve linear system
        ret = solver->solve();
        //! [SolveLinearSystem end] 
      } // end of TimeMonitor

      // Get the number of iterations for this solve.
      if (comm->getRank() == 0)
        std::cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

      // Compute actual residuals.
      int numrhs = 1;
      std::vector<double> actual_resids( numrhs ); //TODO: double?
      std::vector<double> rhs_norm( numrhs );
      RCP<MultiVector> resid = MultiVectorFactory::Build(map, numrhs);

      typedef Belos::OperatorTraits<SC, MV, OP>  OPT;
      typedef Belos::MultiVecTraits<SC, MV>     MVT;

      OPT::Apply( *belosOp, *X, *resid );
      MVT::MvAddMv( -1.0, *resid, 1.0, *B, *resid );
      MVT::MvNorm( *resid, actual_resids );
      MVT::MvNorm( *B, rhs_norm );
      *out<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i = 0; i<numrhs; i++) {
        double actRes = actual_resids[i]/rhs_norm[i];
        *out<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        //if (actRes > tol) { badRes = true; }
      }

    } //try

    catch(...) {
      if (comm->getRank() == 0)
        std::cout << std::endl << "ERROR:  Belos threw an error! " << std::endl;
    }

    //! [CheckConvergence begin]
    // Check convergence
    if (ret != Belos::Converged) {
      if (comm->getRank() == 0) std::cout << std::endl << "ERROR:  Belos did not converge! " << std::endl;
    } else {
      if (comm->getRank() == 0) std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
    }
    //! [CheckConvergence end]
    tm = Teuchos::null;

  } //if (optPrecond)

#endif // HAVE_MUELU_BELOS

  //
  // Timer final summaries
  //

  globalTimeMonitor = Teuchos::null; // stop this timer before summary

  if (optTimings)
    TimeMonitor::summarize();

  //

  return EXIT_SUCCESS;
} // main_

//-----------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}


// TODO: add warning if:
// DEBUG_MODE, LONG_LONG or KLU

/* test direct solve */
/*    if (optMaxLevels == 1) {
      Teuchos::ParameterList amesosList;
      amesosList.set("PrintTiming", true);
      smootherPrototype = rcp(new DirectSolver("", amesosList));
      }
*/

// TODO: option one level
