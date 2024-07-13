// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * Navier2D_epetra.cpp
 *
 *  Created on: Mar 26, 2011
 *      Author: wiesner
 */

#include <unistd.h>
#include <iostream>
#include <fstream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Epetra
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_EpetraMap.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_StridedMap.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_BlockedRAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_BlockedGaussSeidelSmoother.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

#include "Navier2D_Helpers.h"

/*!
 *  2d Navier Stokes example (for Epetra)
 *
 *  using block matrices
 */

int main(int argc, char* argv[]) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Xpetra::EpetraNode Node;
#include "MueLu_UseShortNames.hpp"

  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace MueLuTests;
  using namespace Teuchos;

  oblackholestream blackhole;
  GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  bool success = false;
  bool verbose = true;
  try {
    // default parameters
    LO BGS_nSweeps   = 3;
    Scalar BGS_omega = 1.0;

    // Note: use --help to list available options.
    CommandLineProcessor clp(false);
    clp.setOption("BGS_sweeps", &BGS_nSweeps, "number of sweeps with BGS smoother");
    clp.setOption("BGS_omega", &BGS_omega, "scaling factor for BGS smoother");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
    RCP<FancyOStream> out      = fancyOStream(rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);
    *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

    // Timing
    Time myTime("global");
    TimeMonitor MM(myTime);

    // custom parameters
    LO maxLevels = 2;  // TODO: singular system if MaxLevels > 2?

    GO maxCoarseSize = 1;  // FIXME clp doesn't like long long int

    int globalNumDofs = 1500;  // used for the maps
    int nDofsPerNode  = 3;     // used for generating the fine level null-space

    // build strided maps
    // striding information: 2 velocity dofs and 1 pressure dof = 3 dofs per node
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(2);
    stridingInfo.push_back(1);

    /////////////////////////////////////// build strided maps
    // build strided maps:
    // xstridedfullmap: full map (velocity and pressure dof gids), continous
    // xstridedvelmap: only velocity dof gid maps (i.e. 0,1,3,4,6,7...)
    // xstridedpremap: only pressure dof gid maps (i.e. 2,5,8,...)
    Xpetra::UnderlyingLib lib             = Xpetra::UseEpetra;
    RCP<const StridedMap> xstridedfullmap = StridedMapFactory::Build(lib, globalNumDofs, 0, stridingInfo, comm, -1);
    RCP<const StridedMap> xstridedvelmap  = StridedMapFactory::Build(xstridedfullmap, 0);
    RCP<const StridedMap> xstridedpremap  = StridedMapFactory::Build(xstridedfullmap, 1);

    /////////////////////////////////////// transform Xpetra::Map objects to Epetra
    // this is needed for AztecOO
    const RCP<const Epetra_Map> fullmap = rcpFromRef(Xpetra::toEpetra(*xstridedfullmap));
    RCP<const Epetra_Map> velmap        = rcpFromRef(Xpetra::toEpetra(*xstridedvelmap));
    RCP<const Epetra_Map> premap        = rcpFromRef(Xpetra::toEpetra(*xstridedpremap));

    /////////////////////////////////////// import problem matrix and RHS from files (-> Epetra)

    // read in problem
    Epetra_CrsMatrix* ptrA    = 0;
    Epetra_Vector* ptrf       = 0;
    Epetra_MultiVector* ptrNS = 0;

    *out << "Reading matrix market file" << std::endl;
    EpetraExt::MatrixMarketFileToCrsMatrix("A_re1000_5932.txt", *fullmap, *fullmap, *fullmap, ptrA);
    EpetraExt::MatrixMarketFileToVector("b_re1000_5932.txt", *fullmap, ptrf);
    RCP<Epetra_CrsMatrix> epA    = rcp(ptrA);
    RCP<Epetra_Vector> epv       = rcp(ptrf);
    RCP<Epetra_MultiVector> epNS = rcp(ptrNS);

    /////////////////////////////////////// split system into 2x2 block system

    *out << "Split matrix into 2x2 block matrix" << std::endl;

    // split fullA into A11,..., A22
    RCP<Epetra_CrsMatrix> A11;
    RCP<Epetra_CrsMatrix> A12;
    RCP<Epetra_CrsMatrix> A21;
    RCP<Epetra_CrsMatrix> A22;

    if (SplitMatrix2x2(epA, *velmap, *premap, A11, A12, A21, A22) == false)
      *out << "Problem with splitting matrix" << std::endl;

    /////////////////////////////////////// transform Epetra objects to Xpetra (needed for MueLu)

    // build Xpetra objects from Epetra_CrsMatrix objects
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA11 = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A11));
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA12 = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A12));
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA21 = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A21));
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA22 = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A22));

    /////////////////////////////////////// generate MapExtractor object

    std::vector<RCP<const Xpetra::Map<LO, GO, Node> > > xmaps;
    xmaps.push_back(xstridedvelmap);
    xmaps.push_back(xstridedpremap);

    RCP<const Xpetra::MapExtractor<Scalar, LO, GO, Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar, LO, GO, Node>::Build(xstridedfullmap, xmaps);

    /////////////////////////////////////// build blocked transfer operator
    // using the map extractor
    RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> > bOp = rcp(new Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>(map_extractor, map_extractor, 10));
    bOp->setMatrix(0, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA11)));
    bOp->setMatrix(0, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA12)));
    bOp->setMatrix(1, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA21)));
    bOp->setMatrix(1, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA22)));

    bOp->fillComplete();

    //////////////////////////////////////////////////// create Hierarchy
    RCP<Hierarchy> H = rcp(new Hierarchy());
    H->setDefaultVerbLevel(VERB_HIGH);
    H->SetMaxCoarseSize(maxCoarseSize);

    //////////////////////////////////////////////////////// finest Level
    RCP<MueLu::Level> Finest = H->GetLevel();
    Finest->setDefaultVerbLevel(VERB_HIGH);
    Finest->Set("A", rcp_dynamic_cast<Matrix>(bOp));

    /////////////////////////////////////////////// define subblocks of A
    // make A11 block and A22 block available as variable "A" generated
    // by A11Fact and A22Fact
    RCP<SubBlockAFactory> A11Fact = rcp(new SubBlockAFactory());
    A11Fact->SetFactory("A", MueLu::NoFactory::getRCP());
    A11Fact->SetParameter("block row", ParameterEntry(0));
    A11Fact->SetParameter("block col", ParameterEntry(0));
    RCP<SubBlockAFactory> A22Fact = rcp(new SubBlockAFactory());
    A22Fact->SetFactory("A", MueLu::NoFactory::getRCP());
    A22Fact->SetParameter("block row", ParameterEntry(1));
    A22Fact->SetParameter("block col", ParameterEntry(1));

    ///////////////////////////////////////////// define smoother for A11
    // define block smoother for the first block matrix row in BlockGaussSeidel Smoother
    std::string ifpack11Type;
    ParameterList ifpack11List;
    ifpack11List.set("relaxation: sweeps", (LO)3);
    ifpack11List.set("relaxation: damping factor", (SC)0.5);
    ifpack11Type = "RELAXATION";
    ifpack11List.set("relaxation: type", "Gauss-Seidel");
    RCP<SmootherPrototype> smoProto11 = rcp(new TrilinosSmoother(ifpack11Type, ifpack11List, 0));
    smoProto11->SetFactory("A", A11Fact);
    RCP<SmootherFactory> Smoo11Fact = rcp(new SmootherFactory(smoProto11));

    ////////////////////////////////////////// prepare null space for A11
    RCP<MultiVector> nullspace11 = MultiVectorFactory::Build(xstridedvelmap, 2);  // this is a 2D standard null space

    for (int i = 0; i < nDofsPerNode - 1; ++i) {
      ArrayRCP<Scalar> nsValues = nullspace11->getDataNonConst(i);
      int numBlocks             = nsValues.size() / (nDofsPerNode - 1);
      for (int j = 0; j < numBlocks; ++j) {
        nsValues[j * (nDofsPerNode - 1) + i] = 1.0;
      }
    }

    Finest->Set("Nullspace1", nullspace11);

    ///////////////////////////////////////// define CoalesceDropFactory and Aggregation for A11
    // set up amalgamation for A11. Note: we're using a default null space factory (null)
    RCP<AmalgamationFactory> amalgFact11 = rcp(new AmalgamationFactory());
    amalgFact11->SetFactory("A", A11Fact);
    RCP<CoalesceDropFactory> dropFact11 = rcp(new CoalesceDropFactory());
    dropFact11->SetFactory("A", A11Fact);
    dropFact11->SetFactory("UnAmalgamationInfo", amalgFact11);
    RCP<UncoupledAggregationFactory> CoupledAggFact11 = rcp(new UncoupledAggregationFactory());
    CoupledAggFact11->SetFactory("Graph", dropFact11);
    CoupledAggFact11->SetMinNodesPerAggregate(5);
    CoupledAggFact11->SetMaxNeighAlreadySelected(1);
    CoupledAggFact11->SetOrdering("natural");

    ///////////////////////////////////////// define transfer ops for A11
#if 0
    // use PG-AMG
    RCP<TentativePFactory> P11tentFact = rcp(new TentativePFactory()); // check me

    RCP<PgPFactory> P11Fact = rcp(new PgPFactory());

    RCP<GenericRFactory> R11Fact = rcp(new GenericRFactory());

    RCP<NullspaceFactory> nspFact11 = rcp(new NullspaceFactory("Nullspace1"));
    nspFact11->SetFactory("Nullspace1", P11tentFact);

    RCP<CoarseMapFactory> coarseMapFact11 = rcp(new CoarseMapFactory());
    coarseMapFact11->setStridingData(stridingInfo);
    coarseMapFact11->setStridedBlockId(0);

    //////////////////////////////// define factory manager for (1,1) block
    RCP<FactoryManager> M11 = rcp(new FactoryManager());
    M11->SetFactory("A", A11Fact);
    M11->SetFactory("P", P11Fact);
    M11->SetFactory("R", R11Fact);
    M11->SetFactory("Nullspace", nspFact11);
    M11->SetFactory("Aggregates", CoupledAggFact11);
    M11->SetFactory("UnAmalgamationInfo", amalgFact11);
    M11->SetFactory("Ptent", P11tentFact);
    M11->SetFactory("Smoother", Smoo11Fact);
    M11->SetFactory("CoarseMap", coarseMapFact11);

#else

    RCP<TentativePFactory> P11Fact = rcp(new TentativePFactory());

    RCP<TransPFactory> R11Fact = rcp(new TransPFactory());

    RCP<NullspaceFactory> nspFact11 = rcp(new NullspaceFactory("Nullspace1"));
    nspFact11->SetFactory("Nullspace1", P11Fact);

    RCP<CoarseMapFactory> coarseMapFact11 = rcp(new CoarseMapFactory());
    coarseMapFact11->setStridingData(stridingInfo);
    coarseMapFact11->setStridedBlockId(0);

    //////////////////////////////// define factory manager for (1,1) block
    RCP<FactoryManager> M11 = rcp(new FactoryManager());
    M11->SetFactory("A", A11Fact);
    M11->SetFactory("P", P11Fact);
    M11->SetFactory("R", R11Fact);
    M11->SetFactory("Graph", dropFact11);
    M11->SetFactory("Aggregates", CoupledAggFact11);
    M11->SetFactory("UnAmalgamationInfo", amalgFact11);
    M11->SetFactory("Nullspace", nspFact11);
    M11->SetFactory("Ptent", P11Fact);
    M11->SetFactory("Smoother", Smoo11Fact);
    M11->SetFactory("CoarseMap", coarseMapFact11);

#if OLD
    RCP<TentativePFactory> P11Fact = rcp(new TentativePFactory(CoupledAggFact11, amalgFact11));  // check me

    RCP<TransPFactory> R11Fact = rcp(new TransPFactory(P11Fact));

    RCP<NullspaceFactory> nspFact11 = rcp(new NullspaceFactory("Nullspace1", P11Fact));

    RCP<CoarseMapFactory> coarseMapFact11 = rcp(new CoarseMapFactory(CoupledAggFact11, nspFact11));
    coarseMapFact11->setStridingData(stridingInfo);
    coarseMapFact11->setStridedBlockId(0);

    //////////////////////////////// define factory manager for (1,1) block
    RCP<FactoryManager> M11 = rcp(new FactoryManager());
    M11->SetFactory("A", A11Fact);
    M11->SetFactory("P", P11Fact);
    M11->SetFactory("R", R11Fact);
    M11->SetFactory("Nullspace", nspFact11);
    M11->SetFactory("Ptent", P11Fact);
    M11->SetFactory("Smoother", Smoo11Fact);
    M11->SetFactory("CoarseMap", coarseMapFact11);
#endif  // TODO remove this

#endif
    M11->SetIgnoreUserData(true);  // always use data from factories defined in factory manager

    ///////////////////////////////////////////// define smoother for A22
    std::string ifpack22Type;
    ParameterList ifpack22List;
    ifpack22List.set("relaxation: sweeps", (LO)1);
    ifpack22List.set("relaxation: damping factor", (SC)0.3);
    ifpack22Type = "RELAXATION";
    ifpack22List.set("relaxation: type", "Gauss-Seidel");
    RCP<SmootherPrototype> smoProto22 = rcp(new TrilinosSmoother(ifpack22Type, ifpack22List, 0));
    smoProto22->SetFactory("A", A22Fact);
    RCP<SmootherFactory> Smoo22Fact = rcp(new SmootherFactory(smoProto22));

    ////////////////////////////////////////// prepare null space for A22
    RCP<MultiVector> nullspace22 = MultiVectorFactory::Build(xstridedpremap, 1);  // this is a 2D standard null space
    ArrayRCP<Scalar> nsValues22  = nullspace22->getDataNonConst(0);
    for (int j = 0; j < nsValues22.size(); ++j) {
      nsValues22[j] = 1.0;
    }

    Finest->Set("Nullspace2", nullspace22);

    ///////////////////////////////////////// define transfer ops for A22
#if 0
    // use PGAMG
    RCP<AmalgamationFactory> amalgFact22 = rcp(new AmalgamationFactory(A22Fact));
    RCP<TentativePFactory> P22tentFact = rcp(new TentativePFactory(CoupledAggFact11, amalgFact22)); // check me (fed with A22) wrong column GIDS!!!

    RCP<SaPFactory> P22Fact = rcp(new SaPFactory(P22tentFact));

    //RCP<GenericRFactory> R22Fact = rcp(new GenericRFactory(P22Fact));
    RCP<TransPFactory> R22Fact = rcp(new TransPFactory(P22Fact));

    RCP<NullspaceFactory> nspFact22 = rcp(new NullspaceFactory("Nullspace2",P22tentFact));

    RCP<CoarseMapFactory> coarseMapFact22 = rcp(new CoarseMapFactory(CoupledAggFact11, nspFact22));
    coarseMapFact22->setStridingData(stridingInfo);
    coarseMapFact22->setStridedBlockId(1);

    //////////////////////////////// define factory manager for (2,2) block
    RCP<FactoryManager> M22 = rcp(new FactoryManager());
    M22->SetFactory("A", A22Fact);
    M22->SetFactory("P", P22Fact);
    M22->SetFactory("R", R22Fact);
    M22->SetFactory("UnAmalgamationInfo", amalgFact22);
    M22->SetFactory("Aggregates", AggFact22);
    M22->SetFactory("Nullspace", nspFact22);
    M22->SetFactory("Ptent", P22tentFact);
    M22->SetFactory("Smoother", Smoo22Fact);
    M22->SetFactory("CoarseMap", coarseMapFact22);
    M22->SetIgnoreUserData(true);               // always use data from factories defined in factory manager

#else
    // use TentativePFactory
    RCP<AmalgamationFactory> amalgFact22 = rcp(new AmalgamationFactory());
    RCP<TentativePFactory> P22Fact       = rcp(new TentativePFactory());  // check me (fed with A22) wrong column GIDS!!!

    RCP<TransPFactory> R22Fact = rcp(new TransPFactory());

    RCP<NullspaceFactory> nspFact22 = rcp(new NullspaceFactory("Nullspace2"));
    nspFact22->SetFactory("Nullspace2", P22Fact);

    RCP<CoarseMapFactory> coarseMapFact22 = rcp(new CoarseMapFactory());
    coarseMapFact22->setStridingData(stridingInfo);
    coarseMapFact22->setStridedBlockId(1);

    //////////////////////////////// define factory manager for (2,2) block
    RCP<FactoryManager> M22 = rcp(new FactoryManager());
    M22->SetFactory("A", A22Fact);
    M22->SetFactory("P", P22Fact);
    M22->SetFactory("R", R22Fact);
    M22->SetFactory("UnAmalgamationInfo", amalgFact22);
    M22->SetFactory("Aggregates", CoupledAggFact11);
    M22->SetFactory("Nullspace", nspFact22);
    M22->SetFactory("Ptent", P22Fact);
    M22->SetFactory("Smoother", Smoo22Fact);
    M22->SetFactory("CoarseMap", coarseMapFact22);
    M22->SetIgnoreUserData(true);  // always use data from factories defined in factory manager
#endif

    /////////////////////////////////////////// define blocked transfer ops
    RCP<BlockedPFactory> PFact = rcp(new BlockedPFactory());  // use row map index base from bOp
    PFact->AddFactoryManager(M11);
    PFact->AddFactoryManager(M22);

    RCP<GenericRFactory> RFact = rcp(new GenericRFactory());
    RFact->SetFactory("P", PFact);

    RCP<Factory> AcFact = rcp(new BlockedRAPFactory());
    AcFact->SetFactory("P", PFact);
    AcFact->SetFactory("R", RFact);

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BlockedGaussSeidelSmoother> smootherPrototype = rcp(new BlockedGaussSeidelSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(BGS_nSweeps));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(BGS_omega));
    smootherPrototype->AddFactoryManager(M11, 0);
    smootherPrototype->AddFactoryManager(M22, 1);
    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // Coarse grid correction
    RCP<BlockedGaussSeidelSmoother> coarseSolverPrototype = rcp(new BlockedGaussSeidelSmoother());
    coarseSolverPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(BGS_nSweeps));
    coarseSolverPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(BGS_omega));
    coarseSolverPrototype->AddFactoryManager(M11, 0);
    coarseSolverPrototype->AddFactoryManager(M22, 1);
    RCP<SmootherFactory> coarseSolverFact = rcp(new SmootherFactory(coarseSolverPrototype, null));

    // main factory manager
    FactoryManager M;
    M.SetFactory("A", AcFact);
    M.SetFactory("P", PFact);
    M.SetFactory("R", RFact);
    M.SetFactory("Smoother", smootherFact);      // TODO fix me
    M.SetFactory("PreSmoother", smootherFact);   // TODO fix me
    M.SetFactory("PostSmoother", smootherFact);  // TODO fix me
    M.SetFactory("CoarseSolver", coarseSolverFact);

    //////////////////////////////////// setup multigrid

    H->Setup(M, 0, maxLevels);

    *out << std::endl;
    *out << "print content of multigrid levels:" << std::endl;

    Finest->print(*out);

    RCP<Level> coarseLevel = H->GetLevel(1);
    coarseLevel->print(*out);

    // RCP<Level> coarseLevel2 = H->GetLevel(2);
    // coarseLevel2->print(*out);

    RCP<MultiVector> xLsg = MultiVectorFactory::Build(xstridedfullmap, 1);

    // Use AMG directly as an iterative method
    {
      xLsg->putScalar((SC)0.0);

      // Epetra_Vector -> Xpetra::Vector
      RCP<Vector> xRhs = rcp(new Xpetra::EpetraVectorT<int, Node>(epv));

      // calculate initial (absolute) residual
      Array<ScalarTraits<SC>::magnitudeType> norms(1);
      xRhs->norm2(norms);
      *out << "||x_0|| = " << norms[0] << std::endl;

      // apply ten multigrid iterations
      H->Iterate(*xRhs, *xLsg, 100);

      // calculate and print residual
      RCP<MultiVector> xTmp = MultiVectorFactory::Build(xstridedfullmap, 1);
      bOp->apply(*xLsg, *xTmp, NO_TRANS, (SC)1.0, (SC)0.0);
      xRhs->update((SC)-1.0, *xTmp, (SC)1.0);
      xRhs->norm2(norms);
      *out << "||r|| = " << norms[0] << std::endl;
    }

    // TODO: don't forget to add Aztec as prerequisite in CMakeLists.txt!
    //
    // Solve Ax = b using AMG as a preconditioner in AztecOO
    //
    {
      RCP<Epetra_Vector> X = rcp(new Epetra_Vector(epv->Map()));
      X->PutScalar(0.0);
      Epetra_LinearProblem epetraProblem(epA.get(), X.get(), epv.get());

      AztecOO aztecSolver(epetraProblem);
      aztecSolver.SetAztecOption(AZ_solver, AZ_gmres);

      MueLu::EpetraOperator aztecPrec(H);
      aztecSolver.SetPrecOperator(&aztecPrec);

      int maxIts = 50;
      double tol = 1e-8;

      aztecSolver.Iterate(maxIts, tol);
    }

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
#else
  std::cout << "Epetra (and/or EpetraExt) are not available. Skip test." << std::endl;
  return EXIT_SUCCESS;
#endif  // #if defined(HAVE_MUELU_SERIAL) && defined(HAVE_MUELU_EPETRA)
}
