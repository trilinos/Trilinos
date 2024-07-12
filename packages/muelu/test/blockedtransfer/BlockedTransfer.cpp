// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * BlockedTransfer.cpp
 *
 *  Created on: 01.01.2012
 *      Author: tobias
 */

#include <unistd.h>
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_BlockedRAPFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_BlockedGaussSeidelSmoother.hpp"

// Belos
#ifdef HAVE_MUELU_BELOS
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosXpetraAdapter.hpp"  // this header defines Belos::XpetraOp()
#include "BelosMueLuAdapter.hpp"   // this header defines Belos::MueLuOp()
#endif

/////////////////////////
// helper function

namespace MueLuTests {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > GenerateProblemMatrix(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map, Scalar a = 2.0, Scalar b = -1.0, Scalar c = -1.0) {
#include "MueLu_UseShortNames.hpp"
  Teuchos::RCP<CrsMatrixWrap> mtx = Galeri::Xpetra::MatrixTraits<Map, CrsMatrixWrap>::Build(map, 3);

  LocalOrdinal NumMyElements                               = map->getLocalNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();
  GlobalOrdinal NumGlobalElements                          = map->getGlobalNumElements();
  GlobalOrdinal nIndexBase                                 = map->getIndexBase();

  GlobalOrdinal NumEntries;
  LocalOrdinal nnz = 2;
  std::vector<Scalar> Values(nnz);
  std::vector<GlobalOrdinal> Indices(nnz);

  for (LocalOrdinal i = 0; i < NumMyElements; ++i) {
    if (MyGlobalElements[i] == nIndexBase) {
      // off-diagonal for first row
      Indices[0] = nIndexBase;
      NumEntries = 1;
      Values[0]  = c;
    } else if (MyGlobalElements[i] == nIndexBase + NumGlobalElements - 1) {
      // off-diagonal for last row
      Indices[0] = nIndexBase + NumGlobalElements - 2;
      NumEntries = 1;
      Values[0]  = b;
    } else {
      // off-diagonal for internal row
      Indices[0] = MyGlobalElements[i] - 1;
      Values[1]  = b;
      Indices[1] = MyGlobalElements[i] + 1;
      Values[0]  = c;
      NumEntries = 2;
    }

    // put the off-diagonal entries
    // Xpetra wants ArrayViews (sigh)
    Teuchos::ArrayView<Scalar> av(&Values[0], NumEntries);
    Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], NumEntries);
    mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

    // Put in the diagonal entry
    mtx->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                            Teuchos::tuple<Scalar>(a));

  }  // for (LocalOrdinal i = 0; i < NumMyElements; ++i)

  mtx->fillComplete(map, map);

  return mtx;
}

}  // namespace MueLuTests

/////////////////
// MAIN

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    RCP<Teuchos::FancyOStream> out      = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);
    *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

    //#ifndef HAVE_XPETRA_INT_LONG_LONG
    *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
    //#endif

    /**********************************************************************************/
    /* SET TEST PARAMETERS                                                            */
    /**********************************************************************************/
    Xpetra::Parameters xpetraParameters(clp);  // manage parameters of xpetra

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    xpetraParameters.check();

    if (comm->getRank() == 0) {
      std::cout << xpetraParameters;
      // TODO: print custom parameters // Or use paramList::print()!
    }

    /**********************************************************************************/
    /* CREATE INITIAL MATRIX                                                          */
    /**********************************************************************************/
    RCP<const Map> bigMap;
    RCP<const Map> map1;
    RCP<const Map> map2;
    GO numElements  = 500;
    GO numElements1 = 400;
    GO numElements2 = 100;

    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(1);
    map1 = StridedMapFactory::Build(lib, numElements1, 0, stridingInfo, comm, -1);
    map2 = StridedMapFactory::Build(lib, numElements2, numElements1, stridingInfo, comm, -1);

    std::vector<GlobalOrdinal> localGids;                                               // vector with all local GIDs on cur proc
    Teuchos::ArrayView<const GlobalOrdinal> map1eleList = map1->getLocalElementList();  // append all local gids from map1 and map2
    localGids.insert(localGids.end(), map1eleList.begin(), map1eleList.end());
    Teuchos::ArrayView<const GlobalOrdinal> map2eleList = map2->getLocalElementList();
    localGids.insert(localGids.end(), map2eleList.begin(), map2eleList.end());
    Teuchos::ArrayView<GlobalOrdinal> eleList(&localGids[0], localGids.size());
    bigMap = StridedMapFactory::Build(lib, numElements, eleList, 0, stridingInfo, comm);  // create full big map (concatenation of map1 and map2)

    std::vector<Teuchos::RCP<const Map> > maps;
    maps.push_back(map1);
    maps.push_back(map2);

    Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO, Node> > mapExtractor = Xpetra::MapExtractorFactory<Scalar, LO, GO, Node>::Build(bigMap, maps);

    RCP<CrsMatrixWrap> Op11 = MueLuTests::GenerateProblemMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map1, 2, -1, -1);
    RCP<CrsMatrixWrap> Op22 = MueLuTests::GenerateProblemMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map2, 3, -2, -1);

    /*Op11->describe(*out,Teuchos::VERB_EXTREME);
      Op22->describe(*out,Teuchos::VERB_EXTREME);*/

    // build blocked operator
    Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>(mapExtractor, mapExtractor, 10));

    bOp->setMatrix(0, 0, Op11);
    bOp->setMatrix(1, 1, Op22);

    bOp->fillComplete();

    // build hierarchy
    Hierarchy H;
    H.SetMaxCoarseSize(50);
    RCP<Level> levelOne = H.GetLevel();
    levelOne->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp));  // set blocked operator

    RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory());
    A11Fact->SetFactory("A", MueLu::NoFactory::getRCP());
    A11Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
    A11Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
    RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory());
    A22Fact->SetFactory("A", MueLu::NoFactory::getRCP());
    A22Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
    A22Fact->SetParameter("block col", Teuchos::ParameterEntry(1));

    RCP<TentativePFactory> P11Fact = rcp(new TentativePFactory());
    RCP<TransPFactory> R11Fact     = rcp(new TransPFactory());

    RCP<TentativePFactory> P22TentFact = rcp(new TentativePFactory());
    RCP<PgPFactory> P22Fact            = rcp(new PgPFactory());
    RCP<GenericRFactory> R22Fact       = rcp(new GenericRFactory());

    std::string ifpackType;
    Teuchos::ParameterList ifpackList;
    ifpackList.set("relaxation: sweeps", (LO)5);
    ifpackList.set("relaxation: damping factor", Teuchos::ScalarTraits<SC>::one());
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
    RCP<SmootherPrototype> smoProto11 = rcp(new TrilinosSmoother(ifpackType, ifpackList, 0));
    smoProto11->SetFactory("A", A11Fact);
    RCP<SmootherPrototype> smoProto22 = rcp(new TrilinosSmoother(ifpackType, ifpackList, 0));
    smoProto22->SetFactory("A", A22Fact);

    // RCP<SmootherPrototype> smoProto11     = rcp( new DirectSolver("", Teuchos::ParameterList(), A11Fact) );
    // RCP<SmootherPrototype> smoProto22     = rcp( new DirectSolver("", Teuchos::ParameterList(), A22Fact) );

    RCP<SmootherFactory> Smoo11Fact = rcp(new SmootherFactory(smoProto11));
    RCP<SmootherFactory> Smoo22Fact = rcp(new SmootherFactory(smoProto22));

    RCP<FactoryManager> M11 = rcp(new FactoryManager());
    M11->SetKokkosRefactor(false);
    M11->SetFactory("A", A11Fact);
    M11->SetFactory("P", P11Fact);
    M11->SetFactory("Ptent", P11Fact);  // for Nullspace
    M11->SetFactory("R", R11Fact);
    M11->SetFactory("Smoother", Smoo11Fact);
    M11->SetIgnoreUserData(true);

    RCP<FactoryManager> M22 = rcp(new FactoryManager());
    M22->SetKokkosRefactor(false);
    M22->SetFactory("A", A22Fact);
    M22->SetFactory("P", P22Fact);
    M22->SetFactory("R", R22Fact);
    M22->SetFactory("Ptent", P22TentFact);  // for both P22 and Nullspace
    M22->SetFactory("Smoother", Smoo22Fact);
    M22->SetIgnoreUserData(true);

    RCP<BlockedPFactory> PFact = rcp(new BlockedPFactory());
    PFact->AddFactoryManager(M11);
    PFact->AddFactoryManager(M22);

    RCP<GenericRFactory> RFact = rcp(new GenericRFactory());

    RCP<Factory> AcFact = rcp(new BlockedRAPFactory());

    // Smoothers
    RCP<BlockedGaussSeidelSmoother> smootherPrototype = rcp(new BlockedGaussSeidelSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(2));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::ScalarTraits<Scalar>::one()));
    smootherPrototype->AddFactoryManager(M11, 0);
    smootherPrototype->AddFactoryManager(M22, 1);
    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // Coarse grid correction
    RCP<BlockedGaussSeidelSmoother> coarseSolverPrototype = rcp(new BlockedGaussSeidelSmoother());
    coarseSolverPrototype->AddFactoryManager(M11, 0);
    coarseSolverPrototype->AddFactoryManager(M22, 1);
    RCP<SmootherFactory> coarseSolverFact = rcp(new SmootherFactory(coarseSolverPrototype, Teuchos::null));

    // main factory manager
    FactoryManager M;
    M.SetFactory("A", AcFact);
    M.SetFactory("P", PFact);
    M.SetFactory("R", RFact);
    M.SetFactory("Smoother", smootherFact);  // TODO fix me
    M.SetFactory("CoarseSolver", coarseSolverFact);

    H.SetVerbLevel(MueLu::Test);
    H.Setup(M);

    RCP<Level> l0 = H.GetLevel(0);
    RCP<Level> l1 = H.GetLevel(1);
    RCP<Level> l2 = H.GetLevel(2);

    l0->print(*out, Teuchos::VERB_EXTREME);
    l1->print(*out, Teuchos::VERB_EXTREME);
    l2->print(*out, Teuchos::VERB_EXTREME);

    // Define B
    RCP<Vector> X = VectorFactory::Build(bigMap, 1);
    RCP<Vector> B = VectorFactory::Build(bigMap, 1);
    X->setSeed(846930886);
    X->randomize();
    bOp->apply(*X, *B, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // X = 0
    X->putScalar((SC)0.0);

    LO nIts = 9;
    H.Iterate(*B, *X, nIts);

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
