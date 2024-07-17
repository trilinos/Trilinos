// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_TestHelpersSmoothers.hpp>

#include <MueLu_HierarchyUtils.hpp>
#include <MueLu_SubBlockAFactory.hpp>
#include <MueLu_ReorderBlockAFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_BlockedGaussSeidelSmoother.hpp>
#include <MueLu_BlockedJacobiSmoother.hpp>
#include <MueLu_BraessSarazinSmoother.hpp>
#include <MueLu_SchurComplementFactory.hpp>
#include <MueLu_SimpleSmoother.hpp>
#include <MueLu_UzawaSmoother.hpp>
#include <MueLu_IndefBlockedDiagonalSmoother.hpp>
#include <MueLu_InverseApproximationFactory.hpp>
#include <MueLu_Utilities.hpp>

#include <Xpetra_BlockReorderManager.hpp>
#include <Xpetra_ReorderedBlockedCrsMatrix.hpp>
#include <Xpetra_ReorderedBlockedMultiVector.hpp>

namespace MueLuTests {

// this namespace already has:  #include "MueLu_UseShortNames.hpp"
using namespace TestHelpers::Smoothers;

// Xpetra version of CreateMap
template <class LocalOrdinal, class GlobalOrdinal, class Node, class MapType>
Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > CreateMap(const std::set<GlobalOrdinal>& gids, const Teuchos::Comm<int>& comm) {
  Teuchos::Array<GlobalOrdinal> mapvec;
  mapvec.reserve(gids.size());
  mapvec.assign(gids.begin(), gids.end());
  GlobalOrdinal count = Teuchos::as<GlobalOrdinal>(mapvec.size());
  GlobalOrdinal gcount;
  Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, count, Teuchos::outArg(gcount));

  Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map =
      Teuchos::rcp(new MapType(gcount,
                               mapvec(),
                               0,
                               Teuchos::rcpFromRef(comm)));
  mapvec.clear();
  return map;
}

// Xpetra version of SplitMap
template <class LocalOrdinal, class GlobalOrdinal, class Node, class MapType>
Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > SplitMap(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& Amap, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& Agiven) {
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Amap.getComm();

  GlobalOrdinal count = 0;
  Teuchos::Array<GlobalOrdinal> myaugids(Amap.getLocalNumElements());
  for (size_t i = 0; i < Amap.getLocalNumElements(); ++i) {
    const GlobalOrdinal gid = Amap.getGlobalElement(i);
    if (Agiven.isNodeGlobalElement(gid)) continue;
    myaugids[Teuchos::as<GlobalOrdinal>(count)] = gid;
    ++count;
  }
  myaugids.resize(count);
  GlobalOrdinal gcount;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &count, &gcount);
  return Teuchos::rcp(new MapType(gcount, myaugids(), 0, comm));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class MapType>
Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateBlockDiagonalExampleMatrix(int noBlocks, const Teuchos::Comm<int>& comm) {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix;
  typedef Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixFactory;
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractor;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedCrsMatrix;
  typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixWrap;

  GlobalOrdinal nOverallDOFGidsPerProc = Teuchos::as<GlobalOrdinal>(Teuchos::ScalarTraits<GlobalOrdinal>::pow(2, noBlocks - 2)) * 10;

  GlobalOrdinal procOffset = comm.getRank() * nOverallDOFGidsPerProc;

  std::set<GlobalOrdinal> myDOFGids;
  for (GlobalOrdinal i = 0; i < nOverallDOFGidsPerProc; i++)
    myDOFGids.insert(i + procOffset);

  Teuchos::RCP<Map> fullmap = CreateMap<LocalOrdinal, GlobalOrdinal, Node, MapType>(myDOFGids, comm);

  std::vector<Teuchos::RCP<const Map> > maps(noBlocks, Teuchos::null);
  GlobalOrdinal nPartGIDs            = nOverallDOFGidsPerProc;
  Teuchos::RCP<Map> remainingpartmap = fullmap;
  for (int it = 0; it < noBlocks; it++) {
    if (it == noBlocks - 1) {
      maps[0] = remainingpartmap;
      break;
    }
    // collect first half of GIDs
    nPartGIDs = nPartGIDs / 2;
    std::set<GlobalOrdinal> myHalfGIDs;
    for (GlobalOrdinal j = 0; j < nPartGIDs; j++)
      myHalfGIDs.insert(j + procOffset);

    Teuchos::RCP<Map> halfmap = CreateMap<LocalOrdinal, GlobalOrdinal, Node, MapType>(myHalfGIDs, comm);

    Teuchos::RCP<Map> secondmap = SplitMap<LocalOrdinal, GlobalOrdinal, Node, MapType>(*remainingpartmap, *halfmap);
    remainingpartmap            = halfmap;

    maps[noBlocks - 1 - it] = secondmap;
  }

  // create diagonal blocks
  std::vector<Teuchos::RCP<CrsMatrix> > blocks(noBlocks, Teuchos::null);
  for (int it = 0; it < noBlocks; it++) {
    // std::cout << it << " " << maps[it]->getMinAllGlobalIndex() << " - " << maps[it]->getMaxAllGlobalIndex() << std::endl;
    blocks[it] = CrsMatrixFactory::Build(maps[it], 1);

    LocalOrdinal NumMyElements                               = maps[it]->getLocalNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = maps[it]->getLocalElementList();

    for (LocalOrdinal i = 0; i < NumMyElements; i++)
      blocks[it]->insertGlobalValues(MyGlobalElements[i],
                                     Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                                     Teuchos::tuple<Scalar>(it + 1));
    blocks[it]->fillComplete();
  }

  // create map extractor
  Teuchos::RCP<const MapExtractor> rgMapExtractor = Teuchos::rcp(new MapExtractor(fullmap, maps, false));
  Teuchos::RCP<const MapExtractor> doMapExtractor = Teuchos::rcp(new MapExtractor(fullmap, maps, false));

  // build blocked operator
  Teuchos::RCP<BlockedCrsMatrix> bop = Teuchos::rcp(new BlockedCrsMatrix(rgMapExtractor, doMapExtractor, 1));

  for (int it = 0; it < noBlocks; it++) {
    Teuchos::RCP<CrsMatrixWrap> csrwrap =
        Teuchos::rcp(new CrsMatrixWrap(blocks[it]));
    bop->setMatrix(Teuchos::as<size_t>(it), Teuchos::as<size_t>(it), csrwrap);
  }
  bop->fillComplete();
  return bop;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class MapType>
Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateBlockDiagonalExampleMatrixThyra(int noBlocks, const Teuchos::Comm<int>& comm) {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactory;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix;
  typedef Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixFactory;
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractor;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedCrsMatrix;
  typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixWrap;

  std::vector<Teuchos::RCP<const Map> > maps(noBlocks, Teuchos::null);

  MapType tm(1, 0, Teuchos::rcpFromRef(comm));
  Xpetra::UnderlyingLib lib = tm.lib();

  maps[0] = MapFactory::Build(lib, comm.getSize() * 5, 5, 0, Teuchos::rcpFromRef(comm));
  for (int it = 1; it < noBlocks; it++) {
    GlobalOrdinal localDofs = Teuchos::as<GlobalOrdinal>(Teuchos::ScalarTraits<GlobalOrdinal>::pow(2, it - 1) * 5);
    maps[it]                = MapFactory::Build(lib, comm.getSize() * localDofs, localDofs, 0, Teuchos::rcpFromRef(comm));
  }

  // create diagonal blocks
  std::vector<Teuchos::RCP<CrsMatrix> > blocks(noBlocks, Teuchos::null);
  for (int it = 0; it < noBlocks; it++) {
    // std::cout << it << " " << maps[it]->getMinAllGlobalIndex() << " - " << maps[it]->getMaxAllGlobalIndex() << std::endl;
    blocks[it] = CrsMatrixFactory::Build(maps[it], 1);

    LocalOrdinal NumMyElements                               = maps[it]->getLocalNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = maps[it]->getLocalElementList();

    for (LocalOrdinal i = 0; i < NumMyElements; i++)
      blocks[it]->insertGlobalValues(MyGlobalElements[i],
                                     Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                                     Teuchos::tuple<Scalar>(it + 1));
    blocks[it]->fillComplete();
  }

  // create map extractor
  // To generate the Thyra style map extractor we do not need a full map but only the
  // information about the Map details (i.e. lib and indexBase). We can extract this
  // information from maps[0]
  Teuchos::RCP<const MapExtractor> rgMapExtractor =
      Teuchos::rcp(new MapExtractor(maps[0], maps, true));
  Teuchos::RCP<const MapExtractor> doMapExtractor =
      Teuchos::rcp(new MapExtractor(maps[0], maps, true));
  // build blocked operator
  Teuchos::RCP<BlockedCrsMatrix> bop = Teuchos::rcp(new BlockedCrsMatrix(rgMapExtractor, doMapExtractor, 1));

  for (int it = 0; it < noBlocks; it++) {
    Teuchos::RCP<CrsMatrixWrap> csrwrap = Teuchos::rcp(new CrsMatrixWrap(blocks[it]));
    bop->setMatrix(Teuchos::as<size_t>(it), Teuchos::as<size_t>(it), csrwrap);
  }
  bop->fillComplete();
  return bop;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateBlockedMatrixThyra(const Teuchos::Comm<int>& comm, Xpetra::UnderlyingLib lib) {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractor;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedCrsMatrix;

  std::vector<RCP<const Map> > maps = std::vector<RCP<const Map> >(3, Teuchos::null);
  maps[0]                           = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
  maps[1]                           = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
  maps[2]                           = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
  RCP<Matrix> A00                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[0], 4.0, -1.0, -1.0, lib);
  RCP<Matrix> A01                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[0], -1.0, 0.0, 0.0, lib);
  RCP<Matrix> A10                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], -1.0, 0.0, 0.0, lib);
  RCP<Matrix> A11                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], 4.0, -1.0, -1.0, lib);
  RCP<Matrix> A12                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], -1.0, 0.0, 0.0, lib);
  RCP<Matrix> A21                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[2], -1.0, 0.0, 0.0, lib);
  RCP<Matrix> A22                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[2], 4.0, -1.0, -1.0, lib);

  // create map extractor
  // To generate the Thyra style map extractor we do not need a full map but only the
  // information about the Map details (i.e. lib and indexBase). We can extract this
  // information from maps[0]
  Teuchos::RCP<const MapExtractor> rgMapExtractor =
      Teuchos::rcp(new MapExtractor(maps[0], maps, true));
  Teuchos::RCP<const MapExtractor> doMapExtractor =
      Teuchos::rcp(new MapExtractor(maps[0], maps, true));
  // build blocked operator
  Teuchos::RCP<BlockedCrsMatrix> bop = Teuchos::rcp(new BlockedCrsMatrix(rgMapExtractor, doMapExtractor, 5));
  bop->setMatrix(Teuchos::as<size_t>(0), Teuchos::as<size_t>(0), A00);
  bop->setMatrix(Teuchos::as<size_t>(0), Teuchos::as<size_t>(1), A01);
  bop->setMatrix(Teuchos::as<size_t>(1), Teuchos::as<size_t>(0), A10);
  bop->setMatrix(Teuchos::as<size_t>(1), Teuchos::as<size_t>(1), A11);
  bop->setMatrix(Teuchos::as<size_t>(1), Teuchos::as<size_t>(2), A12);
  bop->setMatrix(Teuchos::as<size_t>(2), Teuchos::as<size_t>(1), A21);
  bop->setMatrix(Teuchos::as<size_t>(2), Teuchos::as<size_t>(2), A22);
  bop->fillComplete();
  return bop;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, Jacobi_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    int noBlocks                             = 4;
    Teuchos::RCP<const BlockedCrsMatrix> bop = CreateBlockDiagonalExampleMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, TpetraMap>(noBlocks, *comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BlockedJacobiSmoother> smootherPrototype = rcp(new BlockedJacobiSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(1));

    std::vector<RCP<SubBlockAFactory> > sA(noBlocks, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(noBlocks, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(noBlocks, Teuchos::null);
    for (int k = 0; k < noBlocks; k++) {
      std::string strInfo = std::string("{ 1 }");
      sA[k]               = rcp(new SubBlockAFactory());
      sA[k]->SetFactory("A", MueLu::NoFactory::getRCP());
      sA[k]->SetParameter("block row", Teuchos::ParameterEntry(k));
      sA[k]->SetParameter("block col", Teuchos::ParameterEntry(k));
      sA[k]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      sA[k]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

      RCP<SmootherPrototype> smoProto = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      smoProto->SetFactory("A", sA[k]);
      sF[k] = rcp(new SmootherFactory(smoProto));

      sM[k] = rcp(new FactoryManager());
      sM[k]->SetFactory("A", sA[k]);
      sM[k]->SetFactory("Smoother", sF[k]);
      sM[k]->SetIgnoreUserData(true);

      smootherPrototype->AddFactoryManager(sM[k], k);
    }

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request Jacobi smoother (and all dependencies) on level
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> jacSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    A->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    jacSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();
    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    out << "solve with zero initial guess, and unreliable nonzeroed vector X" << std::endl;
    X->randomize();
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    jacSmoother->Apply(*X, *RHS, true);  // zero initial guess with nonzero X

    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm2 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm2[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm2[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    out << "solve with random initial guess" << std::endl;
    X->randomize();
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    jacSmoother->Apply(*X, *RHS, false);  // nonzero initial guess

    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm3 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm3[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm3[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    if (comm->getSize() == 1) {
      TEST_EQUALITY(residualNorm1[0] == residualNorm2[0], true);
      TEST_EQUALITY(residualNorm1[0] != residualNorm3[0], true);
    } else {
      out << "Pass/Fail is only checked in serial." << std::endl;
    }
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, BGS_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    int noBlocks                             = 4;
    Teuchos::RCP<const BlockedCrsMatrix> bop = CreateBlockDiagonalExampleMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, TpetraMap>(noBlocks, *comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BlockedGaussSeidelSmoother> smootherPrototype = rcp(new BlockedGaussSeidelSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(1));

    std::vector<RCP<SubBlockAFactory> > sA(noBlocks, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(noBlocks, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(noBlocks, Teuchos::null);
    for (int k = 0; k < noBlocks; k++) {
      std::string strInfo = std::string("{ 1 }");
      sA[k]               = rcp(new SubBlockAFactory());
      sA[k]->SetFactory("A", MueLu::NoFactory::getRCP());
      sA[k]->SetParameter("block row", Teuchos::ParameterEntry(k));
      sA[k]->SetParameter("block col", Teuchos::ParameterEntry(k));
      sA[k]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      sA[k]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

      RCP<SmootherPrototype> smoProto = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      smoProto->SetFactory("A", sA[k]);
      sF[k] = rcp(new SmootherFactory(smoProto));

      sM[k] = rcp(new FactoryManager());
      sM[k]->SetFactory("A", sA[k]);
      sM[k]->SetFactory("Smoother", sF[k]);
      sM[k]->SetIgnoreUserData(true);

      smootherPrototype->AddFactoryManager(sM[k], k);
    }

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> bgsSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    A->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    bgsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    out << "solve with zero initial guess, and unreliable nonzeroed vector X" << std::endl;
    X->randomize();
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    bgsSmoother->Apply(*X, *RHS, true);  // zero initial guess with nonzero X

    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm2 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm2[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm2[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    out << "solve with random initial guess" << std::endl;
    X->randomize();
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    bgsSmoother->Apply(*X, *RHS, false);  // nonzero initial guess

    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm3 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm3[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm3[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    if (comm->getSize() == 1) {
      TEST_EQUALITY(residualNorm1[0] == residualNorm2[0], true);
      TEST_EQUALITY(residualNorm1[0] != residualNorm3[0], true);
    } else {
      out << "Pass/Fail is only checked in serial." << std::endl;
    }
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, Reordered_BGS_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    int noBlocks                             = 4;
    Teuchos::RCP<const BlockedCrsMatrix> bop = CreateBlockDiagonalExampleMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, TpetraMap>(noBlocks, *comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 3 1 2 ]")));

    int noBlocks2 = 3;

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BlockedGaussSeidelSmoother> smootherPrototype = rcp(new BlockedGaussSeidelSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smootherPrototype->SetFactory("A", rAFact);

    std::vector<RCP<SubBlockAFactory> > sA(noBlocks2, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(noBlocks2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(noBlocks2, Teuchos::null);
    for (int k = 0; k < noBlocks2; k++) {
      std::string strInfo = std::string("{ 1 }");
      sA[k]               = rcp(new SubBlockAFactory());
      sA[k]->SetFactory("A", rAFact);
      sA[k]->SetParameter("block row", Teuchos::ParameterEntry(k));
      sA[k]->SetParameter("block col", Teuchos::ParameterEntry(k));
      sA[k]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      sA[k]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

      RCP<SmootherPrototype> smoProto = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      smoProto->SetFactory("A", sA[k]);
      sF[k] = rcp(new SmootherFactory(smoProto));

      sM[k] = rcp(new FactoryManager());
      sM[k]->SetFactory("A", sA[k]);
      sM[k]->SetFactory("Smoother", sF[k]);
      sM[k]->SetIgnoreUserData(true);

      smootherPrototype->AddFactoryManager(sM[k], k);
    }

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("A", rAFact);
    M.SetFactory("Smoother", smootherFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> bgsSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // Random X
    X->putScalar(0.0);
    RHS->putScalar(1.0);

    bgsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);

    Teuchos::RCP<const MapExtractor> doMapExtractor = reorderedbA->getDomainMapExtractor();

    Teuchos::RCP<MultiVector> v0 = doMapExtractor->ExtractVector(X, 0);
    Teuchos::RCP<MultiVector> v1 = doMapExtractor->ExtractVector(X, 1);
    Teuchos::RCP<MultiVector> v2 = doMapExtractor->ExtractVector(X, 2);

    TEST_EQUALITY((v0->getData(0))[0], Teuchos::as<Scalar>(0.25));
    TEST_EQUALITY((v1->getData(0))[0], Teuchos::as<Scalar>(0.5));
    TEST_EQUALITY((v2->getData(0))[0], Teuchos::as<Scalar>(1.0 / 3.0));

    Teuchos::Array<magnitude_type> n0(1);
    v0->norm1(n0);
    Teuchos::Array<magnitude_type> n1(1);
    v1->norm1(n1);
    Teuchos::Array<magnitude_type> n2(1);
    v2->norm1(n2);

    TEST_EQUALITY(n0[0], Teuchos::as<Scalar>(v0->getGlobalLength() * 0.25));
    TEST_EQUALITY(n1[0], Teuchos::as<Scalar>(v1->getGlobalLength() * 0.5));

    TEST_EQUALITY(X->getMap()->getLocalElement(v0->getMap()->getMinGlobalIndex()), 0);
    TEST_EQUALITY(X->getMap()->getLocalElement(v0->getMap()->getMaxGlobalIndex()), 19);
    TEST_EQUALITY(v0->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v0->getMap()->getMaxLocalIndex(), 19);
    TEST_EQUALITY(v0->getMap()->getMinAllGlobalIndex(), 20);
    TEST_EQUALITY(v0->getMap()->getMinGlobalIndex(), comm->getRank() * 40 + 20);
    TEST_EQUALITY(v0->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 39);
    TEST_EQUALITY(v0->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 1);

    TEST_EQUALITY(X->getMap()->getLocalElement(v1->getMap()->getMinGlobalIndex()), 20);
    TEST_EQUALITY(X->getMap()->getLocalElement(v1->getMap()->getMaxGlobalIndex()), 24);
    TEST_EQUALITY(v1->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v1->getMap()->getMaxLocalIndex(), 4);
    TEST_EQUALITY(v1->getMap()->getMinAllGlobalIndex(), 5);
    TEST_EQUALITY(v1->getMap()->getMinGlobalIndex(), comm->getRank() * 40 + 5);
    TEST_EQUALITY(v1->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 9);
    TEST_EQUALITY(v1->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 31);

    TEST_EQUALITY(X->getMap()->getLocalElement(v2->getMap()->getMinGlobalIndex()), 25);
    TEST_EQUALITY(X->getMap()->getLocalElement(v2->getMap()->getMaxGlobalIndex()), 34);
    TEST_EQUALITY(v2->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v2->getMap()->getMaxLocalIndex(), 9);
    TEST_EQUALITY(v2->getMap()->getMinAllGlobalIndex(), 10);
    TEST_EQUALITY(v2->getMap()->getMinGlobalIndex(), comm->getRank() * 40 + 10);
    TEST_EQUALITY(v2->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 19);
    TEST_EQUALITY(v2->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 21);

  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedII30II12II_BGS_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    int noBlocks                             = 4;
    Teuchos::RCP<const BlockedCrsMatrix> bop = CreateBlockDiagonalExampleMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, TpetraMap>(noBlocks, *comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ [ 3 0 ] [1 2] ]")));

    //////////////////////////////////////////////////////////////////////
    // global level
    RCP<BlockedGaussSeidelSmoother> smootherPrototype = rcp(new BlockedGaussSeidelSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smootherPrototype->SetFactory("A", rAFact);

    for (int k = 0; k < 2; k++) {
      Teuchos::RCP<SubBlockAFactory> sA = Teuchos::rcp(new SubBlockAFactory());
      sA->SetFactory("A", rAFact);
      sA->SetParameter("block row", Teuchos::ParameterEntry(k));  // local block indices relative to size of blocked operator
      sA->SetParameter("block col", Teuchos::ParameterEntry(k));

      Teuchos::RCP<BlockedGaussSeidelSmoother> sP = Teuchos::rcp(new BlockedGaussSeidelSmoother());
      sP->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
      sP->SetFactory("A", sA);

      for (int l = 0; l < 2; l++) {
        std::string strInfo                = std::string("{ 1 }");
        Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
        ssA->SetFactory("A", sA);
        ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
        ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
        ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
        ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
        RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
        ssP->SetFactory("A", ssA);
        Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
        Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
        ssM->SetFactory("A", ssA);
        ssM->SetFactory("Smoother", ssF);
        ssM->SetIgnoreUserData(true);
        sP->AddFactoryManager(ssM, l);
      }

      Teuchos::RCP<SmootherFactory> sF = Teuchos::rcp(new SmootherFactory(sP));
      Teuchos::RCP<FactoryManager> sM  = Teuchos::rcp(new FactoryManager());
      sM->SetFactory("A", sA);
      sM->SetFactory("Smoother", sF);
      sM->SetIgnoreUserData(true);
      smootherPrototype->AddFactoryManager(sM, k);
    }

    // create master smoother factory
    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("A", rAFact);
    M.SetFactory("Smoother", smootherFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    smootherFact->Build(level);

    RCP<SmootherBase> bgsSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // Random X
    X->putScalar(0.0);
    RHS->putScalar(1.0);

    bgsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);

    Teuchos::RCP<const MapExtractor> doMapExtractor = reorderedbA->getDomainMapExtractor();

    Teuchos::RCP<MultiVector> v0 = doMapExtractor->ExtractVector(X, 0);
    Teuchos::RCP<MultiVector> v1 = doMapExtractor->ExtractVector(X, 1);

    Teuchos::RCP<BlockedMultiVector> bv0 = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(v0);
    Teuchos::RCP<BlockedMultiVector> bv1 = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(v1);
    TEST_EQUALITY(bv0.is_null(), false);
    TEST_EQUALITY(bv1.is_null(), false);

    Teuchos::RCP<MultiVector> bv00 = bv0->getMultiVector(0, false);
    Teuchos::RCP<MultiVector> bv01 = bv0->getMultiVector(1, false);
    Teuchos::RCP<MultiVector> bv10 = bv1->getMultiVector(0, false);
    Teuchos::RCP<MultiVector> bv11 = bv1->getMultiVector(1, false);

    TEST_EQUALITY((bv00->getData(0))[0], Teuchos::as<Scalar>(0.25));
    TEST_EQUALITY((bv01->getData(0))[0], Teuchos::as<Scalar>(1.0));
    TEST_EQUALITY((bv10->getData(0))[0], Teuchos::as<Scalar>(0.5));

    Teuchos::Array<magnitude_type> n0(1);
    v0->norm1(n0);
    Teuchos::Array<magnitude_type> n1(1);
    v1->norm1(n1);

    TEST_EQUALITY(n0[0], Teuchos::as<Scalar>(comm->getSize() * 10.0));
    TEUCHOS_TEST_COMPARE(n1[0], <, comm->getSize() * 5.84, out, success);
    TEUCHOS_TEST_COMPARE(n1[0], >, comm->getSize() * 5.83, out, success);

    TEST_EQUALITY(v0->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v0->getMap()->getMaxLocalIndex(), 24);
    TEST_EQUALITY(v0->getMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(v0->getMap()->getMinGlobalIndex(), comm->getRank() * 40 + 0);
    TEST_EQUALITY(v0->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 39);
    TEST_EQUALITY(v0->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 1);

    TEST_EQUALITY(v1->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v1->getMap()->getMaxLocalIndex(), 14);
    TEST_EQUALITY(v1->getMap()->getMinAllGlobalIndex(), 5);
    TEST_EQUALITY(v1->getMap()->getMinGlobalIndex(), comm->getRank() * 40 + 5);
    TEST_EQUALITY(v1->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 19);
    TEST_EQUALITY(v1->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 21);

    Teuchos::RCP<BlockedCrsMatrix> b00 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedbA->getMatrix(0, 0));
    TEST_EQUALITY(b00->Rows(), 2);
    TEST_EQUALITY(b00->Cols(), 2);

    Teuchos::RCP<const MapExtractor> me00 = b00->getDomainMapExtractor();
    Teuchos::RCP<MultiVector> v00         = me00->ExtractVector(v0, 0);
    Teuchos::RCP<MultiVector> v01         = me00->ExtractVector(v0, 1);
    TEST_EQUALITY((v00->getData(0))[0], Teuchos::as<Scalar>(0.25));
    TEST_EQUALITY((v01->getData(0))[0], Teuchos::as<Scalar>(1.0));
    TEST_EQUALITY(v00->getLocalLength(), 20);
    TEST_EQUALITY(v01->getLocalLength(), 5);
    TEST_EQUALITY(v00->getGlobalLength(), Teuchos::as<size_t>(comm->getSize() * 20));
    TEST_EQUALITY(v01->getGlobalLength(), Teuchos::as<size_t>(comm->getSize() * 5));

    TEST_EQUALITY(v00->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v00->getMap()->getMaxLocalIndex(), 19);
    TEST_EQUALITY(v00->getMap()->getMinAllGlobalIndex(), 20);
    TEST_EQUALITY(v00->getMap()->getMinGlobalIndex(), comm->getRank() * 40 + 20);
    TEST_EQUALITY(v00->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 39);
    TEST_EQUALITY(v00->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 1);

    TEST_EQUALITY(v01->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v01->getMap()->getMaxLocalIndex(), 4);
    TEST_EQUALITY(v01->getMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(v01->getMap()->getMinGlobalIndex(), comm->getRank() * 40 + 0);
    TEST_EQUALITY(v01->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 4);
    TEST_EQUALITY(v01->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 36);

    Teuchos::RCP<BlockedCrsMatrix> b11 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedbA->getMatrix(1, 1));
    TEST_EQUALITY(b11->Rows(), 2);
    TEST_EQUALITY(b11->Cols(), 2);

    Teuchos::RCP<const MapExtractor> me11 = b11->getDomainMapExtractor();
    Teuchos::RCP<MultiVector> v10         = me11->ExtractVector(v1, 0);
    Teuchos::RCP<MultiVector> v11         = me11->ExtractVector(v1, 1);
    TEST_EQUALITY((v10->getData(0))[0], Teuchos::as<Scalar>(0.5));
    TEST_EQUALITY(v10->getLocalLength(), 5);
    TEST_EQUALITY(v11->getLocalLength(), 10);
    TEST_EQUALITY(v10->getGlobalLength(), Teuchos::as<size_t>(comm->getSize() * 5));
    TEST_EQUALITY(v11->getGlobalLength(), Teuchos::as<size_t>(comm->getSize() * 10));

    TEST_EQUALITY(v10->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v10->getMap()->getMaxLocalIndex(), 4);
    TEST_EQUALITY(v10->getMap()->getMinAllGlobalIndex(), 5);
    TEST_EQUALITY(v10->getMap()->getMinGlobalIndex(), comm->getRank() * 40 + 5);
    TEST_EQUALITY(v10->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 9);
    TEST_EQUALITY(v10->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 31);

    TEST_EQUALITY(v11->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v11->getMap()->getMaxLocalIndex(), 9);
    TEST_EQUALITY(v11->getMap()->getMinAllGlobalIndex(), 10);
    TEST_EQUALITY(v11->getMap()->getMinGlobalIndex(), comm->getRank() * 40 + 10);
    TEST_EQUALITY(v11->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 19);
    TEST_EQUALITY(v11->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 21);

    Teuchos::Array<magnitude_type> n00(1);
    v00->norm1(n00);
    Teuchos::Array<magnitude_type> n11(1);
    v01->norm1(n11);
    Teuchos::Array<magnitude_type> n22(1);
    v10->norm1(n22);
    Teuchos::Array<magnitude_type> n33(1);
    v11->norm1(n33);

    TEST_EQUALITY(n00[0], Teuchos::as<Scalar>(v00->getGlobalLength() * 0.25));
    TEST_EQUALITY(n11[0], Teuchos::as<Scalar>(v01->getGlobalLength() * 1.0));
    TEST_EQUALITY(n22[0], Teuchos::as<Scalar>(v10->getGlobalLength() * 0.5));
    TEST_FLOATING_EQUALITY(n33[0], v11->getGlobalLength() / 3.0, Teuchos::ScalarTraits<magnitude_type>::eps());

  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedII30II12II_BGS_Setup_Apply2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    int noBlocks                             = 4;
    Teuchos::RCP<const BlockedCrsMatrix> bop = CreateBlockDiagonalExampleMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, TpetraMap>(noBlocks, *comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ [ 3 0 ] [1 2] ]")));

    //////////////////////////////////////////////////////////////////////
    // global level
    RCP<BlockedGaussSeidelSmoother> smootherPrototype = rcp(new BlockedGaussSeidelSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smootherPrototype->SetFactory("A", rAFact);

    for (int k = 0; k < 2; k++) {
      Teuchos::RCP<SubBlockAFactory> sA = Teuchos::rcp(new SubBlockAFactory());
      sA->SetFactory("A", rAFact);
      sA->SetParameter("block row", Teuchos::ParameterEntry(k));  // local block indices relative to size of blocked operator
      sA->SetParameter("block col", Teuchos::ParameterEntry(k));

      Teuchos::RCP<BlockedGaussSeidelSmoother> sP = Teuchos::rcp(new BlockedGaussSeidelSmoother());
      sP->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
      sP->SetFactory("A", sA);

      for (int l = 0; l < 2; l++) {
        std::string strInfo                = std::string("{ 1 }");
        Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
        ssA->SetFactory("A", sA);
        ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
        ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
        ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
        ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
        RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
        ssP->SetFactory("A", ssA);
        Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
        Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
        ssM->SetFactory("A", ssA);
        ssM->SetFactory("Smoother", ssF);
        ssM->SetIgnoreUserData(true);
        sP->AddFactoryManager(ssM, l);
      }

      Teuchos::RCP<SmootherFactory> sF = Teuchos::rcp(new SmootherFactory(sP));
      Teuchos::RCP<FactoryManager> sM  = Teuchos::rcp(new FactoryManager());
      sM->SetFactory("A", sA);
      sM->SetFactory("Smoother", sF);
      sM->SetIgnoreUserData(true);
      smootherPrototype->AddFactoryManager(sM, k);
    }

    // create master smoother factory
    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("A", rAFact);
    M.SetFactory("Smoother", smootherFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    smootherFact->Build(level);

    RCP<SmootherBase> bgsSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA                     = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<ReorderedBlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<ReorderedBlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<BlockedMultiVector> bX   = Teuchos::rcp(new BlockedMultiVector(bop->getBlockedDomainMap(), 1, true));
    RCP<BlockedMultiVector> bRHS = Teuchos::rcp(new BlockedMultiVector(bop->getBlockedRangeMap(), 1, true));

    RCP<MultiVector> X   = bX->Merge();
    RCP<MultiVector> RHS = bRHS->Merge();

    // Random X
    X->putScalar(0.0);
    RHS->putScalar(1.0);

    bgsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);

    Teuchos::RCP<const MapExtractor> doMapExtractor = reorderedbA->getDomainMapExtractor();

    {
      RCP<BlockedMultiVector> test = Teuchos::rcp(new BlockedMultiVector(bop->getBlockedDomainMap(), X));
      bX.swap(test);

      Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ [0 1] [2 3] ]");
      test                                                = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(buildReorderedBlockedMultiVector(brm, bX));
      bX.swap(test);
    }

    Teuchos::RCP<MultiVector> v0 = doMapExtractor->ExtractVector(bX, 0);
    Teuchos::RCP<MultiVector> v1 = doMapExtractor->ExtractVector(bX, 1);

    Teuchos::RCP<BlockedMultiVector> bv0 = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(v0);
    Teuchos::RCP<BlockedMultiVector> bv1 = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(v1);
    TEST_EQUALITY(bv0.is_null(), false);
    TEST_EQUALITY(bv1.is_null(), false);

    Teuchos::RCP<MultiVector> bv00 = bv0->getMultiVector(0, false);
    Teuchos::RCP<MultiVector> bv01 = bv0->getMultiVector(1, false);
    Teuchos::RCP<MultiVector> bv10 = bv1->getMultiVector(0, false);
    Teuchos::RCP<MultiVector> bv11 = bv1->getMultiVector(1, false);

    TEST_EQUALITY((bv00->getData(0))[0], Teuchos::as<Scalar>(0.25));
    TEST_EQUALITY((bv01->getData(0))[0], Teuchos::as<Scalar>(0.25));
    TEST_EQUALITY((bv10->getData(0))[0], Teuchos::as<Scalar>(0.25));

    Teuchos::Array<magnitude_type> n0(1);
    v0->norm1(n0);
    Teuchos::Array<magnitude_type> n1(1);
    v1->norm1(n1);

    TEST_EQUALITY(n0[0], Teuchos::as<Scalar>(comm->getSize() * 2.5));
    TEUCHOS_TEST_COMPARE(n1[0], <, comm->getSize() * 13.34, out, success);
    TEUCHOS_TEST_COMPARE(n1[0], >, comm->getSize() * 13.33, out, success);

    TEST_EQUALITY(v0->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v0->getMap()->getMaxLocalIndex(), 9);
    TEST_EQUALITY(v0->getMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(v0->getMap()->getMinGlobalIndex(), comm->getRank() * 40);
    TEST_EQUALITY(v0->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 9);
    TEST_EQUALITY(v0->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 31);

    TEST_EQUALITY(v1->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v1->getMap()->getMaxLocalIndex(), 29);
    TEST_EQUALITY(v1->getMap()->getMinAllGlobalIndex(), 10);
    TEST_EQUALITY(v1->getMap()->getMinGlobalIndex(), comm->getRank() * 40 + 10);
    TEST_EQUALITY(v1->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 39);
    TEST_EQUALITY(v1->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 1);

    Teuchos::RCP<BlockedCrsMatrix> b00 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedbA->getMatrix(0, 0));
    TEST_EQUALITY(b00->Rows(), 2);
    TEST_EQUALITY(b00->Cols(), 2);

    Teuchos::RCP<const MapExtractor> me00 = b00->getDomainMapExtractor();
    Teuchos::RCP<MultiVector> v00         = me00->ExtractVector(v0, 0);
    Teuchos::RCP<MultiVector> v01         = me00->ExtractVector(v0, 1);
    TEST_EQUALITY((v00->getData(0))[0], Teuchos::as<Scalar>(0.25));
    TEST_EQUALITY((v01->getData(0))[0], Teuchos::as<Scalar>(0.25));
    TEST_EQUALITY(v00->getLocalLength(), 5);
    TEST_EQUALITY(v01->getLocalLength(), 5);
    TEST_EQUALITY(v00->getGlobalLength(), Teuchos::as<size_t>(comm->getSize() * 5));
    TEST_EQUALITY(v01->getGlobalLength(), Teuchos::as<size_t>(comm->getSize() * 5));

    TEST_EQUALITY(v00->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v00->getMap()->getMaxLocalIndex(), 4);
    TEST_EQUALITY(v00->getMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(v00->getMap()->getMinGlobalIndex(), comm->getRank() * 40);
    TEST_EQUALITY(v00->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 4);
    TEST_EQUALITY(v00->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 36);

    TEST_EQUALITY(v01->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v01->getMap()->getMaxLocalIndex(), 4);
    TEST_EQUALITY(v01->getMap()->getMinAllGlobalIndex(), 5);
    TEST_EQUALITY(v01->getMap()->getMinGlobalIndex(), comm->getRank() * 40 + 5);
    TEST_EQUALITY(v01->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 9);
    TEST_EQUALITY(v01->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 31);

    Teuchos::RCP<BlockedCrsMatrix> b11 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedbA->getMatrix(1, 1));
    TEST_EQUALITY(b11->Rows(), 2);
    TEST_EQUALITY(b11->Cols(), 2);

    Teuchos::RCP<const MapExtractor> me11 = b11->getDomainMapExtractor();
    Teuchos::RCP<MultiVector> v10         = me11->ExtractVector(v1, 0);
    Teuchos::RCP<MultiVector> v11         = me11->ExtractVector(v1, 1);
    TEST_EQUALITY((v10->getData(0))[0], Teuchos::as<Scalar>(0.25));
    TEST_EQUALITY(v10->getLocalLength(), 10);
    TEST_EQUALITY(v11->getLocalLength(), 20);
    TEST_EQUALITY(v10->getGlobalLength(), Teuchos::as<size_t>(comm->getSize() * 10));
    TEST_EQUALITY(v11->getGlobalLength(), Teuchos::as<size_t>(comm->getSize() * 20));

    TEST_EQUALITY(v10->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v10->getMap()->getMaxLocalIndex(), 9);
    TEST_EQUALITY(v10->getMap()->getMinAllGlobalIndex(), 10);
    TEST_EQUALITY(v10->getMap()->getMinGlobalIndex(), comm->getRank() * 40 + 10);
    TEST_EQUALITY(v10->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 19);
    TEST_EQUALITY(v10->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 21);

    TEST_EQUALITY(v11->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v11->getMap()->getMaxLocalIndex(), 19);
    TEST_EQUALITY(v11->getMap()->getMinAllGlobalIndex(), 20);
    TEST_EQUALITY(v11->getMap()->getMinGlobalIndex(), comm->getRank() * 40 + 20);
    TEST_EQUALITY(v11->getMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 39);
    TEST_EQUALITY(v11->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 1);

    Teuchos::Array<magnitude_type> n00(1);
    v00->norm1(n00);
    Teuchos::Array<magnitude_type> n11(1);
    v01->norm1(n11);
    Teuchos::Array<magnitude_type> n22(1);
    v10->norm1(n22);
    Teuchos::Array<magnitude_type> n33(1);
    v11->norm1(n33);

    TEST_EQUALITY(n00[0], Teuchos::as<Scalar>(v00->getGlobalLength() * 0.25));
    TEST_EQUALITY(n11[0], Teuchos::as<Scalar>(v01->getGlobalLength() * 0.25));
    TEST_EQUALITY(n22[0], Teuchos::as<Scalar>(v10->getGlobalLength() * 0.25));
    TEUCHOS_TEST_COMPARE(n33[0], <, v11->getGlobalLength() * 0.54167, out, success);
    TEUCHOS_TEST_COMPARE(n33[0], >, v11->getGlobalLength() * 0.54166, out, success);

  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, Thyra_BGS_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    RCP<BlockedCrsMatrix> A = CreateBlockedMatrixThyra<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*comm, Xpetra::UseTpetra);

    TEST_EQUALITY(A->Rows(), 3);
    TEST_EQUALITY(A->Cols(), 3);
    TEST_EQUALITY(A->getRangeMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(A->getRangeMap()->getMaxAllGlobalIndex(), comm->getSize() * 300 - 1);
    TEST_EQUALITY(A->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 100);
    TEST_EQUALITY(A->getRangeMap()->getMaxGlobalIndex(), comm->getSize() * 200 + comm->getRank() * 100 + 99);

    TEST_EQUALITY(A->getRangeMap(0, false)->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(A->getRangeMap(0, false)->getMaxAllGlobalIndex(), comm->getSize() * 100 - 1);
    TEST_EQUALITY(A->getRangeMap(0, false)->getMinGlobalIndex(), comm->getRank() * 100);
    TEST_EQUALITY(A->getRangeMap(0, false)->getMaxGlobalIndex(), comm->getRank() * 100 + 99);

    TEST_EQUALITY(A->getRangeMap(1, false)->getMinAllGlobalIndex(), comm->getSize() * 100);
    TEST_EQUALITY(A->getRangeMap(1, false)->getMaxAllGlobalIndex(), comm->getSize() * 100 * 2 - 1);
    TEST_EQUALITY(A->getRangeMap(1, false)->getMinGlobalIndex(), comm->getSize() * 100 + comm->getRank() * 100);
    TEST_EQUALITY(A->getRangeMap(1, false)->getMaxGlobalIndex(), comm->getSize() * 100 + comm->getRank() * 100 + 99);

    TEST_EQUALITY(A->getRangeMap(2, false)->getMinAllGlobalIndex(), comm->getSize() * 200);
    TEST_EQUALITY(A->getRangeMap(2, false)->getMaxAllGlobalIndex(), comm->getSize() * 300 - 1);
    TEST_EQUALITY(A->getRangeMap(2, false)->getMinGlobalIndex(), comm->getSize() * 200 + comm->getRank() * 100);
    TEST_EQUALITY(A->getRangeMap(2, false)->getMaxGlobalIndex(), comm->getSize() * 200 + comm->getRank() * 100 + 99);

    TEST_EQUALITY(A->getRangeMap(0, true)->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(A->getRangeMap(0, true)->getMaxAllGlobalIndex(), comm->getSize() * 100 - 1);
    TEST_EQUALITY(A->getRangeMap(0, true)->getMinGlobalIndex(), comm->getRank() * 100);
    TEST_EQUALITY(A->getRangeMap(0, true)->getMaxGlobalIndex(), comm->getRank() * 100 + 99);

    TEST_EQUALITY(A->getRangeMap(1, true)->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(A->getRangeMap(1, true)->getMaxAllGlobalIndex(), comm->getSize() * 100 - 1);
    TEST_EQUALITY(A->getRangeMap(1, true)->getMinGlobalIndex(), comm->getRank() * 100);
    TEST_EQUALITY(A->getRangeMap(1, true)->getMaxGlobalIndex(), comm->getRank() * 100 + 99);

    TEST_EQUALITY(A->getRangeMap(1, true)->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(A->getRangeMap(1, true)->getMaxAllGlobalIndex(), comm->getSize() * 100 - 1);
    TEST_EQUALITY(A->getRangeMap(1, true)->getMinGlobalIndex(), comm->getRank() * 100);
    TEST_EQUALITY(A->getRangeMap(1, true)->getMaxGlobalIndex(), comm->getRank() * 100 + 99);

    TEST_EQUALITY(A->getMatrix(0, 0)->getRangeMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(A->getMatrix(0, 0)->getRangeMap()->getMaxAllGlobalIndex(), comm->getSize() * 100 - 1);
    TEST_EQUALITY(A->getMatrix(0, 0)->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 100);
    TEST_EQUALITY(A->getMatrix(0, 0)->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 100 + 99);

    TEST_EQUALITY(A->getMatrix(0, 0)->getDomainMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(A->getMatrix(0, 0)->getDomainMap()->getMaxAllGlobalIndex(), comm->getSize() * 100 - 1);
    TEST_EQUALITY(A->getMatrix(0, 0)->getDomainMap()->getMinGlobalIndex(), comm->getRank() * 100);
    TEST_EQUALITY(A->getMatrix(0, 0)->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * 100 + 99);

    TEST_EQUALITY(A->getMatrix(2, 1)->getDomainMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(A->getMatrix(2, 1)->getDomainMap()->getMaxAllGlobalIndex(), comm->getSize() * 100 - 1);
    TEST_EQUALITY(A->getMatrix(2, 1)->getDomainMap()->getMinGlobalIndex(), comm->getRank() * 100);
    TEST_EQUALITY(A->getMatrix(2, 1)->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * 100 + 99);

    // build gloabl vector with one entries
    Teuchos::RCP<Vector> ones = VectorFactory::Build(A->getRangeMap(), true);
    Teuchos::RCP<Vector> res  = VectorFactory::Build(A->getRangeMap(), true);
    ones->putScalar(Teuchos::ScalarTraits<Scalar>::one());

    A->apply(*ones, *res);

    TEST_EQUALITY(res->norm1(), comm->getSize() * 100 * 2 + 6);
    TEST_EQUALITY(res->normInf(), 2);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(A));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BlockedGaussSeidelSmoother> smootherPrototype = rcp(new BlockedGaussSeidelSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(10));

    int noBlocks = Teuchos::as<int>(A->Rows());
    std::vector<RCP<SubBlockAFactory> > sA(noBlocks, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(noBlocks, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(noBlocks, Teuchos::null);
    for (int k = 0; k < noBlocks; k++) {
      std::string strInfo = std::string("{ 1 }");
      sA[k]               = rcp(new SubBlockAFactory());
      sA[k]->SetFactory("A", MueLu::NoFactory::getRCP());
      sA[k]->SetParameter("block row", Teuchos::ParameterEntry(k));
      sA[k]->SetParameter("block col", Teuchos::ParameterEntry(k));
      sA[k]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      sA[k]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

      Teuchos::ParameterList paramList;
      paramList.set("relaxation: sweeps", Teuchos::as<LocalOrdinal>(10));
      paramList.set("relaxation: damping factor", Teuchos::as<Scalar>(0.9));

      RCP<SmootherPrototype> smoProto = rcp(new Ifpack2Smoother(std::string("RELAXATION"), paramList, 0));
      smoProto->SetFactory("A", sA[k]);
      sF[k] = rcp(new SmootherFactory(smoProto));

      sM[k] = rcp(new FactoryManager());
      sM[k]->SetFactory("A", sA[k]);
      sM[k]->SetFactory("Smoother", sF[k]);
      sM[k]->SetIgnoreUserData(true);

      smootherPrototype->AddFactoryManager(sM[k], k);
    }

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> bgsSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    A->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    bgsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 5e-4, out, success);

    out << "solve with random initial guess" << std::endl;
    X->randomize();
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    bgsSmoother->Apply(*X, *RHS, false);  // nonzero initial guess

    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm2 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm2[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm2[0], <, 1e-2, out, success);

    if (comm->getSize() == 1) {
      TEST_EQUALITY(residualNorm1[0] != residualNorm2[0], true);
    } else {
      out << "Pass/Fail is only checked in serial." << std::endl;
    }
  }  // end Tpetra only
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, Thyra_Nested_BGS_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    RCP<BlockedCrsMatrix> A = CreateBlockedMatrixThyra<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*comm, Xpetra::UseTpetra);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(A));

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 0 [1 2] ]")));

    //////////////////////////////////////////////////////////////////////
    // global level
    RCP<BlockedGaussSeidelSmoother> smootherPrototype = rcp(new BlockedGaussSeidelSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(30));
    smootherPrototype->SetFactory("A", rAFact);

    std::string strInfo = std::string("{ 1 }");

    Teuchos::ParameterList paramList;
    paramList.set("relaxation: sweeps", Teuchos::as<LocalOrdinal>(30));
    paramList.set("relaxation: damping factor", Teuchos::as<Scalar>(0.9));

    Teuchos::RCP<SubBlockAFactory> sA00 = Teuchos::rcp(new SubBlockAFactory());
    sA00->SetFactory("A", rAFact);
    sA00->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA00->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA00->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA00->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> sP00 = rcp(new Ifpack2Smoother(std::string("RELAXATION"), paramList, 0));
    sP00->SetFactory("A", sA00);
    Teuchos::RCP<SmootherFactory> sF00 = Teuchos::rcp(new SmootherFactory(sP00));
    Teuchos::RCP<FactoryManager> sM00  = Teuchos::rcp(new FactoryManager());
    sM00->SetFactory("A", sA00);
    sM00->SetFactory("Smoother", sF00);
    sM00->SetIgnoreUserData(true);
    smootherPrototype->AddFactoryManager(sM00, 0);

    // 2x2 blocked operator
    Teuchos::RCP<SubBlockAFactory> sA11 = Teuchos::rcp(new SubBlockAFactory());
    sA11->SetFactory("A", rAFact);
    sA11->SetParameter("block row", Teuchos::ParameterEntry(1));
    sA11->SetParameter("block col", Teuchos::ParameterEntry(1));

    Teuchos::RCP<BlockedGaussSeidelSmoother> sP11 = Teuchos::rcp(new BlockedGaussSeidelSmoother());
    sP11->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    sP11->SetFactory("A", sA11);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", sA11);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), paramList, 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      sP11->AddFactoryManager(ssM, l);
    }

    Teuchos::RCP<SmootherFactory> sF11 = Teuchos::rcp(new SmootherFactory(sP11));
    Teuchos::RCP<FactoryManager> sM11  = Teuchos::rcp(new FactoryManager());
    sM11->SetFactory("A", sA11);
    sM11->SetFactory("Smoother", sF11);
    sM11->SetIgnoreUserData(true);
    smootherPrototype->AddFactoryManager(sM11, 1);

    // create master smoother factory
    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("A", rAFact);
    M.SetFactory("Smoother", smootherFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    smootherFact->Build(level);

    RCP<SmootherBase> bgsSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // Random X
    X->putScalar(0.0);
    RHS->putScalar(1.0);

    bgsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 10 * std::sqrt(Teuchos::ScalarTraits<Scalar>::eps());

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
  }  // end Tpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, Thyra_Nested_BGS_Setup_Apply2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    int noBlocks                             = 4;
    Teuchos::RCP<const BlockedCrsMatrix> bop = CreateBlockDiagonalExampleMatrixThyra<Scalar, LocalOrdinal, GlobalOrdinal, Node, TpetraMap>(noBlocks, *comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ [ 3 0 ] [1 2] ]")));

    //////////////////////////////////////////////////////////////////////
    // global level
    RCP<BlockedGaussSeidelSmoother> smootherPrototype = rcp(new BlockedGaussSeidelSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smootherPrototype->SetFactory("A", rAFact);

    for (int k = 0; k < 2; k++) {
      Teuchos::RCP<SubBlockAFactory> sA = Teuchos::rcp(new SubBlockAFactory());
      sA->SetFactory("A", rAFact);
      sA->SetParameter("block row", Teuchos::ParameterEntry(k));  // local block indices relative to size of blocked operator
      sA->SetParameter("block col", Teuchos::ParameterEntry(k));

      Teuchos::RCP<BlockedGaussSeidelSmoother> sP = Teuchos::rcp(new BlockedGaussSeidelSmoother());
      sP->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
      sP->SetFactory("A", sA);

      for (int l = 0; l < 2; l++) {
        std::string strInfo                = std::string("{ 1 }");
        Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
        ssA->SetFactory("A", sA);
        ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
        ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
        ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
        ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
        RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
        ssP->SetFactory("A", ssA);
        Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
        Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
        ssM->SetFactory("A", ssA);
        ssM->SetFactory("Smoother", ssF);
        ssM->SetIgnoreUserData(true);
        sP->AddFactoryManager(ssM, l);
      }

      Teuchos::RCP<SmootherFactory> sF = Teuchos::rcp(new SmootherFactory(sP));
      Teuchos::RCP<FactoryManager> sM  = Teuchos::rcp(new FactoryManager());
      sM->SetFactory("A", sA);
      sM->SetFactory("Smoother", sF);
      sM->SetIgnoreUserData(true);
      smootherPrototype->AddFactoryManager(sM, k);
    }

    // create master smoother factory
    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("A", rAFact);
    M.SetFactory("Smoother", smootherFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    smootherFact->Build(level);

    RCP<SmootherBase> bgsSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // Random X
    X->putScalar(0.0);
    RHS->putScalar(1.0);

    bgsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);

    Teuchos::RCP<const MapExtractor> doMapExtractor = reorderedbA->getDomainMapExtractor();

    // We should remove the boolean from ExtractVector. It is not necessary
#ifdef HAVE_XPETRA_DEBUG
    TEST_THROW(doMapExtractor->ExtractVector(X, 0, false), Xpetra::Exceptions::RuntimeError);
    TEST_THROW(doMapExtractor->ExtractVector(X, 1, false), Xpetra::Exceptions::RuntimeError);
#endif
    Teuchos::RCP<MultiVector> v0 = doMapExtractor->ExtractVector(X, 0, true);
    Teuchos::RCP<MultiVector> v1 = doMapExtractor->ExtractVector(X, 1, true);

    Teuchos::RCP<BlockedMultiVector> bv0 = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(v0);
    Teuchos::RCP<BlockedMultiVector> bv1 = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(v1);
    TEST_EQUALITY((bv0->getMultiVector(0, true)->getData(0))[0], Teuchos::as<Scalar>(0.25));
    TEST_EQUALITY((bv0->getMultiVector(0, true)->getData(0))[19], Teuchos::as<Scalar>(0.25));
    TEST_EQUALITY((bv1->getMultiVector(0, true)->getData(0))[0], Teuchos::as<Scalar>(0.5));
    TEST_EQUALITY((bv1->getMultiVector(0, true)->getData(0))[4], Teuchos::as<Scalar>(0.5));
    TEST_EQUALITY((bv0->getMultiVector(1, true)->getData(0))[0], Teuchos::as<Scalar>(1.0));
    TEST_EQUALITY((bv0->getMultiVector(1, true)->getData(0))[4], Teuchos::as<Scalar>(1.0));

    Teuchos::Array<magnitude_type> n0(1);
    v0->norm1(n0);
    Teuchos::Array<magnitude_type> n1(1);
    v1->norm1(n1);

    TEST_EQUALITY(n0[0], Teuchos::as<Scalar>(comm->getSize() * 10.0));
    TEUCHOS_TEST_COMPARE(n1[0], <, comm->getSize() * 5.84, out, success);
    TEUCHOS_TEST_COMPARE(n1[0], >, comm->getSize() * 5.83, out, success);

    TEST_EQUALITY(v0->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v0->getMap()->getMaxLocalIndex(), 24);
    TEST_EQUALITY(v0->getMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(v0->getMap()->getMinGlobalIndex(), comm->getRank() * 20);
    TEST_EQUALITY(v0->getMap()->getMaxGlobalIndex(), comm->getSize() * 20 + comm->getRank() * 5 + 4);
    TEST_EQUALITY(v0->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 25 - 1);

    TEST_EQUALITY(v1->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v1->getMap()->getMaxLocalIndex(), 14);
    TEST_EQUALITY(v1->getMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(v1->getMap()->getMinGlobalIndex(), comm->getRank() * 5);
    TEST_EQUALITY(v1->getMap()->getMaxGlobalIndex(), comm->getSize() * 5 + comm->getRank() * 10 + 9);
    TEST_EQUALITY(v1->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 15 - 1);

    Teuchos::RCP<BlockedCrsMatrix> b00 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedbA->getMatrix(0, 0));
    TEST_EQUALITY(b00->Rows(), 2);
    TEST_EQUALITY(b00->Cols(), 2);

    Teuchos::RCP<const MapExtractor> me00 = b00->getDomainMapExtractor();
    Teuchos::RCP<MultiVector> v00         = me00->ExtractVector(v0, 0, true);
    Teuchos::RCP<MultiVector> v01         = me00->ExtractVector(v0, 1, true);
    TEST_EQUALITY((v00->getData(0))[0], Teuchos::as<Scalar>(0.25));
    TEST_EQUALITY((v01->getData(0))[0], Teuchos::as<Scalar>(1.0));
    TEST_EQUALITY(v00->getLocalLength(), 20);
    TEST_EQUALITY(v01->getLocalLength(), 5);
    TEST_EQUALITY(v00->getGlobalLength(), Teuchos::as<size_t>(comm->getSize() * 20));
    TEST_EQUALITY(v01->getGlobalLength(), Teuchos::as<size_t>(comm->getSize() * 5));

    TEST_EQUALITY(v00->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v00->getMap()->getMaxLocalIndex(), 19);
    TEST_EQUALITY(v00->getMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(v00->getMap()->getMinGlobalIndex(), comm->getRank() * 20);
    TEST_EQUALITY(v00->getMap()->getMaxGlobalIndex(), comm->getRank() * 20 + 19);
    TEST_EQUALITY(v00->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 20 - 1);

    TEST_EQUALITY(v01->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v01->getMap()->getMaxLocalIndex(), 4);
    TEST_EQUALITY(v01->getMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(v01->getMap()->getMinGlobalIndex(), comm->getRank() * 5);
    TEST_EQUALITY(v01->getMap()->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
    TEST_EQUALITY(v01->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 5 - 1);

    Teuchos::RCP<BlockedCrsMatrix> b11 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedbA->getMatrix(1, 1));
    TEST_EQUALITY(b11->Rows(), 2);
    TEST_EQUALITY(b11->Cols(), 2);

    Teuchos::RCP<const MapExtractor> me11 = b11->getDomainMapExtractor();
    Teuchos::RCP<MultiVector> v10         = me11->ExtractVector(v1, 0, true);
    Teuchos::RCP<MultiVector> v11         = me11->ExtractVector(v1, 1, true);
    TEST_EQUALITY((v10->getData(0))[0], Teuchos::as<Scalar>(0.5));
    TEST_EQUALITY(v10->getLocalLength(), 5);
    TEST_EQUALITY(v11->getLocalLength(), 10);
    TEST_EQUALITY(v10->getGlobalLength(), Teuchos::as<size_t>(comm->getSize() * 5));
    TEST_EQUALITY(v11->getGlobalLength(), Teuchos::as<size_t>(comm->getSize() * 10));

    TEST_EQUALITY(v10->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v10->getMap()->getMaxLocalIndex(), 4);
    TEST_EQUALITY(v10->getMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(v10->getMap()->getMinGlobalIndex(), comm->getRank() * 5);
    TEST_EQUALITY(v10->getMap()->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
    TEST_EQUALITY(v10->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 5 - 1);

    TEST_EQUALITY(v11->getMap()->getMinLocalIndex(), 0);
    TEST_EQUALITY(v11->getMap()->getMaxLocalIndex(), 9);
    TEST_EQUALITY(v11->getMap()->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(v11->getMap()->getMinGlobalIndex(), comm->getRank() * 10);
    TEST_EQUALITY(v11->getMap()->getMaxGlobalIndex(), comm->getRank() * 10 + 9);
    TEST_EQUALITY(v11->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 10 - 1);

    Teuchos::Array<magnitude_type> n00(1);
    v00->norm1(n00);
    Teuchos::Array<magnitude_type> n11(1);
    v01->norm1(n11);
    Teuchos::Array<magnitude_type> n22(1);
    v10->norm1(n22);
    Teuchos::Array<magnitude_type> n33(1);
    v11->norm1(n33);

    TEST_EQUALITY(n00[0], Teuchos::as<Scalar>(v00->getGlobalLength() * 0.25));
    TEST_EQUALITY(n11[0], Teuchos::as<Scalar>(v01->getGlobalLength() * 1.0));
    TEST_EQUALITY(n22[0], Teuchos::as<Scalar>(v10->getGlobalLength() * 0.5));
    TEST_FLOATING_EQUALITY(n33[0], v11->getGlobalLength() / 3.0, Teuchos::ScalarTraits<magnitude_type>::eps());
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, SIMPLE_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    int noBlocks                             = 2;
    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, noBlocks, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<SimpleSmoother> smootherPrototype = rcp(new SimpleSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smootherPrototype->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", MueLu::NoFactory::getRCP());
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> smoProtoPredict = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoPredict->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoPredict));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->SetVelocityPredictionFactoryManager(sM[0]);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", MueLu::NoFactory::getRCP());

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", SFact);
    sF[1] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->SetSchurCompFactoryManager(sM[1]);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    A->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    out << "solve with random initial guess" << std::endl;
    X->randomize();
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, false);  // nonzero initial guess

    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm2 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm2[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm2[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    if (comm->getSize() == 1) {
      TEST_EQUALITY(residualNorm1[0] != residualNorm2[0], true);
    } else {
      out << "Pass/Fail is only checked in serial." << std::endl;
    }
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedII20I1I_SIMPLE_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ [2 0] 1]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<SimpleSmoother> smootherPrototype = rcp(new SimpleSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smootherPrototype->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    // create a 2x2 SIMPLE for the prediction eq.
    RCP<SimpleSmoother> smoProtoPredict = Teuchos::rcp(new SimpleSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoPredict->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));
    smoProtoPredict->SetFactory("A", sA[0]);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", sA[0]);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      if (l == 0)
        smoProtoPredict->SetVelocityPredictionFactoryManager(ssM);
      else
        smoProtoPredict->SetSchurCompFactoryManager(ssM);
    }

    sF[0] = rcp(new SmootherFactory(smoProtoPredict));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->SetVelocityPredictionFactoryManager(sM[0]);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", SFact);
    sF[1] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->SetSchurCompFactoryManager(sM[1]);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    TEST_EQUALITY(bX.is_null(), false);
    RCP<MultiVector> XX                   = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 10) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 10 && i < 15) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 15 && i < 20) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI2I01II_SIMPLE_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 2 [0 1]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<SimpleSmoother> smootherPrototype = rcp(new SimpleSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smootherPrototype->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->SetVelocityPredictionFactoryManager(sM[0]);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    // create a 2x2 SIMPLE for the prediction eq.
    RCP<SimpleSmoother> smoProtoPredict = Teuchos::rcp(new SimpleSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoPredict->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));
    smoProtoPredict->SetFactory("A", SFact);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      if (l == 0)
        smoProtoPredict->SetVelocityPredictionFactoryManager(ssM);
      else
        smoProtoPredict->SetSchurCompFactoryManager(ssM);
    }

    sF[1] = rcp(new SmootherFactory(smoProtoPredict));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->SetSchurCompFactoryManager(sM[1]);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess
    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    TEST_EQUALITY(bX.is_null(), false);
    RCP<MultiVector> XX = bX->Merge();

    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 10) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 10 && i < 15) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 15 && i < 20) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->describe(out, Teuchos::VERB_EXTREME);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI2I01II_SIMPLE_Setup_Apply2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 2 [0 1]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<SimpleSmoother> smootherPrototype = rcp(new SimpleSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smootherPrototype->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->SetVelocityPredictionFactoryManager(sM[0]);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    // create a 2x2 SIMPLE for the prediction eq.
    RCP<SimpleSmoother> smoProtoPredict = Teuchos::rcp(new SimpleSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoPredict->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));
    smoProtoPredict->SetFactory("A", SFact);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      if (l == 0)
        smoProtoPredict->SetVelocityPredictionFactoryManager(ssM);
      else
        smoProtoPredict->SetSchurCompFactoryManager(ssM);
    }

    sF[1] = rcp(new SmootherFactory(smoProtoPredict));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->SetSchurCompFactoryManager(sM[1]);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA                     = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<ReorderedBlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<ReorderedBlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<BlockedMultiVector> bX   = Teuchos::rcp(new BlockedMultiVector(bop->getBlockedDomainMap(), 1, true));
    RCP<BlockedMultiVector> bRHS = Teuchos::rcp(new BlockedMultiVector(bop->getBlockedRangeMap(), 1, true));

    RCP<MultiVector> X   = bX->Merge();
    RCP<MultiVector> RHS = bRHS->Merge();

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::ArrayRCP<const Scalar> xdata = X->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < X->getLocalLength(); i++) {
      if (i < 5) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 5 && i < 10) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
      if (i >= 10 && i < 20) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    {
      Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 2 [0 1]]");
      RCP<BlockedMultiVector> test                        = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(buildReorderedBlockedMultiVector(brm, bX));
      bX.swap(test);
      test = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(buildReorderedBlockedMultiVector(brm, bRHS));
      bRHS.swap(test);
    }

    bRHS->putScalar((SC)1.0);
    bX->setSeed(846930886);
    bX->randomize();

    RHS = bRHS->Merge();
    X   = bX->Merge();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    bop->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedII20I1I_Thyra_SIMPLE_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrixThyra(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ [2 0] 1]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<SimpleSmoother> smootherPrototype = rcp(new SimpleSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smootherPrototype->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    // create a 2x2 SIMPLE for the prediction eq.
    RCP<SimpleSmoother> smoProtoPredict = Teuchos::rcp(new SimpleSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoPredict->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));
    smoProtoPredict->SetFactory("A", sA[0]);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", sA[0]);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      if (l == 0)
        smoProtoPredict->SetVelocityPredictionFactoryManager(ssM);
      else
        smoProtoPredict->SetSchurCompFactoryManager(ssM);
    }

    sF[0] = rcp(new SmootherFactory(smoProtoPredict));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->SetVelocityPredictionFactoryManager(sM[0]);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", SFact);
    sF[1] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->SetSchurCompFactoryManager(sM[1]);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess
    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    TEST_EQUALITY(bX.is_null(), false);
    RCP<MultiVector> XX                   = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 10) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 10 && i < 15) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 15 && i < 20) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedbA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedbA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI2I01II_Thyra_SIMPLE_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrixThyra(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 2 [0 1]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<SimpleSmoother> smootherPrototype = rcp(new SimpleSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smootherPrototype->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->SetVelocityPredictionFactoryManager(sM[0]);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    // create a 2x2 SIMPLE for the prediction eq.
    RCP<SimpleSmoother> smoProtoPredict = Teuchos::rcp(new SimpleSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoPredict->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));
    smoProtoPredict->SetFactory("A", SFact);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      if (l == 0)
        smoProtoPredict->SetVelocityPredictionFactoryManager(ssM);
      else
        smoProtoPredict->SetSchurCompFactoryManager(ssM);
    }

    sF[1] = rcp(new SmootherFactory(smoProtoPredict));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->SetSchurCompFactoryManager(sM[1]);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess
    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    TEST_EQUALITY(bX.is_null(), false);
    RCP<MultiVector> XX                   = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 10) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 10 && i < 15) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 15 && i < 20) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedbA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedbA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedII01I2I_Thyra_SIMPLE_Setup_Apply2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlocked3x3MatrixThyra(*comm, lib);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ [0 1] 2]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<SimpleSmoother> smootherPrototype = rcp(new SimpleSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(15)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(0.9)));
    smootherPrototype->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    // create a 2x2 BGS for the prediction eq.
    RCP<BlockedGaussSeidelSmoother> smoProtoPredict = Teuchos::rcp(new BlockedGaussSeidelSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetFactory("A", sA[0]);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", sA[0]);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      smoProtoPredict->AddFactoryManager(ssM, l);
    }

    sF[0] = rcp(new SmootherFactory(smoProtoPredict));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->SetVelocityPredictionFactoryManager(sM[0]);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    Teuchos::ParameterList paramList;
    paramList.set("relaxation: sweeps", Teuchos::as<LocalOrdinal>(30));
    paramList.set("relaxation: damping factor", Teuchos::as<Scalar>(0.9));
    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), paramList, 0));
    smoProtoCorrect->SetFactory("A", SFact);
    sF[1] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->SetSchurCompFactoryManager(sM[1]);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 4e-3, out, success);
    TEUCHOS_TEST_COMPARE(residualNorm1[0], >, 2e-3, out, success);
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI0I21II_Thyra_SIMPLE_Setup_Apply3, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlocked3x3MatrixThyra(*comm, lib);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 0 [ 2 1 ]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<SimpleSmoother> smootherPrototype = rcp(new SimpleSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(15)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(0.9)));
    smootherPrototype->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> smoProtoPredict = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoPredict->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoPredict));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->SetVelocityPredictionFactoryManager(sM[0]);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    // create a 2x2 Simple for the prediction eq.
    RCP<SimpleSmoother> smoProtoCorrect = Teuchos::rcp(new SimpleSmoother());
    smoProtoCorrect->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoCorrect->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoCorrect->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));
    smoProtoCorrect->SetFactory("A", SFact);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      if (l == 0)
        smoProtoCorrect->SetVelocityPredictionFactoryManager(ssM);
      else
        smoProtoCorrect->SetSchurCompFactoryManager(ssM);
    }

    sF[1] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->SetSchurCompFactoryManager(sM[1]);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 60e-4, out, success);
    TEUCHOS_TEST_COMPARE(residualNorm1[0], >, 25e-4, out, success);
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, BS_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    int noBlocks                             = 2;
    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, noBlocks, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BraessSarazinSmoother> smootherPrototype = rcp(new BraessSarazinSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

    std::vector<RCP<SmootherFactory> > sF(1, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(1, Teuchos::null);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", MueLu::NoFactory::getRCP());

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", SFact);
    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", SFact);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    A->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    out << "solve with random initial guess" << std::endl;
    X->randomize();
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, false);  // nonzero initial guess

    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm2 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm2[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm2[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    if (comm->getSize() == 1) {
      TEST_EQUALITY(residualNorm1[0] != residualNorm2[0], true);
    } else {
      out << "Pass/Fail is only checked in serial." << std::endl;
    }
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedII20I1I_BS_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ [2 0] 1]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BraessSarazinSmoother> smootherPrototype = rcp(new BraessSarazinSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

    std::vector<RCP<SmootherFactory> > sF(1, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(1, Teuchos::null);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", SFact);
    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", SFact);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess
    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    TEST_EQUALITY(bX.is_null(), false);
    RCP<MultiVector> XX                   = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 10) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 10 && i < 15) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 15 && i < 20) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedbA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedbA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI2I10II_BS_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 2 [ 1 0 ]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BraessSarazinSmoother> smootherPrototype = rcp(new BraessSarazinSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

    std::vector<RCP<SmootherFactory> > sF(1, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(1, Teuchos::null);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    RCP<BraessSarazinSmoother> smoProtoCorrect = Teuchos::rcp(new BraessSarazinSmoother());
    smoProtoCorrect->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoCorrect->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoCorrect->SetFactory("A", SFact);

    std::string strInfo                = std::string("{ 1 }");
    Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
    ssA->SetFactory("A", SFact);
    ssA->SetParameter("block row", Teuchos::ParameterEntry(1));  // local block indices relative to size of blocked operator
    ssA->SetParameter("block col", Teuchos::ParameterEntry(1));
    ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
    RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    ssP->SetFactory("A", ssA);
    Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
    Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
    ssM->SetFactory("A", ssA);
    ssM->SetFactory("Smoother", ssF);
    ssM->SetIgnoreUserData(true);
    smoProtoCorrect->AddFactoryManager(ssM, 0);

    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", SFact);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    smootherFact->Build(level);

    RCP<SmootherBase> bsSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    bsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    TEST_EQUALITY(bX.is_null(), false);
    RCP<MultiVector> XX                   = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 10) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 10 && i < 15) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
      if (i >= 15 && i < 20) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedbA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    bsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedbA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI2I10II_BS_Setup_Apply2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 2 [ 1 0 ]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BraessSarazinSmoother> smootherPrototype = rcp(new BraessSarazinSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

    std::vector<RCP<SmootherFactory> > sF(1, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(1, Teuchos::null);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    RCP<BraessSarazinSmoother> smoProtoCorrect = Teuchos::rcp(new BraessSarazinSmoother());
    smoProtoCorrect->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoCorrect->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoCorrect->SetFactory("A", SFact);

    std::string strInfo                = std::string("{ 1 }");
    Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
    ssA->SetFactory("A", SFact);
    ssA->SetParameter("block row", Teuchos::ParameterEntry(1));  // local block indices relative to size of blocked operator
    ssA->SetParameter("block col", Teuchos::ParameterEntry(1));
    ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
    RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    ssP->SetFactory("A", ssA);
    Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
    Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
    ssM->SetFactory("A", ssA);
    ssM->SetFactory("Smoother", ssF);
    ssM->SetIgnoreUserData(true);
    smoProtoCorrect->AddFactoryManager(ssM, 0);

    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", SFact);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    smootherFact->Build(level);

    RCP<SmootherBase> bsSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA                     = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<ReorderedBlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<ReorderedBlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<BlockedMultiVector> bX   = Teuchos::rcp(new BlockedMultiVector(bop->getBlockedDomainMap(), 1, true));
    RCP<BlockedMultiVector> bRHS = Teuchos::rcp(new BlockedMultiVector(bop->getBlockedRangeMap(), 1, true));

    RCP<MultiVector> X   = bX->Merge();
    RCP<MultiVector> RHS = bRHS->Merge();

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    bsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::ArrayRCP<const Scalar> xdata = X->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < X->getLocalLength(); i++) {
      if (i < 5) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 5 && i < 10) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
      if (i >= 10 && i < 20) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    {
      Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 2 [0 1]]");
      RCP<BlockedMultiVector> test                        = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(buildReorderedBlockedMultiVector(brm, bX));
      bX.swap(test);
      test = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(buildReorderedBlockedMultiVector(brm, bRHS));
      bRHS.swap(test);
    }

    bRHS->putScalar((SC)1.0);
    bX->setSeed(846930886);
    bX->randomize();

    RHS = bRHS->Merge();
    X   = bX->Merge();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    bop->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    bsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedbA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedII02I1I_Thyra_BS_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrixThyra(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ [0 2] 1]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BraessSarazinSmoother> smootherPrototype = rcp(new BraessSarazinSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

    std::vector<RCP<SmootherFactory> > sF(1, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(1, Teuchos::null);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", SFact);
    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", SFact);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess
    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    TEST_EQUALITY(bX.is_null(), false);
    RCP<MultiVector> XX                   = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 5) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 5 && i < 15) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 15 && i < 20) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedbA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedbA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI2I01II_Thyra_BS_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrixThyra(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 2 [ 0 1 ]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BraessSarazinSmoother> smootherPrototype = rcp(new BraessSarazinSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

    std::vector<RCP<SmootherFactory> > sF(1, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(1, Teuchos::null);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    RCP<SimpleSmoother> smoProtoCorrect = Teuchos::rcp(new SimpleSmoother());
    smoProtoCorrect->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoCorrect->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoCorrect->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));
    smoProtoCorrect->SetFactory("A", SFact);

    std::string strInfo = std::string("{ 1 }");
    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      if (l == 0)
        smoProtoCorrect->SetVelocityPredictionFactoryManager(ssM);
      else
        smoProtoCorrect->SetSchurCompFactoryManager(ssM);
    }

    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", SFact);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    smootherFact->Build(level);

    RCP<SmootherBase> bsSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    bsSmoother->Apply(*X, *RHS, true);  // zero initial guess
    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    TEST_EQUALITY(bX.is_null(), false);
    RCP<MultiVector> XX                   = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 10) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 10 && i < 15) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 15 && i < 20) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedbA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    bsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedbA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI2I10II_Thyra_BS_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrixThyra(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 2 [ 1 0 ]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BraessSarazinSmoother> smootherPrototype = rcp(new BraessSarazinSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

    std::vector<RCP<SmootherFactory> > sF(1, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(1, Teuchos::null);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    RCP<BraessSarazinSmoother> smoProtoCorrect = Teuchos::rcp(new BraessSarazinSmoother());
    smoProtoCorrect->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoCorrect->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoCorrect->SetFactory("A", SFact);

    std::string strInfo                = std::string("{ 1 }");
    Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
    ssA->SetFactory("A", SFact);
    ssA->SetParameter("block row", Teuchos::ParameterEntry(1));  // local block indices relative to size of blocked operator
    ssA->SetParameter("block col", Teuchos::ParameterEntry(1));
    ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
    RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    ssP->SetFactory("A", ssA);
    Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
    Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
    ssM->SetFactory("A", ssA);
    ssM->SetFactory("Smoother", ssF);
    ssM->SetIgnoreUserData(true);
    smoProtoCorrect->AddFactoryManager(ssM, 0);

    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", SFact);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    smootherFact->Build(level);

    RCP<SmootherBase> bsSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    bsSmoother->Apply(*X, *RHS, true);  // zero initial guess
    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    TEST_EQUALITY(bX.is_null(), false);
    RCP<MultiVector> XX                   = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 10) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 10 && i < 15) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
      if (i >= 15 && i < 20) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedbA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    bsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedbA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedII01I2I_Thyra_BS_Setup_Apply2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlocked3x3MatrixThyra(*comm, lib);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ [0 1] 2]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BraessSarazinSmoother> smootherPrototype = rcp(new BraessSarazinSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(30)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.4)));

    std::vector<RCP<SmootherFactory> > sF(1, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(1, Teuchos::null);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    Teuchos::ParameterList paramList;
    paramList.set("relaxation: sweeps", Teuchos::as<LocalOrdinal>(30));
    paramList.set("relaxation: damping factor", Teuchos::as<Scalar>(0.9));
    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), paramList, 0));
    smoProtoCorrect->SetFactory("A", SFact);
    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", SFact);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> bsSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    bsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 3e-3, out, success);
    TEUCHOS_TEST_COMPARE(residualNorm1[0], >, 1e-3, out, success);
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI0I21II_Thyra_BS_Setup_Apply3, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlocked3x3MatrixThyra(*comm, lib);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 0 [ 2 1 ]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<BraessSarazinSmoother> smootherPrototype = rcp(new BraessSarazinSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(15)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(0.9)));

    std::vector<RCP<SmootherFactory> > sF(1, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(1, Teuchos::null);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    // create a 2x2 Simple for the prediction eq.
    std::string strInfo                 = std::string("{ 1 }");
    RCP<SimpleSmoother> smoProtoCorrect = Teuchos::rcp(new SimpleSmoother());
    smoProtoCorrect->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoCorrect->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoCorrect->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(false));
    smoProtoCorrect->SetFactory("A", SFact);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      if (l == 0)
        smoProtoCorrect->SetVelocityPredictionFactoryManager(ssM);
      else
        smoProtoCorrect->SetSchurCompFactoryManager(ssM);
    }

    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", SFact);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> bsSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    bsSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 44e-3, out, success);
    TEUCHOS_TEST_COMPARE(residualNorm1[0], >, 15e-3, out, success);
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, Uzawa_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    int noBlocks                             = 2;
    Teuchos::RCP<const BlockedCrsMatrix> bop = CreateBlockDiagonalExampleMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, TpetraMap>(noBlocks, *comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<UzawaSmoother> smootherPrototype = rcp(new UzawaSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(1));

    std::vector<RCP<SubBlockAFactory> > sA(noBlocks, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(noBlocks, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(noBlocks, Teuchos::null);
    for (int k = 0; k < noBlocks; k++) {
      std::string strInfo = std::string("{ 1 }");
      sA[k]               = rcp(new SubBlockAFactory());
      sA[k]->SetFactory("A", MueLu::NoFactory::getRCP());
      sA[k]->SetParameter("block row", Teuchos::ParameterEntry(k));
      sA[k]->SetParameter("block col", Teuchos::ParameterEntry(k));
      sA[k]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      sA[k]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

      RCP<SmootherPrototype> smoProto = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      smoProto->SetFactory("A", sA[k]);
      sF[k] = rcp(new SmootherFactory(smoProto));

      sM[k] = rcp(new FactoryManager());
      sM[k]->SetFactory("A", sA[k]);
      sM[k]->SetFactory("Smoother", sF[k]);
      sM[k]->SetIgnoreUserData(true);

      smootherPrototype->AddFactoryManager(sM[k], k);
    }

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> uzSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    A->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    uzSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    out << "solve with random initial guess" << std::endl;
    X->randomize();
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    uzSmoother->Apply(*X, *RHS, false);  // nonzero initial guess

    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm2 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm2[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm2[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    if (comm->getSize() == 1) {
      TEST_EQUALITY(residualNorm1[0] != residualNorm2[0], true);
    } else {
      out << "Pass/Fail is only checked in serial." << std::endl;
    }
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI2I01II_Uzawa_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 2 [0 1]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<UzawaSmoother> smootherPrototype = rcp(new UzawaSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    // create a 2x2 SIMPLE for the prediction eq.
    RCP<UzawaSmoother> smoProtoPredict = Teuchos::rcp(new UzawaSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoPredict->SetFactory("A", SFact);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      smoProtoPredict->AddFactoryManager(ssM, l);
    }

    sF[1] = rcp(new SmootherFactory(smoProtoPredict));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[1], 1);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess
    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    TEST_EQUALITY(bX.is_null(), false);
    RCP<MultiVector> XX                   = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 10) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 10 && i < 15) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 15 && i < 20) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI2I01II_Uzawa_Setup_Apply2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 2 [0 1]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<UzawaSmoother> smootherPrototype = rcp(new UzawaSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    // create a 2x2 SIMPLE for the prediction eq.
    RCP<UzawaSmoother> smoProtoPredict = Teuchos::rcp(new UzawaSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoPredict->SetFactory("A", SFact);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      smoProtoPredict->AddFactoryManager(ssM, l);
    }

    sF[1] = rcp(new SmootherFactory(smoProtoPredict));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[1], 1);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA                     = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<ReorderedBlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<ReorderedBlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<BlockedMultiVector> bX   = Teuchos::rcp(new BlockedMultiVector(bop->getBlockedDomainMap(), 1, true));
    RCP<BlockedMultiVector> bRHS = Teuchos::rcp(new BlockedMultiVector(bop->getBlockedRangeMap(), 1, true));

    RCP<MultiVector> X   = bX->Merge();
    RCP<MultiVector> RHS = bRHS->Merge();

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::ArrayRCP<const Scalar> xdata = X->getData(0);
    bool bCheck                           = true;
    X->describe(out, Teuchos::VERB_EXTREME);
    for (size_t i = 0; i < X->getLocalLength(); i++) {
      if (i < 5) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 5 && i < 10) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
      if (i >= 10 && i < 20) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    {
      Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 2 [0 1]]");
      RCP<BlockedMultiVector> test                        = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(buildReorderedBlockedMultiVector(brm, bX));
      bX.swap(test);
      test = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(buildReorderedBlockedMultiVector(brm, bRHS));
      bRHS.swap(test);
    }

    bRHS->putScalar((SC)1.0);
    bX->setSeed(846930886);
    bX->randomize();

    RHS = bRHS->Merge();
    X   = bX->Merge();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    bop->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI2I01II_Thyra_Uzawa_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrixThyra(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 2 [0 1]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<UzawaSmoother> smootherPrototype = rcp(new UzawaSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    // create a 2x2 SIMPLE for the prediction eq.
    RCP<UzawaSmoother> smoProtoPredict = Teuchos::rcp(new UzawaSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoPredict->SetFactory("A", SFact);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      smoProtoPredict->AddFactoryManager(ssM, l);
    }

    sF[1] = rcp(new SmootherFactory(smoProtoPredict));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[1], 1);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess
    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    TEST_EQUALITY(bX.is_null(), false);
    RCP<MultiVector> XX                   = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 10) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 10 && i < 15) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 15 && i < 20) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedII01I2I_Thyra_Uzawa_Setup_Apply2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlocked3x3MatrixThyra(*comm, lib);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ [0 1] 2]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<UzawaSmoother> smootherPrototype = rcp(new UzawaSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(15)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(0.9)));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    // create a 2x2 BGS for the prediction eq.
    RCP<BlockedGaussSeidelSmoother> smoProtoPredict = Teuchos::rcp(new BlockedGaussSeidelSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetFactory("A", sA[0]);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", sA[0]);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      smoProtoPredict->AddFactoryManager(ssM, l);
    }

    sF[0] = rcp(new SmootherFactory(smoProtoPredict));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    Teuchos::ParameterList paramList;
    paramList.set("relaxation: sweeps", Teuchos::as<LocalOrdinal>(30));
    paramList.set("relaxation: damping factor", Teuchos::as<Scalar>(0.9));
    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), paramList, 0));
    smoProtoCorrect->SetFactory("A", SFact);
    sF[1] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[1], 1);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 8e-3, out, success);
    TEUCHOS_TEST_COMPARE(residualNorm1[0], >, 5e-3, out, success);
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI0I21II_Thyra_Uzawa_Setup_Apply3, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlocked3x3MatrixThyra(*comm, lib);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 0 [ 2 1 ]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<UzawaSmoother> smootherPrototype = rcp(new UzawaSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(15)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(0.9)));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> smoProtoPredict = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoPredict->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoPredict));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    // create a 2x2 Simple for the prediction eq.
    RCP<UzawaSmoother> smoProtoCorrect = Teuchos::rcp(new UzawaSmoother());
    smoProtoCorrect->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoCorrect->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoCorrect->SetFactory("A", SFact);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      smoProtoCorrect->AddFactoryManager(ssM, l);
    }

    sF[1] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[1], 1);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> simpleSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    simpleSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 11e-3, out, success);
    TEUCHOS_TEST_COMPARE(residualNorm1[0], >, 5e-3, out, success);
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, Indef_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    int noBlocks                             = 2;
    Teuchos::RCP<const BlockedCrsMatrix> bop = CreateBlockDiagonalExampleMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, TpetraMap>(noBlocks, *comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<IndefBlockedDiagonalSmoother> smootherPrototype = rcp(new IndefBlockedDiagonalSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(1));

    std::vector<RCP<SubBlockAFactory> > sA(noBlocks, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(noBlocks, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(noBlocks, Teuchos::null);
    for (int k = 0; k < noBlocks; k++) {
      std::string strInfo = std::string("{ 1 }");
      sA[k]               = rcp(new SubBlockAFactory());
      sA[k]->SetFactory("A", MueLu::NoFactory::getRCP());
      sA[k]->SetParameter("block row", Teuchos::ParameterEntry(k));
      sA[k]->SetParameter("block col", Teuchos::ParameterEntry(k));
      sA[k]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      sA[k]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

      RCP<SmootherPrototype> smoProto = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      smoProto->SetFactory("A", sA[k]);
      sF[k] = rcp(new SmootherFactory(smoProto));

      sM[k] = rcp(new FactoryManager());
      sM[k]->SetFactory("A", sA[k]);
      sM[k]->SetFactory("Smoother", sF[k]);
      sM[k]->SetIgnoreUserData(true);

      smootherPrototype->AddFactoryManager(sM[k], k);
    }

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> inSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    A->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    inSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    out << "solve with random initial guess" << std::endl;
    X->randomize();
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    inSmoother->Apply(*X, *RHS, false);  // nonzero initial guess

    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm2 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm2[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm2[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);

    if (comm->getSize() == 1) {
      TEST_EQUALITY(residualNorm1[0] != residualNorm2[0], true);
    } else {
      out << "Pass/Fail is only checked in serial." << std::endl;
    }
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI2I01II_Indef_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 2 [0 1]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<IndefBlockedDiagonalSmoother> smootherPrototype = rcp(new IndefBlockedDiagonalSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    // create a 2x2 block smoother for the prediction eq.
    RCP<IndefBlockedDiagonalSmoother> smoProtoPredict = Teuchos::rcp(new IndefBlockedDiagonalSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoPredict->SetFactory("A", SFact);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      smoProtoPredict->AddFactoryManager(ssM, l);
    }

    sF[1] = rcp(new SmootherFactory(smoProtoPredict));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[1], 1);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request block smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> inSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    inSmoother->Apply(*X, *RHS, true);  // zero initial guess
    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    TEST_EQUALITY(bX.is_null(), false);
    RCP<MultiVector> XX                   = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 10) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 10 && i < 15) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 15 && i < 20) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    inSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI2I01II_Indef_Setup_Apply2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 2 [0 1]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<IndefBlockedDiagonalSmoother> smootherPrototype = rcp(new IndefBlockedDiagonalSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    // create a 2x2 block smoother for the prediction eq.
    RCP<IndefBlockedDiagonalSmoother> smoProtoPredict = Teuchos::rcp(new IndefBlockedDiagonalSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoPredict->SetFactory("A", SFact);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      smoProtoPredict->AddFactoryManager(ssM, l);
    }

    sF[1] = rcp(new SmootherFactory(smoProtoPredict));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[1], 1);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request block smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> inSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA                     = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<ReorderedBlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<ReorderedBlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<BlockedMultiVector> bX   = Teuchos::rcp(new BlockedMultiVector(bop->getBlockedDomainMap(), 1, true));
    RCP<BlockedMultiVector> bRHS = Teuchos::rcp(new BlockedMultiVector(bop->getBlockedRangeMap(), 1, true));

    RCP<MultiVector> X   = bX->Merge();
    RCP<MultiVector> RHS = bRHS->Merge();

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    inSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::ArrayRCP<const Scalar> xdata = X->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < X->getLocalLength(); i++) {
      if (i < 5) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 5 && i < 10) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
      if (i >= 10 && i < 20) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    {
      Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 2 [0 1]]");
      RCP<BlockedMultiVector> test                        = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(buildReorderedBlockedMultiVector(brm, bX));
      bX.swap(test);
      test = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(buildReorderedBlockedMultiVector(brm, bRHS));
      bRHS.swap(test);
    }

    bRHS->putScalar((SC)1.0);
    bX->setSeed(846930886);
    bX->randomize();

    RHS = bRHS->Merge();
    X   = bX->Merge();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    bop->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    inSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI2I01II_Thyra_Indef_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrixThyra(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 2 [0 1]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<IndefBlockedDiagonalSmoother> smootherPrototype = rcp(new IndefBlockedDiagonalSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(1)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoCorrect->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    // create a 2x2 SIMPLE for the prediction eq.
    RCP<IndefBlockedDiagonalSmoother> smoProtoPredict = Teuchos::rcp(new IndefBlockedDiagonalSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoPredict->SetFactory("A", SFact);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      smoProtoPredict->AddFactoryManager(ssM, l);
    }

    sF[1] = rcp(new SmootherFactory(smoProtoPredict));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[1], 1);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> inSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    inSmoother->Apply(*X, *RHS, true);  // zero initial guess
    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    TEST_EQUALITY(bX.is_null(), false);
    RCP<MultiVector> XX                   = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 10) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 10 && i < 15) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 15 && i < 20) {
        if (xdata[i] != (SC)0.5) bCheck = false;
      }
    }
    TEST_EQUALITY(bCheck, true);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    inSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    magnitude_type tol = 50. * Teuchos::ScalarTraits<Scalar>::eps();

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, tol, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, tol, out, success);
  }  // end useTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedII01I2I_Thyra_Indef_Setup_Apply2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlocked3x3MatrixThyra(*comm, lib);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ [0 1] 2]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<IndefBlockedDiagonalSmoother> smootherPrototype = rcp(new IndefBlockedDiagonalSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(15)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(0.9)));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    // create a 2x2 BGS for the prediction eq.
    RCP<BlockedGaussSeidelSmoother> smoProtoPredict = Teuchos::rcp(new BlockedGaussSeidelSmoother());
    smoProtoPredict->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoPredict->SetFactory("A", sA[0]);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", sA[0]);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      smoProtoPredict->AddFactoryManager(ssM, l);
    }

    sF[0] = rcp(new SmootherFactory(smoProtoPredict));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    Teuchos::ParameterList paramList;
    paramList.set("relaxation: sweeps", Teuchos::as<LocalOrdinal>(30));
    paramList.set("relaxation: damping factor", Teuchos::as<Scalar>(0.9));
    RCP<SmootherPrototype> smoProtoCorrect = rcp(new Ifpack2Smoother(std::string("RELAXATION"), paramList, 0));
    smoProtoCorrect->SetFactory("A", SFact);
    sF[1] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[1], 1);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> inSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    inSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 11e-3, out, success);
    TEUCHOS_TEST_COMPARE(residualNorm1[0], >, 7e-3, out, success);
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, NestedI0I21II_Thyra_Indef_Setup_Apply3, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlocked3x3MatrixThyra(*comm, lib);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 0 [ 2 1 ]]")));

    //////////////////////////////////////////////////////////////////////
    // Smoothers
    RCP<IndefBlockedDiagonalSmoother> smootherPrototype = rcp(new IndefBlockedDiagonalSmoother());
    smootherPrototype->SetFactory("A", rAFact);
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(Teuchos::as<LocalOrdinal>(15)));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(0.9)));

    std::vector<RCP<SubBlockAFactory> > sA(1, Teuchos::null);
    std::vector<RCP<SmootherFactory> > sF(2, Teuchos::null);
    std::vector<RCP<FactoryManager> > sM(2, Teuchos::null);

    // prediction
    std::string strInfo = std::string("{ 1 }");
    sA[0]               = rcp(new SubBlockAFactory());
    sA[0]->SetFactory("A", rAFact);
    sA[0]->SetParameter("block row", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("block col", Teuchos::ParameterEntry(0));
    sA[0]->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
    sA[0]->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));

    RCP<SmootherPrototype> smoProtoPredict = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
    smoProtoPredict->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoPredict));

    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[0], 0);

    // correction
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    // create a 2x2 Simple for the prediction eq.
    RCP<IndefBlockedDiagonalSmoother> smoProtoCorrect = Teuchos::rcp(new IndefBlockedDiagonalSmoother());
    smoProtoCorrect->SetParameter("Sweeps", Teuchos::ParameterEntry(1));
    smoProtoCorrect->SetParameter("Damping factor", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    smoProtoCorrect->SetFactory("A", SFact);

    for (int l = 0; l < 2; l++) {
      Teuchos::RCP<SubBlockAFactory> ssA = rcp(new SubBlockAFactory());
      ssA->SetFactory("A", SFact);
      ssA->SetParameter("block row", Teuchos::ParameterEntry(l));  // local block indices relative to size of blocked operator
      ssA->SetParameter("block col", Teuchos::ParameterEntry(l));
      ssA->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(strInfo));
      ssA->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(strInfo));
      RCP<SmootherPrototype> ssP = rcp(new Ifpack2Smoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
      ssP->SetFactory("A", ssA);
      Teuchos::RCP<SmootherFactory> ssF = Teuchos::rcp(new SmootherFactory(ssP));
      Teuchos::RCP<FactoryManager> ssM  = Teuchos::rcp(new FactoryManager());
      ssM->SetFactory("A", ssA);
      ssM->SetFactory("Smoother", ssF);
      ssM->SetIgnoreUserData(true);
      smoProtoCorrect->AddFactoryManager(ssM, l);
    }

    sF[1] = rcp(new SmootherFactory(smoProtoCorrect));

    sM[1] = rcp(new FactoryManager());
    sM[1]->SetFactory("A", SFact);
    sM[1]->SetFactory("Smoother", sF[1]);
    sM[1]->SetIgnoreUserData(true);

    smootherPrototype->AddFactoryManager(sM[1], 1);

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("A", rAFact);

    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    // smootherFact->DeclareInput(level);
    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> inSmoother = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    RCP<MultiVector> X   = MultiVectorFactory::Build(reorderedA->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(reorderedA->getRangeMap(), 1);

    // Random X
    X->setSeed(846930886);
    X->randomize();

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    // Normalize X
    Array<magnitude_type> norms(1);
    X->norm2(norms);
    X->scale(1 / norms[0]);

    // Compute RHS corresponding to X
    reorderedA->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    inSmoother->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*reorderedA, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 25e-3, out, success);
    TEUCHOS_TEST_COMPARE(residualNorm1[0], >, 9e-3, out, success);
  }  // end UseTpetra
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedSmoother, SplitReorder, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  // TODO test only Tpetra because of Ifpack2 smoother!
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const Matrix> Aconst = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build2DPoisson(comm->getSize() * 3, -1, lib);
    Teuchos::RCP<Matrix> A            = Teuchos::rcp_const_cast<Matrix>(Aconst);

    // split local parts of row map
    Teuchos::RCP<const Map> map = A->getRowMap();

    Teuchos::Array<GlobalOrdinal> myGids1;
    Teuchos::Array<GlobalOrdinal> myGids2;
    GlobalOrdinal count1 = 0;
    GlobalOrdinal count2 = 0;
    for (size_t i = 0; i < map->getLocalNumElements(); ++i) {
      const GlobalOrdinal gid = map->getGlobalElement(i);
      if (gid % 2 == 0) {
        myGids1.push_back(gid);
        count1++;
      } else {
        myGids2.push_back(gid);
        count2++;
      }
    }
    GlobalOrdinal gcount1 = 0;
    GlobalOrdinal gcount2 = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &count1, &gcount1);
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &count2, &gcount2);
    Teuchos::RCP<const Map> map1 = MapFactory::Build(lib, gcount1, myGids1(), 0, comm);
    Teuchos::RCP<const Map> map2 = MapFactory::Build(lib, gcount2, myGids2(), 0, comm);

    // I don't use the testApply infrastructure because it has no provision for an initial guess.
    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);
    level.Set("Map1", map1);
    level.Set("Map2", map2);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetFactory("Map1", MueLu::NoFactory::getRCP());
    rAFact->SetFactory("Map2", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ 0 1 ]")));

    // request BGS smoother (and all dependencies) on level
    level.Request("A", rAFact.get());

    // smootherFact->DeclareInput(level);
    rAFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 2);
    TEST_EQUALITY(reorderedbA->Cols(), 2);

    TEST_EQUALITY(reorderedbA->getRangeMapExtractor()->getFullMap()->isSameAs(*(A->getRowMap())), false);
    TEST_EQUALITY(reorderedbA->getDomainMapExtractor()->getFullMap()->isSameAs(*(A->getDomainMap())), false);
    TEUCHOS_TEST_COMPARE(std::abs(A->getFrobeniusNorm() - reorderedA->getFrobeniusNorm()), <, 1e4 * Teuchos::ScalarTraits<Scalar>::eps(), out, success);
    TEUCHOS_TEST_COMPARE(std::abs(A->getFrobeniusNorm() - reorderedbA->getFrobeniusNorm()), <, 1e4 * Teuchos::ScalarTraits<Scalar>::eps(), out, success);
  }  // end UseTpetra
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, Jacobi_Setup_Apply, SC, LO, GO, NO)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, BGS_Setup_Apply, SC, LO, GO, NO)                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, Reordered_BGS_Setup_Apply, SC, LO, GO, NO)               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedII30II12II_BGS_Setup_Apply, SC, LO, GO, NO)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedII30II12II_BGS_Setup_Apply2, SC, LO, GO, NO)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, Thyra_BGS_Setup_Apply, SC, LO, GO, NO)                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, Thyra_Nested_BGS_Setup_Apply, SC, LO, GO, NO)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, Thyra_Nested_BGS_Setup_Apply2, SC, LO, GO, NO)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, SIMPLE_Setup_Apply, SC, LO, GO, NO)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedII20I1I_SIMPLE_Setup_Apply, SC, LO, GO, NO)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI2I01II_SIMPLE_Setup_Apply, SC, LO, GO, NO)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI2I01II_SIMPLE_Setup_Apply2, SC, LO, GO, NO)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedII20I1I_Thyra_SIMPLE_Setup_Apply, SC, LO, GO, NO)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI2I01II_Thyra_SIMPLE_Setup_Apply, SC, LO, GO, NO)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedII01I2I_Thyra_SIMPLE_Setup_Apply2, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI0I21II_Thyra_SIMPLE_Setup_Apply3, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, BS_Setup_Apply, SC, LO, GO, NO)                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedII20I1I_BS_Setup_Apply, SC, LO, GO, NO)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI2I10II_BS_Setup_Apply, SC, LO, GO, NO)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI2I10II_BS_Setup_Apply2, SC, LO, GO, NO)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedII02I1I_Thyra_BS_Setup_Apply, SC, LO, GO, NO)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI2I01II_Thyra_BS_Setup_Apply, SC, LO, GO, NO)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI2I10II_Thyra_BS_Setup_Apply, SC, LO, GO, NO)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedII01I2I_Thyra_BS_Setup_Apply2, SC, LO, GO, NO)     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI0I21II_Thyra_BS_Setup_Apply3, SC, LO, GO, NO)     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, Uzawa_Setup_Apply, SC, LO, GO, NO)                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI2I01II_Uzawa_Setup_Apply, SC, LO, GO, NO)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI2I01II_Uzawa_Setup_Apply2, SC, LO, GO, NO)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI2I01II_Thyra_Uzawa_Setup_Apply, SC, LO, GO, NO)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedII01I2I_Thyra_Uzawa_Setup_Apply2, SC, LO, GO, NO)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI0I21II_Thyra_Uzawa_Setup_Apply3, SC, LO, GO, NO)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, Indef_Setup_Apply, SC, LO, GO, NO)                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI2I01II_Indef_Setup_Apply, SC, LO, GO, NO)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI2I01II_Indef_Setup_Apply2, SC, LO, GO, NO)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI2I01II_Thyra_Indef_Setup_Apply, SC, LO, GO, NO)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedII01I2I_Thyra_Indef_Setup_Apply2, SC, LO, GO, NO)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, NestedI0I21II_Thyra_Indef_Setup_Apply3, SC, LO, GO, NO)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedSmoother, SplitReorder, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
