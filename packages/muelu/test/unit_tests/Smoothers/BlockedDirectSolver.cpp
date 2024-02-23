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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_TestHelpersSmoothers.hpp>

#include <MueLu_HierarchyUtils.hpp>
#include <MueLu_SubBlockAFactory.hpp>
#include <MueLu_ReorderBlockAFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_SchurComplementFactory.hpp>
#include <MueLu_InverseApproximationFactory.hpp>
#include <MueLu_SimpleSmoother.hpp>
#include <MueLu_BlockedDirectSolver.hpp>
#include <MueLu_Utilities.hpp>

namespace MueLuTests {

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

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedDirectSolver, BlockedDirectSolver_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  std::array<std::string, 3> solver_types{"", "Klu", "Superlu"};

  for (const auto& solver_type : solver_types) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 5, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);
    level.setlib(lib);

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

    // Smoothers
    RCP<BlockedDirectSolver> smootherPrototype = rcp(new BlockedDirectSolver(solver_type));
    RCP<SmootherFactory> smootherFact          = rcp(new SmootherFactory(smootherPrototype));

    // main factory manager
    FactoryManager M;
    M.SetFactory("Smoother", smootherFact);
    MueLu::SetFactoryManager SFM(Teuchos::rcpFromRef(level), Teuchos::rcpFromRef(M));

    // request smoother (and all dependencies) on level
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

    smootherFact->Build(level);

    RCP<SmootherBase> solver = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

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

    solver->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 5e-15, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <,
                         5e-15, out, success);

    out << "solve with random initial guess" << std::endl;
    X->randomize();
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    solver->Apply(*X, *RHS, false);  // nonzero initial guess

    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm2 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm2[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm2[0], <, 5e-15, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <,
                         5e-15, out, success);

    if (comm->getSize() == 1) {
      TEST_EQUALITY(residualNorm1[0] == residualNorm2[0], true);
    } else {
      out << "Pass/Fail is only checked in serial." << std::endl;
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedDirectSolver, NestedI53I42II01II_BlockedDirectSolver_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  std::array<std::string, 3> solver_types{"", "Klu", "Superlu"};

  for (const auto& solver_type : solver_types) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 6, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);
    level.setlib(lib);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[5 3 [4 2] [0 1]]")));

    // Smoothers
    RCP<BlockedDirectSolver> smootherPrototype = rcp(new BlockedDirectSolver(solver_type));
    RCP<SmootherFactory> smootherFact          = rcp(new SmootherFactory(smootherPrototype));

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

    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> solver = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 4);
    TEST_EQUALITY(reorderedbA->Cols(), 4);

    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    solver->Apply(*X, *RHS, true);  // zero initial guess
    Teuchos::RCP<BlockedMultiVector> bX   = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    Teuchos::RCP<MultiVector> XX          = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 5) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 5 && i < 10) {
        if (xdata[i] != (SC)(1.0 / 2.0)) bCheck = false;
      }
      if (i >= 10 && i < 20) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 20 && i < 40) {
        if (xdata[i] != (SC)(1.0 / 4.0)) bCheck = false;
      }
      if (i >= 40 && i < 80) {
        if (xdata[i] != (SC)(1.0 / 5.0)) bCheck = false;
      }
      if (i >= 80 && i < 160) {
        if (xdata[i] != (SC)(1.0 / 6.0)) bCheck = false;
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
    A->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    solver->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 5e-15, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, 5e-15,
                         out, success);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedDirectSolver, NestedI53I42II01II_Thyra_BlockedDirectSolver_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  std::array<std::string, 3> solver_types{"", "Klu", "Superlu"};

  for (const auto& solve_type : solver_types) {
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrixThyra(lib, 6, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);
    level.setlib(lib);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[5 3 [4 2] [0 1]]")));

    // Smoothers
    RCP<BlockedDirectSolver> smootherPrototype = rcp(new BlockedDirectSolver(solve_type));
    RCP<SmootherFactory> smootherFact          = rcp(new SmootherFactory(smootherPrototype));

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

    smootherFact->Build(level);

    level.print(std::cout, Teuchos::VERB_EXTREME);

    RCP<SmootherBase> solver = level.Get<RCP<SmootherBase> >("PreSmoother", smootherFact.get());

    RCP<Matrix> reorderedA            = level.Get<RCP<Matrix> >("A", rAFact.get());
    RCP<BlockedCrsMatrix> reorderedbA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(reorderedA);

    TEST_EQUALITY(reorderedbA->Rows(), 4);
    TEST_EQUALITY(reorderedbA->Cols(), 4);

    RCP<MultiVector> X   = MultiVectorFactory::Build(A->getDomainMap(), 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRangeMap(), 1);

    // apply simple smoother
    RHS->putScalar((SC)1.0);
    X->putScalar((SC)0.0);

    // solve system
    solver->Apply(*X, *RHS, true);  // zero initial guess
    Teuchos::RCP<BlockedMultiVector> bX   = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    Teuchos::RCP<MultiVector> XX          = bX->Merge();
    Teuchos::ArrayRCP<const Scalar> xdata = XX->getData(0);
    bool bCheck                           = true;
    for (size_t i = 0; i < XX->getLocalLength(); i++) {
      if (i < 5) {
        if (xdata[i] != (SC)1.0) bCheck = false;
      }
      if (i >= 5 && i < 10) {
        if (xdata[i] != (SC)(1.0 / 2.0)) bCheck = false;
      }
      if (i >= 10 && i < 20) {
        if (xdata[i] != (SC)(1.0 / 3.0)) bCheck = false;
      }
      if (i >= 20 && i < 40) {
        if (xdata[i] != (SC)(1.0 / 4.0)) bCheck = false;
      }
      if (i >= 40 && i < 80) {
        if (xdata[i] != (SC)(1.0 / 5.0)) bCheck = false;
      }
      if (i >= 80 && i < 160) {
        if (xdata[i] != (SC)(1.0 / 6.0)) bCheck = false;
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
    A->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Reset X to 0
    X->putScalar((SC)0.0);

    RHS->norm2(norms);
    out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    out << "solve with zero initial guess" << std::endl;
    Teuchos::Array<magnitude_type> initialNorms(1);
    X->norm2(initialNorms);
    out << "  ||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

    solver->Apply(*X, *RHS, true);  // zero initial guess

    Teuchos::Array<magnitude_type> finalNorms(1);
    X->norm2(finalNorms);
    Teuchos::Array<magnitude_type> residualNorm1 = Utilities::ResidualNorm(*A, *X, *RHS);
    out << "  ||Residual_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorm1[0] << std::endl;
    out << "  ||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << finalNorms[0] << std::endl;

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 5e-15, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <, 5e-15,
                         out, success);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedDirectSolver, NestedII20I1I_BlockedDirectSolver_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  std::array<std::string, 3> solver_types{"", "Klu", "Superlu"};

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  if (comm->getSize() > 1) {
    out << "Skipping test for " << comm->getSize()
        << " processors as Amesos2 cannot deal with non-standard maps in the non-serial case." << std::endl;
    return;
  }

  for (const auto& solver_type : solver_types) {
    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);
    level.setlib(lib);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ [2 0] 1]")));

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
    RCP<BlockedDirectSolver> smoProtoPredict = Teuchos::rcp(new BlockedDirectSolver(solver_type));
    smoProtoPredict->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoPredict));
    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->SetVelocityPredictionFactoryManager(sM[0]);

    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new TrilinosSmoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
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

    // request smoother (and all dependencies) on level
    level.Request("A", rAFact.get());
    level.Request("Smoother", smootherFact.get());
    level.Request("PreSmoother", smootherFact.get());
    level.Request("PostSmoother", smootherFact.get());

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
    Teuchos::RCP<BlockedMultiVector> bX   = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    Teuchos::RCP<MultiVector> XX          = bX->Merge();
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

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 5e-15, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <,
                         5e-15, out, success);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedDirectSolver, NestedII20I1I_Thyra_BlockedDirectSolver_Setup_Apply, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  std::array<std::string, 3> solver_types{"", "Klu", "Superlu"};

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  if (comm->getSize() > 1) {
    out << "Skipping test for " << comm->getSize()
        << " processors as Amesos2 cannot deal with non-standard maps in the non-serial case." << std::endl;
    return;
  }

  for (const auto& solver_type : solver_types) {
    Teuchos::RCP<const BlockedCrsMatrix> bop = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrixThyra(lib, 3, comm);
    Teuchos::RCP<const Matrix> Aconst        = Teuchos::rcp_dynamic_cast<const Matrix>(bop);
    Teuchos::RCP<Matrix> A                   = Teuchos::rcp_const_cast<Matrix>(Aconst);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);
    level.setlib(lib);

    // Test ReorderBlockAFactory
    Teuchos::RCP<ReorderBlockAFactory> rAFact = Teuchos::rcp(new ReorderBlockAFactory());
    rAFact->SetFactory("A", MueLu::NoFactory::getRCP());
    rAFact->SetParameter(std::string("Reorder Type"), Teuchos::ParameterEntry(std::string("[ [2 0] 1]")));

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
    RCP<BlockedDirectSolver> smoProtoPredict = Teuchos::rcp(new BlockedDirectSolver(solver_type));
    smoProtoPredict->SetFactory("A", sA[0]);
    sF[0] = rcp(new SmootherFactory(smoProtoPredict));
    sM[0] = rcp(new FactoryManager());
    sM[0]->SetFactory("A", sA[0]);
    sM[0]->SetFactory("Smoother", sF[0]);
    sM[0]->SetIgnoreUserData(true);

    smootherPrototype->SetVelocityPredictionFactoryManager(sM[0]);

    RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", rAFact);

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(Teuchos::as<Scalar>(1.0)));
    SFact->SetFactory("A", rAFact);
    SFact->SetFactory("Ainv", AinvFact);

    RCP<SmootherPrototype> smoProtoCorrect = rcp(new TrilinosSmoother(std::string("RELAXATION"), Teuchos::ParameterList(), 0));
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
    Teuchos::RCP<BlockedMultiVector> bX   = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(X);
    Teuchos::RCP<MultiVector> XX          = bX->Merge();
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

    TEUCHOS_TEST_COMPARE(residualNorm1[0], <, 5e-15, out, success);
    TEUCHOS_TEST_COMPARE(finalNorms[0] - Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::one()), <,
                         5e-15, out, success);
  }
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                                                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedDirectSolver, BlockedDirectSolver_Setup_Apply, SC, LO, GO, NO)                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedDirectSolver, NestedI53I42II01II_BlockedDirectSolver_Setup_Apply, SC, LO, GO, NO)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedDirectSolver, NestedI53I42II01II_Thyra_BlockedDirectSolver_Setup_Apply, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedDirectSolver, NestedII20I1I_BlockedDirectSolver_Setup_Apply, SC, LO, GO, NO)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedDirectSolver, NestedII20I1I_Thyra_BlockedDirectSolver_Setup_Apply, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>
}  // namespace MueLuTests
