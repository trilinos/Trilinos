// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Xpetra_UnitTestHelpers.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_as.hpp>

#include <Xpetra_Map.hpp>

#include <Xpetra_MapUtils.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_BlockedMultiVector.hpp>
#include <Xpetra_ReorderedBlockedMultiVector.hpp>
#include <Xpetra_Exceptions.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include <set>

namespace XpetraBlockMatrixTests {

double errorTolSlack = 1e+1;

Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm() {
  return Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
}

/////////////////////////////////////////////////////

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "error-tol-slack", &errorTolSlack,
      "Slack off of machine epsilon used to check test results");
}

//
// Helper routines
//
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
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateMultiVector(int noBlocks, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
  typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;

  GlobalOrdinal nOverallDOFGidsPerProc = Teuchos::as<GlobalOrdinal>(Teuchos::ScalarTraits<GlobalOrdinal>::pow(2, noBlocks - 2)) * 10;

  GlobalOrdinal procOffset = comm->getRank() * nOverallDOFGidsPerProc;

  std::set<GlobalOrdinal> myDOFGids;
  for (GlobalOrdinal i = 0; i < nOverallDOFGidsPerProc; i++)
    myDOFGids.insert(i + procOffset);

  Teuchos::RCP<Map> fullmap = CreateMap<LocalOrdinal, GlobalOrdinal, Node, MapType>(myDOFGids, *comm);

  // create Multivector
  Teuchos::RCP<MultiVector> vv = MultiVectorFactory::Build(fullmap, 2, true);

  // fill multivector data (first multivector contains the GID, the second the LID as scalar)
  Teuchos::ArrayRCP<Scalar> vv1 = vv->getDataNonConst(0);
  Teuchos::ArrayRCP<Scalar> vv2 = vv->getDataNonConst(1);
  for (LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(vv->getLocalLength()); ++i) {
    vv1[i] = Teuchos::as<Scalar>(fullmap->getGlobalElement(i));
    vv2[i] = Teuchos::as<Scalar>(i);
  }

  return vv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class MapType>
Teuchos::RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateBlockedMultiVector(int noBlocks, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedMultiVector;
  typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractor;

  GlobalOrdinal nOverallDOFGidsPerProc = Teuchos::as<GlobalOrdinal>(Teuchos::ScalarTraits<GlobalOrdinal>::pow(2, noBlocks - 2)) * 10;

  GlobalOrdinal procOffset = comm->getRank() * nOverallDOFGidsPerProc;

  std::set<GlobalOrdinal> myDOFGids;
  for (GlobalOrdinal i = 0; i < nOverallDOFGidsPerProc; i++)
    myDOFGids.insert(i + procOffset);

  Teuchos::RCP<Map> fullmap = CreateMap<LocalOrdinal, GlobalOrdinal, Node, MapType>(myDOFGids, *comm);

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

    Teuchos::RCP<Map> halfmap = CreateMap<LocalOrdinal, GlobalOrdinal, Node, MapType>(myHalfGIDs, *comm);

    Teuchos::RCP<Map> secondmap = SplitMap<LocalOrdinal, GlobalOrdinal, Node, MapType>(*remainingpartmap, *halfmap);
    remainingpartmap            = halfmap;

    maps[noBlocks - 1 - it] = secondmap;
  }

  // create map extractor (Xpetra mode)
  Teuchos::RCP<const MapExtractor> xpMapExtractor = Teuchos::rcp(new MapExtractor(fullmap, maps, false));

  // create Multivector
  Teuchos::RCP<MultiVector> vv = MultiVectorFactory::Build(fullmap, 2, true);

  // fill multivector data (first multivector contains the GID, the second the LID as scalar)
  Teuchos::ArrayRCP<Scalar> vv1 = vv->getDataNonConst(0);
  Teuchos::ArrayRCP<Scalar> vv2 = vv->getDataNonConst(1);
  for (LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(vv->getLocalLength()); ++i) {
    vv1[i] = Teuchos::as<Scalar>(fullmap->getGlobalElement(i));
    vv2[i] = Teuchos::as<Scalar>(i);
  }

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = Teuchos::rcp(new BlockedMultiVector(xpMapExtractor, vv));

  return bvv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class MapType>
Teuchos::RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateBlockedMapBlockedMultiVector(int noBlocks, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  // typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedMultiVector;
  // typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;
  typedef Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node> BlockedMap;

  GlobalOrdinal nOverallDOFGidsPerProc = Teuchos::as<GlobalOrdinal>(Teuchos::ScalarTraits<GlobalOrdinal>::pow(2, noBlocks - 2)) * 10;

  GlobalOrdinal procOffset = comm->getRank() * nOverallDOFGidsPerProc;

  std::set<GlobalOrdinal> myDOFGids;
  for (GlobalOrdinal i = 0; i < nOverallDOFGidsPerProc; i++)
    myDOFGids.insert(i + procOffset);

  Teuchos::RCP<Map> fullmap = CreateMap<LocalOrdinal, GlobalOrdinal, Node, MapType>(myDOFGids, *comm);

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

    Teuchos::RCP<Map> halfmap = CreateMap<LocalOrdinal, GlobalOrdinal, Node, MapType>(myHalfGIDs, *comm);

    Teuchos::RCP<Map> secondmap = SplitMap<LocalOrdinal, GlobalOrdinal, Node, MapType>(*remainingpartmap, *halfmap);
    remainingpartmap            = halfmap;

    maps[noBlocks - 1 - it] = secondmap;
  }

  // create map extractor (Xpetra mode)
  Teuchos::RCP<const BlockedMap> xpBmap = Teuchos::rcp(new BlockedMap(fullmap, maps, false));

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = Teuchos::rcp(new BlockedMultiVector(xpBmap, 2, true));

  return bvv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class MapType>
Teuchos::RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateBlockedMultiVectorThyra(int noBlocks, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactory;
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedMultiVector;
  typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractor;

  MapType testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  std::vector<Teuchos::RCP<const Map> > maps(noBlocks, Teuchos::null);
  maps[0] = MapFactory::Build(lib, comm->getSize() * 5, 5, 0, comm);
  for (int it = 1; it < noBlocks; it++) {
    GlobalOrdinal localDofs = Teuchos::as<GlobalOrdinal>(Teuchos::ScalarTraits<GlobalOrdinal>::pow(2, it - 1) * 5);
    maps[it]                = MapFactory::Build(lib, comm->getSize() * localDofs, localDofs, 0, comm);
  }

  // create map extractor
  // To generate the Thyra style map extractor we do not need a full map but only the
  // information about the Map details (i.e. lib and indexBase). We can extract this
  // information from maps[0]
  Teuchos::RCP<const MapExtractor> me =
      Teuchos::rcp(new MapExtractor(maps[0], maps, true));

  // create Multivector
  Teuchos::RCP<MultiVector> vv = MultiVectorFactory::Build(me->getFullMap(), 2, true);

  // fill multivector data (first multivector contains the GID, the second the LID as scalar)
  Teuchos::ArrayRCP<Scalar> vv1 = vv->getDataNonConst(0);
  Teuchos::ArrayRCP<Scalar> vv2 = vv->getDataNonConst(1);
  for (LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(vv->getLocalLength()); ++i) {
    vv1[i] = Teuchos::as<Scalar>(vv->getMap()->getGlobalElement(i));
    vv2[i] = Teuchos::as<Scalar>(i);
  }

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = Teuchos::rcp(new BlockedMultiVector(me, vv));

  return bvv;
}

//
// UNIT TESTS
//

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, Constructor, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_EQUALITY(bvv->getBlockedMap()->getNumMaps(), Teuchos::as<size_t>(noBlocks));
  for (size_t r = 0; r < bvv->getBlockedMap()->getNumMaps(); ++r) {
    Teuchos::RCP<MultiVector> bvvi = bvv->getMultiVector(r);
    TEST_EQUALITY(bvvi->getMap()->isSameAs(*(bvv->getBlockedMap()->getMap(r))), true);
    Teuchos::ArrayRCP<const Scalar> bvvi1 = bvvi->getData(0);
    Teuchos::ArrayRCP<const Scalar> bvvi2 = bvvi->getData(1);
    for (LO l = 0; l < Teuchos::as<LO>(bvvi->getLocalLength()); ++l) {
      TEST_EQUALITY(bvvi1[l], Teuchos::as<Scalar>(bvvi->getMap()->getGlobalElement(l)));
    }
  }

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms(bvv->getNumVectors());
  Teuchos::Array<Magnitude> fnorms(vv->getNumVectors());

  TEST_NOTHROW(vv->norm1(fnorms));
  TEST_NOTHROW(bvv->norm1(bnorms));
  TEST_COMPARE_FLOATING_ARRAYS(fnorms, bnorms, Teuchos::ScalarTraits<Magnitude>::zero());
  TEST_NOTHROW(vv->norm2(fnorms));
  TEST_NOTHROW(bvv->norm2(bnorms));
  TEST_COMPARE_FLOATING_ARRAYS(fnorms, bnorms, Teuchos::ScalarTraits<Magnitude>::zero());
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, Constructor2, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::MultiVectorFactory<Scalar, LO, GO, Node> MultiVectorFactory;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMapBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_EQUALITY(bvv->getBlockedMap()->getNumMaps(), Teuchos::as<size_t>(noBlocks));
  for (size_t r = 0; r < bvv->getBlockedMap()->getNumMaps(); ++r) {
    Teuchos::RCP<MultiVector> vvv = MultiVectorFactory::Build(bvv->getBlockedMap()->getMap(r), 2, true);

    // fill multivector data (first multivector contains the GID, the second the LID as scalar)
    Teuchos::ArrayRCP<Scalar> vv1 = vvv->getDataNonConst(0);
    Teuchos::ArrayRCP<Scalar> vv2 = vvv->getDataNonConst(1);
    for (LO i = 0; i < Teuchos::as<LO>(vvv->getLocalLength()); ++i) {
      vv1[i] = Teuchos::as<Scalar>(vvv->getMap()->getGlobalElement(i));
      vv2[i] = Teuchos::as<Scalar>(i);
    }
    TEST_NOTHROW(bvv->setMultiVector(r, vvv, false));
  }

  for (size_t r = 0; r < bvv->getBlockedMap()->getNumMaps(); ++r) {
    Teuchos::RCP<MultiVector> bvvi = bvv->getMultiVector(r);
    TEST_EQUALITY(bvvi->getMap()->isSameAs(*(bvv->getBlockedMap()->getMap(r))), true);
    Teuchos::ArrayRCP<const Scalar> bvvi1 = bvvi->getData(0);
    Teuchos::ArrayRCP<const Scalar> bvvi2 = bvvi->getData(1);
    for (LO l = 0; l < Teuchos::as<LO>(bvvi->getLocalLength()); ++l) {
      TEST_EQUALITY(bvvi1[l], Teuchos::as<Scalar>(bvvi->getMap()->getGlobalElement(l)));
      TEST_EQUALITY(bvvi2[l], Teuchos::as<Scalar>(l));
    }
  }

  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms(bvv->getNumVectors());
  Teuchos::Array<Magnitude> fnorms(vv->getNumVectors());

  TEST_NOTHROW(vv->norm1(fnorms));
  TEST_NOTHROW(bvv->norm1(bnorms));
  TEST_EQUALITY(fnorms[0], bnorms[0]);
  TEST_INEQUALITY(fnorms[1], bnorms[1]);
  TEST_NOTHROW(vv->norm2(fnorms));
  TEST_NOTHROW(bvv->norm2(bnorms));
  TEST_EQUALITY(fnorms[0], bnorms[0]);
  TEST_INEQUALITY(fnorms[1], bnorms[1]);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, Norm1, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 2;

  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms(bvv->getNumVectors());
  Teuchos::Array<Magnitude> fnorms(vv->getNumVectors());

  TEST_NOTHROW(vv->norm1(fnorms));
  TEST_NOTHROW(bvv->norm1(bnorms));
  TEST_COMPARE_FLOATING_ARRAYS(fnorms, bnorms, Teuchos::ScalarTraits<Magnitude>::zero());
  Magnitude result = Teuchos::ScalarTraits<Magnitude>::zero();
  for (GO gg = 0; gg < Teuchos::as<GO>(vv->getMap()->getGlobalNumElements()); gg++)
    result += Teuchos::as<Magnitude>(gg);
  TEST_EQUALITY(bnorms[0], result);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv2 = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::Array<Magnitude> bnorms2(bvv2->getNumVectors());
  TEST_NOTHROW(bvv2->norm1(bnorms2));
  TEST_EQUALITY(bnorms2[0], result);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, Norm2, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 2;

  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms(bvv->getNumVectors());
  Teuchos::Array<Magnitude> fnorms(vv->getNumVectors());

  TEST_NOTHROW(vv->norm2(fnorms));
  TEST_NOTHROW(bvv->norm2(bnorms));
  TEST_COMPARE_FLOATING_ARRAYS(fnorms, bnorms, Teuchos::ScalarTraits<Magnitude>::zero());
  Magnitude result = Teuchos::ScalarTraits<Magnitude>::zero();
  for (GO gg = 0; gg < Teuchos::as<GO>(vv->getMap()->getGlobalNumElements()); gg++)
    result += Teuchos::as<Magnitude>(gg) * Teuchos::as<Magnitude>(gg);
  result = Teuchos::ScalarTraits<Magnitude>::squareroot(result);
  TEST_EQUALITY(bnorms[0], result);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv2 = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::Array<Magnitude> bnorms2(bvv2->getNumVectors());
  TEST_NOTHROW(bvv2->norm2(bnorms2));
  TEST_COMPARE(bnorms2[0] - result, <, 1e-10);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, NormInf, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 2;

  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms(bvv->getNumVectors());
  Teuchos::Array<Magnitude> fnorms(vv->getNumVectors());

  TEST_NOTHROW(vv->normInf(fnorms));
  TEST_NOTHROW(bvv->normInf(bnorms));
  TEST_COMPARE_FLOATING_ARRAYS(fnorms, bnorms, Teuchos::ScalarTraits<Magnitude>::zero());
  Magnitude result = Teuchos::ScalarTraits<Magnitude>::zero();
  for (GO gg = 0; gg < Teuchos::as<GO>(vv->getMap()->getGlobalNumElements()); gg++)
    result = std::max(result, Teuchos::as<Magnitude>(gg));
  TEST_EQUALITY(bnorms[0], result);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv2 = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::Array<Magnitude> bnorms2(bvv2->getNumVectors());
  TEST_NOTHROW(bvv2->normInf(bnorms2));
  TEST_COMPARE(bnorms2[0] - result, <, 1e-10);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, Scale, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 2;

  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  typedef typename STS::magnitudeType Magnitude;
  Teuchos::Array<Magnitude> bnorms(bvv->getNumVectors());
  Teuchos::Array<Magnitude> fnorms(vv->getNumVectors());

  TEST_NOTHROW(vv->normInf(fnorms));
  TEST_NOTHROW(bvv->normInf(bnorms));
  TEST_COMPARE_FLOATING_ARRAYS(fnorms, bnorms, Teuchos::ScalarTraits<Magnitude>::zero());
  Magnitude myresult = Teuchos::ScalarTraits<Magnitude>::zero();
  for (GO gg = 0; gg < Teuchos::as<GO>(vv->getMap()->getGlobalNumElements()); gg++)
    myresult = std::max(myresult, Teuchos::as<Magnitude>(gg));
  TEST_EQUALITY(bnorms[0], myresult);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv2 = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::Array<Magnitude> bnorms2(bvv2->getNumVectors());
  TEST_NOTHROW(bvv2->normInf(bnorms2));
  TEST_COMPARE(bnorms2[0] - myresult, <, 1e-10);

  bvv->scale(Teuchos::as<Scalar>(2.0));
  vv->scale(Teuchos::as<Scalar>(2.0));
  Teuchos::Array<Magnitude> scaled_bnorms(bvv->getNumVectors());
  Teuchos::Array<Magnitude> scaled_fnorms(vv->getNumVectors());
  TEST_NOTHROW(vv->normInf(scaled_fnorms));
  TEST_NOTHROW(bvv->normInf(scaled_bnorms));
  TEST_COMPARE_FLOATING_ARRAYS(scaled_fnorms, scaled_bnorms, Teuchos::ScalarTraits<Magnitude>::zero());
  myresult = Teuchos::ScalarTraits<Magnitude>::zero();
  for (GO gg = 0; gg < Teuchos::as<GO>(vv->getMap()->getGlobalNumElements()); gg++)
    myresult = std::max(myresult, Teuchos::as<Magnitude>(gg));
  TEST_EQUALITY(scaled_bnorms[0], Teuchos::as<Magnitude>(2.0) * myresult);

  // create BlockedMultiVector
  bvv2 = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);
  bvv2->scale(Teuchos::as<Scalar>(2.0));
  Teuchos::Array<Magnitude> scaled_bnorms2(bvv2->getNumVectors());
  TEST_NOTHROW(bvv2->normInf(scaled_bnorms2));
  TEST_COMPARE(scaled_bnorms2[0] - Teuchos::as<Magnitude>(2.0) * myresult, <, 1e-10);
}

#if 0  // functionality not required any more
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, ExtractVector, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::Map<LO, GO, Node> Map;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractor;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // extract map extractor
  Teuchos::RCP<const MapExtractor> me  = bvv->getMapExtractor();

  for(size_t r = 0; r < me->NumMaps(); ++r) {
    Teuchos::RCP<MultiVector> partB = me->ExtractVector(bvv,r);
    Teuchos::RCP<MultiVector> partV = me->ExtractVector(vv,r);
    TEST_EQUALITY(partB->getMap()->isSameAs(*(partV->getMap())),true);
    Teuchos::ArrayRCP<const Scalar > partBd = partB->getData(0);
    Teuchos::ArrayRCP<const Scalar > partVd = partV->getData(0);
    TEST_COMPARE_FLOATING_ARRAYS(partBd,partVd,Teuchos::ScalarTraits<Magnitude>::zero());
#ifdef HAVE_XPETRA_DEBUG
    TEST_THROW(partB = me->ExtractVector(bvv,r,true),Xpetra::Exceptions::RuntimeError);
#endif
    TEST_THROW(partV = me->ExtractVector(vv,r,true),Xpetra::Exceptions::RuntimeError);
  }

  // create a new faulty MapExtractor
  std::vector<Teuchos::RCP<const Map> > maps(noBlocks, Teuchos::null);
  for(size_t r = 0; r < me->NumMaps(); ++r) {
    maps[me->NumMaps() - 1 - r] = me->getMap(r);
  }
  Teuchos::RCP<const MapExtractor> meFaulty = Teuchos::rcp(new MapExtractor(me->getFullMap(), maps, false));

  // the faulty map extractor reverses the ordering of maps. We have five partial maps (0-4)
  // Therefore, partial map with index 2 are the same in the faulty and the original map extractor and
  // the ExtractVector call then succeeds
  for(size_t r = 0; r < me->NumMaps(); ++r) {
    TEST_NOTHROW(Teuchos::RCP<MultiVector> partV = meFaulty->ExtractVector(vv,r));
    if(r!=2)      {
      TEST_NOTHROW(Teuchos::RCP<MultiVector> partB = meFaulty->ExtractVector(bvv,r));
      Teuchos::RCP<MultiVector> partV = meFaulty->ExtractVector(vv,r);
      Teuchos::RCP<MultiVector> partB = meFaulty->ExtractVector(bvv,r);
      TEST_EQUALITY(partB->getMap()->isSameAs(*(partV->getMap())),false);
    }
    else if(r==2) {
      TEST_NOTHROW(Teuchos::RCP<MultiVector> partB = meFaulty->ExtractVector(bvv,r));
      Teuchos::RCP<MultiVector> partV = meFaulty->ExtractVector(vv,r);
      Teuchos::RCP<MultiVector> partB = meFaulty->ExtractVector(bvv,r);
      TEST_EQUALITY(partB->getMap()->isSameAs(*(partV->getMap())),true);
      Teuchos::ArrayRCP<const Scalar > partBd = partB->getData(0);
      Teuchos::ArrayRCP<const Scalar > partVd = partV->getData(0);
      TEST_COMPARE_FLOATING_ARRAYS(partBd,partVd,Teuchos::ScalarTraits<Magnitude>::zero());
    }
  }

}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, ExtractVectorThyra, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractor;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create full vector
  Teuchos::RCP<MultiVector>         vv = bvv->Merge();

  // extract map extractor
  Teuchos::RCP<const MapExtractor> me  = bvv->getMapExtractor();

  for(size_t r = 0; r < me->NumMaps(); ++r) {
    Teuchos::RCP<MultiVector> partB = me->ExtractVector(bvv,r,true);
    Teuchos::RCP<MultiVector> partV = me->ExtractVector(vv,r,false);
    TEST_EQUALITY(partB->getMap()->isSameAs(*(partV->getMap())) == false || r==0,true);
    TEST_EQUALITY(partB->getMap()->getMinAllGlobalIndex(),0);
    partV = me->ExtractVector(vv,r,true);
    TEST_EQUALITY(partB->getMap()->isSameAs(*(partV->getMap())),true);
    Teuchos::ArrayRCP<const Scalar > partBd = partB->getData(0);
    Teuchos::ArrayRCP<const Scalar > partVd = partV->getData(0);
    TEST_COMPARE_FLOATING_ARRAYS(partBd,partVd,Teuchos::ScalarTraits<Magnitude>::zero());
  }

}


TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, InsertVector, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractor;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_EQUALITY(bvv->getNumVectors(), 2);

  // extract map extractor
  Teuchos::RCP<const MapExtractor> me  = bvv->getMapExtractor();


  for(size_t r = 0; r < me->NumMaps(); ++r) {
    Teuchos::RCP<MultiVector> part = me->getVector(r,bvv->getNumVectors(),false);
    part->putScalar(STS::one());
    me->InsertVector(part,r,bvv);
  }

  for(size_t r = 0; r < me->NumMaps(); ++r) {
    Teuchos::RCP<const MultiVector> part = me->ExtractVector(bvv,r);
    Teuchos::ArrayRCP<const Scalar > partd1 = part->getData(0);
    Teuchos::ArrayRCP<const Scalar > partd2 = part->getData(1);
    for(LO l = 0; l < Teuchos::as<LO>(part->getLocalLength()); l++)
      TEST_EQUALITY(partd1[l], STS::one());
    TEST_COMPARE_FLOATING_ARRAYS(partd1,partd2,Teuchos::ScalarTraits<Magnitude>::zero());
  }

#ifdef HAVE_XPETRA_DEBUG
  // create malicious multivector
  Teuchos::RCP<MultiVector> part1 = me->getVector(0,23,false,true);
  TEST_THROW(me->InsertVector(part1,0,bvv),Xpetra::Exceptions::RuntimeError);
  Teuchos::RCP<MultiVector> part2 = me->getVector(0,2,false,true);
  TEST_THROW(me->InsertVector(part2,1,bvv),Xpetra::Exceptions::RuntimeError);
  TEST_THROW(Teuchos::RCP<MultiVector> part3 = me->getVector(1,2,true,true),Xpetra::Exceptions::RuntimeError);
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( BlockedMultiVector, InsertVectorThyra, M, MA, Scalar, LO, GO, Node )
{
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Xpetra::MapExtractor<Scalar,LO,GO,Node> MapExtractor;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector>         vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVectorThyra<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_EQUALITY(bvv->getNumVectors(), 2);

  // extract map extractor
  Teuchos::RCP<const MapExtractor> me  = bvv->getMapExtractor();
  TEST_EQUALITY(me->getThyraMode(),true);

  for(size_t r = 0; r < me->NumMaps(); ++r) {
    Teuchos::RCP<MultiVector> part = me->getVector(r,bvv->getNumVectors(),true);
    TEST_EQUALITY(part->getMap()->getMinAllGlobalIndex(),0);
    TEST_NOTHROW(part->putScalar(STS::one()));
    TEST_NOTHROW(me->InsertVector(part,r,bvv,me->getThyraMode()));
  }

  for(size_t r = 0; r < me->NumMaps(); ++r) {
    Teuchos::RCP<const MultiVector> part = me->ExtractVector(bvv,r,me->getThyraMode());
    Teuchos::ArrayRCP<const Scalar > partd1 = part->getData(0);
    Teuchos::ArrayRCP<const Scalar > partd2 = part->getData(1);
    for(LO l = 0; l < Teuchos::as<LO>(part->getLocalLength()); l++)
      TEST_EQUALITY(partd1[l], STS::one());
    TEST_COMPARE_FLOATING_ARRAYS(partd1,partd2,Teuchos::ScalarTraits<Magnitude>::zero());
  }

#ifdef HAVE_XPETRA_DEBUG
  // create malicious multivector
  Teuchos::RCP<MultiVector> part1 = me->getVector(0,23,true,true);
  TEST_THROW(me->InsertVector(part1,0,bvv),Xpetra::Exceptions::RuntimeError);
  // unfortunately, in Thyra mode there is no error thrown since the vectors in
  // block 0 and 1 have the same length (and the same GIDs)
  Teuchos::RCP<MultiVector> part2 = me->getVector(0,2,true,true);
  TEST_NOTHROW(me->InsertVector(part2,1,bvv,me->getThyraMode()));
  // This should throw, thought
  Teuchos::RCP<MultiVector> part3 = me->getVector(0,2,true,true);
  TEST_THROW(me->InsertVector(part2,2,bvv,me->getThyraMode()),Xpetra::Exceptions::RuntimeError);
#endif
}

#endif

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, UpdateVector1, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  using STS      = Teuchos::ScalarTraits<Scalar>;
  using mag_type = typename STS::magnitudeType;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv1 = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::RCP<BlockedMultiVector> bvv2 = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_EQUALITY(size_t(bvv1->getNumVectors()), size_t(2));
  TEST_EQUALITY(size_t(bvv2->getNumVectors()), size_t(2));

  // NOTE (mfh 23 Feb 2020) You can't always multiply double by
  // Scalar, e.g., if Scalar=complex<float>.
  TEST_NOTHROW(bvv1->update(mag_type(-0.35) * STS::one(), *bvv2, mag_type(0.7) * STS::one()));
  TEST_NOTHROW(bvv1->update(mag_type(-0.35) * STS::one(), *bvv2, STS::one()));

  Teuchos::Array<mag_type> bnorms1(bvv1->getNumVectors());
  Teuchos::Array<mag_type> bnorms2(vv->getNumVectors());
  TEST_NOTHROW(bvv1->norm1(bnorms1));
  // TAW: CUDA produces a "dirty zero" (not exactly zero)
  // this might be numerical effects caused by the ordering of calculations
  //
  // FIXME (mfh 23 Feb 2020) Tolerances should depend on Scalar.
  TEST_COMPARE(bnorms1[0], <, 3e-12);
  TEST_COMPARE(bnorms1[1], <, 3e-12);
  // TEST_EQUALITY( bnorms1[0], STS::zero());
  // TEST_EQUALITY( bnorms1[1], Teuchos::ScalarTraits<mag_type>STS::zero());
  TEST_NOTHROW(vv->norm1(bnorms1));
  TEST_NOTHROW(bvv2->norm1(bnorms2));
  TEST_EQUALITY(bnorms1[0], bnorms2[0]);
  TEST_EQUALITY(bnorms1[1], bnorms2[1]);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, UpdateVector1b, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  using STS      = Teuchos::ScalarTraits<Scalar>;
  using mag_type = typename STS::magnitudeType;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv1 = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  TEST_EQUALITY(bvv1->getNumVectors(), 2);

  Teuchos::Array<mag_type> bnorms1(bvv1->getNumVectors());
  Teuchos::Array<mag_type> bnorms2(vv->getNumVectors());

  // NOTE (mfh 23 Feb 2020) You can't always multiply double by
  // Scalar, e.g., if Scalar=complex<float>.
  TEST_NOTHROW(bvv1->update(mag_type(-0.5) * STS::one(), *vv, STS::one()));
  TEST_NOTHROW(bvv1->norm1(bnorms1));
  TEST_NOTHROW(vv->norm1(bnorms2));
  TEST_EQUALITY(bnorms1[0], mag_type(0.5) * bnorms2[0]);
  TEST_EQUALITY(bnorms1[1], mag_type(0.5) * bnorms2[1]);

#ifdef HAVE_XPETRA_DEBUG
  // create faulty multivector
  typedef Xpetra::MultiVectorFactory<Scalar, LO, GO, Node> MultiVectorFactory;
  Teuchos::RCP<MultiVector> vvx = MultiVectorFactory::Build(bvv1->getMap(), 1, true);
  TEST_THROW(bvv1->update(STS::one(), *vvx, STS::one()), Xpetra::Exceptions::RuntimeError);
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, UpdateVector2, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType mag_type;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv1 = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::RCP<BlockedMultiVector> bvv2 = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  Teuchos::RCP<BlockedMultiVector> bvv3 = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_EQUALITY(bvv1->getNumVectors(), 2);
  TEST_EQUALITY(bvv2->getNumVectors(), 2);
  TEST_EQUALITY(bvv3->getNumVectors(), 2);

  // NOTE (mfh 23 Feb 2020) You can't always multiply double by
  // Scalar, e.g., if Scalar=complex<float>.
  TEST_NOTHROW(bvv1->update(mag_type(-0.25) * STS::one(), *bvv2, mag_type(-0.25) * STS::one(), *bvv3, mag_type(0.5) * STS::one()));

  Teuchos::Array<mag_type> bnorms1(bvv1->getNumVectors());
  Teuchos::Array<mag_type> bnorms2(vv->getNumVectors());
  TEST_NOTHROW(bvv1->norm1(bnorms1));
  TEST_EQUALITY(bnorms1[0], STS::zero());
  TEST_EQUALITY(bnorms1[1], STS::zero());
  TEST_NOTHROW(vv->norm1(bnorms1));
  TEST_NOTHROW(bvv2->norm1(bnorms2));
  TEST_EQUALITY(bnorms1[0], bnorms2[0]);
  TEST_EQUALITY(bnorms1[1], bnorms2[1]);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, PutScalar, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  using STS      = Teuchos::ScalarTraits<Scalar>;
  using mag_type = typename STS::magnitudeType;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  TEST_NOTHROW(bvv->putScalar(STS::one()));

  Teuchos::Array<mag_type> bnorms(bvv->getNumVectors());
  TEST_NOTHROW(bvv->norm1(bnorms));
  TEST_EQUALITY(bnorms[0], Teuchos::as<mag_type>(bvv->getBlockedMap()->getFullMap()->getGlobalNumElements()));
  TEST_EQUALITY(bnorms[1], Teuchos::as<mag_type>(bvv->getBlockedMap()->getFullMap()->getGlobalNumElements()));

  // NOTE (mfh 23 Feb 2020) You can't always multiply double by
  // Scalar, e.g., if Scalar=complex<float>.
  TEST_NOTHROW(bvv->putScalar(mag_type(3.0) * STS::one()));

  for (size_t r = 0; r < bvv->getBlockedMap()->getNumMaps(); ++r) {
    Teuchos::RCP<const MultiVector> part   = bvv->getMultiVector(r);
    Teuchos::ArrayRCP<const Scalar> partd1 = part->getData(0);
    Teuchos::ArrayRCP<const Scalar> partd2 = part->getData(1);
    for (LO l = 0; l < Teuchos::as<LO>(part->getLocalLength()); l++) {
      TEST_EQUALITY(partd1[l], mag_type(3.0) * STS::one());
      TEST_EQUALITY(partd1[l], partd2[l]);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, MultiVectorFactory, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::MultiVectorFactory<Scalar, LO, GO, Node> MultiVectorFactory;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Xpetra::BlockedMap<LO, GO, Node> BlockedMap;
  typedef Xpetra::Map<LO, GO, Node> Map;
  // typedef Teuchos::ScalarTraits<Scalar> STS;
  // typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // extract BlockedMap
  Teuchos::RCP<const Map> vvm = bvv->getMap();
  TEST_EQUALITY(Teuchos::rcp_dynamic_cast<const BlockedMap>(vvm).is_null(), false);
  TEST_EQUALITY(vvm->getMap()->isSameAs(*(vv->getMap())), true);

  Teuchos::RCP<const BlockedMap> bvvm = bvv->getBlockedMap();
  TEST_EQUALITY(bvvm.is_null(), false);
  TEST_EQUALITY(bvvm->getMap()->isSameAs(*(vv->getMap())), true);
  TEST_EQUALITY(bvvm->getNumMaps(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bvvm->getThyraMode(), false);

  // create new BlockedMultiVector based on bvvm
  Teuchos::RCP<MultiVector> vv2 = MultiVectorFactory::Build(bvvm, 3, true);
  TEST_EQUALITY(vv2.is_null(), false);
  TEST_EQUALITY(Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vv2).is_null(), false);
  Teuchos::RCP<BlockedMultiVector> bvv2 = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vv2);
  TEST_EQUALITY(bvv2->getNumVectors(), 3);
  TEST_EQUALITY(bvv2->getBlockedMap()->getNumMaps(), Teuchos::as<size_t>(noBlocks));
#ifdef HAVE_XPETRA_DEBUG
  TEST_THROW(bvv2->setMultiVector(0, bvv->getMultiVector(0), false), Xpetra::Exceptions::RuntimeError);
#endif

  // create a new standard multivector (with an underlying BlockedMap)
  Teuchos::RCP<MultiVector> vv3 = MultiVectorFactory::Build(vvm, 3, true);
  TEST_EQUALITY(vv3.is_null(), false);
  TEST_EQUALITY(Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vv3).is_null(), false);
  Teuchos::RCP<BlockedMultiVector> bvv3 = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vv3);
  TEST_EQUALITY(bvv3->getNumVectors(), 3);
  TEST_EQUALITY(bvv3->getBlockedMap()->getNumMaps(), Teuchos::as<size_t>(noBlocks));
#ifdef HAVE_XPETRA_DEBUG
  TEST_THROW(bvv3->setMultiVector(0, bvv->getMultiVector(0), false), Xpetra::Exceptions::RuntimeError);
#endif

  // create a new standard multivector
  Teuchos::RCP<MultiVector> vv4 = MultiVectorFactory::Build(vv->getMap(), 3, true);
  TEST_EQUALITY(vv4.is_null(), false);
  TEST_EQUALITY(Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vv4).is_null(), true);
  TEST_EQUALITY(vv4->getMap()->isSameAs(*(vv3->getMap())), true);
  TEST_EQUALITY(vv4->getMap()->isSameAs(*(vv2->getMap())), true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, Merge, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  // typedef Xpetra::Map<LO, GO, Node> Map;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType Magnitude;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 5;

  // create full vector
  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // merge BlockedMultiVector into MultiVector
  Teuchos::RCP<MultiVector> vv2 = bvv->Merge();

  // test merged multivector
  vv2->update(STS::one(), *vv, -STS::one());

  Teuchos::Array<Magnitude> bnorms(vv2->getNumVectors());
  TEST_NOTHROW(vv2->norm1(bnorms));
  TEST_EQUALITY(bnorms[0], STS::zero());
  TEST_EQUALITY(bnorms[1], STS::zero());
  TEST_NOTHROW(vv2->norm2(bnorms));
  TEST_EQUALITY(bnorms[0], STS::zero());
  TEST_EQUALITY(bnorms[1], STS::zero());
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, ConstructorReordered, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Map<LO, GO, Node> Map;
  typedef Xpetra::BlockedMap<LO, GO, Node> BlockedMap;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 6;

  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<const BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  TEST_EQUALITY(bvv.is_null(), false);

  Teuchos::ArrayRCP<const Scalar> vData = bvv->getMultiVector(0)->getData(0);
  for (size_t i = 0; i < bvv->getBlockedMap()->getMap(0, false)->getLocalNumElements(); i++) {
    TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bvv->getBlockedMap()->getMap(0, false)->getGlobalElement(i)));
  }

  vData = bvv->getMultiVector(1)->getData(0);
  for (size_t i = 0; i < bvv->getBlockedMap()->getMap(1, false)->getLocalNumElements(); i++) {
    TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bvv->getBlockedMap()->getMap(1, false)->getGlobalElement(i)));
  }

  vData = bvv->getMultiVector(2)->getData(0);
  for (size_t i = 0; i < bvv->getBlockedMap()->getMap(2, false)->getLocalNumElements(); i++) {
    TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bvv->getBlockedMap()->getMap(2, false)->getGlobalElement(i)));
  }

  // first reordered multivector

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 0 [ 1 2 ] [ [ 3 4 ] 5 ] ]");

  Teuchos::RCP<const MultiVector> bmv = buildReorderedBlockedMultiVector(brm, bvv);
  TEST_EQUALITY(bmv.is_null(), false);
  {
    Teuchos::RCP<const BlockedMultiVector> bbmv = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(bmv);
    TEST_EQUALITY(bbmv.is_null(), false);
    vData = bbmv->getMultiVector(0)->getData(0);
    for (size_t i = 0; i < bbmv->getBlockedMap()->getMap(0, false)->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bbmv->getBlockedMap()->getMap(0, false)->getGlobalElement(i)));
    }
  }

  Teuchos::RCP<const Map> fullmap = bmv->getMap();
  TEST_EQUALITY(fullmap.is_null(), false);
  Teuchos::RCP<const BlockedMap> fullBlockedMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(fullmap);
  TEST_EQUALITY(fullBlockedMap.is_null(), false);
  TEST_EQUALITY(fullBlockedMap->getNumMaps(), 3);
  TEST_EQUALITY(fullBlockedMap->getThyraMode(), false);
  TEST_EQUALITY(fullBlockedMap->getMap(1, false)->getMinAllGlobalIndex(), 5);
  TEST_EQUALITY(fullBlockedMap->getMap(1, false)->getMaxAllGlobalIndex(), (comm->getSize() - 1) * 160 + 19);
  Teuchos::RCP<const Map> map0         = fullBlockedMap->getMap(0, false);
  Teuchos::RCP<const BlockedMap> bmap0 = Teuchos::rcp_dynamic_cast<const BlockedMap>(map0);
  TEST_EQUALITY(bmap0.is_null(), true);
  Teuchos::RCP<const Map> map1         = fullBlockedMap->getMap(1, false);
  Teuchos::RCP<const BlockedMap> bmap1 = Teuchos::rcp_dynamic_cast<const BlockedMap>(map1);
  TEST_EQUALITY(bmap1.is_null(), false);
  TEST_EQUALITY(fullBlockedMap->getMap(2, false)->getMinAllGlobalIndex(), 20);
  TEST_EQUALITY(fullBlockedMap->getMap(2, false)->getMaxAllGlobalIndex(), comm->getSize() * 160 - 1);
  Teuchos::RCP<const Map> map2         = fullBlockedMap->getMap(2, false);
  Teuchos::RCP<const BlockedMap> bmap2 = Teuchos::rcp_dynamic_cast<const BlockedMap>(map2);
  TEST_EQUALITY(bmap2.is_null(), false);

  Teuchos::RCP<const BlockedMultiVector> bbmv = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(bmv);
  TEST_EQUALITY(bbmv.is_null(), false);
  TEST_EQUALITY(bbmv->getNumVectors(), 2);
  TEST_EQUALITY(bbmv->getBlockedMap()->getNumMaps(), 3);
  Teuchos::RCP<const MultiVector> bmv1 = bbmv->getMultiVector(1, false);
  TEST_EQUALITY(bmv1.is_null(), false);
  Teuchos::RCP<const BlockedMultiVector> bmv11 = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(bmv1);
  TEST_EQUALITY(bmv11.is_null(), false);
  TEST_EQUALITY(bmv11->getBlockedMap()->getNumMaps(), 2);
  TEST_EQUALITY(bmv11->getBlockedMap()->getMap(0, false)->getMinAllGlobalIndex(), 5);
  TEST_EQUALITY(bmv11->getBlockedMap()->getMap(0, false)->getMaxAllGlobalIndex(), (comm->getSize() - 1) * 160 + 9);
  TEST_EQUALITY(bmv11->getBlockedMap()->getMap(1, false)->getMinAllGlobalIndex(), 10);
  TEST_EQUALITY(bmv11->getBlockedMap()->getMap(1, false)->getMaxAllGlobalIndex(), (comm->getSize() - 1) * 160 + 19);
  {
    vData = bmv11->getMultiVector(0)->getData(0);
    for (size_t i = 0; i < bmv11->getBlockedMap()->getMap(0, false)->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bmv11->getBlockedMap()->getMap(0, false)->getGlobalElement(i)));
    }
  }
  {
    vData = bmv11->getMultiVector(1)->getData(0);
    for (size_t i = 0; i < bmv11->getBlockedMap()->getMap(1, false)->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bmv11->getBlockedMap()->getMap(1, false)->getGlobalElement(i)));
    }
  }
  Teuchos::RCP<const MultiVector> bmv2 = bbmv->getMultiVector(2, false);
  TEST_EQUALITY(bmv2.is_null(), false);
  Teuchos::RCP<const BlockedMultiVector> bmv21 = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(bmv2);
  TEST_EQUALITY(bmv21.is_null(), false);
  TEST_EQUALITY(bmv21->getBlockedMap()->getNumMaps(), 2);
  TEST_EQUALITY(bmv21->getBlockedMap()->getMap(0, false)->getMinAllGlobalIndex(), 20);
  TEST_EQUALITY(bmv21->getBlockedMap()->getMap(0, false)->getMaxAllGlobalIndex(), (comm->getSize() - 1) * 160 + 79);
  TEST_EQUALITY(bmv21->getBlockedMap()->getMap(1, false)->getMinAllGlobalIndex(), 80);
  TEST_EQUALITY(bmv21->getBlockedMap()->getMap(1, false)->getMaxAllGlobalIndex(), comm->getSize() * 160 - 1);

  Teuchos::RCP<const BlockedMultiVector> bmv211 = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(bmv21->getMultiVector(0, false));
  TEST_EQUALITY(bmv211.is_null(), false);
  TEST_EQUALITY(bmv211->getBlockedMap()->getNumMaps(), 2);
  TEST_EQUALITY(bmv211->getBlockedMap()->getMap(0, false)->getMinAllGlobalIndex(), 20);
  TEST_EQUALITY(bmv211->getBlockedMap()->getMap(0, false)->getMaxAllGlobalIndex(), (comm->getSize() - 1) * 160 + 39);
  TEST_EQUALITY(bmv211->getBlockedMap()->getMap(1, false)->getMinAllGlobalIndex(), 40);
  TEST_EQUALITY(bmv211->getBlockedMap()->getMap(1, false)->getMaxAllGlobalIndex(), (comm->getSize() - 1) * 160 + 79);

  Teuchos::RCP<const MultiVector> bmv212 = bmv21->getMultiVector(1, false);
  TEST_EQUALITY(bmv212.is_null(), false);
  Teuchos::RCP<const BlockedMultiVector> bmv212t = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(bmv212);
  TEST_EQUALITY(bmv212t.is_null(), false);
  TEST_EQUALITY(bmv212t->getBlockedMap()->getNumMaps(), 1);
  TEST_EQUALITY(bmv212t->getBlockedMap()->getMap(0, false)->getMinAllGlobalIndex(), 80);
  TEST_EQUALITY(bmv212t->getBlockedMap()->getMap(0, false)->getMaxAllGlobalIndex(), comm->getSize() * 160 - 1);
  TEST_EQUALITY(bmv212t->getMap()->getMinAllGlobalIndex(), 80);
  TEST_EQUALITY(bmv212t->getMap()->getMaxAllGlobalIndex(), comm->getSize() * 160 - 1);

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm2 = Xpetra::blockedReorderFromString("[ 4 1 ]");

  Teuchos::RCP<const MultiVector> bvvv = buildReorderedBlockedMultiVector(brm2, bvv);
  TEST_EQUALITY(bvvv.is_null(), false);
  Teuchos::RCP<const BlockedMultiVector> brvvv = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(bvvv);
  TEST_EQUALITY(brvvv.is_null(), false);
  TEST_EQUALITY(brvvv->getBlockedMap()->getNumMaps(), 2);
  TEST_EQUALITY(brvvv->getBlockedMap()->getMap(0, false)->getMinAllGlobalIndex(), 40);
  TEST_EQUALITY(brvvv->getBlockedMap()->getMap(0, false)->getMaxAllGlobalIndex(), (comm->getSize() - 1) * 160 + 79);
  TEST_EQUALITY(brvvv->getBlockedMap()->getMap(1, false)->getMinAllGlobalIndex(), 5);
  TEST_EQUALITY(brvvv->getBlockedMap()->getMap(1, false)->getMaxAllGlobalIndex(), (comm->getSize() - 1) * 160 + 9);
  TEST_EQUALITY(brvvv->getMap()->getMinAllGlobalIndex(), 5);
  TEST_EQUALITY(brvvv->getMap()->getMaxAllGlobalIndex(), (comm->getSize() - 1) * 160 + 79);

  Teuchos::RCP<const MultiVector> bm = brvvv->Merge();
  {
    vData = bm->getData(0);
    for (size_t i = 0; i < bm->getLocalLength(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bm->getMap()->getGlobalElement(i)));
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, ConstructorReorderedSmall, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 3;

  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<const BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  TEST_EQUALITY(bvv.is_null(), false);

  Teuchos::ArrayRCP<const Scalar> vData  = bvv->getMultiVector(0)->getData(0);
  Teuchos::ArrayRCP<const Scalar> vData2 = bvv->getMultiVector(0)->getData(1);
  for (size_t i = 0; i < bvv->getBlockedMap()->getMap(0, false)->getLocalNumElements(); i++) {
    TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bvv->getBlockedMap()->getMap(0, false)->getGlobalElement(i)));
    TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(i));
  }

  vData  = bvv->getMultiVector(1)->getData(0);
  vData2 = bvv->getMultiVector(1)->getData(1);
  for (size_t i = 0; i < bvv->getBlockedMap()->getMap(1, false)->getLocalNumElements(); i++) {
    TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bvv->getBlockedMap()->getMap(1, false)->getGlobalElement(i)));
    TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(5 + i));
  }

  vData  = bvv->getMultiVector(2)->getData(0);
  vData2 = bvv->getMultiVector(2)->getData(1);
  for (size_t i = 0; i < bvv->getBlockedMap()->getMap(2, false)->getLocalNumElements(); i++) {
    TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bvv->getBlockedMap()->getMap(2, false)->getGlobalElement(i)));
    TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(10 + i));
  }

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 0 [ 1 2 ] ]");

  Teuchos::RCP<const MultiVector> bmv = buildReorderedBlockedMultiVector(brm, bvv);
  TEST_EQUALITY(bmv.is_null(), false);
  Teuchos::RCP<const BlockedMultiVector> bbmv = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(bmv);
  TEST_EQUALITY(bbmv.is_null(), false);
  Teuchos::RCP<const BlockedMultiVector> bbmv1 = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(bbmv->getMultiVector(1));
  TEST_EQUALITY(bbmv1.is_null(), false);

  {
    vData  = bbmv->getMultiVector(0)->getData(0);
    vData2 = bbmv->getMultiVector(0)->getData(1);
    for (size_t i = 0; i < bbmv->getBlockedMap()->getMap(0, false)->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bbmv->getBlockedMap()->getMap(0, false)->getGlobalElement(i)));
      TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(i));
    }

    vData  = bbmv1->getMultiVector(0, false)->getData(0);
    vData2 = bbmv1->getMultiVector(0, false)->getData(1);
    for (size_t i = 0; i < bbmv1->getBlockedMap()->getMap(0, false)->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bbmv1->getBlockedMap()->getMap(0, false)->getGlobalElement(i)));
      TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(5 + i));
    }

    vData  = bbmv1->getMultiVector(1, false)->getData(0);
    vData2 = bbmv1->getMultiVector(1, false)->getData(1);
    for (size_t i = 0; i < bbmv1->getBlockedMap()->getMap(1, false)->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bbmv1->getBlockedMap()->getMap(1, false)->getGlobalElement(i)));
      TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(10 + i));
    }
  }

  Teuchos::RCP<const MultiVector> mmv = bbmv->Merge();
  TEST_EQUALITY(mmv.is_null(), false);

  {
    vData  = mmv->getData(0);
    vData2 = mmv->getData(1);
    for (size_t i = 0; i < mmv->getMap()->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(mmv->getMap()->getGlobalElement(i)));
      TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(i));
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, ConstructorReorderedSmall2, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 3;

  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<const BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  TEST_EQUALITY(bvv.is_null(), false);

  Teuchos::ArrayRCP<const Scalar> vData  = bvv->getMultiVector(0)->getData(0);
  Teuchos::ArrayRCP<const Scalar> vData2 = bvv->getMultiVector(0)->getData(1);
  for (size_t i = 0; i < bvv->getBlockedMap()->getMap(0, false)->getLocalNumElements(); i++) {
    TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bvv->getBlockedMap()->getMap(0, false)->getGlobalElement(i)));
    TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(i));
  }

  vData  = bvv->getMultiVector(1)->getData(0);
  vData2 = bvv->getMultiVector(1)->getData(1);
  for (size_t i = 0; i < bvv->getBlockedMap()->getMap(1, false)->getLocalNumElements(); i++) {
    TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bvv->getBlockedMap()->getMap(1, false)->getGlobalElement(i)));
    TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(5 + i));
  }

  vData  = bvv->getMultiVector(2)->getData(0);
  vData2 = bvv->getMultiVector(2)->getData(1);
  for (size_t i = 0; i < bvv->getBlockedMap()->getMap(2, false)->getLocalNumElements(); i++) {
    TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bvv->getBlockedMap()->getMap(2, false)->getGlobalElement(i)));
    TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(10 + i));
  }

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ [ 2 0 ] 1 ]");

  Teuchos::RCP<const MultiVector> bmv = buildReorderedBlockedMultiVector(brm, bvv);
  TEST_EQUALITY(bmv.is_null(), false);
  Teuchos::RCP<const BlockedMultiVector> bbmv = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(bmv);
  TEST_EQUALITY(bbmv.is_null(), false);
  Teuchos::RCP<const BlockedMultiVector> bbmv0 = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(bbmv->getMultiVector(0));
  TEST_EQUALITY(bbmv0.is_null(), false);

  {
    vData  = bbmv->getMultiVector(1)->getData(0);
    vData2 = bbmv->getMultiVector(1)->getData(1);
    for (size_t i = 0; i < bbmv->getBlockedMap()->getMap(1, false)->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bbmv->getBlockedMap()->getMap(1, false)->getGlobalElement(i)));
      TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(5 + i));
    }

    vData  = bbmv0->getMultiVector(0, false)->getData(0);
    vData2 = bbmv0->getMultiVector(0, false)->getData(1);
    for (size_t i = 0; i < bbmv0->getBlockedMap()->getMap(0, false)->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bbmv0->getBlockedMap()->getMap(0, false)->getGlobalElement(i)));
      TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(10 + i));
    }

    vData  = bbmv0->getMultiVector(1, false)->getData(0);
    vData2 = bbmv0->getMultiVector(1, false)->getData(1);
    for (size_t i = 0; i < bbmv0->getBlockedMap()->getMap(1, false)->getLocalNumElements(); i++) {
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(bbmv0->getBlockedMap()->getMap(1, false)->getGlobalElement(i)));
      TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(i));
    }
  }

  Teuchos::RCP<const MultiVector> mmv = bbmv->Merge();
  TEST_EQUALITY(mmv.is_null(), false);

  {
    vData  = mmv->getData(0);
    vData2 = mmv->getData(1);
    for (size_t i = 0; i < mmv->getMap()->getLocalNumElements(); i++) {
      GO expected = 42, expected2 = 43;
      if (i < 10) expected = comm->getRank() * 20 + 10 + i;
      if (i >= 10 && i < 15) expected = comm->getRank() * 20 + i - 10;
      if (i >= 15 && i < 20) expected = comm->getRank() * 20 + 5 + i - 15;
      if (i < 10) expected2 = 10 + i;
      if (i >= 10 && i < 15) expected2 = i - 10;
      if (i >= 15 && i < 20) expected2 = 5 + i - 15;
      TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(expected));
      TEST_EQUALITY(vData2[i], Teuchos::as<Scalar>(expected2));
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, BlockedMapDeepCopy, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::BlockedMap<LO, GO, Node> BlockedMap;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 3;

  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<const BlockedMultiVector> bvv = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  TEST_EQUALITY(bvv.is_null(), false);

  Teuchos::RCP<const BlockedMap> ppbm = bvv->getBlockedMap();
  TEST_EQUALITY(ppbm.is_null(), false);
  TEST_EQUALITY(ppbm->getThyraMode(), false);

  Teuchos::RCP<const BlockedMap> ppbm2 = Teuchos::rcp(new BlockedMap(*ppbm));
  TEST_EQUALITY(ppbm.is_null(), false);

  TEST_EQUALITY(ppbm->isSameAs(*ppbm2), true);

  TEST_EQUALITY(ppbm->getMap(0, false)->isSameAs(*(ppbm2->getMap(0, false))), true);
  TEST_EQUALITY(ppbm->getMap(1, false)->isSameAs(*(ppbm2->getMap(1, false))), true);
#ifdef HAVE_XPETRA_DEBUG
  TEST_THROW(ppbm->getMap(0, true)->isSameAs(*(ppbm2->getMap(0, true))), Xpetra::Exceptions::RuntimeError);
  TEST_THROW(ppbm->getMap(1, true)->isSameAs(*(ppbm2->getMap(1, true))), Xpetra::Exceptions::RuntimeError);
#endif
  TEST_EQUALITY(ppbm->getMap(1, false)->getMinGlobalIndex(), ppbm2->getMap(1, false)->getMinGlobalIndex());

  ppbm = Teuchos::null;

  TEST_EQUALITY(ppbm2.is_null(), false);
  TEST_EQUALITY(ppbm2->getThyraMode(), false);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedMultiVector, BlockedVectorDeepCopy, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  int noBlocks = 3;

  Teuchos::RCP<MultiVector> vv = CreateMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);

  // create BlockedMultiVector
  Teuchos::RCP<const BlockedMultiVector> bvec = CreateBlockedMultiVector<Scalar, LO, GO, Node, M>(noBlocks, comm);
  TEST_EQUALITY(bvec.is_null(), false);

  //
  Teuchos::RCP<BlockedMultiVector> bvec2 = Teuchos::rcp(new BlockedMultiVector(bvec->getBlockedMap(), 22));
  *bvec2                                 = *bvec;  // deep copy
  TEST_EQUALITY(bvec2.is_null(), false);
  TEST_EQUALITY(bvec2->getBlockedMap()->isSameAs(*(bvec->getBlockedMap())), true);
  TEST_EQUALITY(bvec2->getNumVectors(), bvec->getNumVectors());

  Teuchos::Array<typename STS::magnitudeType> nn(bvec->getNumVectors());
  Teuchos::Array<typename STS::magnitudeType> nn2(bvec2->getNumVectors());
  TEST_NOTHROW(bvec->norm1(nn));

  bvec = Teuchos::null;

  TEST_NOTHROW(bvec2->norm1(nn2));
  for (size_t t = 0; t < bvec2->getNumVectors(); t++) {
    TEST_EQUALITY(nn[t], nn2[t]);
  }
}

//
// INSTANTIATIONS
//
#ifdef HAVE_XPETRA_TPETRA

#define XPETRA_TPETRA_TYPES(S, LO, GO, N)                     \
  typedef typename Xpetra::TpetraMap<LO, GO, N> M##LO##GO##N; \
  typedef typename Xpetra::TpetraMultiVector<S, LO, GO, N> MV##S##LO##GO##N;

#endif

#ifdef HAVE_XPETRA_EPETRA

#define XPETRA_EPETRA_TYPES(S, LO, GO, N)                  \
  typedef typename Xpetra::EpetraMapT<GO, N> M##LO##GO##N; \
  typedef typename Xpetra::EpetraMultiVectorT<GO, N> MV##S##LO##GO##N;

#endif

#define XP_BLOCKEDMULTIVECTOR_INSTANT(S, LO, GO, N)                                                                                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, Constructor, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, Constructor2, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)               \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, Norm1, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, Norm2, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, NormInf, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, Scale, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, UpdateVector1, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, UpdateVector1b, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)             \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, UpdateVector2, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, PutScalar, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, MultiVectorFactory, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, Merge, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, ConstructorReordered, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, ConstructorReorderedSmall, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, ConstructorReorderedSmall2, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, BlockedMapDeepCopy, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedMultiVector, BlockedVectorDeepCopy, M##LO##GO##N, MV##S##LO##GO##N, S, LO, GO, N)

// TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, ExtractVector, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N )
// TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, ExtractVectorThyra, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N )
// TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, InsertVector, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N )
// TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( BlockedMultiVector, InsertVectorThyra, M##LO##GO##N , MV##S##LO##GO##N, S, LO, GO, N )

// List of tests which run only with Tpetra
#define XP_TPETRA_BLOCKEDMULTIVECTOR_INSTANT(S, LO, GO, N)

// List of tests which run only with Epetra
#define XP_EPETRA_BLOCKEDMULTIVECTOR_INSTANT(S, LO, GO, N)

#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XPETRA_TPETRA_TYPES)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_TPETRA_BLOCKEDMULTIVECTOR_INSTANT)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_BLOCKEDMULTIVECTOR_INSTANT)

#endif

#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp"  // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double, int, int, EpetraNode)
XP_EPETRA_BLOCKEDMULTIVECTOR_INSTANT(double, int, int, EpetraNode)
XP_BLOCKEDMULTIVECTOR_INSTANT(double, int, int, EpetraNode)
#endif
// EpetraExt routines are not working with 64 bit
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double, int, LongLong, EpetraNode)
XP_EPETRA_BLOCKEDMULTIVECTOR_INSTANT(double, int, LongLong, EpetraNode)
XP_EPETRA_BLOCKEDMULTIVECTOR_INSTANT(double, int, LongLong, EpetraNode)
#endif

#endif

}  // namespace XpetraBlockMatrixTests
