// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iterator>
#include <set>
#include <Epetra_CrsMatrix.h>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_CrsMatrixFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

namespace XpetraBlockMatrixTests {

//////////////////////////////////////////////////////////////////////////
// EPETRA helper functions
Teuchos::RCP<Epetra_Map> SplitMap(const Epetra_Map& Amap, const Epetra_Map& Agiven) {
  const Epetra_Comm& Comm = Amap.Comm();
  const Epetra_Map& Ag    = Agiven;

  int count = 0;
  std::vector<int> myaugids(Amap.NumMyElements());
  for (int i = 0; i < Amap.NumMyElements(); ++i) {
    const int gid = Amap.GID(i);
    if (Ag.MyGID(gid)) continue;
    myaugids[count] = gid;
    ++count;
  }
  myaugids.resize(count);
  int gcount;
  Comm.SumAll(&count, &gcount, 1);
  Teuchos::RCP<Epetra_Map> Aunknown = Teuchos::rcp(new Epetra_Map(gcount, count, &myaugids[0], 0, Comm));

  return Aunknown;
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

Teuchos::RCP<Epetra_Map> CreateMap(const std::set<int>& gids, const Epetra_Comm& comm) {
  std::vector<int> mapvec;
  mapvec.reserve(gids.size());
  mapvec.assign(gids.begin(), gids.end());
  Teuchos::RCP<Epetra_Map> map =
      Teuchos::rcp(new Epetra_Map(-1,
                                  mapvec.size(),
                                  &mapvec[0],
                                  0,
                                  comm));
  mapvec.clear();
  return map;
}

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

Teuchos::RCP<Epetra_Map> MergeMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps) {
  if (maps.size() == 0)
    std::cout << "no maps to merge" << std::endl;
  for (unsigned i = 0; i < maps.size(); ++i) {
    if (maps[i] == Teuchos::null)
      std::cout << "can not merge extractor with null maps" << std::endl;
    if (maps[i]->UniqueGIDs() == false)
      std::cout << "map " << i << " not unique" << std::endl;
  }
  std::set<int> mapentries;
  for (unsigned i = 0; i < maps.size(); ++i) {
    const Epetra_Map& map = *maps[i];
    std::copy(map.MyGlobalElements(),
              map.MyGlobalElements() + map.NumMyElements(),
              std::inserter(mapentries, mapentries.begin()));
  }
  return CreateMap(mapentries, maps[0]->Comm());
}

}  // end namespace XpetraBlockMatrixTests
