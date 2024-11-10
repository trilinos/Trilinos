// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>

#include "Teuchos_Array.hpp"
#include "Teuchos_as.hpp"

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMap.hpp"
#include "Tpetra_Details_Behavior.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMap.hpp"
#endif

namespace {

using Teuchos::Array;
using Teuchos::RCP;

bool testMpi = true;

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used.");
}

RCP<const Teuchos::Comm<int> > getDefaultComm() {
  if (testMpi) {
    return Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
  }
  return Teuchos::rcp(new Teuchos::SerialComm<int>());
}

//
// UNIT TESTS
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MapFactory, ContigUniformMapFact, M, LO, GO, N) {
  // create a comm
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  const Xpetra::global_size_t numGlobalElements = 20;

  // Get underlying lib in a very crummy way
  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // Create the map via a MapFactory
  using MapFactory   = Xpetra::MapFactory<LO, GO, N>;
  using Map          = Xpetra::Map<LO, GO, N>;
  RCP<const Map> map = MapFactory::Build(lib, numGlobalElements, Teuchos::ScalarTraits<GO>::zero(), comm);

  TEST_INEQUALITY(map, Teuchos::null);
  TEST_EQUALITY_CONST(map->getGlobalNumElements(), numGlobalElements);

  // All procs fail if any proc fails
  int globalSuccess_int = -1;
  reduceAll(*comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int));
  TEST_EQUALITY_CONST(globalSuccess_int, 0);
}  // ContigUniformMapFact

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MapFactory, ContigUserMapFact, M, LO, GO, N) {
  // create a comm
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  const size_t numRanks                         = comm->getSize();
  const size_t numLocalElements                 = 5;
  const Xpetra::global_size_t numGlobalElements = Teuchos::as<Xpetra::global_size_t>(numLocalElements * numRanks);

  // Get underlying lib in a very crummy way
  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // Create the map via a MapFactory
  using MapFactory   = Xpetra::MapFactory<LO, GO, N>;
  using Map          = Xpetra::Map<LO, GO, N>;
  RCP<const Map> map = MapFactory::Build(lib, numGlobalElements, numLocalElements, Teuchos::ScalarTraits<GO>::zero(), comm);

  TEST_INEQUALITY(map, Teuchos::null);
  TEST_EQUALITY_CONST(map->getLocalNumElements(), numLocalElements);
  TEST_EQUALITY_CONST(map->getGlobalNumElements(), numGlobalElements);

  // All procs fail if any proc fails
  int globalSuccess_int = -1;
  reduceAll(*comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int));
  TEST_EQUALITY_CONST(globalSuccess_int, 0);
}  // ContigUserMapFact

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MapFactory, NonContigUserMapFact, M, LO, GO, N) {
  // create a comm
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  const int myRank                              = comm->getRank();
  const size_t numRanks                         = comm->getSize();
  const size_t numLocalElements                 = 5;
  const Xpetra::global_size_t numGlobalElements = Teuchos::as<Xpetra::global_size_t>(numLocalElements * numRanks);

  // Get underlying lib in a very crummy way
  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  const GO offsetPerRank = 2 * myRank * numLocalElements;
  Array<GO> gidList(numLocalElements);
  for (size_t i = 0; i < numLocalElements; ++i)
    gidList[i] = i + offsetPerRank;

  // Create the map via a MapFactory
  using MapFactory   = Xpetra::MapFactory<LO, GO, N>;
  using Map          = Xpetra::Map<LO, GO, N>;
  RCP<const Map> map = MapFactory::Build(lib, numGlobalElements, gidList, Teuchos::ScalarTraits<GO>::zero(), comm);

  TEST_INEQUALITY(map, Teuchos::null);
  TEST_EQUALITY_CONST(map->getLocalNumElements(), numLocalElements);
  TEST_EQUALITY_CONST(map->getGlobalNumElements(), numGlobalElements);
  TEST_EQUALITY_CONST(map->getMinGlobalIndex(), offsetPerRank);
  TEST_EQUALITY(map->getMaxGlobalIndex(), Teuchos::as<GO>(offsetPerRank + numLocalElements - 1));

  // All procs fail if any proc fails
  int globalSuccess_int = -1;
  reduceAll(*comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int));
  TEST_EQUALITY_CONST(globalSuccess_int, 0);
}  // NonContigUserMapFact

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MapFactory, TransformNumDofsPerNode, M, LO, GO, N) {
  // create a comm
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  const Xpetra::global_size_t numNodes = 10;
  const LO numDofsPerNode              = 2;

  // Get underlying lib in a very crummy way
  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // Create the map via a MapFactory
  using MapFactory       = Xpetra::MapFactory<LO, GO, N>;
  using Map              = Xpetra::Map<LO, GO, N>;
  RCP<const Map> nodeMap = MapFactory::Build(lib, numNodes, Teuchos::ScalarTraits<GO>::zero(), comm);

  TEST_INEQUALITY(nodeMap, Teuchos::null);
  TEST_EQUALITY_CONST(nodeMap->getGlobalNumElements(), numNodes);

  RCP<const Map> dofMap = MapFactory::Build(nodeMap, numDofsPerNode);

  TEST_INEQUALITY(dofMap, Teuchos::null);
  TEST_EQUALITY(dofMap->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(numNodes * numDofsPerNode));
  TEST_EQUALITY(dofMap->getMinAllGlobalIndex(), Teuchos::ScalarTraits<GO>::zero());
  TEST_EQUALITY(dofMap->getMaxAllGlobalIndex(), Teuchos::as<GO>(numNodes * numDofsPerNode - 1));

  // All procs fail if any proc fails
  int globalSuccess_int = -1;
  reduceAll(*comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int));
  TEST_EQUALITY_CONST(globalSuccess_int, 0);
}  // TransformNumDofsPerNode

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MapFactory, TransformNumDofsPerNodeWithOffset, M, LO, GO, N) {
  // create a comm
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  const GO numNodes           = 10;
  const LO numDofsPerNode     = 2;
  const GO gidOffsetForDofMap = 3;

  // Get underlying lib in a very crummy way
  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // Create the map via a MapFactory
  using MapFactory       = Xpetra::MapFactory<LO, GO, N>;
  using Map              = Xpetra::Map<LO, GO, N>;
  RCP<const Map> nodeMap = MapFactory::Build(lib, numNodes, Teuchos::ScalarTraits<GO>::zero(), comm);

  TEST_INEQUALITY(nodeMap, Teuchos::null);
  TEST_EQUALITY_CONST(nodeMap->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(numNodes));

  RCP<const Map> dofMap = MapFactory::Build(nodeMap, numDofsPerNode, gidOffsetForDofMap);

  TEST_INEQUALITY(dofMap, Teuchos::null);
  TEST_EQUALITY(dofMap->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(numNodes * numDofsPerNode));
  TEST_EQUALITY_CONST(dofMap->getMinAllGlobalIndex(), Teuchos::as<GO>(gidOffsetForDofMap));
  TEST_EQUALITY(dofMap->getMaxAllGlobalIndex(), Teuchos::as<GO>(gidOffsetForDofMap + numNodes * numDofsPerNode - 1));

  // All procs fail if any proc fails
  int globalSuccess_int = -1;
  reduceAll(*comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int));
  TEST_EQUALITY_CONST(globalSuccess_int, 0);
}  // TransformNumDofsPerNodeWithOffset

//
// INSTANTIATIONS
//
#ifdef HAVE_XPETRA_TPETRA

#define XPETRA_TPETRA_TYPES(LO, GO, N) \
  typedef typename Xpetra::TpetraMap<LO, GO, N> M##LO##GO##N;

#endif

#ifdef HAVE_XPETRA_EPETRA

#define XPETRA_EPETRA_TYPES(LO, GO, N) \
  typedef typename Xpetra::EpetraMapT<GO, N> M##LO##GO##N;

#endif

// List of tests (which run both on Epetra and Tpetra)
#define XP_MAP_INSTANT(LO, GO, N)                                                                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MapFactory, ContigUniformMapFact, M##LO##GO##N, LO, GO, N)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MapFactory, ContigUserMapFact, M##LO##GO##N, LO, GO, N)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MapFactory, NonContigUserMapFact, M##LO##GO##N, LO, GO, N)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MapFactory, TransformNumDofsPerNode, M##LO##GO##N, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MapFactory, TransformNumDofsPerNodeWithOffset, M##LO##GO##N, LO, GO, N)

#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
// no ordinal types as scalar for testing as some tests use ScalarTraits::eps...
TPETRA_INSTANTIATE_LGN(XPETRA_TPETRA_TYPES)
TPETRA_INSTANTIATE_LGN(XP_MAP_INSTANT)

#endif

#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp"  // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(int, int, EpetraNode)
XP_MAP_INSTANT(int, int, EpetraNode)
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(int, LongLong, EpetraNode)
XP_MAP_INSTANT(int, LongLong, EpetraNode)
#endif
#endif

}  // namespace
