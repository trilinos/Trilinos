// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMap.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMap.hpp"
#endif

#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_StridedMap.hpp"

namespace {
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;
using Xpetra::DefaultPlatform;
using Xpetra::global_size_t;

bool testMpi         = true;
double errorTolSlack = 1e+1;

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used.");
  clp.setOption(
      "error-tol-slack", &errorTolSlack,
      "Slack off of machine epsilon used to check test results");
}

RCP<const Comm<int> > getDefaultComm() {
  if (testMpi) {
    return DefaultPlatform::getDefaultPlatform().getComm();
  }
  return rcp(new Teuchos::SerialComm<int>());
}

//
// UNIT TESTS
//

////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StridedMapFactory, CreateStridedMap1, M, LO, GO, N) {
  typedef Xpetra::StridedMap<LO, GO, N> SM;

  // create a comm
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
  const int numImages                          = comm->getSize();

  // test constructor for Xpetra::StridedMaps
  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  for (int indexBase = 0; indexBase < 2; indexBase++) {
    GO offset = 111;

    // constructor calls: (num global elements, index base)
    global_size_t numGlobalElements = 120 * numImages;
    size_t numLocalElements         = 120;
    std::vector<size_t> stridedInfo;
    stridedInfo.push_back(3);
    stridedInfo.push_back(4);
    stridedInfo.push_back(5);

    Teuchos::RCP<SM> map = Xpetra::StridedMapFactory<LO, GO, N>::Build(lib, numGlobalElements, indexBase, stridedInfo, comm, -1, offset);

    TEST_EQUALITY_CONST(map->getFixedBlockSize(), 12);
    TEST_EQUALITY_CONST(map->isStrided(), true);
    TEST_EQUALITY_CONST(map->isBlocked(), true);
    TEST_EQUALITY_CONST(map->getMinAllGlobalIndex(), indexBase + offset);
    TEST_EQUALITY_CONST(map->getMaxAllGlobalIndex(), indexBase + offset + Teuchos::as<GO>(numGlobalElements) - 1);
    TEST_EQUALITY_CONST(map->isContiguous(), false);
    TEST_EQUALITY_CONST(map->getLocalNumElements() % 12, 0);

    // Teuchos::RCP<SM> emap2 = Teuchos::null;
    for (size_t k = 0; k < stridedInfo.size(); k++) {
      Teuchos::RCP<SM> map2 = Xpetra::StridedMapFactory<LO, GO, N>::Build(lib, numGlobalElements, indexBase, stridedInfo, comm, k, offset);
      TEST_EQUALITY_CONST(map2->getFixedBlockSize(), 12);
      TEST_EQUALITY_CONST(map2->isStrided(), true);
      TEST_EQUALITY_CONST(map2->isBlocked(), true);
      TEST_EQUALITY_CONST(map2->isContiguous(), false);
      TEST_EQUALITY_CONST(map2->getLocalNumElements() % stridedInfo[k], 0);
      TEST_EQUALITY_CONST(map2->getLocalNumElements(), numLocalElements / map2->getFixedBlockSize() * stridedInfo[k]);
    }
  }
}
// TODO add test routines for remaining constructors of StridedMapFactory

////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StridedMapFactory, CreateStridedMap2, M, LO, GO, N) {
  typedef Xpetra::StridedMap<LO, GO, N> SM;

  // create a comm
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
  const int numImages                          = comm->getSize();

  // test constructor for Xpetra::StridedMaps
  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  GO offset = 111;

  // constructor calls: (num global elements, index base)
  global_size_t numGlobalElements = 120 * numImages;
  size_t numLocalElements         = 120;
  std::vector<size_t> stridedInfo;
  stridedInfo.push_back(3);
  stridedInfo.push_back(4);
  stridedInfo.push_back(5);

  Teuchos::RCP<SM> map = Xpetra::StridedMapFactory<LO, GO, N>::Build(lib, numGlobalElements, 0, stridedInfo, comm, -1, offset);

  TEST_EQUALITY_CONST(map->getFixedBlockSize(), 12);
  TEST_EQUALITY_CONST(map->isStrided(), true);
  TEST_EQUALITY_CONST(map->isBlocked(), true);
  TEST_EQUALITY_CONST(map->getMinAllGlobalIndex(), offset);
  TEST_EQUALITY_CONST(map->getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 1);
  TEST_EQUALITY_CONST(map->isContiguous(), false);
  TEST_EQUALITY_CONST(map->getLocalNumElements() % 12, 0);

  Teuchos::RCP<SM> map2 = Xpetra::StridedMapFactory<LO, GO, N>::Build(map, 0);
  TEST_EQUALITY_CONST(map2->getFixedBlockSize(), 12);
  TEST_EQUALITY_CONST(map2->isStrided(), true);
  TEST_EQUALITY_CONST(map2->isBlocked(), true);
  TEST_EQUALITY_CONST(map2->getMinAllGlobalIndex(), offset);
  TEST_EQUALITY_CONST(map2->getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 10);
  TEST_EQUALITY_CONST(map2->isContiguous(), false);
  TEST_EQUALITY_CONST(map2->getLocalNumElements() % 3, 0);
  TEST_EQUALITY_CONST(map2->getLocalNumElements(), numLocalElements / map2->getFixedBlockSize() * stridedInfo[0]);

  Teuchos::RCP<SM> map3 = Xpetra::StridedMapFactory<LO, GO, N>::Build(map, 1);
  TEST_EQUALITY_CONST(map3->getFixedBlockSize(), 12);
  TEST_EQUALITY_CONST(map3->isStrided(), true);
  TEST_EQUALITY_CONST(map3->isBlocked(), true);
  TEST_EQUALITY_CONST(map3->getMinAllGlobalIndex(), offset + 3);
  TEST_EQUALITY_CONST(map3->getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 6);
  TEST_EQUALITY_CONST(map3->isContiguous(), false);
  TEST_EQUALITY_CONST(map3->getLocalNumElements() % 4, 0);
  TEST_EQUALITY_CONST(map3->getLocalNumElements(), numLocalElements / map3->getFixedBlockSize() * stridedInfo[1]);

  Teuchos::RCP<SM> map4 = Xpetra::StridedMapFactory<LO, GO, N>::Build(map, 2);
  TEST_EQUALITY_CONST(map4->getFixedBlockSize(), 12);
  TEST_EQUALITY_CONST(map4->isStrided(), true);
  TEST_EQUALITY_CONST(map4->isBlocked(), true);
  TEST_EQUALITY_CONST(map4->getMinAllGlobalIndex(), offset + 7);
  TEST_EQUALITY_CONST(map4->getMaxAllGlobalIndex(), offset + Teuchos::as<GO>(numGlobalElements) - 1);
  TEST_EQUALITY_CONST(map4->isContiguous(), false);
  TEST_EQUALITY_CONST(map4->getLocalNumElements() % 5, 0);
  TEST_EQUALITY_CONST(map4->getLocalNumElements(), numLocalElements / map4->getFixedBlockSize() * stridedInfo[2]);
}

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
#define XP_MAP_INSTANT(LO, GO, N)                                                                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StridedMapFactory, CreateStridedMap1, M##LO##GO##N, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StridedMapFactory, CreateStridedMap2, M##LO##GO##N, LO, GO, N)

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

}  // end namespace
