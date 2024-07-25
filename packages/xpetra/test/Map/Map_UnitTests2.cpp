// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

#include "Xpetra_DefaultPlatform.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMap.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMap.hpp"
#endif

// This file regroups tests that are specific to Xpetra.

namespace {
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::Comm;
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

RCP<const Comm<int> > getDefaultComm() {
  using Xpetra::DefaultPlatform;

  if (testMpi) {
    return DefaultPlatform::getDefaultPlatform().getComm();
  }
  return rcp(new Teuchos::SerialComm<int>());
}

// Test of correctness for the return value of the function: Map::getRemoteIndexList()
// (getRemoteIndexList() uses Xpetra::LookupStatus toXpetra(int))
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Map, getRemoteIndexList, M, LO, GO, N) {
  typedef typename Teuchos::ArrayView<int>::const_iterator IntConstIt;
  typedef typename ArrayRCP<GO>::size_type GOSize;

  RCP<const Comm<int> > comm    = getDefaultComm();
  const int numImages           = comm->getSize();
  const size_t numGlobalEntries = numImages * 4;

  const M map(numGlobalEntries, 0, comm);

  size_t numID = numGlobalEntries;
  ArrayRCP<GO> GIDList(numID);
  ArrayRCP<int> nodeIDList(GIDList.size());
  Xpetra::LookupStatus ls;

  // Test 1 (AllIDsPresent)
  for (GOSize i = 0; i < GIDList.size(); i++) {
    GIDList[i] = i;
  }
  ls = map.getRemoteIndexList(GIDList(), nodeIDList());
  TEST_EQUALITY(ls, Xpetra::AllIDsPresent);
  out << "nodeIDList(): ";
  for (IntConstIt it = nodeIDList.begin(); it != nodeIDList.end(); ++it) {
    out << *it << " ";
  }

  // Test 2 (IDNotPresent)
  for (GOSize i = 0; i < GIDList.size(); i++) {
    GIDList[i] = i;
  }
  GIDList[GIDList.size() - 1] = numGlobalEntries + 1;  // IDNotPresent
  ls                          = map.getRemoteIndexList(GIDList(), nodeIDList());
  TEST_EQUALITY(ls, Xpetra::IDNotPresent);
  out << "nodeIDList(): ";
  for (IntConstIt it = nodeIDList.begin(); it != nodeIDList.end(); ++it) {
    out << *it << " ";
  }
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
#define XP_MAP_INSTANT(LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Map, getRemoteIndexList, M##LO##GO##N, LO, GO, N)

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
