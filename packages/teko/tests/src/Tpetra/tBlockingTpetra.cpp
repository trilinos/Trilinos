// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Thyra testing tools
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"

// Thyra includes
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "tBlockingTpetra.hpp"

#include "Teko_BlockingTpetra.hpp"

using namespace Teko::TpetraHelpers;

namespace Teko {
namespace Test {

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Thyra::createMember;
using Thyra::LinearOpBase;
using Thyra::LinearOpTester;
using Thyra::VectorBase;

void tBlockingTpetra::initializeTest() { tolerance_ = 1e-14; }

void buildGIDs(std::vector<std::vector<GO> >& gids, const Tpetra::Map<LO, GO, NT>& map) {
  LO numLocal = map.getLocalNumElements();
  LO numHalf  = numLocal / 2;
  numHalf += ((numHalf % 2 == 0) ? 0 : 1);

  gids.clear();
  gids.resize(3);

  std::vector<GO>& blk0 = gids[0];
  std::vector<GO>& blk1 = gids[1];
  std::vector<GO>& blk2 = gids[2];

  // loop over global IDs: treat first block as strided
  GO gid = -1;
  for (LO i = 0; i < numHalf; i += 2) {
    gid = map.getGlobalElement(i);
    blk0.push_back(gid);

    gid = map.getGlobalElement(i + 1);
    blk1.push_back(gid);
  }

  // loop over global IDs: treat remainder as contiguous
  for (LO i = numHalf; i < numLocal; i++) {
    gid = map.getGlobalElement(i);
    blk2.push_back(gid);
  }
}

int tBlockingTpetra::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                             int& totalrun) {
  bool allTests = true;
  bool status   = true;
  int failcount = 0;

  failstrm << "tBlockingTpetra";

  status = test_buildMaps(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"buildMaps\" ... PASSED", "   \"buildMaps\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_one2many(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"one2many\" ... PASSED", "   \"one2many\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_many2one(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"many2one\" ... PASSED", "   \"many2one\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_buildSubBlock(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"buildSubBlock\" ... PASSED",
                       "   \"buildSubBlock\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tBlockingTpetra...PASSED", "tBlockingTpetra...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tBlockingTpetra...FAILED");
  }

  return failcount;
}

bool tBlockingTpetra::test_buildMaps(int verbosity, std::ostream& os) {
  using Teuchos::as;

  bool status    = false;
  bool allPassed = true;

  int size = 30;

  TEST_MSG("\n   Builing Tpetra::Map");
  RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(size, 0, GetComm_tpetra()));

  // build gids
  std::vector<std::vector<GO> > gids;
  buildGIDs(gids, *map);

  std::vector<GO>& gid0 = gids[0];
  std::vector<GO>& gid1 = gids[1];
  std::vector<GO>& gid2 = gids[2];

  TEST_MSG("\n   Building sub maps");
  Blocking::MapPair map0 = Blocking::buildSubMap(gid0, *GetComm_tpetra());
  Blocking::MapPair map1 = Blocking::buildSubMap(gid1, *GetComm_tpetra());
  Blocking::MapPair map2 = Blocking::buildSubMap(gid2, *GetComm_tpetra());

  TEST_ASSERT(map0.first->getLocalNumElements() == gid0.size() &&
                  map0.second->getLocalNumElements() == gid0.size(),
              "   tBlockingTpetra::test_buildMaps ("
                  << Teko::Test::toString(status) << "): "
                  << " Checking map size: first=" << map0.first->getLocalNumElements()
                  << ", second=" << map0.second->getLocalNumElements() << ", gid=" << gid0.size());
  TEST_ASSERT(map1.first->getLocalNumElements() == gid1.size() &&
                  map0.second->getLocalNumElements() == gid0.size(),
              "   tBlockingTpetra::test_buildMaps ("
                  << Teko::Test::toString(status) << "): "
                  << " Checking map size: first=" << map1.first->getLocalNumElements()
                  << ", second=" << map1.second->getLocalNumElements() << ", gid=" << gid1.size());
  TEST_ASSERT(map2.first->getLocalNumElements() == gid2.size() &&
                  map0.second->getLocalNumElements() == gid0.size(),
              "   tBlockingTpetra::test_buildMaps ("
                  << Teko::Test::toString(status) << "): "
                  << " Checking map size: first=" << map2.first->getLocalNumElements()
                  << ", second=" << map2.second->getLocalNumElements() << ", gid=" << gid2.size());

  std::vector<Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > globalMaps(3);
  std::vector<Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > contigMaps(3);

  // get sub maps for convenient use and access
  globalMaps[0] = map0.first;
  globalMaps[1] = map1.first;
  globalMaps[2] = map2.first;

  contigMaps[0] = map0.second;
  contigMaps[1] = map1.second;
  contigMaps[2] = map2.second;

  // test that the extra data is attached
  TEST_ASSERT(contigMaps[0] != Teuchos::null,
              "   tBlockingTpetra::test_buildMaps (" << Teko::Test::toString(status) << ")");
  TEST_ASSERT(contigMaps[1] != Teuchos::null,
              "   tBlockingTpetra::test_buildMaps (" << Teko::Test::toString(status) << ")");
  TEST_ASSERT(contigMaps[2] != Teuchos::null,
              "   tBlockingTpetra::test_buildMaps (" << Teko::Test::toString(status) << ")");
  TEST_MSG("   tBlockingTpetra::test_buildMaps: extracted \"contigMaps\"");

  bool test;

  // check contiguous and global maps
  test = true;
  for (size_t i = 0; i < globalMaps[0]->getLocalNumElements(); i++) {
    GO gid = globalMaps[0]->getGlobalElement(i);
    GO cid = contigMaps[0]->getGlobalElement(i);

    test &= gid == gid0[i];
    test &= cid == (GO)(i + contigMaps[0]->getMinGlobalIndex());
  }
  TEST_ASSERT(test, "   tBlockingTpetra::test_buildMaps ("
                        << Teko::Test::toString(status) << "): "
                        << "checked that block maps were internally consitent");

  test = true;
  for (size_t i = 0; i < globalMaps[1]->getLocalNumElements(); i++) {
    GO gid = globalMaps[1]->getGlobalElement(i);
    GO cid = contigMaps[1]->getGlobalElement(i);

    test &= gid == gid1[i];
    test &= cid == (GO)(i + contigMaps[1]->getMinGlobalIndex());
  }
  TEST_ASSERT(test, "   tBlockingTpetra::test_buildMaps ("
                        << Teko::Test::toString(status) << "): "
                        << "checked that block maps were internally consitent");

  test = true;
  for (size_t i = 0; i < globalMaps[2]->getLocalNumElements(); i++) {
    GO gid = globalMaps[2]->getGlobalElement(i);
    GO cid = contigMaps[2]->getGlobalElement(i);

    test &= gid == gid2[i];
    test &= cid == (GO)(i + contigMaps[2]->getMinGlobalIndex());
  }
  TEST_ASSERT(test, "   tBlockingTpetra::test_buildMaps ("
                        << Teko::Test::toString(status) << "): "
                        << "checked that block maps were internally consitent");

  return allPassed;
}

bool tBlockingTpetra::test_one2many(int verbosity, std::ostream& os) {
  bool allPassed = true;

  GO size = 3 * 1000;
  TEST_MSG("\n   tBlockingTpetra::test_one2many: Builing Epetra_Map and source vector");
  RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(size, 0, GetComm_tpetra()));
  RCP<Tpetra::MultiVector<ST, LO, GO, NT> > v =
      rcp(new Tpetra::MultiVector<ST, LO, GO, NT>(map, 1));
  v->randomize();

  TEST_MSG("\n   Building sub maps");
  // build gids
  std::vector<std::vector<GO> > gids;
  buildGIDs(gids, *map);

  std::vector<GO>& gid0 = gids[0];
  std::vector<GO>& gid1 = gids[1];
  std::vector<GO>& gid2 = gids[2];

  std::vector<Blocking::MapPair> maps(3);
  maps[0] = Blocking::buildSubMap(gid0, *GetComm_tpetra());
  maps[1] = Blocking::buildSubMap(gid1, *GetComm_tpetra());
  maps[2] = Blocking::buildSubMap(gid2, *GetComm_tpetra());

  TEST_MSG("\n   tBlockingTpetra::test_one2many: Building Export/Import");
  std::vector<RCP<Tpetra::Import<LO, GO, NT> > > subImport(3);
  std::vector<RCP<Tpetra::Export<LO, GO, NT> > > subExport(3);
  for (int i = 0; i < 3; i++) {
    Blocking::ImExPair imex = Blocking::buildExportImport(*map, maps[i]);
    subImport[i]            = imex.first;
    subExport[i]            = imex.second;
  }

  TEST_MSG("\n   tBlockingTpetra::test_one2many: Building sub vectors");
  std::vector<RCP<Tpetra::MultiVector<ST, LO, GO, NT> > > subVectors;
  Blocking::buildSubVectors(maps, subVectors, 1);

  TEST_MSG("\n   tBlockingTpetra::test_one2many: Performing data copy");
  Blocking::one2many(subVectors, *v, subImport);

  // just assume it works! :)

  return allPassed;
}

bool tBlockingTpetra::test_many2one(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  GO size = 3 * 1000;
  TEST_MSG("\n   tBlockingTpetra::test_many2one: Building Tpetra::Map and source vector");
  RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(size, 0, GetComm_tpetra()));
  RCP<Tpetra::MultiVector<ST, LO, GO, NT> > v =
      rcp(new Tpetra::MultiVector<ST, LO, GO, NT>(map, 4));
  v->randomize();

  TEST_MSG("\n   Building sub maps");
  // build gids
  std::vector<std::vector<GO> > gids;
  buildGIDs(gids, *map);

  std::vector<GO>& gid0 = gids[0];
  std::vector<GO>& gid1 = gids[1];
  std::vector<GO>& gid2 = gids[2];

  std::vector<Blocking::MapPair> maps(3);
  maps[0] = Blocking::buildSubMap(gid0, *GetComm_tpetra());
  maps[1] = Blocking::buildSubMap(gid1, *GetComm_tpetra());
  maps[2] = Blocking::buildSubMap(gid2, *GetComm_tpetra());

  TEST_MSG("\n   tBlockingTpetra::test_many2one: Building Export/Import");
  std::vector<RCP<Tpetra::Import<LO, GO, NT> > > subImport(3);
  std::vector<RCP<Tpetra::Export<LO, GO, NT> > > subExport(3);
  for (int i = 0; i < 3; i++) {
    Blocking::ImExPair imex = Blocking::buildExportImport(*map, maps[i]);
    subImport[i]            = imex.first;
    subExport[i]            = imex.second;
  }

  TEST_MSG("\n   tBlockingTpetra::test_many2one: Building sub vectors");
  std::vector<RCP<Tpetra::MultiVector<ST, LO, GO, NT> > > subVectors;
  Blocking::buildSubVectors(maps, subVectors, 4);

  TEST_MSG("\n   tBlockingTpetra::test_many2one: Performing one2many");
  Blocking::one2many(subVectors, *v, subImport);

  std::vector<RCP<const Tpetra::MultiVector<ST, LO, GO, NT> > > cSubVectors;
  std::vector<RCP<Tpetra::MultiVector<ST, LO, GO, NT> > >::const_iterator itr;
  for (itr = subVectors.begin(); itr != subVectors.end(); ++itr) cSubVectors.push_back(*itr);

  TEST_MSG("\n   tBlockingTpetra::test_many2one: Performing many2one");
  RCP<Tpetra::MultiVector<ST, LO, GO, NT> > one =
      rcp(new Tpetra::MultiVector<ST, LO, GO, NT>(map, 4));
  one->putScalar(0.0);
  Blocking::many2one(*one, cSubVectors, subExport);
  ST norm2[4] = {0, 0, 0, 0};
  one->norm2(Teuchos::ArrayView<ST>(norm2, 4));

  one->update(1.0, *v, -1.0);

  ST diff[4] = {0, 0, 0, 0};
  ST max = -1.0, maxn = -1, maxn2 = -1;
  ST norm[4] = {0, 0, 0, 0};
  one->norm2(Teuchos::ArrayView<ST>(diff, 4));
  v->norm2(Teuchos::ArrayView<ST>(norm, 4));
  for (int i = 0; i < 4; i++) {
    max   = max > diff[i] / norm[i] ? max : diff[i] / norm[i];
    maxn  = maxn > norm[i] ? maxn : norm[i];
    maxn2 = maxn2 > norm2[i] ? maxn2 : norm2[i];
  }
  TEST_ASSERT(maxn > 0.0, "   tBlockingTpetra::test_many2one maxn>0? maxn = " << maxn);
  TEST_ASSERT(max <= tolerance_, "   tBlockingTpetra::test_many2one ("
                                     << Teko::Test::toString(status) << "): "
                                     << "norm must be better than the tolerance ( " << max
                                     << " <=? " << tolerance_ << " maxn = " << maxn
                                     << ", maxn2 = " << maxn2 << " )");

  return allPassed;
}

bool tBlockingTpetra::test_buildSubBlock(int verbosity, std::ostream& os) {
  bool allPassed = true;

  int numProc   = GetComm_tpetra()->getSize();
  LO mypid      = GetComm_tpetra()->getRank();
  GO myElmts[2] = {2 * mypid + 0, 2 * mypid + 1};
  int grid      = -1;

  // build epetra operator
  /////////////////////////////////////////
  Tpetra::Map<LO, GO, NT> map(-1, Teuchos::ArrayView<GO>(myElmts, 2), 0, GetComm_tpetra());
  std::vector<GO> indices(numProc * 2);
  std::vector<ST> values(numProc * 2);
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > A =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(rcpFromRef(map), numProc * 2);
  for (std::size_t i = 0; i < indices.size(); i++) indices[i] = i;

  // build local row 0
  for (std::size_t i = 0; i < indices.size() / 2; i++) {
    values[2 * i + 0] = (mypid + 1.0) * 1.0 + i;
    values[2 * i + 1] = (mypid + 1.0) * 2.0 + i;
  }
  grid = A->getRowMap()->getGlobalElement(0);
  A->replaceGlobalValues(grid, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(values));

  // build local row 1
  for (std::size_t i = 0; i < indices.size() / 2; i++) {
    values[2 * i + 0] = (mypid + 1.0) * 3.0 + i;
    values[2 * i + 1] = (mypid + 1.0) * 4.0 + i;
  }
  grid = A->getRowMap()->getGlobalElement(1);
  A->replaceGlobalValues(grid, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(values));
  A->fillComplete();

  // build sub maps
  /////////////////////////////////////////
  std::vector<std::vector<GO> > v;
  v.resize(2);
  v[0].push_back(myElmts[0]);
  v[1].push_back(myElmts[1]);
  std::vector<Blocking::MapPair> mapPairs;
  mapPairs.push_back(Blocking::buildSubMap(v[0], *GetComm_tpetra()));
  mapPairs.push_back(Blocking::buildSubMap(v[1], *GetComm_tpetra()));

  return allPassed;
}

}  // namespace Test
}  // namespace Teko
