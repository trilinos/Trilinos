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
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraOperatorWrapper.hpp"
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

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

#include "tBlockingEpetra.hpp"

#include "Teko_BlockingEpetra.hpp"

using namespace Teko::Epetra;

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

void tBlockingEpetra::initializeTest() { tolerance_ = 1e-14; }

void buildGIDs(std::vector<std::vector<int> >& gids, const Epetra_Map& map) {
  int numLocal = map.NumMyElements();
  int numHalf  = numLocal / 2;
  numHalf += ((numHalf % 2 == 0) ? 0 : 1);

  gids.clear();
  gids.resize(3);

  std::vector<int>& blk0 = gids[0];
  std::vector<int>& blk1 = gids[1];
  std::vector<int>& blk2 = gids[2];

  // loop over global IDs: treat first block as strided
  int gid = -1;
  for (int i = 0; i < numHalf; i += 2) {
    gid = map.GID(i);
    blk0.push_back(gid);

    gid = map.GID(i + 1);
    blk1.push_back(gid);
  }

  // loop over global IDs: treat remainder as contiguous
  for (int i = numHalf; i < numLocal; i++) {
    gid = map.GID(i);
    blk2.push_back(gid);
  }
}

int tBlockingEpetra::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                             int& totalrun) {
  bool allTests = true;
  bool status   = true;
  int failcount = 0;

  failstrm << "tBlockingEpetra";

  status = test_buildMaps(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"buildMaps\" ... PASSED", "   \"buildMaps\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_one2many(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"one2many\" ... PASSED", "   \"one2many\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_many2one(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"many2one\" ... PASSED", "   \"many2one\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_buildSubBlock(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"buildSubBlock\" ... PASSED", "   \"buildSubBlock\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tBlockingEpetra...PASSED", "tBlockingEpetra...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tBlockingEpetra...FAILED");
  }

  return failcount;
}

bool tBlockingEpetra::test_buildMaps(int verbosity, std::ostream& os) {
  using Teuchos::as;

  bool status    = false;
  bool allPassed = true;

  int size = 30;

  TEST_MSG("\n   Builing Epetra_Map");
  RCP<Epetra_Map> map = rcp(new Epetra_Map(size, 0, *GetComm()));

  // build gids
  std::vector<std::vector<int> > gids;
  buildGIDs(gids, *map);

  std::vector<int>& gid0 = gids[0];
  std::vector<int>& gid1 = gids[1];
  std::vector<int>& gid2 = gids[2];

  TEST_MSG("\n   Building sub maps");
  Blocking::MapPair map0 = Blocking::buildSubMap(gid0, *GetComm());
  Blocking::MapPair map1 = Blocking::buildSubMap(gid1, *GetComm());
  Blocking::MapPair map2 = Blocking::buildSubMap(gid2, *GetComm());

  TEST_ASSERT(map0.first->NumMyElements() == as<int>(gid0.size()) &&
                  map0.second->NumMyElements() == as<int>(gid0.size()),
              "   tBlockingEpetra::test_buildMaps ("
                  << Teko::Test::toString(status) << "): "
                  << " Checking map size: first=" << map0.first->NumMyElements()
                  << ", second=" << map0.second->NumMyElements() << ", gid=" << gid0.size());
  TEST_ASSERT(map1.first->NumMyElements() == as<int>(gid1.size()) &&
                  map0.second->NumMyElements() == as<int>(gid0.size()),
              "   tBlockingEpetra::test_buildMaps ("
                  << Teko::Test::toString(status) << "): "
                  << " Checking map size: first=" << map1.first->NumMyElements()
                  << ", second=" << map1.second->NumMyElements() << ", gid=" << gid1.size());
  TEST_ASSERT(map2.first->NumMyElements() == as<int>(gid2.size()) &&
                  map0.second->NumMyElements() == as<int>(gid0.size()),
              "   tBlockingEpetra::test_buildMaps ("
                  << Teko::Test::toString(status) << "): "
                  << " Checking map size: first=" << map2.first->NumMyElements()
                  << ", second=" << map2.second->NumMyElements() << ", gid=" << gid2.size());

  std::vector<Teuchos::RCP<Epetra_Map> > globalMaps(3);
  std::vector<Teuchos::RCP<Epetra_Map> > contigMaps(3);

  // get sub maps for convenient use and access
  globalMaps[0] = map0.first;
  globalMaps[1] = map1.first;
  globalMaps[2] = map2.first;

  contigMaps[0] = map0.second;
  contigMaps[1] = map1.second;
  contigMaps[2] = map2.second;

  // test that the extra data is attached
  TEST_ASSERT(contigMaps[0] != Teuchos::null,
              "   tBlockingEpetra::test_buildMaps (" << Teko::Test::toString(status) << ")");
  TEST_ASSERT(contigMaps[1] != Teuchos::null,
              "   tBlockingEpetra::test_buildMaps (" << Teko::Test::toString(status) << ")");
  TEST_ASSERT(contigMaps[2] != Teuchos::null,
              "   tBlockingEpetra::test_buildMaps (" << Teko::Test::toString(status) << ")");
  TEST_MSG("   tBlockingEpetra::test_buildMaps: extracted \"contigMaps\"");

  bool test;

  // check contiguous and global maps
  test = true;
  for (int i = 0; i < globalMaps[0]->NumMyElements(); i++) {
    int gid = globalMaps[0]->GID(i);
    int cid = contigMaps[0]->GID(i);

    test &= gid == gid0[i];
    test &= cid == (i + contigMaps[0]->MinMyGID());
  }
  TEST_ASSERT(test, "   tBlockingEpetra::test_buildMaps ("
                        << Teko::Test::toString(status) << "): "
                        << "checked that block maps were internally consitent");

  test = true;
  for (int i = 0; i < globalMaps[1]->NumMyElements(); i++) {
    int gid = globalMaps[1]->GID(i);
    int cid = contigMaps[1]->GID(i);

    test &= gid == gid1[i];
    test &= cid == (i + contigMaps[1]->MinMyGID());
  }
  TEST_ASSERT(test, "   tBlockingEpetra::test_buildMaps ("
                        << Teko::Test::toString(status) << "): "
                        << "checked that block maps were internally consitent");

  test = true;
  for (int i = 0; i < globalMaps[2]->NumMyElements(); i++) {
    int gid = globalMaps[2]->GID(i);
    int cid = contigMaps[2]->GID(i);

    test &= gid == gid2[i];
    test &= cid == (i + contigMaps[2]->MinMyGID());
  }
  TEST_ASSERT(test, "   tBlockingEpetra::test_buildMaps ("
                        << Teko::Test::toString(status) << "): "
                        << "checked that block maps were internally consitent");

  return allPassed;
}

bool tBlockingEpetra::test_one2many(int verbosity, std::ostream& os) {
  bool allPassed = true;

  int size = 3 * 1000;
  TEST_MSG("\n   tBlockingEpetra::test_one2many: Builing Epetra_Map and source vector");
  RCP<Epetra_Map> map       = rcp(new Epetra_Map(size, 0, *GetComm()));
  RCP<Epetra_MultiVector> v = rcp(new Epetra_MultiVector(*map, 1));
  v->Random();

  TEST_MSG("\n   Building sub maps");
  // build gids
  std::vector<std::vector<int> > gids;
  buildGIDs(gids, *map);

  std::vector<int>& gid0 = gids[0];
  std::vector<int>& gid1 = gids[1];
  std::vector<int>& gid2 = gids[2];

  std::vector<Blocking::MapPair> maps(3);
  maps[0] = Blocking::buildSubMap(gid0, *GetComm());
  maps[1] = Blocking::buildSubMap(gid1, *GetComm());
  maps[2] = Blocking::buildSubMap(gid2, *GetComm());

  TEST_MSG("\n   tBlockingEpetra::test_one2many: Building Export/Import");
  std::vector<RCP<Epetra_Import> > subImport(3);
  std::vector<RCP<Epetra_Export> > subExport(3);
  for (int i = 0; i < 3; i++) {
    Blocking::ImExPair imex = Blocking::buildExportImport(*map, maps[i]);
    subImport[i]            = imex.first;
    subExport[i]            = imex.second;
  }

  TEST_MSG("\n   tBlockingEpetra::test_one2many: Building sub vectors");
  std::vector<RCP<Epetra_MultiVector> > subVectors;
  Blocking::buildSubVectors(maps, subVectors, 1);

  TEST_MSG("\n   tBlockingEpetra::test_one2many: Performing data copy");
  Blocking::one2many(subVectors, *v, subImport);

  // just assume it works! :)

  return allPassed;
}

bool tBlockingEpetra::test_many2one(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  int size = 3 * 1000;
  TEST_MSG("\n   tBlockingEpetra::test_many2one: Builing Epetra_Map and source vector");
  RCP<Epetra_Map> map       = rcp(new Epetra_Map(size, 0, *GetComm()));
  RCP<Epetra_MultiVector> v = rcp(new Epetra_MultiVector(*map, 4));
  v->Random();

  TEST_MSG("\n   Building sub maps");
  // build gids
  std::vector<std::vector<int> > gids;
  buildGIDs(gids, *map);

  std::vector<int>& gid0 = gids[0];
  std::vector<int>& gid1 = gids[1];
  std::vector<int>& gid2 = gids[2];

  std::vector<Blocking::MapPair> maps(3);
  maps[0] = Blocking::buildSubMap(gid0, *GetComm());
  maps[1] = Blocking::buildSubMap(gid1, *GetComm());
  maps[2] = Blocking::buildSubMap(gid2, *GetComm());

  TEST_MSG("\n   tBlockingEpetra::test_many2one: Building Export/Import");
  std::vector<RCP<Epetra_Import> > subImport(3);
  std::vector<RCP<Epetra_Export> > subExport(3);
  for (int i = 0; i < 3; i++) {
    Blocking::ImExPair imex = Blocking::buildExportImport(*map, maps[i]);
    subImport[i]            = imex.first;
    subExport[i]            = imex.second;
  }

  TEST_MSG("\n   tBlockingEpetra::test_many2one: Building sub vectors");
  std::vector<RCP<Epetra_MultiVector> > subVectors;
  Blocking::buildSubVectors(maps, subVectors, 4);

  TEST_MSG("\n   tBlockingEpetra::test_many2one: Performing one2many");
  Blocking::one2many(subVectors, *v, subImport);

  std::vector<RCP<const Epetra_MultiVector> > cSubVectors;
  std::vector<RCP<Epetra_MultiVector> >::const_iterator itr;
  for (itr = subVectors.begin(); itr != subVectors.end(); ++itr) cSubVectors.push_back(*itr);

  TEST_MSG("\n   tBlockingEpetra::test_many2one: Performing many2one");
  RCP<Epetra_MultiVector> one = rcp(new Epetra_MultiVector(*map, 4));
  one->PutScalar(0.0);
  Blocking::many2one(*one, cSubVectors, subExport);
  double norm2[4] = {0, 0, 0, 0};
  one->Norm2(norm2);

  one->Update(1.0, *v, -1.0);

  double diff[4] = {0, 0, 0, 0};
  double max = -1.0, maxn = -1, maxn2 = -1;
  double norm[4] = {0, 0, 0, 0};
  one->Norm2(diff);
  v->Norm2(norm);
  for (int i = 0; i < 4; i++) {
    max   = max > diff[i] / norm[i] ? max : diff[i] / norm[i];
    maxn  = maxn > norm[i] ? maxn : norm[i];
    maxn2 = maxn2 > norm2[i] ? maxn2 : norm2[i];
  }
  TEST_ASSERT(maxn > 0.0, "   tBlockingEpetra::test_many2one maxn>0? maxn = " << maxn);
  TEST_ASSERT(max <= tolerance_, "   tBlockingEpetra::test_many2one ("
                                     << Teko::Test::toString(status) << "): "
                                     << "norm must be better than the tolerance ( " << max
                                     << " <=? " << tolerance_ << " maxn = " << maxn
                                     << ", maxn2 = " << maxn2 << " )");

  return allPassed;
}

bool tBlockingEpetra::test_buildSubBlock(int verbosity, std::ostream& os) {
  bool allPassed = true;

  int numProc    = GetComm()->NumProc();
  int mypid      = GetComm()->MyPID();
  int myElmts[2] = {2 * mypid + 0, 2 * mypid + 1};
  int grid       = -1;

  // build epetra operator
  /////////////////////////////////////////
  Epetra_Map map(-1, 2, myElmts, 0, *GetComm());
  std::vector<int> indices(numProc * 2);
  std::vector<double> values(numProc * 2);
  Epetra_CrsMatrix A(Copy, map, numProc * 2, false);
  for (std::size_t i = 0; i < indices.size(); i++) indices[i] = i;

  // build local row 0
  for (std::size_t i = 0; i < indices.size() / 2; i++) {
    values[2 * i + 0] = (mypid + 1.0) * 1.0 + i;
    values[2 * i + 1] = (mypid + 1.0) * 2.0 + i;
  }
  grid = A.GRID(0);
  A.InsertGlobalValues(grid, values.size(), &values[0], &indices[0]);

  // build local row 1
  for (std::size_t i = 0; i < indices.size() / 2; i++) {
    values[2 * i + 0] = (mypid + 1.0) * 3.0 + i;
    values[2 * i + 1] = (mypid + 1.0) * 4.0 + i;
  }
  grid = A.GRID(1);
  A.InsertGlobalValues(grid, values.size(), &values[0], &indices[0]);
  A.FillComplete();

  // build sub maps
  /////////////////////////////////////////
  std::vector<std::vector<int> > v;
  v.resize(2);
  v[0].push_back(myElmts[0]);
  v[1].push_back(myElmts[1]);
  std::vector<Blocking::MapPair> mapPairs;
  mapPairs.push_back(Blocking::buildSubMap(v[0], *GetComm()));
  mapPairs.push_back(Blocking::buildSubMap(v[1], *GetComm()));

  return allPassed;
}

}  // namespace Test
}  // namespace Teko
