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

#include "tInterlacedTpetra.hpp"

#include "Teko_InterlacedTpetra.hpp"

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

void tInterlacedTpetra::initializeTest() { tolerance_ = 1e-14; }

int tInterlacedTpetra::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                               int& totalrun) {
  bool allTests = true;
  bool status   = true;
  int failcount = 0;

  failstrm << "tInterlacedTpetra";

  status = test_buildSubMaps_num(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"buildSubMaps_num\" ... PASSED",
                       "   \"buildSubMaps_num\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_buildSubMaps_vec(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"buildSubMaps_vec\" ... PASSED",
                       "   \"buildSubMaps_vec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

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

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tInterlacedTpetra...PASSED", "tInterlacedTpetra...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tInterlacedTpetra...FAILED");
  }

  return failcount;
}

bool tInterlacedTpetra::test_buildSubMaps_num(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  const Teuchos::Comm<int>& comm = *GetComm_tpetra();

  try {
    std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > > subMaps;
    GO globals  = 10;
    int numVars = 3;

    // build a set of submaps: this should fail
    Strided::buildSubMaps(globals, numVars, comm, subMaps);

    TEST_ASSERT(false,
                "\n   tInterlacedTpetra::test_buildSubMaps_num: "
                    << toString(status)
                    << ": "
                       "buildSubMaps(int,vector<pair<int,RCP<Tpetra::Map> > >) did not throw "
                       "with incorrect parameters");
  } catch (...) {
    TEST_MSG("\n   tInterlacedTpetra::test_buildSubMaps_num: "
             << "correctly threw an exception on incorrect parameters");
  }

  try {
    std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > > subMaps;
    GO globals  = 9;
    int numVars = 3;

    // build a set of submaps: this should fail
    Strided::buildSubMaps(globals, numVars, comm, subMaps);

    TEST_EQUALITY(subMaps.size(), 3,
                  "\n   tInterlacedTpetra::test_buildSubMaps_num: "
                      << toString(status) << ": "
                      << "testing number of maps built ( " << subMaps.size() << " == " << 3
                      << "? ) ");

    bool cur = true;
    for (int i = 0; i < 3; i++) cur &= (subMaps[i].first == 1);
    TEST_ASSERT(cur,
                "\n   tInterlacedTpetra::test_buildSubMaps_num: "
                    << toString(status) << ": "
                    << "testing that maps are associated with the correct numbers of variables");

    Teuchos::ArrayView<const GO> gids;

    // test the first of three used maps
    TEST_EQUALITY(subMaps[0].second->getGlobalNumElements(), 3,
                  "\n   tInterlacedTpetra::test_buildSubMaps_num: "
                      << toString(status) << ": "
                      << "testing that first map has correct number of global elements ( "
                      << subMaps[0].second->getGlobalNumElements() << " == " << 3 << " ?)");
    gids = subMaps[0].second->getLocalElementList();
    cur  = (gids[0] == 0 && gids[1] == 3 && gids[2] == 6);
    TEST_ASSERT(cur, "\n   tInterlacedTpetra::test_buildSubMaps_num: "
                         << toString(status) << ": "
                         << "testing that first map is created correctly");

    // test the second of three used maps
    TEST_EQUALITY(subMaps[1].second->getGlobalNumElements(), 3,
                  "\n   tInterlacedTpetra::test_buildSubMaps_num: "
                      << toString(status) << ": "
                      << "testing that second map has correct number of global elements ( "
                      << subMaps[1].second->getGlobalNumElements() << " == " << 3 << " ?)");
    gids = subMaps[1].second->getLocalElementList();
    cur  = (gids[0] == 1 && gids[1] == 4 && gids[2] == 7);
    TEST_ASSERT(cur, "\n   tInterlacedTpetra::test_buildSubMaps_num: "
                         << toString(status) << ": "
                         << "testing that second map is created correctly");

    // test the first of three used maps
    TEST_EQUALITY(subMaps[2].second->getGlobalNumElements(), 3,
                  "\n   tInterlacedTpetra::test_buildSubMaps_num: "
                      << toString(status) << ": "
                      << "testing that the third map has correct number of global elements ( "
                      << subMaps[2].second->getGlobalNumElements() << " == " << 3 << " ?)");
    gids = subMaps[2].second->getLocalElementList();
    cur  = (gids[0] == 2 && gids[1] == 5 && gids[2] == 8);
    TEST_ASSERT(cur, "\n   tInterlacedTpetra::test_buildSubMaps_num: "
                         << toString(status) << ": "
                         << "testing that third map is created correctly");
  } catch (...) {
    TEST_ASSERT(false,
                "\n   tInterlacedTpetra::test_buildSubMaps_num: "
                    << toString(status)
                    << ": "
                       "buildSubMaps(int,vector<pair<int,RCP<Epetra_Map> > >) threw an unexpected "
                       "exception");
  }

  return allPassed;
}

bool tInterlacedTpetra::test_buildSubMaps_vec(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  const Teuchos::Comm<int>& comm = *GetComm_tpetra();

  try {
    std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > > subMaps;
    GO globals = 15;

    std::vector<int> vars(3);
    vars[0] = 2;
    vars[1] = 1;
    vars[2] = 3;

    // build a set of submaps: this should fail
    Strided::buildSubMaps(globals, vars, comm, subMaps);

    TEST_ASSERT(false,
                "\n   tInterlacedTpetra::test_buildSubMaps_vec: "
                    << toString(status)
                    << ": "
                       "buildSubMaps(int,vector<pair<int,RCP<Tpetra::Map> > >) did not throw "
                       "with incorrect parameters");
  } catch (...) {
    TEST_MSG("\n   tInterlacedTpetra::test_buildSubMaps_vec: "
             << "correctly threw an exception on incorrect parameters");
  }

  try {
    std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > > subMaps;
    GO globals = 18;

    std::vector<int> vars(3);
    vars[0] = 2;
    vars[1] = 1;
    vars[2] = 3;

    // build a set of submaps: this should fail
    Strided::buildSubMaps(globals, vars, comm, subMaps);

    TEST_EQUALITY(subMaps.size(), 3,
                  "\n   tInterlacedTpetra::test_buildSubMaps_vec: "
                      << toString(status) << ": "
                      << "testing number of maps built ( " << subMaps.size() << " == " << 3
                      << "? ) ");

    bool cur = true;
    for (int i = 0; i < 3; i++) cur &= (subMaps[i].first == vars[i]);
    TEST_ASSERT(cur,
                "\n   tInterlacedTpetra::test_buildSubMaps_vec: "
                    << toString(status) << ": "
                    << "testing that maps are associated with the correct numbers of variables");

    Teuchos::ArrayView<const GO> gids;

    // test the first of three used maps
    TEST_EQUALITY(subMaps[0].second->getGlobalNumElements(), 6,
                  "\n   tInterlacedTpetra::test_buildSubMaps_vec: "
                      << toString(status) << ": "
                      << "testing that first map has correct number of global elements ( "
                      << subMaps[0].second->getGlobalNumElements() << " == " << 6 << " ?)");
    gids = subMaps[0].second->getLocalElementList();
    cur  = (gids[0] == 0 && gids[1] == 1 && gids[2] == 6 && gids[3] == 7 && gids[4] == 12 &&
           gids[5] == 13);
    TEST_ASSERT(cur, "\n   tInterlacedTpetra::test_buildSubMaps_vec: "
                         << toString(status) << ": "
                         << "testing that first map is created correctly");

    // test the second of three used maps
    TEST_EQUALITY(subMaps[1].second->getGlobalNumElements(), 3,
                  "\n   tInterlacedTpetra::test_buildSubMaps_vec: "
                      << toString(status) << ": "
                      << "testing that second map has correct number of global elements ( "
                      << subMaps[1].second->getGlobalNumElements() << " == " << 3 << " ?)");
    gids = subMaps[1].second->getLocalElementList();
    cur  = (gids[0] == 2 && gids[1] == 8 && gids[2] == 14);
    TEST_ASSERT(cur, "\n   tInterlacedTpetra::test_buildSubMaps_vec: "
                         << toString(status) << ": "
                         << "testing that second map is created correctly");

    // test the first of three used maps
    TEST_EQUALITY(subMaps[2].second->getGlobalNumElements(), 9,
                  "\n   tInterlacedTpetra::test_buildSubMaps_vec: "
                      << toString(status) << ": "
                      << "testing that the third map has correct number of global elements ( "
                      << subMaps[2].second->getGlobalNumElements() << " == " << 9 << " ?)");
    gids = subMaps[2].second->getLocalElementList();
    cur  = (gids[0] == 3 && gids[1] == 4 && gids[2] == 5 && gids[3] == 9 && gids[4] == 10 &&
           gids[5] == 11 && gids[6] == 15 && gids[7] == 16 && gids[8] == 17);
    TEST_ASSERT(cur, "\n   tInterlacedTpetra::test_buildSubMaps_vec: "
                         << toString(status) << ": "
                         << "testing that third map is created correctly");
  } catch (...) {
    TEST_ASSERT(false,
                "\n   tInterlacedTpetra::test_buildSubMaps_vec: "
                    << toString(status)
                    << ": "
                       "buildSubMaps(int,vector<pair<int,RCP<Epetra_Map> > >) threw an unexpected "
                       "exception");
  }

  return allPassed;
}

bool tInterlacedTpetra::test_buildMaps(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  GO size = 3 * 1000;

  TEST_MSG("\n   Builing Tpetra::Map");
  RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(size, 0, GetComm_tpetra()));

  TEST_MSG("\n   Building sub maps");
  std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > > subMaps;
  std::vector<int> vec(2);
  vec[0] = 2;
  vec[1] = 1;
  Strided::buildSubMaps(size, vec, *GetComm_tpetra(), subMaps);

  TEST_ASSERT(subMaps[0].first == vec[0], "   tInterlacedTpetra::test_buildMaps ("
                                              << Teko::Test::toString(status) << "): "
                                              << "  first map unknowns is " << subMaps[0].first
                                              << " ( should be " << vec[0] << ")");
  TEST_ASSERT(subMaps[1].first == vec[1], "   tInterlacedTpetra::test_buildMaps ("
                                              << Teko::Test::toString(status) << "): "
                                              << "  second map unknowns is " << subMaps[1].first
                                              << " ( should be " << vec[1] << ")");

  std::vector<Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > globalMaps(2);
  std::vector<Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > contigMaps(2);

  // get sub maps for convenient use and access
  globalMaps[0] = subMaps[0].second;
  globalMaps[1] = subMaps[1].second;

  contigMaps[0] =
      Teuchos::get_extra_data<Teuchos::RCP<Tpetra::Map<LO, GO, NT> > >(globalMaps[0], "contigMap");
  contigMaps[1] =
      Teuchos::get_extra_data<Teuchos::RCP<Tpetra::Map<LO, GO, NT> > >(globalMaps[1], "contigMap");

  // test that the extra data is attached
  TEST_ASSERT(contigMaps[0] != Teuchos::null,
              "   tInterlacedTpetra::test_buildMaps (" << Teko::Test::toString(status) << ")");
  TEST_ASSERT(contigMaps[1] != Teuchos::null,
              "   tInterlacedTpetra::test_buildMaps (" << Teko::Test::toString(status) << ")");
  TEST_MSG("   tInterlacedTpetra::test_buildMaps: extracted \"contigMaps\" from RCP");

  // make sure all maps have the correct size
  TEST_ASSERT(globalMaps[0]->getGlobalNumElements() == 2000,
              "   tInterlacedTpetra::test_buildMaps (" << Teko::Test::toString(status) << "): "
                                                       << "checking number of global elements");
  TEST_ASSERT(globalMaps[1]->getGlobalNumElements() == 1000,
              "   tInterlacedTpetra::test_buildMaps (" << Teko::Test::toString(status) << "): "
                                                       << "checking number of global elements");
  TEST_ASSERT(contigMaps[0]->getGlobalNumElements() == 2000,
              "   tInterlacedTpetra::test_buildMaps (" << Teko::Test::toString(status) << "): "
                                                       << "checking number of global elements");
  TEST_ASSERT(contigMaps[1]->getGlobalNumElements() == 1000,
              "   tInterlacedTpetra::test_buildMaps (" << Teko::Test::toString(status) << "): "
                                                       << "checking number of global elements");

  TEST_ASSERT(contigMaps[0]->getLocalNumElements() == globalMaps[0]->getLocalNumElements(),
              "   tInterlacedTpetra::test_buildMaps ("
                  << Teko::Test::toString(status) << "): "
                  << "check for lineup of number of local elements");
  TEST_ASSERT(contigMaps[1]->getLocalNumElements() == globalMaps[1]->getLocalNumElements(),
              "   tInterlacedTpetra::test_buildMaps ("
                  << Teko::Test::toString(status) << "): "
                  << "check for lineup of number of local elements");

  bool test;

  // check contiguous and global maps
  test = true;
  for (size_t i = 0; i < globalMaps[0]->getLocalNumElements(); i += 2) {
    int block;
    GO gid = globalMaps[0]->getGlobalElement(i);
    GO cid = contigMaps[0]->getGlobalElement(i);

    block = gid / 3;
    test &= cid == 2 * block;

    GO gidp1 = globalMaps[0]->getGlobalElement(i + 1);
    GO cidp1 = contigMaps[0]->getGlobalElement(i + 1);

    block = (gidp1 - 1) / 3;
    test &= cidp1 == 2 * block + 1;
  }
  TEST_ASSERT(test, "   tInterlacedTpetra::test_buildMaps ("
                        << Teko::Test::toString(status) << "): "
                        << "checked that block maps were internally consitent");

  test = true;
  for (size_t i = 0; i < globalMaps[1]->getLocalNumElements(); i++) {
    int block;
    GO gid = globalMaps[1]->getGlobalElement(i);
    GO cid = contigMaps[1]->getGlobalElement(i);

    block = (gid - 2) / 3;
    test &= cid == block;
  }
  TEST_ASSERT(test, "   tInterlacedTpetra::test_buildMaps ("
                        << Teko::Test::toString(status) << "): "
                        << "checked that block maps were internally consitent");

  return allPassed;
}

bool tInterlacedTpetra::test_one2many(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  GO size = 3 * 1000;
  TEST_MSG("\n   tInterlacedTpetra::test_one2many: Builing Epetra_Map and source vector");
  RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(size, 0, GetComm_tpetra()));
  RCP<Tpetra::MultiVector<ST, LO, GO, NT> > v =
      rcp(new Tpetra::MultiVector<ST, LO, GO, NT>(map, 1));
  v->randomize();

  TEST_MSG("\n   tInterlacedTpetra::test_one2many: Building sub maps");
  std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > > subMaps;
  std::vector<int> vec(2);
  vec[0] = 2;
  vec[1] = 1;
  Strided::buildSubMaps(size, vec, *GetComm_tpetra(), subMaps);

  TEST_MSG("\n   tInterlacedTpetra::test_one2many: Building Export/Import");
  std::vector<RCP<Tpetra::Export<LO, GO, NT> > > subExport;
  std::vector<RCP<Tpetra::Import<LO, GO, NT> > > subImport;
  Strided::buildExportImport(*map, subMaps, subExport, subImport);

  TEST_MSG("\n   tInterlacedTpetra::test_one2many: Building sub vectors");
  std::vector<RCP<Tpetra::MultiVector<ST, LO, GO, NT> > > subVectors;
  Strided::buildSubVectors(subMaps, subVectors, 1);

  std::vector<RCP<Tpetra::MultiVector<ST, LO, GO, NT> > >::const_iterator itr;
  for (itr = subVectors.begin(); itr != subVectors.end(); ++itr) {
    RCP<const Tpetra::Map<LO, GO, NT> > lm =
        Teuchos::get_extra_data<RCP<Tpetra::Map<LO, GO, NT> > >(*itr, "globalMap");
    TEST_ASSERT(lm != Teuchos::null, "   tInterlacedTpetra::test_buildMaps ("
                                         << Teko::Test::toString(status) << "): "
                                         << "check that vector contains \"globalMap\" in RCP");
  }

  TEST_MSG("\n   tInterlacedTpetra::test_one2many: Performing data copy");
  Strided::one2many(subVectors, *v, subImport);

  // just assume it works! :)

  return allPassed;
}

bool tInterlacedTpetra::test_many2one(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  GO size = 3 * 1000;
  TEST_MSG("\n   tInterlacedTpetra::test_one2many: Builing Tpetra_Map and source vector");
  RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(size, 0, GetComm_tpetra()));
  RCP<Tpetra::MultiVector<ST, LO, GO, NT> > v =
      rcp(new Tpetra::MultiVector<ST, LO, GO, NT>(map, 4));
  v->randomize();

  TEST_MSG("\n   tInterlacedTpetra::test_one2many: Building sub maps");
  std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > > subMaps;
  std::vector<int> vec(2);
  vec[0] = 2;
  vec[1] = 1;
  Strided::buildSubMaps(size, vec, *GetComm_tpetra(), subMaps);

  TEST_MSG("\n   tInterlacedTpetra::test_one2many: Building Export/Import");
  std::vector<RCP<Tpetra::Export<LO, GO, NT> > > subExport;
  std::vector<RCP<Tpetra::Import<LO, GO, NT> > > subImport;
  Strided::buildExportImport(*map, subMaps, subExport, subImport);

  TEST_MSG("\n   tInterlacedTpetra::test_one2many: Building sub vectors");
  std::vector<RCP<Tpetra::MultiVector<ST, LO, GO, NT> > > subVectors;
  Strided::buildSubVectors(subMaps, subVectors, 4);

  std::vector<RCP<Tpetra::MultiVector<ST, LO, GO, NT> > >::const_iterator itr;
  for (itr = subVectors.begin(); itr != subVectors.end(); ++itr) {
    RCP<const Tpetra::Map<LO, GO, NT> > lm =
        Teuchos::get_extra_data<RCP<Tpetra::Map<LO, GO, NT> > >(*itr, "globalMap");
    TEST_ASSERT(lm != Teuchos::null, "   tInterlacedTpetra::test_buildMaps ("
                                         << Teko::Test::toString(status) << "): "
                                         << "check that vector contains \"globalMap\" in RCP");
  }

  TEST_MSG("\n   tInterlacedTpetra::test_one2many: Performing one2many");
  Strided::one2many(subVectors, *v, subImport);

  std::vector<RCP<const Tpetra::MultiVector<ST, LO, GO, NT> > > cSubVectors;
  for (itr = subVectors.begin(); itr != subVectors.end(); ++itr) cSubVectors.push_back(*itr);

  TEST_MSG("\n   tInterlacedTpetra::test_one2many: Performing many2one");
  RCP<Tpetra::MultiVector<ST, LO, GO, NT> > one =
      rcp(new Tpetra::MultiVector<ST, LO, GO, NT>(map, 4));
  Strided::many2one(*one, cSubVectors, subExport);

  one->update(1.0, *v, -1.0);

  ST diff[4] = {0, 0, 0, 0};
  ST max = 0.0, maxn = 0;
  ST norm[4] = {0, 0, 0, 0};
  one->norm2(Teuchos::ArrayView<ST>(diff, 4));
  v->norm2(Teuchos::ArrayView<ST>(norm, 4));
  for (int i = 0; i < 4; i++) {
    max  = max > diff[i] / norm[i] ? max : diff[i] / norm[i];
    maxn = maxn > norm[i] ? maxn : norm[i];
  }
  TEST_ASSERT(max <= tolerance_, "   tInterlacedTpetra::test_buildMaps ("
                                     << Teko::Test::toString(status) << "): "
                                     << "norm must be better than the tolerance ( " << max
                                     << " <=? " << tolerance_ << " maxn = " << maxn << " )");

  return allPassed;
}

}  // namespace Test
}  // namespace Teko
