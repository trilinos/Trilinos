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

#include "tInterlacedEpetra.hpp"

#include "Teko_InterlacedEpetra.hpp"

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

void tInterlacedEpetra::initializeTest() { tolerance_ = 1e-14; }

int tInterlacedEpetra::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                               int& totalrun) {
  bool allTests = true;
  bool status   = true;
  int failcount = 0;

  failstrm << "tInterlacedEpetra";

  status = test_buildSubMaps_num(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"buildSubMaps_num\" ... PASSED",
                "   \"buildSubMaps_num\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_buildSubMaps_vec(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"buildSubMaps_vec\" ... PASSED",
                "   \"buildSubMaps_vec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

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

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tInterlacedEpetra...PASSED", "tInterlacedEpetra...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tInterlacedEpetra...FAILED");
  }

  return failcount;
}

bool tInterlacedEpetra::test_buildSubMaps_num(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  const Epetra_Comm& comm = *GetComm();

  try {
    std::vector<std::pair<int, RCP<Epetra_Map> > > subMaps;
    int globals = 10;
    int numVars = 3;

    // build a set of submaps: this should fail
    Strided::buildSubMaps(globals, numVars, comm, subMaps);

    TEST_ASSERT(false, "\n   tInerlacedEpetra::test_buildSubMaps_num: "
                           << toString(status)
                           << ": "
                              "buildSubMaps(int,vector<pair<int,RCP<Epetra_Map> > >) did not throw "
                              "with incorrect parameters");
  } catch (...) {
    TEST_MSG("\n   tInerlacedEpetra::test_buildSubMaps_num: "
             << "correctly threw an exception on incorrect parameters");
  }

  try {
    std::vector<std::pair<int, RCP<Epetra_Map> > > subMaps;
    int globals = 9;
    int numVars = 3;

    // build a set of submaps: this should fail
    Strided::buildSubMaps(globals, numVars, comm, subMaps);

    TEST_EQUALITY(subMaps.size(), 3,
                  "\n   tInerlacedEpetra::test_buildSubMaps_num: "
                      << toString(status) << ": "
                      << "testing number of maps built ( " << subMaps.size() << " == " << 3
                      << "? ) ");

    bool cur = true;
    for (int i = 0; i < 3; i++) cur &= (subMaps[i].first == 1);
    TEST_ASSERT(cur,
                "\n   tInerlacedEpetra::test_buildSubMaps_num: "
                    << toString(status) << ": "
                    << "testing that maps are associated with the correct numbers of variables");

    int* gids;

    // test the first of three used maps
    TEST_EQUALITY(subMaps[0].second->NumGlobalElements(), 3,
                  "\n   tInerlacedEpetra::test_buildSubMaps_num: "
                      << toString(status) << ": "
                      << "testing that first map has correct number of global elements ( "
                      << subMaps[0].second->NumGlobalElements() << " == " << 3 << " ?)");
    gids = subMaps[0].second->MyGlobalElements();
    cur  = (gids[0] == 0 && gids[1] == 3 && gids[2] == 6);
    TEST_ASSERT(cur, "\n   tInerlacedEpetra::test_buildSubMaps_num: "
                         << toString(status) << ": "
                         << "testing that first map is created correctly");

    // test the second of three used maps
    TEST_EQUALITY(subMaps[1].second->NumGlobalElements(), 3,
                  "\n   tInerlacedEpetra::test_buildSubMaps_num: "
                      << toString(status) << ": "
                      << "testing that second map has correct number of global elements ( "
                      << subMaps[1].second->NumGlobalElements() << " == " << 3 << " ?)");
    gids = subMaps[1].second->MyGlobalElements();
    cur  = (gids[0] == 1 && gids[1] == 4 && gids[2] == 7);
    TEST_ASSERT(cur, "\n   tInerlacedEpetra::test_buildSubMaps_num: "
                         << toString(status) << ": "
                         << "testing that second map is created correctly");

    // test the first of three used maps
    TEST_EQUALITY(subMaps[2].second->NumGlobalElements(), 3,
                  "\n   tInerlacedEpetra::test_buildSubMaps_num: "
                      << toString(status) << ": "
                      << "testing that the third map has correct number of global elements ( "
                      << subMaps[2].second->NumGlobalElements() << " == " << 3 << " ?)");
    gids = subMaps[2].second->MyGlobalElements();
    cur  = (gids[0] == 2 && gids[1] == 5 && gids[2] == 8);
    TEST_ASSERT(cur, "\n   tInerlacedEpetra::test_buildSubMaps_num: "
                         << toString(status) << ": "
                         << "testing that third map is created correctly");
  } catch (...) {
    TEST_ASSERT(false,
                "\n   tInerlacedEpetra::test_buildSubMaps_num: "
                    << toString(status)
                    << ": "
                       "buildSubMaps(int,vector<pair<int,RCP<Epetra_Map> > >) threw an unexpected "
                       "exception");
  }

  return allPassed;
}

bool tInterlacedEpetra::test_buildSubMaps_vec(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  const Epetra_Comm& comm = *GetComm();

  try {
    std::vector<std::pair<int, RCP<Epetra_Map> > > subMaps;
    int globals = 15;

    std::vector<int> vars(3);
    vars[0] = 2;
    vars[1] = 1;
    vars[2] = 3;

    // build a set of submaps: this should fail
    Strided::buildSubMaps(globals, vars, comm, subMaps);

    TEST_ASSERT(false, "\n   tInerlacedEpetra::test_buildSubMaps_vec: "
                           << toString(status)
                           << ": "
                              "buildSubMaps(int,vector<pair<int,RCP<Epetra_Map> > >) did not throw "
                              "with incorrect parameters");
  } catch (...) {
    TEST_MSG("\n   tInerlacedEpetra::test_buildSubMaps_vec: "
             << "correctly threw an exception on incorrect parameters");
  }

  try {
    std::vector<std::pair<int, RCP<Epetra_Map> > > subMaps;
    int globals = 18;

    std::vector<int> vars(3);
    vars[0] = 2;
    vars[1] = 1;
    vars[2] = 3;

    // build a set of submaps: this should fail
    Strided::buildSubMaps(globals, vars, comm, subMaps);

    TEST_EQUALITY(subMaps.size(), 3,
                  "\n   tInerlacedEpetra::test_buildSubMaps_vec: "
                      << toString(status) << ": "
                      << "testing number of maps built ( " << subMaps.size() << " == " << 3
                      << "? ) ");

    bool cur = true;
    for (int i = 0; i < 3; i++) cur &= (subMaps[i].first == vars[i]);
    TEST_ASSERT(cur,
                "\n   tInerlacedEpetra::test_buildSubMaps_vec: "
                    << toString(status) << ": "
                    << "testing that maps are associated with the correct numbers of variables");

    int* gids;

    // test the first of three used maps
    TEST_EQUALITY(subMaps[0].second->NumGlobalElements(), 6,
                  "\n   tInerlacedEpetra::test_buildSubMaps_vec: "
                      << toString(status) << ": "
                      << "testing that first map has correct number of global elements ( "
                      << subMaps[0].second->NumGlobalElements() << " == " << 6 << " ?)");
    gids = subMaps[0].second->MyGlobalElements();
    cur  = (gids[0] == 0 && gids[1] == 1 && gids[2] == 6 && gids[3] == 7 && gids[4] == 12 &&
           gids[5] == 13);
    TEST_ASSERT(cur, "\n   tInerlacedEpetra::test_buildSubMaps_vec: "
                         << toString(status) << ": "
                         << "testing that first map is created correctly");

    // test the second of three used maps
    TEST_EQUALITY(subMaps[1].second->NumGlobalElements(), 3,
                  "\n   tInerlacedEpetra::test_buildSubMaps_vec: "
                      << toString(status) << ": "
                      << "testing that second map has correct number of global elements ( "
                      << subMaps[1].second->NumGlobalElements() << " == " << 3 << " ?)");
    gids = subMaps[1].second->MyGlobalElements();
    cur  = (gids[0] == 2 && gids[1] == 8 && gids[2] == 14);
    TEST_ASSERT(cur, "\n   tInerlacedEpetra::test_buildSubMaps_vec: "
                         << toString(status) << ": "
                         << "testing that second map is created correctly");

    // test the first of three used maps
    TEST_EQUALITY(subMaps[2].second->NumGlobalElements(), 9,
                  "\n   tInerlacedEpetra::test_buildSubMaps_vec: "
                      << toString(status) << ": "
                      << "testing that the third map has correct number of global elements ( "
                      << subMaps[2].second->NumGlobalElements() << " == " << 9 << " ?)");
    gids = subMaps[2].second->MyGlobalElements();
    cur  = (gids[0] == 3 && gids[1] == 4 && gids[2] == 5 && gids[3] == 9 && gids[4] == 10 &&
           gids[5] == 11 && gids[6] == 15 && gids[7] == 16 && gids[8] == 17);
    TEST_ASSERT(cur, "\n   tInerlacedEpetra::test_buildSubMaps_vec: "
                         << toString(status) << ": "
                         << "testing that third map is created correctly");
  } catch (...) {
    TEST_ASSERT(false,
                "\n   tInerlacedEpetra::test_buildSubMaps_vec: "
                    << toString(status)
                    << ": "
                       "buildSubMaps(int,vector<pair<int,RCP<Epetra_Map> > >) threw an unexpected "
                       "exception");
  }

  return allPassed;
}

bool tInterlacedEpetra::test_buildMaps(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  int size = 3 * 1000;

  TEST_MSG("\n   Builing Epetra_Map");
  RCP<Epetra_Map> map = rcp(new Epetra_Map(size, 0, *GetComm()));

  TEST_MSG("\n   Building sub maps");
  std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > > subMaps;
  std::vector<int> vec(2);
  vec[0] = 2;
  vec[1] = 1;
  Strided::buildSubMaps(size, vec, *GetComm(), subMaps);

  TEST_ASSERT(subMaps[0].first == vec[0], "   tInterlacedEpetra::test_buildMaps ("
                                              << Teko::Test::toString(status) << "): "
                                              << "  first map unknowns is " << subMaps[0].first
                                              << " ( should be " << vec[0] << ")");
  TEST_ASSERT(subMaps[1].first == vec[1], "   tInterlacedEpetra::test_buildMaps ("
                                              << Teko::Test::toString(status) << "): "
                                              << "  second map unknowns is " << subMaps[1].first
                                              << " ( should be " << vec[1] << ")");

  std::vector<Teuchos::RCP<Epetra_Map> > globalMaps(2);
  std::vector<Teuchos::RCP<Epetra_Map> > contigMaps(2);

  // get sub maps for convenient use and access
  globalMaps[0] = subMaps[0].second;
  globalMaps[1] = subMaps[1].second;

  contigMaps[0] = Teuchos::get_extra_data<Teuchos::RCP<Epetra_Map> >(globalMaps[0], "contigMap");
  contigMaps[1] = Teuchos::get_extra_data<Teuchos::RCP<Epetra_Map> >(globalMaps[1], "contigMap");

  // test that the extra data is attached
  TEST_ASSERT(contigMaps[0] != Teuchos::null,
              "   tInterlacedEpetra::test_buildMaps (" << Teko::Test::toString(status) << ")");
  TEST_ASSERT(contigMaps[1] != Teuchos::null,
              "   tInterlacedEpetra::test_buildMaps (" << Teko::Test::toString(status) << ")");
  TEST_MSG("   tInterlacedEpetra::test_buildMaps: extracted \"contigMaps\" from RCP");

  // make sure all maps have the correct size
  TEST_ASSERT(globalMaps[0]->NumGlobalElements() == 2000,
              "   tInterlacedEpetra::test_buildMaps (" << Teko::Test::toString(status) << "): "
                                                       << "checking number of global elements");
  TEST_ASSERT(globalMaps[1]->NumGlobalElements() == 1000,
              "   tInterlacedEpetra::test_buildMaps (" << Teko::Test::toString(status) << "): "
                                                       << "checking number of global elements");
  TEST_ASSERT(contigMaps[0]->NumGlobalElements() == 2000,
              "   tInterlacedEpetra::test_buildMaps (" << Teko::Test::toString(status) << "): "
                                                       << "checking number of global elements");
  TEST_ASSERT(contigMaps[1]->NumGlobalElements() == 1000,
              "   tInterlacedEpetra::test_buildMaps (" << Teko::Test::toString(status) << "): "
                                                       << "checking number of global elements");

  TEST_ASSERT(contigMaps[0]->NumMyElements() == globalMaps[0]->NumMyElements(),
              "   tInterlacedEpetra::test_buildMaps ("
                  << Teko::Test::toString(status) << "): "
                  << "check for lineup of number of local elements");
  TEST_ASSERT(contigMaps[1]->NumMyElements() == globalMaps[1]->NumMyElements(),
              "   tInterlacedEpetra::test_buildMaps ("
                  << Teko::Test::toString(status) << "): "
                  << "check for lineup of number of local elements");

  bool test;

  // check contiguous and global maps
  test = true;
  for (int i = 0; i < globalMaps[0]->NumMyElements(); i += 2) {
    int block;
    int gid = globalMaps[0]->GID(i);
    int cid = contigMaps[0]->GID(i);

    block = gid / 3;
    test &= cid == 2 * block;

    int gidp1 = globalMaps[0]->GID(i + 1);
    int cidp1 = contigMaps[0]->GID(i + 1);

    block = (gidp1 - 1) / 3;
    test &= cidp1 == 2 * block + 1;
  }
  TEST_ASSERT(test, "   tInterlacedEpetra::test_buildMaps ("
                        << Teko::Test::toString(status) << "): "
                        << "checked that block maps were internally consitent");

  test = true;
  for (int i = 0; i < globalMaps[1]->NumMyElements(); i++) {
    int block;
    int gid = globalMaps[1]->GID(i);
    int cid = contigMaps[1]->GID(i);

    block = (gid - 2) / 3;
    test &= cid == block;
  }
  TEST_ASSERT(test, "   tInterlacedEpetra::test_buildMaps ("
                        << Teko::Test::toString(status) << "): "
                        << "checked that block maps were internally consitent");

  return allPassed;
}

bool tInterlacedEpetra::test_one2many(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  int size = 3 * 1000;
  TEST_MSG("\n   tInterlacedEpetra::test_one2many: Builing Epetra_Map and source vector");
  RCP<Epetra_Map> map       = rcp(new Epetra_Map(size, 0, *GetComm()));
  RCP<Epetra_MultiVector> v = rcp(new Epetra_MultiVector(*map, 1));
  v->Random();

  TEST_MSG("\n   tInterlacedEpetra::test_one2many: Building sub maps");
  std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > > subMaps;
  std::vector<int> vec(2);
  vec[0] = 2;
  vec[1] = 1;
  Strided::buildSubMaps(size, vec, *GetComm(), subMaps);

  TEST_MSG("\n   tInterlacedEpetra::test_one2many: Building Export/Import");
  std::vector<RCP<Epetra_Export> > subExport;
  std::vector<RCP<Epetra_Import> > subImport;
  Strided::buildExportImport(*map, subMaps, subExport, subImport);

  TEST_MSG("\n   tInterlacedEpetra::test_one2many: Building sub vectors");
  std::vector<RCP<Epetra_MultiVector> > subVectors;
  Strided::buildSubVectors(subMaps, subVectors, 1);

  std::vector<RCP<Epetra_MultiVector> >::const_iterator itr;
  for (itr = subVectors.begin(); itr != subVectors.end(); ++itr) {
    RCP<const Epetra_Map> lm = Teuchos::get_extra_data<RCP<Epetra_Map> >(*itr, "globalMap");
    TEST_ASSERT(lm != Teuchos::null, "   tInterlacedEpetra::test_buildMaps ("
                                         << Teko::Test::toString(status) << "): "
                                         << "check that vector contains \"globalMap\" in RCP");
  }

  TEST_MSG("\n   tInterlacedEpetra::test_one2many: Performing data copy");
  Strided::one2many(subVectors, *v, subImport);

  // just assume it works! :)

  return allPassed;
}

bool tInterlacedEpetra::test_many2one(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  int size = 3 * 1000;
  TEST_MSG("\n   tInterlacedEpetra::test_one2many: Builing Epetra_Map and source vector");
  RCP<Epetra_Map> map       = rcp(new Epetra_Map(size, 0, *GetComm()));
  RCP<Epetra_MultiVector> v = rcp(new Epetra_MultiVector(*map, 4));
  v->Random();

  TEST_MSG("\n   tInterlacedEpetra::test_one2many: Building sub maps");
  std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > > subMaps;
  std::vector<int> vec(2);
  vec[0] = 2;
  vec[1] = 1;
  Strided::buildSubMaps(size, vec, *GetComm(), subMaps);

  TEST_MSG("\n   tInterlacedEpetra::test_one2many: Building Export/Import");
  std::vector<RCP<Epetra_Export> > subExport;
  std::vector<RCP<Epetra_Import> > subImport;
  Strided::buildExportImport(*map, subMaps, subExport, subImport);

  TEST_MSG("\n   tInterlacedEpetra::test_one2many: Building sub vectors");
  std::vector<RCP<Epetra_MultiVector> > subVectors;
  Strided::buildSubVectors(subMaps, subVectors, 4);

  std::vector<RCP<Epetra_MultiVector> >::const_iterator itr;
  for (itr = subVectors.begin(); itr != subVectors.end(); ++itr) {
    RCP<const Epetra_Map> lm = Teuchos::get_extra_data<RCP<Epetra_Map> >(*itr, "globalMap");
    TEST_ASSERT(lm != Teuchos::null, "   tInterlacedEpetra::test_buildMaps ("
                                         << Teko::Test::toString(status) << "): "
                                         << "check that vector contains \"globalMap\" in RCP");
  }

  TEST_MSG("\n   tInterlacedEpetra::test_one2many: Performing one2many");
  Strided::one2many(subVectors, *v, subImport);

  std::vector<RCP<const Epetra_MultiVector> > cSubVectors;
  for (itr = subVectors.begin(); itr != subVectors.end(); ++itr) cSubVectors.push_back(*itr);

  TEST_MSG("\n   tInterlacedEpetra::test_one2many: Performing many2one");
  RCP<Epetra_MultiVector> one = rcp(new Epetra_MultiVector(*map, 1));
  Strided::many2one(*one, cSubVectors, subExport);

  one->Update(1.0, *v, -1.0);

  double diff[4] = {0, 0, 0, 0};
  double max = 0.0, maxn = 0;
  double norm[4] = {0, 0, 0, 0};
  one->Norm2(diff);
  v->Norm2(norm);
  for (int i = 0; i < 4; i++) {
    max  = max > diff[i] / norm[i] ? max : diff[i] / norm[i];
    maxn = maxn > norm[i] ? maxn : norm[i];
  }
  TEST_ASSERT(max <= tolerance_, "   tInterlacedEpetra::test_buildMaps ("
                                     << Teko::Test::toString(status) << "): "
                                     << "norm must be better than the tolerance ( " << max
                                     << " <=? " << tolerance_ << " maxn = " << maxn << " )");

  return allPassed;
}

}  // namespace Test
}  // namespace Teko
