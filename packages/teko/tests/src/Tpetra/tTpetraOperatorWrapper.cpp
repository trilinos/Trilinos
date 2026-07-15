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
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
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

// Galeri includes
#include "Galeri_XpetraMaps.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraParameters.hpp"

// Teko-Package includes
#include "Teko_ConfigDefs.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_TpetraHelpers.hpp"

#include "tTpetraOperatorWrapper.hpp"

namespace Teko {
namespace Test {

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;
using Thyra::createMember;
using Thyra::LinearOpBase;
using Thyra::LinearOpTester;
using Thyra::MultiVectorBase;
using Thyra::VectorBase;

namespace {

using ST = Teko::ST;
using LO = Teko::LO;
using GO = Teko::GO;
using NT = Teko::NT;

using map_t = Tpetra::Map<LO, GO, NT>;
using crs_t = Tpetra::CrsMatrix<ST, LO, GO, NT>;
using vec_t = Tpetra::Vector<ST, LO, GO, NT>;

RCP<crs_t> buildRecirc2DMatrix(const RCP<const Teuchos::Comm<int>>& comm, GO nx, GO ny) {
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  galeriList.set("mx", comm->getSize());
  galeriList.set("my", 1);

  auto tMap = Galeri::Xpetra::CreateMap<LO, GO, map_t>("Cartesian2D", comm, galeriList);

  auto problem =
      Galeri::Xpetra::BuildProblem<ST, LO, GO, map_t, crs_t, Tpetra::MultiVector<ST, LO, GO, NT>>(
          "Recirc2D", tMap, galeriList);

  return problem->BuildMatrix();
}

RCP<crs_t> buildLaplace2DMatrix(const RCP<const Teuchos::Comm<int>>& comm, GO nx, GO ny) {
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  galeriList.set("mx", comm->getSize());
  galeriList.set("my", 1);

  auto tMap = Galeri::Xpetra::CreateMap<LO, GO, map_t>("Cartesian2D", comm, galeriList);

  auto problem =
      Galeri::Xpetra::BuildProblem<ST, LO, GO, map_t, crs_t, Tpetra::MultiVector<ST, LO, GO, NT>>(
          "Laplace2D", tMap, galeriList);

  return problem->BuildMatrix();
}

RCP<crs_t> buildDiagMatrix(const RCP<const Tpetra::Map<LO, GO, NT>>& rangeMap,
                           const RCP<const Tpetra::Map<LO, GO, NT>>& domainMap, double value) {
  using crs_t = Tpetra::CrsMatrix<ST, LO, GO, NT>;

  auto A = Teuchos::rcp(new crs_t(rangeMap, 1));

  Teuchos::Array<GO> cols(1);
  Teuchos::Array<ST> vals(1);
  vals[0] = value;

  const size_t localNumRows = rangeMap->getLocalNumElements();

  for (size_t localRow = 0; localRow < localNumRows; ++localRow) {
    const GO gid = rangeMap->getGlobalElement(static_cast<LO>(localRow));

    cols[0] = gid;

    A->insertGlobalValues(gid, cols(), vals());
  }

  A->fillComplete(domainMap, rangeMap);

  return A;
}

}  // namespace

void tTpetraOperatorWrapper::initializeTest() {}

int tTpetraOperatorWrapper::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                                    int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tTpetraOperatorWrapper";

  status = test_functionality(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"functionality\" ... PASSED",
                       "   \"functionality\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tTpetraOperatorWrapper...PASSED",
                         "tTpetraOperatorWrapper...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tTpetraOperatorWrapper...FAILED");
  }

  return failcount;
}

bool tTpetraOperatorWrapper::test_functionality(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  RCP<const Teuchos::Comm<int>> comm_tpetra = GetComm_tpetra();

  TEST_MSG("\n   tTpetraOperatorWrapper::test_functionality: "
           << "Running on " << comm_tpetra->getSize() << " processors");

  GO nx = 39;  // essentially random values
  GO ny = 53;

  TEST_MSG("   tTpetraOperatorWrapper::test_functionality: "
           << "Using Galeri to create test matrices");

  // create some big blocks to play with
  RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tpetraF = buildRecirc2DMatrix(comm_tpetra, nx, ny);
  RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tpetraC = buildLaplace2DMatrix(comm_tpetra, nx, ny);

  RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tpetraB =
      buildDiagMatrix(tpetraC->getRangeMap(), tpetraF->getDomainMap(), 5.0);

  RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT>> tpetraBt =
      buildDiagMatrix(tpetraF->getRangeMap(), tpetraC->getDomainMap(), 3.0);

  // load'em up in a thyra operator
  TEST_MSG("   tTpetraOperatorWrapper::test_functionality: "
           << " Building block2x2 Thyra matrix ... wrapping in TpetraOperatorWrapper");
  const RCP<const LinearOpBase<double>> A = Thyra::block2x2<double>(
      Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraF->getDomainMap()),
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraF->getRangeMap()), tpetraF),
      Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraBt->getDomainMap()),
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraBt->getRangeMap()), tpetraBt),
      Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraB->getDomainMap()),
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraB->getRangeMap()), tpetraB),
      Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraC->getDomainMap()),
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraC->getRangeMap()), tpetraC),
      "A");

  // const RCP<Thyra::TpetraOperatorWrapper> epetra_A = rcp(new Thyra::TpetraOperatorWrapper(A));
  const RCP<Teko::TpetraHelpers::TpetraOperatorWrapper> tpetra_A =
      rcp(new Teko::TpetraHelpers::TpetraOperatorWrapper(A));

  // begin the tests!
  const RCP<const Tpetra::Map<LO, GO, NT>>& rangeMap  = tpetra_A->getRangeMap();
  const RCP<const Tpetra::Map<LO, GO, NT>>& domainMap = tpetra_A->getDomainMap();

  // check to see that the number of global elements is correct
  TEST_EQUALITY(rangeMap->getGlobalNumElements(), (Tpetra::global_size_t)2 * nx * ny,
                "   tTpetraOperatorWrapper::test_functionality: "
                    << toString(status) << ": "
                    << "checking rangeMap size "
                    << "( map = " << rangeMap->getGlobalNumElements() << ", true = " << 2 * nx * ny
                    << " )");

  // check to see that the number of global elements is correct
  TEST_EQUALITY(domainMap->getGlobalNumElements(), (Tpetra::global_size_t)2 * nx * ny,
                "   tTpetraOperatorWrapper::test_functionality: "
                    << toString(status) << ": "
                    << "checking domainMap size "
                    << "( map = " << domainMap->getGlobalNumElements() << ", true = " << 2 * nx * ny
                    << " )");

  // largest global ID should be one less then the # of elements
  TEST_EQUALITY(rangeMap->getGlobalNumElements() - 1,
                (Tpetra::global_size_t)rangeMap->getMaxAllGlobalIndex(),
                "   tTpetraOperatorWrapper::test_functionality: "
                    << toString(status) << ": "
                    << " checking largest range element "
                    << "( largest = " << rangeMap->getMaxAllGlobalIndex()
                    << ", true = " << rangeMap->getGlobalNumElements() - 1 << " )");
  TEST_EQUALITY(domainMap->getGlobalNumElements() - 1,
                (Tpetra::global_size_t)domainMap->getMaxAllGlobalIndex(),
                "   tTpetraOperatorWrapper::test_functionality: "
                    << toString(status) << ": "
                    << " checking largest domain element "
                    << "( largest = " << domainMap->getMaxAllGlobalIndex()
                    << ", true = " << domainMap->getGlobalNumElements() - 1 << " )");

  RCP<const Teko::TpetraHelpers::MappingStrategy> ms = tpetra_A->getMapStrategy();

  // create a vector to test: copyThyraIntoTpetra
  //////////////////////////////////////////////////////////////
  {
    const RCP<MultiVectorBase<ST>> tv = Thyra::createMembers(A->domain(), 1);
    Thyra::randomize(-100.0, 100.0, tv.ptr());
    // const Thyra::ConstVector<double> handle_tv(tv);
    const RCP<const MultiVectorBase<ST>> tv_0 =
        Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<ST>>(tv)->getMultiVectorBlock(
            0);
    const RCP<const MultiVectorBase<ST>> tv_1 =
        Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<ST>>(tv)->getMultiVectorBlock(
            1);
    const Thyra::ConstDetachedSpmdVectorView<ST> vv_0(tv_0->col(0));
    const Thyra::ConstDetachedSpmdVectorView<ST> vv_1(tv_1->col(0));

    LO off_0 = vv_0.globalOffset();
    LO off_1 = vv_1.globalOffset();

    // create its Tpetra counter part
    const RCP<Tpetra::Vector<ST, LO, GO, NT>> ev =
        rcp(new Tpetra::Vector<ST, LO, GO, NT>(tpetra_A->getDomainMap()));
    ms->copyThyraIntoTpetra(tv, *ev);

    // compare tv to ev!
    TEST_EQUALITY((Tpetra::global_size_t)tv->range()->dim(), ev->getGlobalLength(),
                  "   tTpetraOperatorWrapper::test_functionality: "
                      << toString(status) << ": "
                      << " checking ThyraIntoTpetra copy "
                      << "( thyra dim = " << tv->range()->dim()
                      << ", global length = " << ev->getGlobalLength() << " )");
    LO numMyElements = domainMap->getLocalNumElements();
    TEST_MSG("domainMap->getLocalNumElements() = " << domainMap->getLocalNumElements());
    bool compareThyraToTpetraValue = true;
    ST tval                        = 0.0;
    for (LO i = 0; i < numMyElements; i++) {
      GO gid = domainMap->getGlobalElement(i);
      if (gid - off_0 < nx * ny) {
        tval = vv_0[gid - off_0];
      } else {
        tval = vv_1[gid - off_1 - nx * ny];
      }
      compareThyraToTpetraValue &= (ev->get1dView()[i] == tval);
    }
    TEST_ASSERT(compareThyraToTpetraValue, "   tTpetraOperatorWrapper::test_functionality: "
                                               << toString(status) << ": "
                                               << " comparing Thyra to Tpetra values");
  }

  // create a vector to test: copyTpetraIntoThyra
  //////////////////////////////////////////////////////////////
  {
    // create an Tpetra vector
    const RCP<Tpetra::Vector<ST, LO, GO, NT>> ev =
        rcp(new Tpetra::Vector<ST, LO, GO, NT>(tpetra_A->getDomainMap()));
    ev->randomize();

    // create its thyra counterpart
    const RCP<MultiVectorBase<ST>> tv = Thyra::createMembers(A->domain(), 1);
    const RCP<const MultiVectorBase<ST>> tv_0 =
        Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<ST>>(tv)->getMultiVectorBlock(
            0);
    const RCP<const MultiVectorBase<ST>> tv_1 =
        Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<ST>>(tv)->getMultiVectorBlock(
            1);
    const Thyra::ConstDetachedSpmdVectorView<ST> vv_0(tv_0->col(0));
    const Thyra::ConstDetachedSpmdVectorView<ST> vv_1(tv_1->col(0));

    LO off_0 = rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<ST>>(tv_0->range())->localOffset();
    LO off_1 = rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<ST>>(tv_1->range())->localOffset();

    ms->copyTpetraIntoThyra(*ev, tv.ptr());

    // compare handle_tv to ev!
    TEST_EQUALITY((Tpetra::global_size_t)tv->range()->dim(), ev->getGlobalLength(),
                  "   tTpetraOperatorWrapper::test_functionality: "
                      << toString(status) << ": "
                      << " checking TpetraIntoThyra copy "
                      << "( thyra dim = " << tv->range()->dim()
                      << ", global length = " << ev->getGlobalLength() << " )");
    LO numMyElements               = domainMap->getLocalNumElements();
    bool compareTpetraToThyraValue = true;
    ST tval                        = 0.0;
    for (LO i = 0; i < numMyElements; i++) {
      GO gid = domainMap->getGlobalElement(i);
      if (gid - off_0 < nx * ny) {
        tval = vv_0[gid - off_0];
      } else {
        tval = vv_1[gid - off_1 - nx * ny];
      }
      compareTpetraToThyraValue &= (ev->get1dView()[i] == tval);
    }
    TEST_ASSERT(compareTpetraToThyraValue, "   tTpetraOperatorWrapper::test_functionality: "
                                               << toString(status) << ": "
                                               << " comparing Thyra to Tpetra values");
  }

  return allPassed;
}

}  // namespace Test
}  // namespace Teko
