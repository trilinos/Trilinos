// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tDiagonalPreconditionerFactory_tpetra.hpp"
#include "Teko_DiagonalPreconditionerFactory.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

// Tpetra includes
#include "Tpetra_Vector.hpp"

// Galeri / Xpetra
#include "Galeri_XpetraMaps.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraParameters.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_TpetraHelpers.hpp"
#include "Thyra_TpetraLinearOp.hpp"

#include <vector>
#include <cmath>

namespace Teko {
namespace Test {

using namespace Teuchos;
using namespace Thyra;

namespace {

using ST = Teko::ST;
using LO = Teko::LO;
using GO = Teko::GO;
using NT = Teko::NT;

using map_t = Tpetra::Map<LO, GO, NT>;
using crs_t = Tpetra::CrsMatrix<ST, LO, GO, NT>;
using vec_t = Tpetra::Vector<ST, LO, GO, NT>;
using mv_t  = Tpetra::MultiVector<ST, LO, GO, NT>;

RCP<crs_t> buildLaplace2DMatrix(const RCP<const Teuchos::Comm<int> >& comm, GO nx, GO ny) {
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  galeriList.set("mx", comm->getSize());
  galeriList.set("my", 1);

  auto tMap = Galeri::Xpetra::CreateMap<LO, GO, map_t>("Cartesian2D", comm, galeriList);

  auto problem =
      Galeri::Xpetra::BuildProblem<ST, LO, GO, map_t, crs_t, mv_t>("Laplace2D", tMap, galeriList);

  return problem->BuildMatrix();
}

}  // namespace

tDiagonalPreconditionerFactory_tpetra::~tDiagonalPreconditionerFactory_tpetra() {
  delete fact;
  delete pstate;
}

void tDiagonalPreconditionerFactory_tpetra::initializeTest() {
  const RCP<const Teuchos::Comm<int> > comm_tpetra = GetComm_tpetra();

  tolerance_ = 1.0e-14;

  GO nx = 39;
  GO ny = 53;

  auto tmpF = buildLaplace2DMatrix(comm_tpetra, nx, ny);
  tpetraF   = tmpF;
  F_        = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraF->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraF->getRangeMap()), tpetraF);
}

void tDiagonalPreconditionerFactory_tpetra::buildParameterList(int blocksize) {
  const crs_t* F = &*tpetraF;
  TEUCHOS_ASSERT(F);

  List_ = Teuchos::ParameterList();

  if (blocksize > 0) {
    const LO Nr = F->getLocalNumRows();
    const int Nb =
        static_cast<int>(std::ceil(static_cast<double>(Nr) / static_cast<double>(blocksize)));

    block_starts.resize(Nb + 1);
    block_gids.resize(Nr);

    block_starts[0] = 0;
    for (int i = 0; i < Nb; i++) block_starts[i + 1] = block_starts[i] + blocksize;
    block_starts[Nb] = Nr;

    for (LO i = 0; i < Nr; i++) block_gids[i] = F->getRowMap()->getGlobalElement(i);

    Teuchos::ParameterList sublist;
    List_.set("Diagonal Type", "BlkDiag");
    List_.set("number of local blocks", Nb);
    List_.set("block start index", block_starts);
    List_.set("block entry gids", block_gids);
    sublist.set("apply mode", "invert");
    List_.set("blockdiagmatrix: list", sublist);
  } else {
    List_.set("Diagonal Type", "Diagonal");
  }
}

int tDiagonalPreconditionerFactory_tpetra::runTest(int verbosity, std::ostream& stdstrm,
                                                   std::ostream& failstrm, int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tDiagonalPreconditionerFactory_tpetra";

  status = test_createPrec(verbosity, failstrm, 2);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"createPrec\" ... PASSED", "   \"createPrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"initializePrec\" ... PASSED",
                       "   \"initializePrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_canApply(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"canApply\" ... PASSED", "   \"canApply\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tDiagonalPreconditionedFactory...PASSED",
                         "tDiagonalPreconditionedFactory...FAILED");
  } else {
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tDiagonalPreconditionedFactory...FAILED");
  }

  return failcount;
}

bool tDiagonalPreconditionerFactory_tpetra::test_initializePrec(int verbosity, std::ostream& os) {
  delete pstate;

  pstate = new DiagonalPrecondState();
  pop    = fact->buildPreconditionerOperator(F_, *pstate);

  // For blocksize==2 this should be an explicit sparse operator, not a diagonal op
  if (rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(pop).is_null()) return false;

  return true;
}

bool tDiagonalPreconditionerFactory_tpetra::test_createPrec(int verbosity, std::ostream& os,
                                                            int blocksize) {
  buildParameterList(blocksize);

  delete fact;
  fact = new DiagonalPreconditionerFactory();
  fact->initializeFromParameterList(List_);
  if (!fact) return false;

  return true;
}

bool tDiagonalPreconditionerFactory_tpetra::test_canApply(int verbosity, std::ostream& os) {
  RCP<const Tpetra::Map<LO, GO, NT> > domain_ = tpetraF->getDomainMap();
  RCP<const Tpetra::Map<LO, GO, NT> > range_  = tpetraF->getRangeMap();

  RCP<Tpetra::Vector<ST, LO, GO, NT> > X = Tpetra::createVector<ST, LO, GO, NT>(domain_);
  RCP<Tpetra::Vector<ST, LO, GO, NT> > Y = Tpetra::createVector<ST, LO, GO, NT>(range_);
  RCP<Tpetra::Vector<ST, LO, GO, NT> > Z = Tpetra::createVector<ST, LO, GO, NT>(range_);
  Y->putScalar(0.0);
  Z->putScalar(1.0);

  // Let X = diag(F). Then applying the preconditioner to X should yield a vector of ones
  tpetraF->getLocalDiagCopy(*X);

  // Build Thyra wrappers
  MultiVector tX =
      Thyra::createVector<ST, LO, GO, NT>(X, Thyra::createVectorSpace<ST, LO, GO, NT>(domain_));
  MultiVector tY =
      Thyra::createVector<ST, LO, GO, NT>(Y, Thyra::createVectorSpace<ST, LO, GO, NT>(range_));

  // Do the apply via thyra
  Teko::applyOp(pop, tX, tY, 1.0, 0.0);

  // Compare solutions
  double znrm, dnrm;
  znrm = Z->norm2();
  Z->update(-1.0, *Y, 1.0);
  dnrm = Z->norm2();

  if (!tpetraF->getComm()->getRank()) std::cout << "||Z-Y||/||Z|| = " << dnrm / znrm << std::endl;
  if (dnrm / znrm > 1e-12)
    return false;
  else
    return true;
}

}  // end namespace Test
}  // end namespace Teko
