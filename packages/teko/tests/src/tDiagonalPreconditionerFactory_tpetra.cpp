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
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraVectorSpace.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

// Tpetra includes
#include "Tpetra_Vector.hpp"

// Galeri
#include "Galeri_XpetraMaps.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraParameters.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_TpetraHelpers.hpp"

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

RCP<crs_t> buildTwoByTwoBlockDiagMatrix(const RCP<const Teuchos::Comm<int>>& comm) {
  auto map = rcp(new const map_t(4, 0, comm));
  auto A   = rcp(new crs_t(map, 2));

  // block 1:
  // [2 1]
  // [1 3]
  {
    GO row = 0;
    std::vector<GO> cols{0, 1};
    std::vector<ST> vals{2.0, 1.0};
    A->insertGlobalValues(row, Teuchos::ArrayView<const GO>(cols),
                          Teuchos::ArrayView<const ST>(vals));
  }
  {
    GO row = 1;
    std::vector<GO> cols{0, 1};
    std::vector<ST> vals{1.0, 3.0};
    A->insertGlobalValues(row, Teuchos::ArrayView<const GO>(cols),
                          Teuchos::ArrayView<const ST>(vals));
  }

  // block 2:
  // [4 1]
  // [1 2]
  {
    GO row = 2;
    std::vector<GO> cols{2, 3};
    std::vector<ST> vals{4.0, 1.0};
    A->insertGlobalValues(row, Teuchos::ArrayView<const GO>(cols),
                          Teuchos::ArrayView<const ST>(vals));
  }
  {
    GO row = 3;
    std::vector<GO> cols{2, 3};
    std::vector<ST> vals{1.0, 2.0};
    A->insertGlobalValues(row, Teuchos::ArrayView<const GO>(cols),
                          Teuchos::ArrayView<const ST>(vals));
  }

  A->fillComplete();
  return A;
}

void fillVector(vec_t& v, const std::vector<ST>& vals) {
  for (size_t i = 0; i < vals.size(); ++i) {
    v.replaceGlobalValue(static_cast<GO>(i), vals[i]);
  }
}

double relativeDifference(const RCP<const vec_t>& a, const RCP<const vec_t>& b) {
  auto diff = rcp(new vec_t(a->getMap()));
  diff->assign(*a);
  diff->update(1.0, *b, -1.0);
  const double dn = diff->norm2();
  const double bn = b->norm2();
  return (bn == 0.0 ? dn : dn / bn);
}

}  // namespace

void tDiagonalPreconditionerFactory_tpetra::initializeTest() {
  const RCP<const Teuchos::Comm<int>> comm_tpetra = GetComm_tpetra();

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

  status = test_createPrec(verbosity, failstrm, 1);
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

  status = test_blkdiag_exact_2x2_blocks(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"blkdiag_exact_2x2_blocks\" ... PASSED",
                       "   \"blkdiag_exact_2x2_blocks\" ... FAILED");
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
  pstate = Teuchos::make_rcp<DiagonalPrecondState>();
  pop    = fact->buildPreconditionerOperator(F_, *pstate);

  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT>> top =
      rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT>>(pop);
  if (top.is_null()) return false;

  return true;
}

bool tDiagonalPreconditionerFactory_tpetra::test_createPrec(int verbosity, std::ostream& os,
                                                            int blocksize) {
  buildParameterList(blocksize);

  fact = Teuchos::make_rcp<DiagonalPreconditionerFactory>();
  fact->initializeFromParameterList(List_);
  if (!fact) return false;

  return true;
}

bool tDiagonalPreconditionerFactory_tpetra::test_canApply(int verbosity, std::ostream& os) {
  RCP<const Tpetra::Map<LO, GO, NT>> domain_ = tpetraF->getDomainMap();
  RCP<const Tpetra::Map<LO, GO, NT>> range_  = tpetraF->getRangeMap();

  RCP<Tpetra::Vector<ST, LO, GO, NT>> X = Tpetra::createVector<ST, LO, GO, NT>(domain_);
  RCP<Tpetra::Vector<ST, LO, GO, NT>> Y = Tpetra::createVector<ST, LO, GO, NT>(range_);
  RCP<Tpetra::Vector<ST, LO, GO, NT>> Z = Tpetra::createVector<ST, LO, GO, NT>(range_);
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

  // Copy Thyra result back to Y
  RCP<const Thyra::TpetraVector<ST, LO, GO, NT>> thyraY =
      Teuchos::rcp_dynamic_cast<const Thyra::TpetraVector<ST, LO, GO, NT>>(tY, true);
  Y->assign(*(thyraY->getConstTpetraVector()));

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

bool tDiagonalPreconditionerFactory_tpetra::test_blkdiag_exact_2x2_blocks(int verbosity,
                                                                          std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  auto comm   = GetComm_tpetra();
  auto smallF = buildTwoByTwoBlockDiagMatrix(comm);

  Teko::LinearOp smallFop = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(smallF->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(smallF->getRangeMap()), smallF);

  Teuchos::ParameterList list;
  list.set("Diagonal Type", "BlkDiag");
  list.set("contiguous block size", 2);

  DiagonalPreconditionerFactory localFact;
  localFact.initializeFromParameterList(list);

  DiagonalPrecondState localState;
  auto H = localFact.buildPreconditionerOperator(smallFop, localState);

  auto Ht   = rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT>>(H, true);
  auto Hcrs = rcp_dynamic_cast<const crs_t>(Ht->getConstTpetraOperator(), true);

  RCP<const Tpetra::Map<LO, GO, NT>> map = smallF->getDomainMap();

  auto x1 = rcp(new vec_t(map));
  auto y1 = rcp(new vec_t(map));
  auto z1 = rcp(new vec_t(map));
  fillVector(*x1, {1.0, 0.0, 1.0, 0.0});
  y1->putScalar(0.0);
  fillVector(*z1, {3.0 / 5.0, -1.0 / 5.0, 2.0 / 7.0, -1.0 / 7.0});

  Hcrs->apply(*x1, *y1);
  double rel       = relativeDifference(y1, z1);
  const double eps = 1e-14;
  if (not(rel <= eps) || verbosity >= 10) {
    os << "blkdiag exact test 1 rel = " << rel << std::endl;
  }
  status = (rel <= eps);
  allPassed &= status;

  auto x2 = rcp(new vec_t(map));
  auto y2 = rcp(new vec_t(map));
  auto z2 = rcp(new vec_t(map));
  fillVector(*x2, {0.0, 1.0, 0.0, 1.0});
  y2->putScalar(0.0);
  fillVector(*z2, {-1.0 / 5.0, 2.0 / 5.0, -1.0 / 7.0, 4.0 / 7.0});

  Hcrs->apply(*x2, *y2);
  rel = relativeDifference(y2, z2);
  if (not(rel <= eps) || verbosity >= 10) {
    os << "blkdiag exact test 2 rel = " << rel << std::endl;
  }
  status = (rel <= eps);
  allPassed &= status;

  auto x3 = rcp(new vec_t(map));
  auto y3 = rcp(new vec_t(map));
  auto z3 = rcp(new vec_t(map));
  fillVector(*x3, {1.0, 2.0, 3.0, 4.0});
  y3->putScalar(0.0);
  fillVector(*z3, {1.0 / 5.0, 3.0 / 5.0, 2.0 / 7.0, 13.0 / 7.0});

  Hcrs->apply(*x3, *y3);
  rel = relativeDifference(y3, z3);
  if (not(rel <= eps) || verbosity >= 10) {
    os << "blkdiag exact test 3 rel = " << rel << std::endl;
  }
  status = (rel <= eps);
  allPassed &= status;

  return allPassed;
}

}  // end namespace Test
}  // end namespace Teko
