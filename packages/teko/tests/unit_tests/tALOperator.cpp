// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * Author: Zhen Wang
 * Email: wangz@ornl.gov
 *        zhen.wang@alum.emory.edu
 */

/*
 * This test constructs a tiny Stokes-like blocked operator directly in-memory
 * and tests Teko::NS::ALOperator on the Tpetra stack by comparing against an
 * explicitly assembled reference augmented-Lagrangian operator and augmented RHS.
 */

#include "Teko_Config.h"

#include <iostream>
#include <vector>
#include <sys/types.h>
#include <unistd.h>

// Teuchos
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"

// Tpetra
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"

// Thyra
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_LinearOpTester.hpp"

// Teko
#include "Teko_ALOperator.hpp"
#include "Teko_BlockedTpetraOperator.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_ConfigDefs.hpp"
#include "Teko_Utilities.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace {

using ST = Teko::ST;
using LO = Teko::LO;
using GO = Teko::GO;
using NT = Teko::NT;

using map_t = Tpetra::Map<LO, GO, NT>;
using crs_t = Tpetra::CrsMatrix<ST, LO, GO, NT>;
using vec_t = Tpetra::Vector<ST, LO, GO, NT>;
using mv_t  = Tpetra::MultiVector<ST, LO, GO, NT>;

RCP<const map_t> buildMap(const RCP<const Teuchos::Comm<int> >& comm, GO n) {
  return rcp(new map_t(n, 0, comm));
}

RCP<crs_t> buildDiagMatrix(const RCP<const map_t>& map, const std::vector<ST>& diag) {
  auto A = rcp(new crs_t(map, 1));
  for (size_t i = 0; i < diag.size(); ++i) {
    GO gid = map->getGlobalElement(static_cast<LO>(i));
    GO col = gid;
    ST val = diag[i];
    A->insertGlobalValues(gid, Teuchos::ArrayView<const GO>(&col, 1),
                          Teuchos::ArrayView<const ST>(&val, 1));
  }
  A->fillComplete();
  return A;
}

RCP<crs_t> buildRectMatrix(const RCP<const map_t>& rowMap, const RCP<const map_t>& colMap,
                           const std::vector<std::vector<std::pair<GO, ST> > >& rows) {
  auto A = rcp(new crs_t(rowMap, 4));
  for (size_t i = 0; i < rows.size(); ++i) {
    GO rowGid = rowMap->getGlobalElement(static_cast<LO>(i));
    std::vector<GO> cols;
    std::vector<ST> vals;
    for (const auto& entry : rows[i]) {
      cols.push_back(entry.first);
      vals.push_back(entry.second);
    }
    if (!cols.empty()) {
      A->insertGlobalValues(rowGid, Teuchos::ArrayView<const GO>(cols),
                            Teuchos::ArrayView<const ST>(vals));
    }
  }
  A->fillComplete(colMap, rowMap);
  return A;
}

RCP<crs_t> buildFullSystem(const RCP<const Teuchos::Comm<int> >& comm,
                           std::vector<std::vector<GO> >& blockedVec,
                           Teko::LinearOp& pressureMassMatrix) {
  const GO nUx = 3, nUy = 3, nP = 2;
  const GO nAll = nUx + nUy + nP;

  auto fullMap = buildMap(comm, nAll);

  blockedVec.resize(3);
  blockedVec[0] = {0, 2, 4};
  blockedVec[1] = {1, 3, 5};
  blockedVec[2] = {6, 7};

  auto pMap          = buildMap(comm, nP);
  auto Mp            = buildDiagMatrix(pMap, {3.0, 4.0});
  pressureMassMatrix = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(Mp->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(Mp->getDomainMap()), Mp);

  std::vector<std::vector<std::pair<GO, ST> > > rows(8);

  rows[0] = {{0, 1.0}, {6, 1.0}};
  rows[1] = {{1, 1.0}, {7, 1.0}};
  rows[2] = {{2, 1.0}, {6, 1.0}};
  rows[3] = {{3, 2.0}, {7, 1.0}};
  rows[4] = {{4, 2.0}, {6, 1.0}};
  rows[5] = {{5, 2.0}, {7, 1.0}};
  rows[6] = {{0, 2.0}, {2, 2.0}, {4, 2.0}, {6, -3.0}};
  rows[7] = {{1, 2.0}, {3, 2.0}, {5, 2.0}, {7, -4.0}};

  auto A = buildRectMatrix(fullMap, fullMap, rows);
  return A;
}

void printVector(const std::string& label, const RCP<const vec_t>& v, Teuchos::FancyOStream& out) {
  auto view = v->getLocalViewHost(Tpetra::Access::ReadOnly);
  out << label << ":" << std::endl;
  for (size_t i = 0; i < v->getLocalLength(); ++i) {
    out << "  [" << i << "] = " << view(i, 0) << std::endl;
  }
}

RCP<const vec_t> applyReferenceBlockedOpByBlockedVec(
    const Teko::LinearOp& op, const std::vector<std::vector<GO> >& blockedVec,
    const RCP<const map_t>& flatMap) {
  auto xThyra = Thyra::createMembers(op->domain(), 1);
  auto yThyra = Thyra::createMembers(op->range(), 1);

  Thyra::assign(xThyra.ptr(), 1.0);
  Thyra::assign(yThyra.ptr(), 0.0);

  op->apply(Thyra::NOTRANS, *xThyra, yThyra.ptr(), 1.0, 0.0);

  auto yProd = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<ST> >(yThyra, true);

  auto yFlat = rcp(new vec_t(flatMap));

  for (size_t b = 0; b < blockedVec.size(); ++b) {
    auto blk = yProd->getMultiVectorBlock(static_cast<int>(b));
    auto col = blk->col(0);
    Thyra::ConstDetachedVectorView<ST> view(*col);

    for (size_t k = 0; k < blockedVec[b].size(); ++k) {
      yFlat->replaceGlobalValue(blockedVec[b][k], view[k]);
    }
  }

  return yFlat;
}

Teko::LinearOp buildReferenceALOperator(const RCP<const crs_t>& fullA,
                                        const std::vector<std::vector<GO> >& blockedVec,
                                        const Teko::LinearOp& pressureMassMatrix, double gamma) {
  Teko::TpetraHelpers::BlockedTpetraOperator blk(blockedVec, fullA);

  auto A00 = blk.GetBlock(0, 0);
  auto A01 = blk.GetBlock(0, 1);
  auto A02 = blk.GetBlock(0, 2);

  auto A10 = blk.GetBlock(1, 0);
  auto A11 = blk.GetBlock(1, 1);
  auto A12 = blk.GetBlock(1, 2);

  auto A20 = blk.GetBlock(2, 0);
  auto A21 = blk.GetBlock(2, 1);
  auto A22 = blk.GetBlock(2, 2);

  Teko::LinearOp tA00 = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A00->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A00->getDomainMap()), A00);
  Teko::LinearOp tA01 = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A01->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A01->getDomainMap()), A01);
  Teko::LinearOp tA02 = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A02->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A02->getDomainMap()), A02);

  Teko::LinearOp tA10 = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A10->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A10->getDomainMap()), A10);
  Teko::LinearOp tA11 = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A11->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A11->getDomainMap()), A11);
  Teko::LinearOp tA12 = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A12->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A12->getDomainMap()), A12);

  Teko::LinearOp tA20 = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A20->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A20->getDomainMap()), A20);
  Teko::LinearOp tA21 = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A21->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A21->getDomainMap()), A21);
  Teko::LinearOp tA22 = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A22->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A22->getDomainMap()), A22);

  Teko::LinearOp Winv = Teko::getInvDiagonalOp(pressureMassMatrix);

  auto ref = Teko::createBlockedOp();
  Teko::beginBlockFill(ref, 3, 3);

  Teko::setBlock(0, 0, ref,
                 Thyra::add(tA00, Thyra::scale(gamma, Thyra::multiply(tA02, Winv, tA20))));
  Teko::setBlock(0, 1, ref,
                 Thyra::add(tA01, Thyra::scale(gamma, Thyra::multiply(tA02, Winv, tA21))));
  Teko::setBlock(0, 2, ref,
                 Thyra::add(tA02, Thyra::scale(-gamma, Thyra::multiply(tA02, Winv, tA22))));

  Teko::setBlock(1, 0, ref,
                 Thyra::add(tA10, Thyra::scale(gamma, Thyra::multiply(tA12, Winv, tA20))));
  Teko::setBlock(1, 1, ref,
                 Thyra::add(tA11, Thyra::scale(gamma, Thyra::multiply(tA12, Winv, tA21))));
  Teko::setBlock(1, 2, ref,
                 Thyra::add(tA12, Thyra::scale(-gamma, Thyra::multiply(tA12, Winv, tA22))));

  Teko::setBlock(2, 0, ref, tA20);
  Teko::setBlock(2, 1, ref, tA21);
  Teko::setBlock(2, 2, ref, tA22);

  Teko::endBlockFill(ref);

  return ref;
}

Teko::LinearOp buildReferenceAugmentedRHSOperator(const RCP<const crs_t>& fullA,
                                                  const std::vector<std::vector<GO> >& blockedVec,
                                                  const Teko::LinearOp& pressureMassMatrix,
                                                  double gamma) {
  Teko::TpetraHelpers::BlockedTpetraOperator blk(blockedVec, fullA);

  auto A02 = blk.GetBlock(0, 2);
  auto A12 = blk.GetBlock(1, 2);
  auto A22 = blk.GetBlock(2, 2);

  Teko::LinearOp tA02 = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A02->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A02->getDomainMap()), A02);
  Teko::LinearOp tA12 = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A12->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A12->getDomainMap()), A12);
  Teko::LinearOp tA22 = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A22->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A22->getDomainMap()), A22);

  Teko::LinearOp Winv = Teko::getInvDiagonalOp(pressureMassMatrix);

  auto rhsOp = Teko::createBlockedOp();
  Teko::beginBlockFill(rhsOp, 3, 3);

  Teko::setBlock(0, 0, rhsOp, Thyra::identity<ST>(tA02->range()));
  Teko::setBlock(1, 1, rhsOp, Thyra::identity<ST>(tA12->range()));
  Teko::setBlock(2, 2, rhsOp, Thyra::identity<ST>(tA22->range()));

  Teko::setBlock(0, 2, rhsOp, Thyra::scale(gamma, Thyra::multiply(tA02, Winv)));
  Teko::setBlock(1, 2, rhsOp, Thyra::scale(gamma, Thyra::multiply(tA12, Winv)));

  Teko::endBlockFill(rhsOp);

  return rhsOp;
}

void printOriginalBlocks(const Teko::TpetraHelpers::BlockedTpetraOperator& blk,
                         Teuchos::FancyOStream& out) {
  out << "Original operator block Frobenius norms:" << std::endl;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      auto op  = blk.GetBlock(i, j);
      auto crs = Teuchos::rcp_dynamic_cast<const crs_t>(op, true);
      out << "  block(" << i << "," << j << ") = " << crs->getFrobeniusNorm() << std::endl;
    }
  }
}

}  // namespace

TEUCHOS_UNIT_TEST(tALOperator, test_tpetra_apply_and_rhs) {
  auto comm = Tpetra::getDefaultComm();

  const int myPID = comm->getRank();
  out << "MPI_PID = " << myPID << ", UNIX_PID = " << getpid() << std::endl;

  std::vector<std::vector<GO> > blockedVec;
  Teko::LinearOp pressureMassMatrix;
  RCP<crs_t> mat = buildFullSystem(comm, blockedVec, pressureMassMatrix);
  TEST_ASSERT(!mat.is_null());
  TEST_ASSERT(pressureMassMatrix != Teuchos::null);

  double gamma = 0.05;

  Teko::NS::ALOperator al(blockedVec, mat, pressureMassMatrix, gamma);
  Teko::TpetraHelpers::BlockedTpetraOperator blk(blockedVec, mat);

  printOriginalBlocks(blk, out);

  Teko::LinearOp refApplyOp = buildReferenceALOperator(mat, blockedVec, pressureMassMatrix, gamma);
  Teko::LinearOp refRhsOp =
      buildReferenceAugmentedRHSOperator(mat, blockedVec, pressureMassMatrix, gamma);

  {
    Thyra::LinearOpTester<ST> tester;
    tester.show_all_tests(true);
    tester.dump_all(true);

    bool result = tester.compare(*al.getThyraOp(), *refApplyOp, Teuchos::ptrFromRef(out));
    if (!result) {
      out << "Blocked Thyra operator comparison FAILED." << std::endl;
    } else {
      out << "Blocked Thyra operator comparison PASSED." << std::endl;
    }
    TEST_ASSERT(result);
  }

  {
    vec_t x(mat->getDomainMap()), y(mat->getRangeMap());
    x.putScalar(1.0);

    al.apply(x, y);

    RCP<const vec_t> yRef =
        applyReferenceBlockedOpByBlockedVec(refApplyOp, blockedVec, mat->getRangeMap());

    auto diff = rcp(new vec_t(y.getMap()));
    diff->assign(y);
    diff->update(1.0, *yRef, -1.0);

    printVector("Computed AL apply", Teuchos::rcpFromRef(y), out);
    printVector("Reference AL apply", yRef, out);
    printVector("Difference", diff, out);

    double dn  = diff->norm2();
    double yn  = yRef->norm2();
    double rel = (yn == 0.0 ? dn : dn / yn);
    if (rel < 1.0e-14) {
      out << "Test:ALOperator(Tpetra Apply): Passed." << std::endl;
    } else {
      out << "Test:ALOperator(Tpetra Apply): Failed. rel = " << rel << std::endl;
    }

    TEST_ASSERT(rel < 1.0e-14);
  }

  {
    vec_t rhs(mat->getRangeMap()), rhsAug(mat->getRangeMap());
    rhs.putScalar(1.0);

    al.augmentRHS(rhs, rhsAug);

    RCP<const vec_t> rhsRef =
        applyReferenceBlockedOpByBlockedVec(refRhsOp, blockedVec, mat->getRangeMap());

    auto diff = rcp(new vec_t(rhsAug.getMap()));
    diff->assign(rhsAug);
    diff->update(1.0, *rhsRef, -1.0);

    printVector("Computed augmented RHS", Teuchos::rcpFromRef(rhsAug), out);
    printVector("Reference augmented RHS", rhsRef, out);
    printVector("Difference RHS", diff, out);

    double dn  = diff->norm2();
    double yn  = rhsRef->norm2();
    double rel = (yn == 0.0 ? dn : dn / yn);
    if (rel < 1.0e-14) {
      out << "Test:ALOperator(Tpetra augmentRHS): Passed." << std::endl;
    } else {
      out << "Test:ALOperator(Tpetra augmentRHS): Failed. rel = " << rel << std::endl;
    }

    TEST_ASSERT(rel < 1.0e-14);
  }
}
