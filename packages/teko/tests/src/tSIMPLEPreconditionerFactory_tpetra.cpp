// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_Config.h"
#include "tSIMPLEPreconditionerFactory_tpetra.hpp"
#include "Teko_SIMPLEPreconditionerFactory.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_InverseLibrary.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Array.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraVectorSpace.hpp"
#include "Teko_Utilities.hpp"

#include <vector>

// This whole test rig is based on inverting the matrix
//
//      [  1  2  1 -1 ]
//  A = [  2  1 -3  1 ]
//      [  1 -3  1  2 ]
//      [ -1  1  1  2 ]
//
// see the matlab file

namespace Teko {
namespace Test {

using namespace Teuchos;
using namespace Thyra;
using namespace Teko::NS;

namespace {

using ST    = Teko::ST;
using LO    = Teko::LO;
using GO    = Teko::GO;
using NT    = Teko::NT;
using vec_t = Tpetra::Vector<ST, LO, GO, NT>;

Teuchos::Array<int> buildBlockStartsArray(const Teuchos::ArrayRCP<GO>& starts) {
  Teuchos::Array<int> arr(starts.size());
  for (Teuchos::ArrayRCP<GO>::size_type i = 0; i < starts.size(); ++i)
    arr[i] = static_cast<int>(starts[i]);
  return arr;
}

Teuchos::Array<GO> buildBlockGidsArray(const Teuchos::ArrayRCP<GO>& gids) {
  Teuchos::Array<GO> arr(gids.size());
  for (Teuchos::ArrayRCP<GO>::size_type i = 0; i < gids.size(); ++i) arr[i] = gids[i];
  return arr;
}

}  // namespace

void tSIMPLEPreconditionerFactory_tpetra::initializeTest() {
  std::vector<GO> indices(2);
  std::vector<ST> row0(2), row1(2);

  tolerance_ = 1.0e-13;

  comm = GetComm_tpetra();

  const RCP<const Tpetra::Map<LO, GO, NT> > map =
      rcp(new const Tpetra::Map<LO, GO, NT>(2, 0, comm));

  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrF =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrB =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrBt =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrC =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);

  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrInvF =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrInvS =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);

  indices[0] = 0;
  indices[1] = 1;

  // build F matrix
  row0[0] = 1.0;
  row0[1] = 2.0;
  row1[0] = 2.0;
  row1[1] = 1.0;
  ptrF->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrF->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrF->fillComplete();
  F_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrF->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrF->getRangeMap()), ptrF);

  // build block info for block diagonal (one block)
  block_starts_    = arcp(new GO[2], 0, 2);
  block_gids_      = arcp(new GO[2], 0, 2);
  GO* block_starts = block_starts_.getRawPtr();
  GO* block_gids   = block_gids_.getRawPtr();
  block_starts[0]  = 0;
  block_starts[1]  = 2;
  block_gids[0]    = ptrF->getRowMap()->getGlobalElement(0);
  block_gids[1]    = ptrF->getRowMap()->getGlobalElement(1);

  // build B matrix
  row0[0] = 1.0;
  row0[1] = -3.0;
  row1[0] = -1.0;
  row1[1] = 1.0;
  ptrB->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrB->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrB->fillComplete();
  B_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrB->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrB->getRangeMap()), ptrB);

  // build Bt matrix
  row0[0] = 1.0;
  row0[1] = -1.0;
  row1[0] = -3.0;
  row1[1] = 1.0;
  ptrBt->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrBt->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrBt->fillComplete();
  Bt_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrBt->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrBt->getRangeMap()), ptrBt);

  // build C matrix
  row0[0] = 1.0;
  row0[1] = 2.0;
  row1[0] = 2.0;
  row1[1] = 1.0;
  ptrC->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrC->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrC->fillComplete();
  C_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrC->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrC->getRangeMap()), ptrC);

  // build inv(F) matrix
  row0[0] = -1.0 / 3.0;
  row0[1] = 2.0 / 3.0;
  row1[0] = 2.0 / 3.0;
  row1[1] = -1.0 / 3.0;
  ptrInvF->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrInvF->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrInvF->fillComplete();
  invF_ = rcp(new StaticOpInverseFactory(Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvF->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvF->getRangeMap()), ptrInvF)));

  // build inv(Pschur) matrix
  row0[0] = 0.037037037037037;
  row0[1] = 0.222222222222222;
  row1[0] = 0.222222222222222;
  row1[1] = 0.333333333333333;
  ptrInvS->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrInvS->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrInvS->fillComplete();
  invS_ = rcp(new StaticOpInverseFactory(Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvS->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvS->getRangeMap()), ptrInvS)));

  A_ = Thyra::block2x2<ST>(F_, Bt_, B_, C_, "A");
}

int tSIMPLEPreconditionerFactory_tpetra::runTest(int verbosity, std::ostream& stdstrm,
                                                 std::ostream& failstrm, int& totalrun) {
  bool allTests = true;
  bool status   = true;
  int failcount = 0;

  failstrm << "tSIMPLEPreconditionerFactory_tpetra";

  status = test_createPrec(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"createPrec\" ... PASSED", "   \"createPrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm, -2);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"initializePrec(lumped)\" ... PASSED",
                       "   \"initializePrec(lumped)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm, -1);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"initializePrec(diag)\" ... PASSED",
                       "   \"initializePrec(diag)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm, 0);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"initializePrec(absrowsum)\" ... PASSED",
                       "   \"initializePrec(absrowsum)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm, 1);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"initializePrec(block-1)\" ... PASSED",
                       "   \"initializePrec(block-1)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm, 2);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"initializePrec(block-2)\" ... PASSED",
                       "   \"initializePrec(block-2)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_uninitializePrec(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"uninitializePrec\" ... PASSED",
                       "   \"uninitializePrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_isCompatable(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"isCompatable\" ... PASSED",
                       "   \"isCompatable\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_iterativeSolves(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"iterativeSolves\" ... PASSED",
                       "   \"iterativeSolves\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_hierarchicalSolves(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"hierarchicalSolves\" ... PASSED",
                       "   \"hierarchicalSolves\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_diagonal(verbosity, failstrm, 0);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"diagonal(diag)\" ... PASSED",
                       "   \"diagonal(diag)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_diagonal(verbosity, failstrm, 1);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"diagonal(block-1)\" ... PASSED",
                       "   \"diagonal(block-1)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_diagonal(verbosity, failstrm, 2);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"diagonal(block-2)\" ... PASSED",
                       "   \"diagonal(block-2)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_result(verbosity, failstrm, 0);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"result(diag)\" ... PASSED",
                       "   \"result(diag)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_result(verbosity, failstrm, 1);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"result(block-1)\" ... PASSED",
                       "   \"result(block-1)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_result(verbosity, failstrm, 2);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"result(block-2)\" ... PASSED",
                       "   \"result(block-2)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tSIMPLEPreconditionedFactory...PASSED",
                         "tSIMPLEPreconditionedFactory...FAILED");
  } else {
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tSIMPLEPreconditionedFactory...FAILED");
  }

  return failcount;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_createPrec(int verbosity, std::ostream& os) {
  const RCP<const Thyra::PreconditionerFactoryBase<double> > fact =
      rcp(new SIMPLEPreconditionerFactory(invF_, 0.9));

  try {
    rcp_dynamic_cast<DefaultPreconditioner<double> >(fact->createPrec(), true);
  } catch (std::exception& e) {
    os << std::endl
       << "   test_createPrec: dynamic cast to \"DefaultPreconditioner\" FAILED" << std::endl;
    os << "   Descriptive exception \"" << e.what() << "\"" << std::endl;
    return false;
  }

  return true;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_initializePrec(int verbosity, std::ostream& os,
                                                              int use_blocking) {
  bool status    = false;
  bool allPassed = true;

  RCP<SIMPLEPreconditionerFactory> sFactory = rcp(new SIMPLEPreconditionerFactory(invF_, 0.9));
  const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory = sFactory;
  RCP<Thyra::PreconditionerBase<double> > prec = precFactory->createPrec();

  if (use_blocking == -2) {
    Teuchos::ParameterList List;
    List.set("Explicit Velocity Inverse Type", "Lumped");
    List.set("Inverse Pressure Type", "Ifpack2");
    List.set("Inverse Velocity Type", "Ifpack2");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == -1) {
    Teuchos::ParameterList List;
    List.set("Explicit Velocity Inverse Type", "Diagonal");
    List.set("Inverse Pressure Type", "Ifpack2");
    List.set("Inverse Velocity Type", "Ifpack2");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == 0) {
    Teuchos::ParameterList List;
    List.set("Explicit Velocity Inverse Type", "AbsRowSum");
    List.set("Inverse Pressure Type", "Ifpack2");
    List.set("Inverse Velocity Type", "Ifpack2");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == 1) {
    Teuchos::ParameterList List, BlkList;
    BlkList.set("number of local blocks", 1);
    BlkList.set("block start index", buildBlockStartsArray(block_starts_));
    BlkList.set("block entry gids", buildBlockGidsArray(block_gids_));
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Ifpack2");
    List.set("Inverse Velocity Type", "Ifpack2");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == 2) {
    Teuchos::ParameterList List, BlkList;
    BlkList.set("contiguous block size", 2);
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Ifpack2");
    List.set("Inverse Velocity Type", "Ifpack2");
    sFactory->initializeFromParameterList(List);
  }

  precFactory->initializePrec(Thyra::defaultLinearOpSource(A_), &*prec);

  RCP<const Thyra::LinearOpBase<double> > op;

  op     = prec->getUnspecifiedPrecOp();
  status = (op != Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory_tpetra::test_initializePrec " << toString(status)
       << std::endl;
    os << "      "
       << "Preconditioner \"getUnspecifiedPrecOp\" is null (it should not be!)" << std::endl;
  }
  allPassed &= status;

  op     = prec->getRightPrecOp();
  status = (op == Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory_tpetra::test_initializePrec " << toString(status)
       << std::endl;
    os << "      "
       << "Preconditioner \"getRightPrecOp\" is not null (it should be!)" << std::endl;
  }
  allPassed &= status;

  op     = prec->getLeftPrecOp();
  status = (op == Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory_tpetra::test_initializePrec " << toString(status)
       << std::endl;
    os << "      "
       << "Preconditioner \"getLeftPrecOp\" is not null (it should be!)" << std::endl;
  }
  allPassed &= status;

  return allPassed;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_uninitializePrec(int verbosity, std::ostream& os) {
  return true;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_isCompatable(int verbosity, std::ostream& os) {
  return true;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_diagonal(int verbosity, std::ostream& os,
                                                        int use_blocking) {
  typedef RCP<const Thyra::LinearOpBase<ST> > LinearOp;

  bool status    = false;
  bool allPassed = true;
  double diff    = 0.0;

  std::vector<GO> indices(2);
  indices[0] = 0;
  indices[1] = 1;

  auto map = rcp(new const Tpetra::Map<LO, GO, NT>(2, 0, comm));

  auto buildDiag = [&](ST a0, ST a1) -> LinearOp {
    auto mat = Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 1);
    std::vector<GO> c0(1), c1(1);
    std::vector<ST> r0(1), r1(1);
    c0[0] = 0;
    c1[0] = 1;
    r0[0] = a0;
    r1[0] = a1;
    mat->insertGlobalValues(0, Teuchos::ArrayView<GO>(c0), Teuchos::ArrayView<ST>(r0));
    mat->insertGlobalValues(1, Teuchos::ArrayView<GO>(c1), Teuchos::ArrayView<ST>(r1));
    mat->fillComplete();
    return Thyra::tpetraLinearOp<ST, LO, GO, NT>(
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat->getDomainMap()),
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat->getRangeMap()), mat);
  };

  LinearOp F  = buildDiag(1.0, 2.0);
  LinearOp G  = buildDiag(7.0, 8.0);
  LinearOp D  = buildDiag(5.0, 6.0);
  LinearOp C  = buildDiag(3.0, 4.0);
  LinearOp iF = buildDiag(1.0, 0.5);
  LinearOp iS = buildDiag(-0.03125, -0.05);

  RCP<Teko::InverseFactory> invF = rcp(new Teko::StaticOpInverseFactory(iF));
  RCP<Teko::InverseFactory> invS = rcp(new Teko::StaticOpInverseFactory(iS));

  LinearOp A                                = Thyra::block2x2(F, G, D, C);
  RCP<SIMPLEPreconditionerFactory> sFactory = rcp(new SIMPLEPreconditionerFactory(invF, invS, 0.9));
  const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory = sFactory;
  RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory, A);

  if (use_blocking == 1) {
    Teuchos::ParameterList List, BlkList;
    BlkList.set("number of local blocks", 1);
    BlkList.set("block start index", buildBlockStartsArray(block_starts_));
    BlkList.set("block entry gids", buildBlockGidsArray(block_gids_));
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Amesos2");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == 2) {
    Teuchos::ParameterList List, BlkList;
    BlkList.set("contiguous block size", 2);
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Amesos2");
    sFactory->initializeFromParameterList(List);
  }

  RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp();

  auto ea = rcp(new vec_t(map));
  auto eb = rcp(new vec_t(map));
  auto ef = rcp(new vec_t(map));
  auto eg = rcp(new vec_t(map));

  const RCP<const Thyra::MultiVectorBase<ST> > x = Teko::Test::BlockVector(*ea, *eb, A->domain());
  const RCP<const Thyra::MultiVectorBase<ST> > z = Teko::Test::BlockVector(*ef, *eg, A->domain());
  const RCP<Thyra::MultiVectorBase<ST> > y       = Thyra::createMembers(A->range(), 1);

  auto setVec = [](const RCP<vec_t>& v, ST a0, ST a1) {
    v->replaceGlobalValue(0, a0);
    v->replaceGlobalValue(1, a1);
  };

  setVec(ea, 0.0, 1.0);
  setVec(eb, 1.0, 3.0);
  setVec(ef, 0.21875, 0.5);
  setVec(eg, -0.028125, 0.0);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory_tpetra::test_diagonal " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    Teko::Test::Print(os, "x", x);
    Teko::Test::Print(os, "y", y);
    Teko::Test::Print(os, "z", z);
  }
  allPassed &= status;

  setVec(ea, -2.0, 4.0);
  setVec(eb, 7.0, 9.0);
  setVec(ef, 1.71875, 1.4);
  setVec(eg, -0.478125, 0.135);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory_tpetra::test_diagonal " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    Teko::Test::Print(os, "x", x);
    Teko::Test::Print(os, "y", y);
    Teko::Test::Print(os, "z", z);
  }
  allPassed &= status;

  setVec(ea, 1.0, 0.0);
  setVec(eb, 0.0, -5.0);
  setVec(ef, -0.09375, -1.0);
  setVec(eg, 0.140625, 0.225);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory_tpetra::test_diagonal " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    Teko::Test::Print(os, "x", x);
    Teko::Test::Print(os, "y", y);
    Teko::Test::Print(os, "z", z);
  }
  allPassed &= status;

  setVec(ea, 4.0, -4.0);
  setVec(eb, 6.0, 12.0);
  setVec(ef, 0.9375, 2.800000000000001);
  setVec(eg, 0.39375, -1.08);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory_tpetra::test_diagonal " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    Teko::Test::Print(os, "x", x);
    Teko::Test::Print(os, "y", y);
    Teko::Test::Print(os, "z", z);
  }
  allPassed &= status;

  return allPassed;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_result(int verbosity, std::ostream& os,
                                                      int use_blocking) {
  bool status    = false;
  bool allPassed = true;
  double diff    = -1000.0;

  RCP<SIMPLEPreconditionerFactory> sFactory =
      rcp(new SIMPLEPreconditionerFactory(invF_, invS_, 0.9));
  const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory = sFactory;
  RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory, A_);

  if (use_blocking == 1) {
    Teuchos::ParameterList List, BlkList;
    BlkList.set("number of local blocks", 1);
    BlkList.set("block start index", buildBlockStartsArray(block_starts_));
    BlkList.set("block entry gids", buildBlockGidsArray(block_gids_));
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Amesos2");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == 2) {
    Teuchos::ParameterList List, BlkList;
    BlkList.set("contiguous block size", 2);
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Amesos2");
    sFactory->initializeFromParameterList(List);
  }

  RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp();

  auto map = rcp(new const Tpetra::Map<LO, GO, NT>(2, 0, comm));

  auto ea = rcp(new vec_t(map));
  auto eb = rcp(new vec_t(map));
  auto ef = rcp(new vec_t(map));
  auto eg = rcp(new vec_t(map));

  const RCP<const Thyra::MultiVectorBase<ST> > x = Teko::Test::BlockVector(*ea, *eb, A_->domain());
  const RCP<const Thyra::MultiVectorBase<ST> > z = Teko::Test::BlockVector(*ef, *eg, A_->domain());
  const RCP<Thyra::MultiVectorBase<ST> > y       = Thyra::createMembers(A_->range(), 1);

  auto setVec = [](const RCP<vec_t>& v, ST a0, ST a1) {
    v->replaceGlobalValue(0, a0);
    v->replaceGlobalValue(1, a1);
  };

  setVec(ea, 0.0, 1.0);
  setVec(eb, 1.0, 3.0);
  setVec(ef, 0.987654320987654, 1.074074074074074);
  setVec(eg, 0.777777777777778, 1.066666666666667);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory_tpetra::test_result " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    Teko::Test::Print(os, "x", x);
    Teko::Test::Print(os, "y", y);
    Teko::Test::Print(os, "z", z);
  }
  allPassed &= status;

  setVec(ea, -2.0, 4.0);
  setVec(eb, 7.0, 9.0);
  setVec(ef, 4.197530864197531, 2.814814814814815);
  setVec(eg, 2.855555555555555, 3.633333333333334);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory_tpetra::test_result " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    Teko::Test::Print(os, "x", x);
    Teko::Test::Print(os, "y", y);
    Teko::Test::Print(os, "z", z);
  }
  allPassed &= status;

  setVec(ea, 1.0, 0.0);
  setVec(eb, 0.0, -5.0);
  setVec(ef, -0.567901234567901, -1.592592592592592);
  setVec(eg, -1.122222222222222, -1.333333333333333);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory_tpetra::test_result " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    Teko::Test::Print(os, "x", x);
    Teko::Test::Print(os, "y", y);
    Teko::Test::Print(os, "z", z);
  }
  allPassed &= status;

  setVec(ea, 4.0, -4.0);
  setVec(eb, 6.0, 12.0);
  setVec(ef, 0.518518518518519, 2.888888888888889);
  setVec(eg, 1.533333333333334, 5.600000000000001);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory_tpetra::test_result " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    Teko::Test::Print(os, "x", x);
    Teko::Test::Print(os, "y", y);
    Teko::Test::Print(os, "z", z);
  }
  allPassed &= status;

  return allPassed;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_iterativeSolves(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  RCP<ParameterList> params = Teuchos::rcp(new ParameterList());
  ParameterList& tekoList   = params->sublist("Preconditioner Types").sublist("Teko");
  tekoList.set("Inverse Type", "SIMPLE");
  ParameterList& ifl = tekoList.sublist("Inverse Factory Library");
  ifl.sublist("SIMPLE").set("Type", "NS SIMPLE");
  ifl.sublist("SIMPLE").set("Inverse Velocity Type", "Belos");
  ifl.sublist("SIMPLE").set("Preconditioner Velocity Type", "Ifpack2");
  ifl.sublist("SIMPLE").set("Inverse Pressure Type", "Belos");
  ifl.sublist("SIMPLE").set("Preconditioner Pressure Type", "Ifpack2");

  RCP<Teko::InverseLibrary> invLib        = Teko::InverseLibrary::buildFromParameterList(ifl);
  RCP<const Teko::InverseFactory> invFact = invLib->getInverseFactory("SIMPLE");

  Teko::ModifiableLinearOp inv = Teko::buildInverse(*invFact, A_);
  TEST_ASSERT(!inv.is_null(), "Constructed preconditioner is null");

  if (!inv.is_null()) {
    Teko::rebuildInverse(*invFact, A_, inv);
    TEST_ASSERT(!inv.is_null(), "Constructed preconditioner is null");
  }

  return true;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_hierarchicalSolves(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  const RCP<Thyra::PhysicallyBlockedLinearOpBase<ST> > blkOp = Thyra::defaultBlockedLinearOp<ST>();
  blkOp->beginBlockFill(3, 3);
  blkOp->setBlock(0, 0, F_);
  blkOp->setBlock(0, 1, Bt_);
  blkOp->setBlock(1, 0, B_);
  blkOp->setBlock(1, 1, C_);
  blkOp->setBlock(1, 2, B_);
  blkOp->setBlock(2, 1, Bt_);
  blkOp->setBlock(2, 2, F_);
  blkOp->endBlockFill();

  RCP<ParameterList> params = Teuchos::rcp(new ParameterList());
  ParameterList& tekoList   = params->sublist("Preconditioner Types").sublist("Teko");
  tekoList.set("Inverse Type", "HBGS");
  ParameterList& ifl = tekoList.sublist("Inverse Factory Library");

  ifl.sublist("HBGS").set("Type", "Hierarchical Block Gauss-Seidel");

  ifl.sublist("HBGS").sublist("Hierarchical Block 1").set("Included Subblocks", "1,2");
  ifl.sublist("HBGS").sublist("Hierarchical Block 1").set("Inverse Type", "Belos");
  ifl.sublist("HBGS").sublist("Hierarchical Block 1").set("Preconditioner Type", "SIMPLE");

  ifl.sublist("HBGS").sublist("Hierarchical Block 2").set("Included Subblocks", "3");
  ifl.sublist("HBGS").sublist("Hierarchical Block 2").set("Inverse Type", "Belos");
  ifl.sublist("HBGS").sublist("Hierarchical Block 2").set("Preconditioner Type", "SINGLE_BLOCK");

  ifl.sublist("SIMPLE").set("Type", "NS SIMPLE");
  ifl.sublist("SIMPLE").set("Inverse Velocity Type", "Belos");
  ifl.sublist("SIMPLE").set("Preconditioner Velocity Type", "Ifpack2");
  ifl.sublist("SIMPLE").set("Inverse Pressure Type", "Belos");
  ifl.sublist("SIMPLE").set("Preconditioner Pressure Type", "Ifpack2");

  ifl.sublist("SINGLE_BLOCK").set("Type", "Ifpack2");

  RCP<Teko::InverseLibrary> invLib        = Teko::InverseLibrary::buildFromParameterList(ifl);
  RCP<const Teko::InverseFactory> invFact = invLib->getInverseFactory("HBGS");

  Teko::ModifiableLinearOp inv = Teko::buildInverse(*invFact, blkOp);
  TEST_ASSERT(!inv.is_null(), "Constructed preconditioner is null");

  if (!inv.is_null()) {
    Teko::rebuildInverse(*invFact, blkOp, inv);
    TEST_ASSERT(!inv.is_null(), "Constructed preconditioner is null");
  }

  return true;
}

}  // end namespace Test
}  // end namespace Teko
