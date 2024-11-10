// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tExplicitOps_tpetra.hpp"

#include <string>

// Tpetra includes
#include "Tpetra_Export.hpp"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_TpetraHelpers.hpp"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

// Test-rig
#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::epetraLinearOp;

void tExplicitOps_tpetra::initializeTest() {
  const Epetra_Comm& comm_epetra                   = *GetComm();
  const RCP<const Teuchos::Comm<int> > comm_tpetra = GetComm_tpetra();

  tolerance_ = 1.0e-4;

  int nx = 39;  // essentially random values
  int ny = 53;

  // create some big blocks to play with
  Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d", comm_epetra, false);
  FGallery.Set("nx", nx);
  FGallery.Set("ny", ny);
  Epetra_CrsMatrix& epetraF = FGallery.GetMatrixRef();
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > tpetraF =
      Teko::TpetraHelpers::nonConstEpetraCrsMatrixToTpetra(rcpFromRef(epetraF), comm_tpetra);
  F_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraF->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraF->getRangeMap()), tpetraF);

  // create some big blocks to play with
  Trilinos_Util::CrsMatrixGallery GGallery("laplace_2d", comm_epetra, false);
  GGallery.Set("nx", nx);
  GGallery.Set("ny", ny);
  Epetra_CrsMatrix& epetraG = GGallery.GetMatrixRef();
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > tpetraG =
      Teko::TpetraHelpers::nonConstEpetraCrsMatrixToTpetra(rcpFromRef(epetraG), comm_tpetra);
  G_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraG->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraG->getRangeMap()), tpetraG);

  // create some big blocks to play with
  Trilinos_Util::CrsMatrixGallery HGallery("eye", comm_epetra, false);
  HGallery.Set("nx", nx * ny);
  // HGallery.Set("ny",ny);
  Epetra_CrsMatrix& epetraH = HGallery.GetMatrixRef();
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > tpetraH =
      Teko::TpetraHelpers::nonConstEpetraCrsMatrixToTpetra(rcpFromRef(epetraH), comm_tpetra);
  H_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraH->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraH->getRangeMap()), tpetraH);

  RCP<Tpetra::Vector<ST, LO, GO, NT> > v =
      rcp(new Tpetra::Vector<ST, LO, GO, NT>(tpetraF->getRangeMap()));
  v->randomize();
  RCP<Thyra::VectorBase<ST> > tV = Thyra::createVector<ST, LO, GO, NT>(
      v, Thyra::createVectorSpace<ST, LO, GO, NT>(tpetraF->getRowMap()));
  D_ = Thyra::diagonal(tV);
}

int tExplicitOps_tpetra::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                                 int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tExplicitOps_tpetra";

  status = test_mult_diagScaleMatProd(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"mult_diagScaleMatProd\" ... PASSED",
                "   \"mult_diagScaleMatProd\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_mult_diagScaling(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"mult_diagScaling\" ... PASSED",
                "   \"mult_diagScaling\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_add(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"add\" ... PASSED", "   \"add\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_mult_modScaleMatProd(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"mult_modScaleMatProd\" ... PASSED",
                "   \"mult_modScaleMatProd\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_add_mod(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"add_mod\" ... PASSED", "   \"add\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tExplicitOps_tpetra...PASSED", "tExplicitOps_tpetra...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tExplicitOps_tpetra...FAILED");
  }

  return failcount;
}

bool tExplicitOps_tpetra::test_mult_diagScaleMatProd(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  Thyra::LinearOpTester<ST> tester;
  tester.set_all_error_tol(1e-10);
  tester.show_all_tests(true);

  RCP<const Thyra::LinearOpBase<ST> > thyOp;
  Teko::LinearOp expOp;

  thyOp = Thyra::multiply(Teko::scale(-4.0, F_), D_, Teko::adjoint(G_));
  expOp = Teko::explicitMultiply(Teko::scale(-4.0, F_), D_, Teko::adjoint(G_));

  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*thyOp, *expOp, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tExplicitOps_tpetra::test_diagScaleMatProd "
                            << ": Testing triple matrix product");
    if (not result || verbosity >= 10) os << ss.str();
  }

  thyOp = Teko::multiply(Teko::scale(-4.0, F_), Teko::adjoint(G_));
  expOp = Teko::explicitMultiply(Teko::scale(-4.0, F_), Teko::adjoint(G_));

  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*thyOp, *expOp, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tExplicitOps_tpetra::test_diagScaleMatProd "
                            << ": Testing triple matrix product");
    if (not result || verbosity >= 10) os << ss.str();
  }

  return allPassed;
}

bool tExplicitOps_tpetra::test_mult_diagScaling(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  Thyra::LinearOpTester<ST> tester;
  tester.set_all_error_tol(1e-10);
  tester.show_all_tests(true);

  RCP<const Thyra::LinearOpBase<ST> > thyOp;
  Teko::LinearOp expOp;

  thyOp = Teko::multiply(Teko::scale(-4.0, F_), D_);
  expOp = Teko::explicitMultiply(Teko::scale(-4.0, F_), D_);

  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*thyOp, *expOp, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tExplicitOps_tpetra::test_diagScaleMatProd "
                            << ": Testing diagonal scaling");
    if (not result || verbosity >= 10) os << ss.str();
  }

  thyOp = Teko::multiply(D_, Teko::scale(-9.0, F_));
  expOp = Teko::explicitMultiply(D_, Teko::scale(-9.0, F_));

  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*thyOp, *expOp, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tExplicitOps_tpetra::test_diagScaleMatProd "
                            << ": Testing diagonal scaling");
    if (not result || verbosity >= 10) os << ss.str();
  }

  return allPassed;
}

bool tExplicitOps_tpetra::test_mult_modScaleMatProd(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  Thyra::LinearOpTester<ST> tester;
  tester.set_all_error_tol(1e-10);
  tester.show_all_tests(true);

  Teko::LinearOp thyOp;
  Teko::ModifiableLinearOp expOp;

  thyOp = Teko::multiply(Teko::scale(-4.0, F_), D_, Teko::adjoint(G_));
  expOp = Teko::explicitMultiply(Teko::scale(-4.0, F_), D_, Teko::adjoint(G_), expOp);

  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp1 =
      Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(expOp, true);
  RCP<const Tpetra::Operator<ST, LO, GO, NT> > eop1 = tOp1->getConstTpetraOperator();

  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*thyOp, *expOp, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tExplicitOps_tpetra::test_modScaleMatProd "
                            << ": Testing triple matrix product");
    if (not result || verbosity >= 10) os << ss.str();
  }

  RCP<Thyra::TpetraLinearOp<ST, LO, GO, NT> > tF =
      Teuchos::rcp_dynamic_cast<Thyra::TpetraLinearOp<ST, LO, GO, NT> >(F_);
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > crsF =
      Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<ST, LO, GO, NT> >(tF->getTpetraOperator(), true);
  crsF->scale(5.0);
  RCP<Thyra::TpetraLinearOp<ST, LO, GO, NT> > tG =
      Teuchos::rcp_dynamic_cast<Thyra::TpetraLinearOp<ST, LO, GO, NT> >(G_);
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > crsG =
      Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<ST, LO, GO, NT> >(tG->getTpetraOperator(), true);
  crsG->scale(2.0);

  // do some random violence (oh my brothers) to one row
  size_t numEntries = crsF->getNumEntriesInLocalRow(3);
  auto indices1 = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_local_inds_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), numEntries);
  auto values1 = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_values_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), numEntries);
  crsF->getLocalRowCopy(3, indices1, values1, numEntries);
  for (size_t i = 0; i < numEntries; i++) values1(i) *= values1(i) * ST(i + 1) * 0.92;
  crsF->replaceLocalValues(3, indices1, values1);

  // do some random violence (oh my brothers) to one row
  numEntries    = crsF->getNumEntriesInLocalRow(7);
  auto indices2 = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_local_inds_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), numEntries);
  auto values2 = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_values_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), numEntries);
  crsF->getLocalRowCopy(7, indices2, values2, numEntries);
  for (size_t i = 0; i < numEntries; i++) values2(i) *= values2(i) * ST(i + 1) * 0.92;
  crsF->replaceLocalValues(7, indices2, values2);

  // perform the next test
  thyOp = Teko::multiply(Teko::scale(-4.0, F_), D_, Teko::adjoint(G_));
  expOp = Teko::explicitMultiply(Teko::scale(-4.0, F_), D_, Teko::adjoint(G_), expOp);

  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp2 =
      Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(expOp, true);
  RCP<const Tpetra::Operator<ST, LO, GO, NT> > eop2 = tOp2->getConstTpetraOperator();

  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*thyOp, *expOp, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tExplicitOps_tpetra::test_modScaleMatProd "
                            << ": Testing triple matrix product");
    if (not result || verbosity >= 10) os << ss.str();
  }

  return allPassed;
}

bool tExplicitOps_tpetra::test_add(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  Thyra::LinearOpTester<ST> tester;
  tester.set_all_error_tol(1e-10);
  tester.show_all_tests(true);

  RCP<const Thyra::LinearOpBase<ST> > thyOp;
  Teko::LinearOp expOp;

  thyOp = Teko::add(Teko::scale(-4.0, F_), Teko::adjoint(G_));
  expOp = Teko::explicitAdd(Teko::scale(-4.0, F_), Teko::adjoint(G_));

  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*thyOp, *expOp, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tExplicitOps_tpetra::test_add "
                            << ": Testing explicit add");
    if (not result || verbosity >= 10) os << ss.str();
  }

  return allPassed;
}

bool tExplicitOps_tpetra::test_add_mod(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  Thyra::LinearOpTester<ST> tester;
  tester.set_all_error_tol(1e-10);
  tester.show_all_tests(true);

  RCP<const Thyra::LinearOpBase<ST> > thyOp;
  Teko::ModifiableLinearOp expOp;

  thyOp = Teko::add(Teko::scale(-4.0, F_), Teko::adjoint(G_));
  expOp = Teko::explicitAdd(Teko::scale(-4.0, F_), Teko::adjoint(G_), expOp);

  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp1 =
      Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(expOp, true);
  RCP<const Tpetra::Operator<ST, LO, GO, NT> > eop1 = tOp1->getConstTpetraOperator();

  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*thyOp, *expOp, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tExplicitOps_tpetra::test_add_mod"
                            << ": Testing explicit add");
    if (not result || verbosity >= 10) os << ss.str();
  }

  RCP<Thyra::TpetraLinearOp<ST, LO, GO, NT> > tF =
      Teuchos::rcp_dynamic_cast<Thyra::TpetraLinearOp<ST, LO, GO, NT> >(F_);
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > crsF =
      Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<ST, LO, GO, NT> >(tF->getTpetraOperator(), true);
  crsF->scale(5.0);
  RCP<Thyra::TpetraLinearOp<ST, LO, GO, NT> > tG =
      Teuchos::rcp_dynamic_cast<Thyra::TpetraLinearOp<ST, LO, GO, NT> >(G_);
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > crsG =
      Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<ST, LO, GO, NT> >(tG->getTpetraOperator(), true);
  crsG->scale(2.0);

  // do some random violence (oh my brothers) to one row
  size_t numEntries = crsF->getNumEntriesInLocalRow(3);
  auto indices1 = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_local_inds_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), numEntries);
  auto values1 = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_values_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), numEntries);
  crsF->getLocalRowCopy(3, indices1, values1, numEntries);
  for (size_t i = 0; i < numEntries; i++) values1(i) *= values1(i) * ST(i + 1) * 0.92;
  crsF->replaceLocalValues(3, indices1, values1);

  // do some random violence (oh my brothers) to one row
  numEntries    = crsF->getNumEntriesInLocalRow(7);
  auto indices2 = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_local_inds_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), numEntries);
  auto values2 = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_values_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), numEntries);
  crsF->getLocalRowCopy(7, indices2, values2, numEntries);
  for (size_t i = 0; i < numEntries; i++) values2(i) *= values2(i) * ST(i + 1) * 0.92;
  crsF->replaceLocalValues(7, indices2, values2);

  // perform the next test
  thyOp = Teko::add(Teko::scale(-4.0, F_), Teko::adjoint(G_));
  expOp = Teko::explicitAdd(Teko::scale(-4.0, F_), Teko::adjoint(G_), expOp);

  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp2 =
      Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(expOp, true);
  RCP<const Tpetra::Operator<ST, LO, GO, NT> > eop2 = tOp2->getConstTpetraOperator();

  // check that underlying pointers are the same
  // since the sparsity pattern of expOp entering the explicitAdd
  // should match the sparsity pattern returned by the explicitAdd
  // the pointer will be reused
  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    TEST_ASSERT(
        eop2.getRawPtr() == eop1.getRawPtr(), std::endl
                                                  << " tExplicitOps_tpetra::test_add_mod"
                                                  << ": Testing matrix addition preserves pointer");
    if (not(eop2.getRawPtr() == eop1.getRawPtr()) || verbosity >= 10) os << ss.str();
  }

  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*thyOp, *expOp, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tExplicitOps_tpetra::test_add_mod"
                            << ": Testing matrix addition");
    if (not result || verbosity >= 10) os << ss.str();
  }

  Teko::ModifiableLinearOp expOp2;

  // Create a target linear operator, with a known sparsity pattern (diagonal)
  expOp2 = Teko::explicitAdd(H_, H_, expOp2);
  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp3 =
      Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(expOp2, true);
  RCP<const Tpetra::Operator<ST, LO, GO, NT> > eop3 = tOp3->getConstTpetraOperator();

  // Create a target linear operator, with a known sparsity pattern (diagonal) for a second time
  expOp2 = Teko::explicitAdd(H_, H_, expOp2);
  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp3b =
      Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(expOp2, true);
  RCP<const Tpetra::Operator<ST, LO, GO, NT> > eop3b = tOp3b->getConstTpetraOperator();

  // check that underlying pointers are the same
  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    TEST_ASSERT(eop3.getRawPtr() == eop3b.getRawPtr(),
                std::endl
                    << " tExplicitOps_tpetra::test_add_mod"
                    << ": Testing matrix addition returns new pointer");
    if (not(eop3.getRawPtr() == eop3b.getRawPtr()) || verbosity >= 10) os << ss.str();
  }

  // Try to add matricies with sparsity patterns that differ from the target operator
  expOp2 = Teko::explicitAdd(Teko::scale(-4.0, F_), Teko::adjoint(G_), expOp2);
  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp4 =
      Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(expOp2, true);
  RCP<const Tpetra::Operator<ST, LO, GO, NT> > eop4 = tOp4->getConstTpetraOperator();

  // check that underlying pointers are different
  // since the sparsity pattern of expOp entering the explicitAdd
  // does not match the sparsity pattern returned by the explicitAdd
  // a new tpetra operator will be created
  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    TEST_ASSERT(eop3.getRawPtr() != eop4.getRawPtr(),
                std::endl
                    << " tExplicitOps_tpetra::test_add_mod"
                    << ": Testing matrix addition returns new pointer");
    if (not(eop3.getRawPtr() == eop4.getRawPtr()) || verbosity >= 10) os << ss.str();
  }
  // check the result
  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*thyOp, *expOp2, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tExplicitOps_tpetra::test_add_mod"
                            << ": Testing matrix addition");
    if (not result || verbosity >= 10) os << ss.str();
  }

  // Create a target linear operator, with a known sparsity pattern (diagonal) and aliased args
  expOp2 = Teko::explicitAdd(H_, H_, expOp2);
  tOp3   = Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(expOp2, true);
  eop3   = tOp3->getConstTpetraOperator();
  auto expOp3 = Teko::explicitAdd(H_, expOp2, expOp2);
  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp5 =
      Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(expOp3, true);
  RCP<const Tpetra::Operator<ST, LO, GO, NT> > eop5 = tOp5->getConstTpetraOperator();

  // check that underlying pointers are different
  // since the sparsity pattern of expOp entering the explicitAdd
  // does not match the sparsity pattern returned by the explicitAdd
  // a new tpetra operator will be created
  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    TEST_ASSERT(eop3.getRawPtr() != eop5.getRawPtr(),
                std::endl
                    << " tExplicitOps_tpetra::test_add_mod2"
                    << ": Testing matrix addition returns new pointer");
    if (not(eop3.getRawPtr() == eop5.getRawPtr()) || verbosity >= 10) os << ss.str();
  }

  // Create a target linear operator, with a known sparsity pattern (diagonal) and aliased args
  expOp2 = Teko::explicitAdd(H_, H_, expOp2);
  tOp3   = Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(expOp2, true);
  eop3   = tOp3->getConstTpetraOperator();
  auto expOp4 = Teko::explicitAdd(expOp2, H_, expOp2);
  RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp6 =
      Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(expOp4, true);
  RCP<const Tpetra::Operator<ST, LO, GO, NT> > eop6 = tOp6->getConstTpetraOperator();

  // check that underlying pointers are different
  // since the sparsity pattern of expOp entering the explicitAdd
  // does not match the sparsity pattern returned by the explicitAdd
  // a new tpetra operator will be created
  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    TEST_ASSERT(eop3.getRawPtr() != eop6.getRawPtr(),
                std::endl
                    << " tExplicitOps_tpetra::test_add_mod3"
                    << ": Testing matrix addition returns new pointer");
    if (not(eop3.getRawPtr() == eop6.getRawPtr()) || verbosity >= 10) os << ss.str();
  }

  return allPassed;
}

}  // namespace Test
}  // end namespace Teko
