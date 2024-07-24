// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

#include <string>
#include <iostream>

#include "Teko_ConfigDefs.hpp"

#ifdef TEKO_HAVE_EPETRA
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Thyra_EpetraLinearOp.hpp"
#endif

#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_LU2x2InverseOp.hpp"
#include "Teko_ReorderedLinearOp.hpp"

// Thyra includes
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

// Test-rig

typedef Teko::ST ST;
typedef Teko::LO LO;
typedef Teko::GO GO;
typedef Teko::NT NT;

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::tpetraLinearOp;

#ifdef TEKO_HAVE_EPETRA
using Thyra::epetraLinearOp;
const RCP<const Thyra::LinearOpBase<double> > build2x2(const Epetra_Comm& comm, double a, double b,
                                                       double c, double d) {
  RCP<Epetra_Map> map = rcp(new Epetra_Map(2, 0, comm));

  int indicies[2];
  double row0[2];
  double row1[2];

  indicies[0] = 0;
  indicies[1] = 1;

  // build a CrsMatrix
  RCP<Epetra_CrsMatrix> blk = rcp(new Epetra_CrsMatrix(Copy, *map, 2));
  row0[0]                   = a;
  row0[1]                   = b;  // do a transpose here!
  row1[0]                   = c;
  row1[1]                   = d;
  blk->InsertGlobalValues(0, 2, &row0[0], &indicies[0]);
  blk->InsertGlobalValues(1, 2, &row1[0], &indicies[0]);
  blk->FillComplete();

  return Thyra::epetraLinearOp(blk);
}
#endif

const RCP<const Thyra::LinearOpBase<ST> > build2x2(
    const Teuchos::RCP<const Teuchos::Comm<int> > comm, ST a, ST b, ST c, ST d) {
  RCP<const Tpetra::Map<LO, GO, NT> > map = rcp(new const Tpetra::Map<LO, GO, NT>(2, 0, comm));

  GO indices[2];
  ST row0[2];
  ST row1[2];

  indices[0] = 0;
  indices[1] = 1;

  // build a CrsMatrix
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > blk = Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);
  row0[0]                                     = a;
  row0[1]                                     = b;  // do a transpose here!
  row1[0]                                     = c;
  row1[1]                                     = d;
  blk->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices, 2), Teuchos::ArrayView<ST>(row0, 2));
  blk->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices, 2), Teuchos::ArrayView<ST>(row1, 2));
  blk->fillComplete();

  return Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(blk->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(blk->getDomainMap()), blk);
}

#ifdef TEKO_HAVE_EPETRA
TEUCHOS_UNIT_TEST(tLU2x2InverseOp, exact_test) {
// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Teko::LinearOp A_00     = build2x2(Comm, 1, 2, 3, 4);
  Teko::LinearOp A_01     = build2x2(Comm, 5, 6, 7, 8);
  Teko::LinearOp A_02     = build2x2(Comm, 9, 10, 11, 12);
  Teko::LinearOp A_10     = build2x2(Comm, 9, 10, 11, 12);
  Teko::LinearOp A_11     = build2x2(Comm, -13, -14, -15, -16);
  Teko::LinearOp A_12     = build2x2(Comm, -1, -4, -5, -6);
  Teko::LinearOp A_20     = build2x2(Comm, -9, -10, -11, -12);
  Teko::LinearOp A_21     = build2x2(Comm, 13, 14, 15, 16);
  Teko::LinearOp A_22     = build2x2(Comm, 1, 4, 5, 6);
  Teko::BlockedLinearOp A = Teko::createBlockedOp();
  Teko::beginBlockFill(A, 3, 3);
  Teko::setBlock(0, 0, A, A_00);
  Teko::setBlock(0, 1, A, A_01);
  Teko::setBlock(0, 2, A, A_02);

  Teko::setBlock(1, 0, A, A_10);
  Teko::setBlock(1, 1, A, A_11);
  Teko::setBlock(1, 2, A, A_12);

  Teko::setBlock(2, 0, A, A_20);
  Teko::setBlock(2, 1, A, A_21);
  Teko::setBlock(2, 2, A, A_22);
  Teko::endBlockFill(A);

  std::string reorderType                           = "[ [0 1] 2]";
  Teuchos::RCP<const Teko::BlockReorderManager> brm = Teko::blockedReorderFromString(reorderType);
  Teko::ModifiableLinearOp re_A =
      Teuchos::rcp_const_cast<Thyra::LinearOpBase<double> >(Teko::buildReorderedLinearOp(*brm, A));
  Teko::ModifiableLinearOp final_A = Teuchos::rcp(new Teko::ReorderedLinearOp(brm, re_A));

  Thyra::LinearOpTester<double> tester;
  tester.dump_all(true);
  tester.show_all_tests(true);

  {
    const bool result = tester.compare(*final_A, *A, Teuchos::ptrFromRef(out));
    if (!result) {
      out << "Apply: FAILURE" << std::endl;
      success = false;
    } else
      out << "Apply: SUCCESS" << std::endl;
  }
}
#endif

TEUCHOS_UNIT_TEST(tLU2x2InverseOp, exact_test_tpetra) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  Teko::LinearOp A_00     = build2x2(Comm, 1, 2, 3, 4);
  Teko::LinearOp A_01     = build2x2(Comm, 5, 6, 7, 8);
  Teko::LinearOp A_02     = build2x2(Comm, 9, 10, 11, 12);
  Teko::LinearOp A_10     = build2x2(Comm, 9, 10, 11, 12);
  Teko::LinearOp A_11     = build2x2(Comm, -13, -14, -15, -16);
  Teko::LinearOp A_12     = build2x2(Comm, -1, -4, -5, -6);
  Teko::LinearOp A_20     = build2x2(Comm, -9, -10, -11, -12);
  Teko::LinearOp A_21     = build2x2(Comm, 13, 14, 15, 16);
  Teko::LinearOp A_22     = build2x2(Comm, 1, 4, 5, 6);
  Teko::BlockedLinearOp A = Teko::createBlockedOp();
  Teko::beginBlockFill(A, 3, 3);
  Teko::setBlock(0, 0, A, A_00);
  Teko::setBlock(0, 1, A, A_01);
  Teko::setBlock(0, 2, A, A_02);

  Teko::setBlock(1, 0, A, A_10);
  Teko::setBlock(1, 1, A, A_11);
  Teko::setBlock(1, 2, A, A_12);

  Teko::setBlock(2, 0, A, A_20);
  Teko::setBlock(2, 1, A, A_21);
  Teko::setBlock(2, 2, A, A_22);
  Teko::endBlockFill(A);

  std::string reorderType                           = "[ [0 1] 2]";
  Teuchos::RCP<const Teko::BlockReorderManager> brm = Teko::blockedReorderFromString(reorderType);
  Teko::ModifiableLinearOp re_A =
      Teuchos::rcp_const_cast<Thyra::LinearOpBase<double> >(Teko::buildReorderedLinearOp(*brm, A));
  Teko::ModifiableLinearOp final_A = Teuchos::rcp(new Teko::ReorderedLinearOp(brm, re_A));

  Thyra::LinearOpTester<double> tester;
  tester.dump_all(true);
  tester.show_all_tests(true);

  {
    const bool result = tester.compare(*final_A, *A, Teuchos::ptrFromRef(out));
    if (!result) {
      out << "Apply: FAILURE" << std::endl;
      success = false;
    } else
      out << "Apply: SUCCESS" << std::endl;
  }
}
