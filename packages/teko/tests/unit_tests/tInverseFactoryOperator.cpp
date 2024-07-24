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

#include "Teko_ConfigDefs.hpp"

#ifdef TEKO_HAVE_EPETRA
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Teko_EpetraOperatorWrapper.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Teko_InverseFactoryOperator.hpp"
#endif

#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_TpetraInverseFactoryOperator.hpp"
#include "Teko_JacobiPreconditionerFactory.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"

#include "Teko_TpetraOperatorWrapper.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

#define SS_ECHO(ops)      \
  {                       \
    std::stringstream ss; \
    ss << ops;            \
    ECHO(ss.str())        \
  };

typedef Teko::ST ST;
typedef Teko::LO LO;
typedef Teko::GO GO;
typedef Teko::NT NT;

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

///////////////////////////////////////////////////////////

#ifdef TEKO_HAVE_EPETRA
const RCP<Epetra_Operator> buildSystem(const Epetra_Comm& comm, int size) {
  Epetra_Map map(size, 0, comm);

  RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy, map, 0));

  double values[] = {-1.0, 2.0, -1.0};
  int iTemp[]     = {-1, 0, 1}, indices[3];
  double* vPtr;
  int* iPtr;

  for (int i = 0; i < map.NumMyElements(); i++) {
    int count = 3;
    int gid   = map.GID(i);

    vPtr = values;
    iPtr = indices;

    indices[0] = gid + iTemp[0];
    indices[1] = gid + iTemp[1];
    indices[2] = gid + iTemp[2];

    if (gid == 0) {
      vPtr  = &values[1];
      iPtr  = &indices[1];
      count = 2;
    } else if (gid == map.MaxAllGID())
      count = 2;

    mat->InsertGlobalValues(gid, count, vPtr, iPtr);
  }

  mat->FillComplete();

  // return Thyra::nonconstEpetraLinearOp(mat);
  return mat;
}
#endif

const RCP<Tpetra::Operator<ST, LO, GO, NT> > buildSystem(
    const Teuchos::RCP<const Teuchos::Comm<int> > comm, GO size) {
  RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(size, 0, comm));

  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > mat = Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 3);

  ST values[] = {-1.0, 2.0, -1.0};
  GO iTemp[]  = {-1, 0, 1}, indices[3];
  ST* vPtr;
  GO* iPtr;
  for (size_t i = 0; i < map->getLocalNumElements(); i++) {
    int count = 3;
    GO gid    = map->getGlobalElement(i);

    vPtr = values;
    iPtr = indices;

    indices[0] = gid + iTemp[0];
    indices[1] = gid + iTemp[1];
    indices[2] = gid + iTemp[2];

    if (gid == 0) {
      vPtr  = &values[1];
      iPtr  = &indices[1];
      count = 2;
    } else if (gid == map->getMaxAllGlobalIndex())
      count = 2;

    mat->insertGlobalValues(gid, Teuchos::ArrayView<GO>(iPtr, count),
                            Teuchos::ArrayView<ST>(vPtr, count));
  }

  mat->fillComplete();

  return mat;
  // return
  // Thyra::tpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(mat->getRangeMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(mat->getDomainMap()),mat);
}

#ifdef TEKO_HAVE_EPETRA
TEUCHOS_UNIT_TEST(tInverseFactoryOperator, test_Direct_Solve) {
// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  Teuchos::RCP<Teko::InverseLibrary> invLib     = Teko::InverseLibrary::buildFromStratimikos();
  Teuchos::RCP<Teko::InverseFactory> invFactory = invLib->getInverseFactory("Amesos");

  Teuchos::RCP<Epetra_Operator> eA = buildSystem(comm, 50);
  Teko::LinearOp A                 = Thyra::epetraLinearOp(eA);
  Teko::ModifiableLinearOp invA    = Teko::buildInverse(*invFactory, A);

  Teko::Epetra::InverseFactoryOperator invFactOp(invFactory);
  invFactOp.buildInverseOperator(eA);

  {
    // because InverseFactoryOperator is a "Preconditioner" then need to
    // call Epetra_Operator::ApplyInverse
    Teko::LinearOp testInvA = Thyra::epetraLinearOp(Teuchos::rcpFromRef(invFactOp), Thyra::NOTRANS,
                                                    Thyra::EPETRA_OP_APPLY_APPLY_INVERSE);

    Thyra::LinearOpTester<double> tester;
    tester.show_all_tests(true);
    tester.set_all_error_tol(1e-14);

    const bool result = tester.compare(*invA, *testInvA, Teuchos::ptrFromRef(out));
    if (!result) {
      out << "Apply 0: FAILURE" << std::endl;
      success = false;
    } else
      out << "Apply 0: SUCCESS" << std::endl;
  }

  invFactOp.rebuildInverseOperator(eA);
  {
    // because InverseFactoryOperator is a "Preconditioner" then need to
    // call Epetra_Operator::ApplyInverse
    Teko::LinearOp testInvA = Thyra::epetraLinearOp(Teuchos::rcpFromRef(invFactOp), Thyra::NOTRANS,
                                                    Thyra::EPETRA_OP_APPLY_APPLY_INVERSE);

    Thyra::LinearOpTester<double> tester;
    tester.show_all_tests(true);
    tester.set_all_error_tol(1e-14);

    const bool result = tester.compare(*invA, *testInvA, Teuchos::ptrFromRef(out));
    if (!result) {
      out << "Apply 0: FAILURE" << std::endl;
      success = false;
    } else
      out << "Apply 0: SUCCESS" << std::endl;
  }
}
#endif

TEUCHOS_UNIT_TEST(tInverseFactoryOperator, test_Direct_Solve_tpetra) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  Teuchos::RCP<Teko::InverseLibrary> invLib     = Teko::InverseLibrary::buildFromStratimikos();
  Teuchos::RCP<Teko::InverseFactory> invFactory = invLib->getInverseFactory("Ifpack2");

  Teuchos::RCP<Tpetra::Operator<ST, LO, GO, NT> > eA = buildSystem(comm, 50);
  Teko::LinearOp A                                   = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(eA->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(eA->getDomainMap()), eA);
  Teko::ModifiableLinearOp invA = Teko::buildInverse(*invFactory, A);

  Teko::TpetraHelpers::InverseFactoryOperator invFactOp(invFactory);
  invFactOp.buildInverseOperator(eA);

  {
    // because InverseFactoryOperator is a "Preconditioner" then need to
    // call Tpetra_Operator::ApplyInverse
    Teko::LinearOp testInvA = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(invFactOp.getRangeMap()),
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(invFactOp.getDomainMap()),
        Teuchos::rcpFromRef(invFactOp));

    Thyra::LinearOpTester<ST> tester;
    tester.show_all_tests(true);
    tester.set_all_error_tol(1e-14);

    const bool result = tester.compare(*invA, *testInvA, Teuchos::ptrFromRef(out));
    if (!result) {
      out << "Apply 0: FAILURE" << std::endl;
      success = false;
    } else
      out << "Apply 0: SUCCESS" << std::endl;
  }

  invFactOp.rebuildInverseOperator(eA);
  {
    // because InverseFactoryOperator is a "Preconditioner" then need to
    // call Tpetra_Operator::ApplyInverse
    Teko::LinearOp testInvA = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(invFactOp.getRangeMap()),
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(invFactOp.getDomainMap()),
        Teuchos::rcpFromRef(invFactOp));

    Thyra::LinearOpTester<ST> tester;
    tester.show_all_tests(true);
    tester.set_all_error_tol(1e-14);

    const bool result = tester.compare(*invA, *testInvA, Teuchos::ptrFromRef(out));
    if (!result) {
      out << "Apply 0: FAILURE" << std::endl;
      success = false;
    } else
      out << "Apply 0: SUCCESS" << std::endl;
  }
}

#ifdef TEKO_HAVE_EPETRA
TEUCHOS_UNIT_TEST(tInverseFactoryOperator, test_Block_Solve) {
// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  Teuchos::RCP<Teko::InverseLibrary> invLib        = Teko::InverseLibrary::buildFromStratimikos();
  Teuchos::RCP<Teko::InverseFactory> amesosFactory = invLib->getInverseFactory("Amesos");

  Teuchos::RCP<Epetra_Operator> eA00 = buildSystem(comm, 50);
  Teko::LinearOp A_00                = Thyra::epetraLinearOp(eA00);
  Teko::ModifiableLinearOp invA_00   = Teko::buildInverse(*amesosFactory, A_00);

  Teko::LinearOp A    = Thyra::block2x2<double>(A_00, Teuchos::null, Teuchos::null, A_00);
  Teko::LinearOp invA = Thyra::block2x2<double>(invA_00, Teuchos::null, Teuchos::null, invA_00);
  Teuchos::RCP<Epetra_Operator> eInvA = Teuchos::rcp(new Teko::Epetra::EpetraOperatorWrapper(invA));
  Teko::LinearOp cmpInvA              = Thyra::epetraLinearOp(eInvA);

  Teuchos::RCP<Teko::PreconditionerFactory> jacFact =
      Teuchos::rcp(new Teko::JacobiPreconditionerFactory(invA_00, invA_00));
  Teuchos::RCP<Teko::InverseFactory> invFactory =
      Teuchos::rcp(new Teko::PreconditionerInverseFactory(jacFact, Teuchos::null));
  Teuchos::RCP<Epetra_Operator> eA = Teuchos::rcp(new Teko::Epetra::EpetraOperatorWrapper(A));

  Teko::Epetra::InverseFactoryOperator invFactOp(invFactory);
  invFactOp.buildInverseOperator(eA);

  {
    // because InverseFactoryOperator is a "Preconditioner" then need to
    // call Epetra_Operator::ApplyInverse
    Teko::LinearOp testInvA = Thyra::epetraLinearOp(Teuchos::rcpFromRef(invFactOp), Thyra::NOTRANS,
                                                    Thyra::EPETRA_OP_APPLY_APPLY_INVERSE);

    Thyra::LinearOpTester<double> tester;
    tester.show_all_tests(true);
    tester.set_all_error_tol(1e-14);

    const bool result = tester.compare(*cmpInvA, *testInvA, Teuchos::ptrFromRef(out));
    if (!result) {
      out << "Apply 0: FAILURE" << std::endl;
      success = false;
    } else
      out << "Apply 0: SUCCESS" << std::endl;
  }

  invFactOp.rebuildInverseOperator(eA);
  {
    // because InverseFactoryOperator is a "Preconditioner" then need to
    // call Epetra_Operator::ApplyInverse
    Teko::LinearOp testInvA = Thyra::epetraLinearOp(Teuchos::rcpFromRef(invFactOp), Thyra::NOTRANS,
                                                    Thyra::EPETRA_OP_APPLY_APPLY_INVERSE);

    Thyra::LinearOpTester<double> tester;
    tester.show_all_tests(true);
    tester.set_all_error_tol(1e-14);

    const bool result = tester.compare(*cmpInvA, *testInvA, Teuchos::ptrFromRef(out));
    if (!result) {
      out << "Apply 0: FAILURE" << std::endl;
      success = false;
    } else
      out << "Apply 0: SUCCESS" << std::endl;
  }
}
#endif

TEUCHOS_UNIT_TEST(tInverseFactoryOperator, test_Block_Solve_tpetra) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  Teuchos::RCP<Teko::InverseLibrary> invLib        = Teko::InverseLibrary::buildFromStratimikos();
  Teuchos::RCP<Teko::InverseFactory> amesosFactory = invLib->getInverseFactory("Ifpack2");

  Teuchos::RCP<Tpetra::Operator<ST, LO, GO, NT> > eA00 = buildSystem(comm, 50);
  Teko::LinearOp A_00                                  = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(eA00->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(eA00->getDomainMap()), eA00);
  Teko::ModifiableLinearOp invA_00 = Teko::buildInverse(*amesosFactory, A_00);

  Teko::LinearOp A    = Thyra::block2x2<ST>(A_00, Teuchos::null, Teuchos::null, A_00);
  Teko::LinearOp invA = Thyra::block2x2<ST>(invA_00, Teuchos::null, Teuchos::null, invA_00);
  Teuchos::RCP<Tpetra::Operator<ST, LO, GO, NT> > eInvA =
      Teuchos::rcp(new Teko::TpetraHelpers::TpetraOperatorWrapper(invA));
  Teko::LinearOp cmpInvA = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(eInvA->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(eInvA->getDomainMap()), eInvA);

  Teuchos::RCP<Teko::PreconditionerFactory> jacFact =
      Teuchos::rcp(new Teko::JacobiPreconditionerFactory(invA_00, invA_00));
  Teuchos::RCP<Teko::InverseFactory> invFactory =
      Teuchos::rcp(new Teko::PreconditionerInverseFactory(jacFact, Teuchos::null));
  Teuchos::RCP<Tpetra::Operator<ST, LO, GO, NT> > eA =
      Teuchos::rcp(new Teko::TpetraHelpers::TpetraOperatorWrapper(A));

  Teko::TpetraHelpers::InverseFactoryOperator invFactOp(invFactory);
  invFactOp.buildInverseOperator(eA);

  {
    // because InverseFactoryOperator is a "Preconditioner" then need to
    // call Tpetra_Operator::ApplyInverse
    Teko::LinearOp testInvA = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(invFactOp.getRangeMap()),
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(invFactOp.getDomainMap()),
        Teuchos::rcpFromRef(invFactOp));

    Thyra::LinearOpTester<ST> tester;
    tester.show_all_tests(true);
    tester.set_all_error_tol(1e-14);

    const bool result = tester.compare(*cmpInvA, *testInvA, Teuchos::ptrFromRef(out));
    if (!result) {
      out << "Apply 0: FAILURE" << std::endl;
      success = false;
    } else
      out << "Apply 0: SUCCESS" << std::endl;
  }

  invFactOp.rebuildInverseOperator(eA);
  {
    // because InverseFactoryOperator is a "Preconditioner" then need to
    // call Tpetra_Operator::ApplyInverse
    Teko::LinearOp testInvA = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(invFactOp.getRangeMap()),
        Thyra::tpetraVectorSpace<ST, LO, GO, NT>(invFactOp.getDomainMap()),
        Teuchos::rcpFromRef(invFactOp));

    Thyra::LinearOpTester<ST> tester;
    tester.show_all_tests(true);
    tester.set_all_error_tol(1e-14);

    const bool result = tester.compare(*cmpInvA, *testInvA, Teuchos::ptrFromRef(out));
    if (!result) {
      out << "Apply 0: FAILURE" << std::endl;
      success = false;
    } else
      out << "Apply 0: SUCCESS" << std::endl;
  }
}
