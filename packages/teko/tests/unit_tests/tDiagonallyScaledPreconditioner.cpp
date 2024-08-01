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
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#include "Teko_Config.h"

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
#include "Teko_ConfigDefs.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_DiagonallyScaledPreconditionerFactory.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_PreconditionerLinearOp.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

// Test-rig

typedef Teko::ST ST;
typedef Teko::LO LO;
typedef Teko::GO GO;
typedef Teko::NT NT;

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

#ifdef TEKO_HAVE_EPETRA
using Thyra::epetraLinearOp;
const RCP<Thyra::LinearOpBase<double> > buildSystem(const Epetra_Comm& comm, int size) {
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

  return Thyra::nonconstEpetraLinearOp(mat);
}
#endif

const RCP<Thyra::LinearOpBase<ST> > buildSystem(const Teuchos::RCP<const Teuchos::Comm<int> > comm,
                                                GO size) {
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

  return Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat->getDomainMap()), mat);
}

#ifdef TEKO_HAVE_EPETRA
TEUCHOS_UNIT_TEST(tDiagonallyScaledPreconditioner, invfactory_test) {
// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::ParameterList pl;
  Teuchos::ParameterList& diagList = pl.sublist("DiagScal");
  diagList.set<std::string>("Type", "Diagonal Scaling");
  diagList.set<std::string>("Inverse Factory", "Amesos");

  RCP<Teko::InverseLibrary> invLib  = Teko::InverseLibrary::buildFromParameterList(pl);
  RCP<Teko::InverseFactory> invFact = invLib->getInverseFactory("DiagScal");
  RCP<Teko::InverseFactory> dirFact = invLib->getInverseFactory("Amesos");

  RCP<Thyra::LinearOpBase<double> > A = buildSystem(Comm, 50);

  Teko::LinearOp invA      = Teko::buildInverse(*invFact, A);
  Teko::LinearOp invExactA = Teko::buildInverse(*dirFact, A);

  Thyra::LinearOpTester<double> tester;
  tester.show_all_tests(true);
  tester.set_all_error_tol(8e-14);

  const bool result = tester.compare(*invA, *invExactA, Teuchos::ptrFromRef(out));
  if (!result) {
    out << "Apply 0: FAILURE" << std::endl;
    success = false;
  } else
    out << "Apply 0: SUCCESS" << std::endl;
}
#endif

TEUCHOS_UNIT_TEST(tDiagonallyScaledPreconditioner, invfactory_test_tpetra) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  Teuchos::ParameterList pl;
  Teuchos::ParameterList& diagList = pl.sublist("DiagScal");
  diagList.set<std::string>("Type", "Diagonal Scaling");
  diagList.set<std::string>("Inverse Factory", "Belos");

  RCP<Teko::InverseLibrary> invLib  = Teko::InverseLibrary::buildFromParameterList(pl);
  RCP<Teko::InverseFactory> invFact = invLib->getInverseFactory("DiagScal");
  RCP<Teko::InverseFactory> dirFact = invLib->getInverseFactory("Ifpack2");

  RCP<Thyra::LinearOpBase<ST> > A = buildSystem(Comm, 50);

  Teko::LinearOp invA      = Teko::buildInverse(*invFact, A);
  Teko::LinearOp invExactA = Teko::buildInverse(*dirFact, A);

  Thyra::LinearOpTester<ST> tester;
  tester.show_all_tests(true);
  tester.set_all_error_tol(8e-14);

  const bool result = tester.compare(*invA, *invExactA, Teuchos::ptrFromRef(out));
  if (!result) {
    out << "Apply 0: FAILURE" << std::endl;
    success = false;
  } else
    out << "Apply 0: SUCCESS" << std::endl;
}

#ifdef TEKO_HAVE_EPETRA
TEUCHOS_UNIT_TEST(tDiagonallyScaledPreconditioner, application_test_row) {
// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // build linear op tester
  bool result;
  Thyra::LinearOpTester<double> tester;
  tester.show_all_tests(true);
  tester.set_all_error_tol(8e-14);

  // build operators and factories
  RCP<Teko::InverseLibrary> invLib      = Teko::InverseLibrary::buildFromStratimikos();
  RCP<Teko::InverseFactory> subsolve    = invLib->getInverseFactory("Amesos");
  RCP<Teko::InverseFactory> subsolve_ml = invLib->getInverseFactory("ML");

  RCP<Thyra::LinearOpBase<double> > A = buildSystem(Comm, 50);

  typedef Teko::DiagonallyScaledPreconditionerFactory DSPF;

  RCP<Teko::PreconditionerFactory> dspf =
      rcp(new Teko::DiagonallyScaledPreconditionerFactory(subsolve, DSPF::ROW_SCALING));
  RCP<Teko::PreconditionerFactory> dspf_ml =
      rcp(new Teko::DiagonallyScaledPreconditionerFactory(subsolve_ml, DSPF::ROW_SCALING));
  RCP<Teko::InverseFactory> invFact =
      rcp(new Teko::PreconditionerInverseFactory(dspf, Teuchos::null));
  RCP<Teko::InverseFactory> invFact_ml =
      rcp(new Teko::PreconditionerInverseFactory(dspf_ml, Teuchos::null));

  // test build inverse capability
  Teko::ModifiableLinearOp invA      = Teko::buildInverse(*invFact, A);
  Teko::ModifiableLinearOp invA_ml   = Teko::buildInverse(*invFact_ml, A);
  Teko::ModifiableLinearOp invExactA = Teko::buildInverse(*subsolve, A);

  result = tester.compare(*invA, *invExactA, Teuchos::ptrFromRef(out));
  if (!result) {
    out << "Apply 0: FAILURE" << std::endl;
    success = false;
  } else
    out << "Apply 0: SUCCESS" << std::endl;

  // test Diagonally scaled rebuild capability
  Teko::rebuildInverse(*invFact, A, invA);
  Teko::rebuildInverse(*invFact_ml, A, invA_ml);  // here we tested repeatability using ML

  result = tester.compare(*invA, *invExactA, Teuchos::ptrFromRef(out));
  if (!result) {
    out << "Apply 0: FAILURE" << std::endl;
    success = false;
  } else
    out << "Apply 0: SUCCESS" << std::endl;
}
#endif

TEUCHOS_UNIT_TEST(tDiagonallyScaledPreconditioner, application_test_row_tpetra) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  // build linear op tester
  bool result;
  Thyra::LinearOpTester<ST> tester;
  tester.show_all_tests(true);
  tester.set_all_error_tol(5.8e-14);

  // build operators and factories
  RCP<Teko::InverseLibrary> invLib      = Teko::InverseLibrary::buildFromStratimikos();
  RCP<Teko::InverseFactory> subsolve    = invLib->getInverseFactory("Belos");
  RCP<Teko::InverseFactory> subsolve_ml = invLib->getInverseFactory("Ifpack2");

  RCP<Thyra::LinearOpBase<ST> > A = buildSystem(Comm, 50);

  typedef Teko::DiagonallyScaledPreconditionerFactory DSPF;

  RCP<Teko::PreconditionerFactory> dspf =
      rcp(new Teko::DiagonallyScaledPreconditionerFactory(subsolve, DSPF::ROW_SCALING));
  RCP<Teko::PreconditionerFactory> dspf_ml =
      rcp(new Teko::DiagonallyScaledPreconditionerFactory(subsolve_ml, DSPF::ROW_SCALING));
  RCP<Teko::InverseFactory> invFact =
      rcp(new Teko::PreconditionerInverseFactory(dspf, Teuchos::null));
  RCP<Teko::InverseFactory> invFact_ml =
      rcp(new Teko::PreconditionerInverseFactory(dspf_ml, Teuchos::null));

  // test build inverse capability
  Teko::ModifiableLinearOp invA      = Teko::buildInverse(*invFact, A);
  Teko::ModifiableLinearOp invA_ml   = Teko::buildInverse(*invFact_ml, A);
  Teko::ModifiableLinearOp invExactA = Teko::buildInverse(*subsolve, A);

  result = tester.compare(*invA, *invExactA, Teuchos::ptrFromRef(out));
  if (!result) {
    out << "Apply 0: FAILURE" << std::endl;
    success = false;
  } else
    out << "Apply 0: SUCCESS" << std::endl;

  // test Diagonally scaled rebuild capability
  Teko::rebuildInverse(*invFact, A, invA);
  Teko::rebuildInverse(*invFact_ml, A, invA_ml);  // here we tested repeatability using ML

  result = tester.compare(*invA, *invExactA, Teuchos::ptrFromRef(out));
  if (!result) {
    out << "Apply 0: FAILURE" << std::endl;
    success = false;
  } else
    out << "Apply 0: SUCCESS" << std::endl;
}

#ifdef TEKO_HAVE_EPETRA
TEUCHOS_UNIT_TEST(tDiagonallyScaledPreconditioner, application_test_column) {
// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // build linear op tester
  bool result;
  Thyra::LinearOpTester<double> tester;
  tester.show_all_tests(true);
  tester.set_all_error_tol(8e-14);

  // build operators and factories
  RCP<Teko::InverseLibrary> invLib      = Teko::InverseLibrary::buildFromStratimikos();
  RCP<Teko::InverseFactory> subsolve    = invLib->getInverseFactory("Amesos");
  RCP<Teko::InverseFactory> subsolve_ml = invLib->getInverseFactory("ML");

  RCP<Thyra::LinearOpBase<double> > A = buildSystem(Comm, 50);

  RCP<Teko::PreconditionerFactory> dspf =
      rcp(new Teko::DiagonallyScaledPreconditionerFactory(subsolve));
  RCP<Teko::PreconditionerFactory> dspf_ml =
      rcp(new Teko::DiagonallyScaledPreconditionerFactory(subsolve_ml));
  RCP<Teko::InverseFactory> invFact =
      rcp(new Teko::PreconditionerInverseFactory(dspf, Teuchos::null));
  RCP<Teko::InverseFactory> invFact_ml =
      rcp(new Teko::PreconditionerInverseFactory(dspf_ml, Teuchos::null));

  // test build inverse capability
  Teko::ModifiableLinearOp invA      = Teko::buildInverse(*invFact, A);
  Teko::ModifiableLinearOp invA_ml   = Teko::buildInverse(*invFact_ml, A);
  Teko::ModifiableLinearOp invExactA = Teko::buildInverse(*subsolve, A);

  result = tester.compare(*invA, *invExactA, Teuchos::ptrFromRef(out));
  if (!result) {
    out << "Apply 0: FAILURE" << std::endl;
    success = false;
  } else
    out << "Apply 0: SUCCESS" << std::endl;

  // test Diagonally scaled rebuild capability
  Teko::rebuildInverse(*invFact, A, invA);
  Teko::rebuildInverse(*invFact_ml, A, invA_ml);  // here we tested repeatability using ML

  result = tester.compare(*invA, *invExactA, Teuchos::ptrFromRef(out));
  if (!result) {
    out << "Apply 0: FAILURE" << std::endl;
    success = false;
  } else
    out << "Apply 0: SUCCESS" << std::endl;
}
#endif

TEUCHOS_UNIT_TEST(tDiagonallyScaledPreconditioner, application_test_column_tpetra) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  // build linear op tester
  bool result;
  Thyra::LinearOpTester<ST> tester;
  tester.show_all_tests(true);
  tester.set_all_error_tol(5e-13);

  // build operators and factories
  RCP<Teko::InverseLibrary> invLib      = Teko::InverseLibrary::buildFromStratimikos();
  RCP<Teko::InverseFactory> subsolve    = invLib->getInverseFactory("Belos");
  RCP<Teko::InverseFactory> subsolve_ml = invLib->getInverseFactory("Ifpack2");

  RCP<Thyra::LinearOpBase<ST> > A = buildSystem(Comm, 50);

  RCP<Teko::PreconditionerFactory> dspf =
      rcp(new Teko::DiagonallyScaledPreconditionerFactory(subsolve));
  RCP<Teko::PreconditionerFactory> dspf_ml =
      rcp(new Teko::DiagonallyScaledPreconditionerFactory(subsolve_ml));
  RCP<Teko::InverseFactory> invFact =
      rcp(new Teko::PreconditionerInverseFactory(dspf, Teuchos::null));
  RCP<Teko::InverseFactory> invFact_ml =
      rcp(new Teko::PreconditionerInverseFactory(dspf_ml, Teuchos::null));

  // test build inverse capability
  Teko::ModifiableLinearOp invA      = Teko::buildInverse(*invFact, A);
  Teko::ModifiableLinearOp invA_ml   = Teko::buildInverse(*invFact_ml, A);
  Teko::ModifiableLinearOp invExactA = Teko::buildInverse(*subsolve, A);

  result = tester.compare(*invA, *invExactA, Teuchos::ptrFromRef(out));
  if (!result) {
    out << "Apply 0: FAILURE" << std::endl;
    success = false;
  } else
    out << "Apply 0: SUCCESS" << std::endl;

  // test Diagonally scaled rebuild capability
  Teko::rebuildInverse(*invFact, A, invA);
  Teko::rebuildInverse(*invFact_ml, A, invA_ml);  // here we tested repeatability using ML

  result = tester.compare(*invA, *invExactA, Teuchos::ptrFromRef(out));
  if (!result) {
    out << "Apply 0: FAILURE" << std::endl;
    success = false;
  } else
    out << "Apply 0: SUCCESS" << std::endl;
}

#ifdef TEKO_HAVE_EPETRA
TEUCHOS_UNIT_TEST(tDiagonalOperator, replaceValues) {
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  RCP<Thyra::LinearOpBase<double> > A = buildSystem(Comm, 50);

  Teko::MultiVector diag = Teko::getDiagonal(A, Teko::AbsRowSum);
  Teko::replaceValue(diag, 0.0, 1.0);
  Teko::LinearOp invDiagOp = Teko::buildInvDiagonal(diag);
}
#endif

TEUCHOS_UNIT_TEST(tDiagonalOperator, replaceValues_tpetra) {
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  RCP<Thyra::LinearOpBase<ST> > A = buildSystem(Comm, 50);

  Teko::MultiVector diag = Teko::getDiagonal(A, Teko::AbsRowSum);
  Teko::replaceValue(diag, 0.0, 1.0);
  Teko::LinearOp invDiagOp = Teko::buildInvDiagonal(diag);
}
