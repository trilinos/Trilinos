// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_Config.h"

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#ifdef TEKO_HAVE_EPETRA

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Thyra_EpetraLinearOp.hpp"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

#endif  // TEKO_HAVE_EPETRA

#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

// Teko-Package includes
#include "Teko_ConfigDefs.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_DiagnosticLinearOp.hpp"
#include "Teko_DiagnosticPreconditionerFactory.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_PreconditionerLinearOp.hpp"

#include "Thyra_TpetraLinearOp.hpp"

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
#endif  // TEKO_HAVE_EPETRA

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
TEUCHOS_UNIT_TEST(tDiagnosticLinearOp, application_test) {
// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  RCP<Teko::InverseLibrary> invLibrary = Teko::InverseLibrary::buildFromStratimikos();
  RCP<Teko::InverseFactory> invFact    = invLibrary->getInverseFactory("Amesos");
  Teko::LinearOp A                     = buildSystem(Comm, 10000);
  Teko::ModifiableLinearOp invA        = Teko::buildInverse(*invFact, A);

  Teuchos::RCP<std::ostream> rcp_out = Teuchos::rcpFromRef(out);
  Teuchos::RCP<Teko::DiagnosticLinearOp> diag_A =
      rcp(new Teko::DiagnosticLinearOp(rcp_out, invA, "descriptive_label"));
  Teko::LinearOp diag_Alo = diag_A;

  Teko::MultiVector x = Thyra::createMember(A->domain());
  Teko::MultiVector y = Thyra::createMember(A->range());
  Thyra::randomize(-1.0, 1.0, x.ptr());

  Teuchos::Time timer("test-time");
  int count = 50;
  for (int i = 0; i < count; i++) {
    Teuchos::TimeMonitor monitor(timer, false);
    Teko::applyOp(diag_Alo, x, y);
  }
  // TEST_FLOATING_EQUALITY(timer.totalElapsedTime(),diag_A->totalTime(),0.05); // within 5% should
  // be good enough
  TEST_EQUALITY(count, diag_A->numApplications());
}
#endif  // TEKO_HAVE_EPETRA

TEUCHOS_UNIT_TEST(tDiagnosticLinearOp, application_test_tpetra) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  RCP<Teko::InverseLibrary> invLibrary = Teko::InverseLibrary::buildFromStratimikos();
  RCP<Teko::InverseFactory> invFact    = invLibrary->getInverseFactory("Ifpack2");
  Teko::LinearOp A                     = buildSystem(Comm, 10000);
  Teko::ModifiableLinearOp invA        = Teko::buildInverse(*invFact, A);

  Teuchos::RCP<std::ostream> rcp_out = Teuchos::rcpFromRef(out);
  Teuchos::RCP<Teko::DiagnosticLinearOp> diag_A =
      rcp(new Teko::DiagnosticLinearOp(rcp_out, invA, "descriptive_label"));
  Teko::LinearOp diag_Alo = diag_A;

  Teko::MultiVector x = Thyra::createMember(A->domain());
  Teko::MultiVector y = Thyra::createMember(A->range());
  Thyra::randomize(-1.0, 1.0, x.ptr());

  Teuchos::Time timer("test-time");
  int count = 50;
  for (int i = 0; i < count; i++) {
    Teuchos::TimeMonitor monitor(timer, false);
    Teko::applyOp(diag_Alo, x, y);
  }
  // TEST_FLOATING_EQUALITY(timer.totalElapsedTime(),diag_A->totalTime(),0.05); // within 5% should
  // be good enough
  TEST_EQUALITY(count, diag_A->numApplications());
}

#ifdef TEKO_HAVE_EPETRA
TEUCHOS_UNIT_TEST(tDiagnosticPreconditionerFactory, inverse_lib_test) {
// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // setup diagnostic inverse parameter list: uses Amesos under neady
  Teuchos::ParameterList pl;
  Teuchos::ParameterList& diagList = pl.sublist("Diagnostic");
  diagList.set<std::string>("Type", "Diagnostic Inverse");
  diagList.set<std::string>("Inverse Factory", "Amesos");
  diagList.set<std::string>("Descriptive Label", "the_descriptive_label");

  // build inverse factory
  RCP<Teko::InverseLibrary> invLibrary = Teko::InverseLibrary::buildFromParameterList(pl);
  RCP<Teko::InverseFactory> invFact    = invLibrary->getInverseFactory("Diagnostic");

  // build inverse operator
  Teko::LinearOp A    = buildSystem(Comm, 2000);
  Teko::LinearOp invA = Teko::buildInverse(*invFact, A);

  // apply inverse
  Teko::MultiVector x = Thyra::createMember(invA->domain());
  Teko::MultiVector y = Thyra::createMember(invA->range());
  Thyra::randomize(-1.0, 1.0, x.ptr());
  Teko::applyOp(invA, x, y);
}
#endif  // TEKO_HAVE_EPETRA

TEUCHOS_UNIT_TEST(tDiagnosticPreconditionerFactory, inverse_lib_test_tpetra) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  // setup diagnostic inverse parameter list: uses Amesos under neady
  Teuchos::ParameterList pl;
  Teuchos::ParameterList& diagList = pl.sublist("Diagnostic");
  diagList.set<std::string>("Type", "Diagnostic Inverse");
  diagList.set<std::string>("Inverse Factory", "Ifpack2");
  diagList.set<std::string>("Descriptive Label", "the_descriptive_label");

  // build inverse factory
  RCP<Teko::InverseLibrary> invLibrary = Teko::InverseLibrary::buildFromParameterList(pl);
  RCP<Teko::InverseFactory> invFact    = invLibrary->getInverseFactory("Diagnostic");

  // build inverse operator
  Teko::LinearOp A    = buildSystem(Comm, 2000);
  Teko::LinearOp invA = Teko::buildInverse(*invFact, A);

  // apply inverse
  Teko::MultiVector x = Thyra::createMember(invA->domain());
  Teko::MultiVector y = Thyra::createMember(invA->range());
  Thyra::randomize(-1.0, 1.0, x.ptr());
  Teko::applyOp(invA, x, y);
}

TEUCHOS_UNIT_TEST(tDiagnosticPreconditionerFactory, inverse_lib_test_tpetra_preconditioned_solve) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  // setup diagnostic inverse parameter list: uses Amesos under neady
  Teuchos::ParameterList pl;
  Teuchos::ParameterList& diagList = pl.sublist("Diagnostic");
  diagList.set<std::string>("Type", "Diagnostic Inverse");
  diagList.set<std::string>("Inverse Factory", "Belos");
  diagList.set<std::string>("Preconditioner Factory", "Ifpack2");
  diagList.set<std::string>("Descriptive Label", "the_descriptive_label");

  // build inverse factory
  RCP<Teko::InverseLibrary> invLibrary = Teko::InverseLibrary::buildFromParameterList(pl);
  RCP<Teko::InverseFactory> invFact    = invLibrary->getInverseFactory("Diagnostic");

  // build inverse operator
  Teko::LinearOp A    = buildSystem(Comm, 2000);
  Teko::LinearOp invA = Teko::buildInverse(*invFact, A);

  // apply inverse
  Teko::MultiVector x = Thyra::createMember(invA->domain());
  Teko::MultiVector y = Thyra::createMember(invA->range());
  Thyra::randomize(-1.0, 1.0, x.ptr());
  Teko::applyOp(invA, x, y);
}

#ifdef TEKO_HAVE_EPETRA
TEUCHOS_UNIT_TEST(tDiagnosticPreconditionerFactory, construction_test) {
// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // build operator and solver
  Teko::LinearOp A                     = buildSystem(Comm, 2000);
  RCP<Teko::InverseLibrary> invLibrary = Teko::InverseLibrary::buildFromStratimikos();
  RCP<Teko::InverseFactory> direct     = invLibrary->getInverseFactory("Ifpack");

  // build diagnostic preconditioner
  Teko::DiagnosticPreconditionerFactory dpf(direct, "Direct Diagnostic", Teuchos::rcpFromRef(out));
  RCP<Teko::InverseFactory> invFact =
      Teuchos::rcp(new Teko::PreconditionerInverseFactory(Teuchos::rcpFromRef(dpf), Teuchos::null));

  // test rebuild functionality
  int count = 10;
  Teuchos::Time buildTime("build-time");
  Teko::ModifiableLinearOp invA;
  for (int i = 0; i < count; i++) {
    // do a timed build of linear operators
    {
      Teuchos::TimeMonitor monitor(buildTime, false);
      invA = Teko::buildInverse(*invFact, A);
    }

    RCP<const Teko::PreconditionerLinearOp<double> > precOp =
        rcp_dynamic_cast<const Teko::PreconditionerLinearOp<double> >(invA, true);
    RCP<const Teko::DiagnosticLinearOp> diagOp =
        rcp_dynamic_cast<const Teko::DiagnosticLinearOp>(precOp->getOperator(), true);
  }
  TEST_EQUALITY(dpf.numInitialBuilds(), buildTime.numCalls());
  // TEST_FLOATING_EQUALITY(dpf.totalInitialBuildTime(),
  //                        buildTime.totalElapsedTime(),0.05);  // within 5% should be good enough

  // test rebuild functionality
  Teuchos::Time rebuildTime("rebuild-time");
  for (int i = 0; i < count; i++) {
    // do a timed build of linear operators
    {
      Teuchos::TimeMonitor monitor(rebuildTime, false);
      Teko::rebuildInverse(*invFact, A, invA);
    }
  }
  TEST_EQUALITY(dpf.numRebuilds(), rebuildTime.numCalls());
  // TEST_FLOATING_EQUALITY(dpf.totalRebuildTime(),
  //                        rebuildTime.totalElapsedTime(),0.05);  // within 5% should be good
  //                        enough
}
#endif  // TEKO_HAVE_EPETRA

TEUCHOS_UNIT_TEST(tDiagnosticPreconditionerFactory, construction_test_tpetra) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  // build operator and solver
  Teko::LinearOp A                     = buildSystem(Comm, 2000);
  RCP<Teko::InverseLibrary> invLibrary = Teko::InverseLibrary::buildFromStratimikos();
  RCP<Teko::InverseFactory> direct     = invLibrary->getInverseFactory("Ifpack2");

  // build diagnostic preconditioner
  Teko::DiagnosticPreconditionerFactory dpf(direct, "Direct Diagnostic", Teuchos::rcpFromRef(out));
  RCP<Teko::InverseFactory> invFact =
      Teuchos::rcp(new Teko::PreconditionerInverseFactory(Teuchos::rcpFromRef(dpf), Teuchos::null));

  // test rebuild functionality
  int count = 10;
  Teuchos::Time buildTime("build-time");
  Teko::ModifiableLinearOp invA;
  for (int i = 0; i < count; i++) {
    // do a timed build of linear operators
    {
      Teuchos::TimeMonitor monitor(buildTime, false);
      invA = Teko::buildInverse(*invFact, A);
    }

    RCP<const Teko::PreconditionerLinearOp<double> > precOp =
        rcp_dynamic_cast<const Teko::PreconditionerLinearOp<double> >(invA, true);
    RCP<const Teko::DiagnosticLinearOp> diagOp =
        rcp_dynamic_cast<const Teko::DiagnosticLinearOp>(precOp->getOperator(), true);
  }
  TEST_EQUALITY(dpf.numInitialBuilds(), buildTime.numCalls());
  // TEST_FLOATING_EQUALITY(dpf.totalInitialBuildTime(),
  //                        buildTime.totalElapsedTime(),0.05);  // within 5% should be good enough

  // test rebuild functionality
  Teuchos::Time rebuildTime("rebuild-time");
  for (int i = 0; i < count; i++) {
    // do a timed build of linear operators
    {
      Teuchos::TimeMonitor monitor(rebuildTime, false);
      Teko::rebuildInverse(*invFact, A, invA);
    }
  }
  TEST_EQUALITY(dpf.numRebuilds(), rebuildTime.numCalls());
  // TEST_FLOATING_EQUALITY(dpf.totalRebuildTime(),
  //                        rebuildTime.totalElapsedTime(),0.05);  // within 5% should be good
  //                        enough
}

#ifdef TEKO_HAVE_EPETRA
TEUCHOS_UNIT_TEST(tDiagnosticLinearOp, residual_test) {
// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  RCP<Teko::InverseLibrary> invLibrary = Teko::InverseLibrary::buildFromStratimikos();
  RCP<Teko::InverseFactory> invFact    = invLibrary->getInverseFactory("Amesos");
  Teko::LinearOp A                     = buildSystem(Comm, 10000);
  Teko::ModifiableLinearOp invA        = Teko::buildInverse(*invFact, A);

  Teuchos::RCP<std::ostream> rcp_out = Teuchos::rcpFromRef(out);
  Teuchos::RCP<Teko::DiagnosticLinearOp> diag_invA =
      rcp(new Teko::DiagnosticLinearOp(rcp_out, A, invA, "descriptive_label"));
  Teko::LinearOp diag_invAlo = diag_invA;

  // simple default value
  {
    Teko::MultiVector x = Thyra::createMember(invA->domain());
    Teko::MultiVector y = Thyra::createMember(invA->range());
    Thyra::randomize(-1.0, 1.0, x.ptr());

    Teko::MultiVector residual = Teko::deepcopy(x);

    Teko::applyOp(diag_invAlo, x, y);
    Teko::applyOp(A, y, residual, -1.0, 1.0);

    double myresid = Teko::norm_2(residual, 0);

    TEST_FLOATING_EQUALITY(myresid, diag_invA->getResidualNorm(), 1e-14);
  }

  // arbitrary alpha and beta
  {
    double alpha        = 3.141;
    double beta         = 1.618;
    Teko::MultiVector x = Thyra::createMember(invA->domain());
    Teko::MultiVector y = Thyra::createMember(invA->range());
    Thyra::randomize(-1.0, 1.0, x.ptr());
    Thyra::randomize(-1.0, 1.0, y.ptr());

    Teko::MultiVector residual = Teko::deepcopy(x);
    Teko::MultiVector z        = Teko::deepcopy(y);

    Teko::applyOp(diag_invAlo, x, z, alpha, beta);

    Teko::applyOp(A, z, residual, -1.0, alpha);
    // alpha x - A z

    Teko::applyOp(A, y, residual, beta, 1.0);
    // alpha (x - A z) - beta A y

    double myresid = Teko::norm_2(residual, 0);

    TEST_FLOATING_EQUALITY(myresid, diag_invA->getResidualNorm(), 1e-14);
  }
}
#endif  // TEKO_HAVE_EPETRA

TEUCHOS_UNIT_TEST(tDiagnosticLinearOp, residual_test_tpetra) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  RCP<Teko::InverseLibrary> invLibrary = Teko::InverseLibrary::buildFromStratimikos();
  RCP<Teko::InverseFactory> invFact    = invLibrary->getInverseFactory("Ifpack2");
  Teko::LinearOp A                     = buildSystem(Comm, 10000);
  Teko::ModifiableLinearOp invA        = Teko::buildInverse(*invFact, A);

  Teuchos::RCP<std::ostream> rcp_out = Teuchos::rcpFromRef(out);
  Teuchos::RCP<Teko::DiagnosticLinearOp> diag_invA =
      rcp(new Teko::DiagnosticLinearOp(rcp_out, A, invA, "descriptive_label"));
  Teko::LinearOp diag_invAlo = diag_invA;

  // simple default value
  {
    Teko::MultiVector x = Thyra::createMember(invA->domain());
    Teko::MultiVector y = Thyra::createMember(invA->range());
    Thyra::randomize(-1.0, 1.0, x.ptr());

    Teko::MultiVector residual = Teko::deepcopy(x);

    Teko::applyOp(diag_invAlo, x, y);
    Teko::applyOp(A, y, residual, -1.0, 1.0);

    double myresid = Teko::norm_2(residual, 0);

    // residual is O(1e-10), so check using rel tolerance of 1e-6
    TEST_FLOATING_EQUALITY(myresid, diag_invA->getResidualNorm(), 1e-6);
  }

  // arbitrary alpha and beta
  {
    double alpha        = 3.141;
    double beta         = 1.618;
    Teko::MultiVector x = Thyra::createMember(invA->domain());
    Teko::MultiVector y = Thyra::createMember(invA->range());
    Thyra::randomize(-1.0, 1.0, x.ptr());
    Thyra::randomize(-1.0, 1.0, y.ptr());

    Teko::MultiVector residual = Teko::deepcopy(x);
    Teko::MultiVector z        = Teko::deepcopy(y);

    Teko::applyOp(diag_invAlo, x, z, alpha, beta);

    Teko::applyOp(A, z, residual, -1.0, alpha);
    // alpha x - A z

    Teko::applyOp(A, y, residual, beta, 1.0);
    // alpha (x - A z) - beta A y

    double myresid = Teko::norm_2(residual, 0);

    // residual is O(1e-10), so check using rel tolerance of 1e-6
    TEST_FLOATING_EQUALITY(myresid, diag_invA->getResidualNorm(), 1e-6);
  }
}
