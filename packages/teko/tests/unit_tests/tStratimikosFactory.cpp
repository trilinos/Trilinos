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

#include "Kokkos_Core.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "EpetraExt_RowMatrixOut.h"
#include "MatrixMarket_Tpetra.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_InverseFactoryOperator.hpp"
#include "Teko_JacobiPreconditionerFactory.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_EpetraOperatorWrapper.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_StratimikosFactory.hpp"
#include "Teko_StridedEpetraOperator.hpp"
#include "Teko_StridedTpetraOperator.hpp"

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

#include "Teuchos_AbstractFactoryStd.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

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

  return mat;
}

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

const RCP<Epetra_Operator> buildStridedSystem(const Epetra_Comm& comm, int size) {
  Epetra_Map map(2 * size, 0, comm);

  RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy, map, 0));

  int numUnks      = 2;
  double valuesA[] = {-1.0, 2.0, 7.0, -1.0};
  double valuesB[] = {-1.0, 2.0, 9.0, -1.0};
  int iTempA[]     = {-numUnks, 0, 1, numUnks}, indices[4];
  int iTempB[]     = {-numUnks, 0, -1, numUnks};
  double* vPtr;
  int* iPtr;

  for (int i = 0; i < map.NumMyElements() / numUnks; i++) {
    int count = 4;
    int gidA  = map.GID(2 * i);
    int gidB  = gidA + 1;

    for (int n = 0; n < numUnks; n++) {
      int* iTemp = (n == 0) ? iTempA : iTempB;
      int gid    = (n == 0) ? gidA : gidB;

      indices[0] = gid + iTemp[0];
      indices[1] = gid + iTemp[1];
      indices[2] = gid + iTemp[2];
      indices[3] = gid + iTemp[3];

      vPtr = (n == 0) ? valuesA : valuesB;
      iPtr = indices;

      if (gid < numUnks) {
        vPtr++;
        iPtr  = &indices[1];
        count = 3;
      } else if (gid > map.MaxAllGID() - numUnks)
        count = 3;

      mat->InsertGlobalValues(gid, count, vPtr, iPtr);
    }
  }

  mat->FillComplete();

  EpetraExt::RowMatrixToMatrixMarketFile("strided.mm", *mat);

  return mat;
}

const RCP<Tpetra::Operator<ST, LO, GO, NT> > buildStridedSystem(
    const Teuchos::RCP<Teuchos::Comm<int> > comm, GO size) {
  RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(2 * size, 0, comm));

  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > mat = Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 4);

  int numUnks  = 2;
  ST valuesA[] = {-1.0, 2.0, 7.0, -1.0};
  ST valuesB[] = {-1.0, 2.0, 9.0, -1.0};
  GO iTempA[]  = {-numUnks, 0, 1, numUnks}, indices[4];
  GO iTempB[]  = {-numUnks, 0, -1, numUnks};
  ST* vPtr;
  GO* iPtr;

  for (size_t i = 0; i < map->getLocalNumElements() / numUnks; i++) {
    int count = 4;
    GO gidA   = map->getGlobalElement(2 * i);
    GO gidB   = gidA + 1;

    for (int n = 0; n < numUnks; n++) {
      GO* iTemp = (n == 0) ? iTempA : iTempB;
      GO gid    = (n == 0) ? gidA : gidB;

      indices[0] = gid + iTemp[0];
      indices[1] = gid + iTemp[1];
      indices[2] = gid + iTemp[2];
      indices[3] = gid + iTemp[3];

      vPtr = (n == 0) ? valuesA : valuesB;
      iPtr = indices;

      if (gid < numUnks) {
        vPtr++;
        iPtr  = &indices[1];
        count = 3;
      } else if (gid > map->getMaxAllGlobalIndex() - numUnks)
        count = 3;

      mat->insertGlobalValues(gid, Teuchos::ArrayView<GO>(iPtr, count),
                              Teuchos::ArrayView<ST>(vPtr, count));
    }
  }

  mat->fillComplete();

  Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<ST, LO, GO, NT> >::writeSparseFile("strided.mm",
                                                                                    mat);

  return mat;
}

TEUCHOS_UNIT_TEST(tStratimikosFactory, test_Defaults) {
  using Teuchos::ParameterList;
  using Teuchos::RCP;

// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // build epetra operator
  RCP<Epetra_Operator> eA              = buildSystem(comm, 5);
  RCP<Thyra::LinearOpBase<double> > tA = Thyra::nonconstEpetraLinearOp(eA);

  // build stratimikos factory, adding Teko's version
  Stratimikos::DefaultLinearSolverBuilder stratFactory;
  stratFactory.setPreconditioningStrategyFactory(
      Teuchos::abstractFactoryStd<Thyra::PreconditionerFactoryBase<double>,
                                  Teko::StratimikosFactory>(),
      "Teko");
  RCP<const ParameterList> validParams = stratFactory.getValidParameters();
  stratFactory.setParameterList(Teuchos::rcp(new Teuchos::ParameterList(*validParams)));

  // print out Teko's parameter list and fail if it doesn't exist!
  TEST_NOTHROW(
      validParams->sublist("Preconditioner Types")
          .sublist("Teko")
          .print(out, ParameterList::PrintOptions().showDoc(true).indent(2).showTypes(true)));

  // build teko preconditioner factory
  RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      stratFactory.createPreconditioningStrategy("Teko");

  // make sure factory is built
  TEST_ASSERT(precFactory != Teuchos::null);

  // build preconditioner
  RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory, tA);
  TEST_ASSERT(prec != Teuchos::null);

  // build an operator to test against
  RCP<const Teko::InverseLibrary> invLib  = Teko::InverseLibrary::buildFromStratimikos();
  RCP<const Teko::InverseFactory> invFact = invLib->getInverseFactory("Amesos");
  Teko::LinearOp testOp                   = Teko::buildInverse(*invFact, tA);

  Teko::LinearOp precOp = prec->getUnspecifiedPrecOp();
  TEST_ASSERT(precOp != Teuchos::null);

  Thyra::LinearOpTester<double> tester;
  tester.show_all_tests(true);
  tester.set_all_error_tol(0);
  TEST_ASSERT(tester.compare(*precOp, *testOp, Teuchos::ptrFromRef(out)));
}

TEUCHOS_UNIT_TEST(tStratimikosFactory, test_MultipleCalls) {
  using Teuchos::ParameterList;
  using Teuchos::RCP;

// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // build epetra operator
  RCP<Epetra_Operator> eA              = buildStridedSystem(comm, 5);
  RCP<Thyra::LinearOpBase<double> > tA = Thyra::nonconstEpetraLinearOp(eA);

  Stratimikos::DefaultLinearSolverBuilder stratFactory;
  RCP<ParameterList> params = Teuchos::rcp(new ParameterList(*stratFactory.getValidParameters()));
  ParameterList& tekoList   = params->sublist("Preconditioner Types").sublist("Teko");
  tekoList.set("Write Block Operator", false);
  tekoList.set("Test Block Operator", false);
  tekoList.set("Strided Blocking", "1 1");
  tekoList.set("Inverse Type", "BGS");
  ParameterList& ifl = tekoList.sublist("Inverse Factory Library");
  ifl.sublist("BGS").set("Type", "Block Gauss-Seidel");
  ifl.sublist("BGS").set("Inverse Type", "Amesos");

  Teko::addTekoToStratimikosBuilder(stratFactory, "Teko");
  stratFactory.setParameterList(params);
  RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      stratFactory.createPreconditioningStrategy("Teko");

  // build teko preconditioner factory
  RCP<Thyra::PreconditionerBase<double> > prec_a = Thyra::prec<double>(*precFactory, tA);
  Teko::LinearOp precOp_a                        = prec_a->getUnspecifiedPrecOp();
  TEST_ASSERT(precOp_a != Teuchos::null);

  // try to do it again
  RCP<Thyra::PreconditionerBase<double> > prec_b = Thyra::prec<double>(*precFactory, tA);
  Teko::LinearOp precOp_b                        = prec_b->getUnspecifiedPrecOp();
  TEST_ASSERT(precOp_b != Teuchos::null);
}

TEUCHOS_UNIT_TEST(tStratimikosFactory, test_BlockGaussSeidel) {
  using Teuchos::ParameterList;
  using Teuchos::RCP;

// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // build epetra operator
  RCP<Epetra_Operator> eA              = buildStridedSystem(comm, 5);
  RCP<Thyra::LinearOpBase<double> > tA = Thyra::nonconstEpetraLinearOp(eA);

  // build stratimikos factory, adding Teko's version
  Stratimikos::DefaultLinearSolverBuilder stratFactory;
  stratFactory.setPreconditioningStrategyFactory(
      Teuchos::abstractFactoryStd<Thyra::PreconditionerFactoryBase<double>,
                                  Teko::StratimikosFactory>(),
      "Teko");
  RCP<ParameterList> params = Teuchos::rcp(new ParameterList(*stratFactory.getValidParameters()));
  ParameterList& tekoList   = params->sublist("Preconditioner Types").sublist("Teko");
  tekoList.set("Write Block Operator", false);
  tekoList.set("Test Block Operator", false);
  tekoList.set("Strided Blocking", "1 1");
  tekoList.set("Inverse Type", "BGS");
  ParameterList& ifl = tekoList.sublist("Inverse Factory Library");
  ifl.sublist("BGS").set("Type", "Block Gauss-Seidel");
  ifl.sublist("BGS").set("Inverse Type", "Amesos");

  // RCP<Thyra::PreconditionerFactoryBase<double> > precFactory
  //       = stratFactory.createPreconditioningStrategy("Teko");

  // build operator to test against
  Teko::LinearOp testOp;
  {
    Teuchos::ParameterList iflCopy(ifl);
    RCP<Epetra_Operator> strided_eA  = Teuchos::rcp(new Teko::Epetra::StridedEpetraOperator(2, eA));
    RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(iflCopy);
    RCP<const Teko::InverseFactory> invFact = invLib->getInverseFactory("BGS");
    RCP<Teko::Epetra::InverseFactoryOperator> invFactOp =
        Teuchos::rcp(new Teko::Epetra::InverseFactoryOperator(invFact));
    invFactOp->buildInverseOperator(strided_eA);

    testOp = Thyra::epetraLinearOp(invFactOp, Thyra::NOTRANS, Thyra::EPETRA_OP_APPLY_APPLY_INVERSE);
  }

  stratFactory.setParameterList(params);
  RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      stratFactory.createPreconditioningStrategy("Teko");

  // build teko preconditioner factory
  RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory, tA);
  Teko::LinearOp precOp                        = prec->getUnspecifiedPrecOp();
  TEST_ASSERT(precOp != Teuchos::null);

  Thyra::LinearOpTester<double> tester;
  tester.show_all_tests(true);
  tester.set_all_error_tol(0);
  TEST_ASSERT(tester.compare(*precOp, *testOp, Teuchos::ptrFromRef(out)));
}

TEUCHOS_UNIT_TEST(tStratimikosFactory, test_RelatedFunctions) {
  using Teuchos::ParameterList;
  using Teuchos::RCP;

// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // build epetra operator
  RCP<Epetra_Operator> eA              = buildStridedSystem(comm, 5);
  RCP<Thyra::LinearOpBase<double> > tA = Thyra::nonconstEpetraLinearOp(eA);

  {
    // build stratimikos factory, adding Teko's version
    Stratimikos::DefaultLinearSolverBuilder stratFactory;

    Teko::addTekoToStratimikosBuilder(stratFactory);
    TEST_THROW(Teko::addTekoToStratimikosBuilder(stratFactory), std::logic_error);
    Teko::addTekoToStratimikosBuilder(stratFactory, "Teko-2");

    TEST_NOTHROW(
        stratFactory.getValidParameters()->sublist("Preconditioner Types").sublist("Teko"));
    TEST_NOTHROW(
        stratFactory.getValidParameters()->sublist("Preconditioner Types").sublist("Teko-2"));
  }

  {
    Teuchos::RCP<Teko::RequestHandler> rh = Teuchos::rcp(new Teko::RequestHandler);

    // build stratimikos factory, adding Teko's version
    Stratimikos::DefaultLinearSolverBuilder stratFactory;

    Teko::addTekoToStratimikosBuilder(stratFactory, rh);
    TEST_THROW(Teko::addTekoToStratimikosBuilder(stratFactory, rh), std::logic_error);
    Teko::addTekoToStratimikosBuilder(stratFactory, rh, "Teko-2");

    TEST_NOTHROW(
        stratFactory.getValidParameters()->sublist("Preconditioner Types").sublist("Teko"));
    TEST_NOTHROW(
        stratFactory.getValidParameters()->sublist("Preconditioner Types").sublist("Teko-2"));

    RCP<ParameterList> params = Teuchos::rcp(new ParameterList(*stratFactory.getValidParameters()));
    ParameterList& tekoList   = params->sublist("Preconditioner Types").sublist("Teko");
    tekoList.set("Write Block Operator", false);
    tekoList.set("Test Block Operator", false);
    tekoList.set("Strided Blocking", "1 1");
    tekoList.set("Inverse Type", "BGS");
    ParameterList& ifl = tekoList.sublist("Inverse Factory Library");
    ifl.sublist("BGS").set("Type", "Block Gauss-Seidel");
    ifl.sublist("BGS").set("Inverse Type", "Amesos");
    stratFactory.setParameterList(params);

    RCP<Teko::StratimikosFactory> precFactory = Teuchos::rcp_dynamic_cast<Teko::StratimikosFactory>(
        stratFactory.createPreconditioningStrategy("Teko-2"));
    TEST_EQUALITY(precFactory->getRequestHandler(), rh);
  }
}

TEUCHOS_UNIT_TEST(tStratimikosFactory, test_multi_use) {
  using Teuchos::ParameterList;
  using Teuchos::RCP;

// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // build epetra operator
  RCP<Epetra_Operator> eA              = buildSystem(comm, 5);
  RCP<Thyra::LinearOpBase<double> > tA = Thyra::nonconstEpetraLinearOp(eA);

  // build stratimikos factory, adding Teko's version
  Stratimikos::DefaultLinearSolverBuilder stratFactory;
  stratFactory.setPreconditioningStrategyFactory(
      Teuchos::abstractFactoryStd<Thyra::PreconditionerFactoryBase<double>,
                                  Teko::StratimikosFactory>(),
      "Teko");
  RCP<const ParameterList> validParams = stratFactory.getValidParameters();
  stratFactory.setParameterList(Teuchos::rcp(new Teuchos::ParameterList(*validParams)));

  // print out Teko's parameter list and fail if it doesn't exist!
  TEST_NOTHROW(
      validParams->sublist("Preconditioner Types")
          .sublist("Teko")
          .print(out, ParameterList::PrintOptions().showDoc(true).indent(2).showTypes(true)));

  // build teko preconditioner factory
  RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      stratFactory.createPreconditioningStrategy("Teko");

  // make sure factory is built
  TEST_ASSERT(precFactory != Teuchos::null);

  // try using a different preconditioner each time
  RCP<Thyra::PreconditionerBase<double> > prec;
  for (int i = 0; i < 2; i++) {
    prec = precFactory->createPrec();

    RCP<const Thyra::LinearOpSourceBase<double> > losb =
        rcp(new Thyra::DefaultLinearOpSource<double>(tA));

    precFactory->initializePrec(losb, prec.get());

    RCP<Teko::StratimikosFactory> stratFact =
        rcp_dynamic_cast<Teko::StratimikosFactory>(precFactory);
    const std::vector<int>& decomp = stratFact->getDecomposition();

    TEST_EQUALITY(decomp.size(), 1);
    TEST_EQUALITY(decomp[0], 1);
  }

  // try using a single preconditioner multiple times
  prec = precFactory->createPrec();
  for (int i = 0; i < 2; i++) {
    RCP<const Thyra::LinearOpSourceBase<double> > losb =
        rcp(new Thyra::DefaultLinearOpSource<double>(tA));

    precFactory->initializePrec(losb, prec.get());

    RCP<Teko::StratimikosFactory> stratFact =
        rcp_dynamic_cast<Teko::StratimikosFactory>(precFactory);
    const std::vector<int>& decomp = stratFact->getDecomposition();

    TEST_EQUALITY(decomp.size(), 1);
    TEST_EQUALITY(decomp[0], 1);
  }
}
