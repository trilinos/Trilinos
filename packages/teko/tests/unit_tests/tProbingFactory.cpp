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

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_Config.h"
#include "Teko_DiagonallyScaledPreconditionerFactory.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_PreconditionerLinearOp.hpp"
#include "Teko_ProbingPreconditionerFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_ConfigDefs.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace {

using ST = Teko::ST;
using LO = Teko::LO;
using GO = Teko::GO;
using NT = Teko::NT;

const RCP<Thyra::LinearOpBase<double> > buildSystem(const RCP<const Teuchos::Comm<int> >& comm,
                                                    GO size) {
  RCP<const Tpetra::Map<LO, GO, NT> > map = rcp(new const Tpetra::Map<LO, GO, NT>(size, 0, comm));

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
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat->getRangeMap()), mat);
}

}  // namespace

TEUCHOS_UNIT_TEST(tProbingFactory, basic_test) {
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  Teko::LinearOp lo = buildSystem(Comm, 10);

  Teuchos::RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
  Teuchos::RCP<Teko::InverseFactory> directSolveFactory = invLib->getInverseFactory("Ifpack2");

  Teuchos::RCP<Teko::ProbingPreconditionerFactory> probeFact =
      rcp(new Teko::ProbingPreconditionerFactory);
  probeFact->setGraphOperator(lo);
  probeFact->setInverseFactory(directSolveFactory);

  RCP<Teko::InverseFactory> invFact =
      Teuchos::rcp(new Teko::PreconditionerInverseFactory(probeFact, Teuchos::null));

  Teko::LinearOp probedInverse = Teko::buildInverse(*invFact, lo);
  Teko::LinearOp invLo         = Teko::buildInverse(*directSolveFactory, lo);

  Thyra::LinearOpTester<double> tester;
  tester.dump_all(true);
  tester.show_all_tests(true);

  {
    const bool result = tester.compare(*probedInverse, *invLo, Teuchos::ptrFromRef(out));
    if (!result) {
      out << "Apply: FAILURE" << std::endl;
      success = false;
    } else
      out << "Apply: SUCCESS" << std::endl;
  }
}

TEUCHOS_UNIT_TEST(tProbingFactory, parameterlist_constr) {
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  Teko::LinearOp lo = buildSystem(Comm, 10);

  Teuchos::RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
  Teuchos::RCP<Teko::InverseFactory> directSolveFactory = invLib->getInverseFactory("Ifpack2");

  {
    Teuchos::ParameterList pl;
    pl.set("Inverse Type", "Ifpack2");
    pl.set("Probing Graph Operator", lo);

    Teuchos::RCP<Teko::ProbingPreconditionerFactory> probeFact =
        rcp(new Teko::ProbingPreconditionerFactory);
    probeFact->initializeFromParameterList(pl);

    RCP<Teko::InverseFactory> invFact =
        Teuchos::rcp(new Teko::PreconditionerInverseFactory(probeFact, Teuchos::null));

    Teko::LinearOp probedInverse = Teko::buildInverse(*invFact, lo);
    Teko::LinearOp invLo         = Teko::buildInverse(*directSolveFactory, lo);

    Thyra::LinearOpTester<double> tester;
    tester.dump_all(true);
    tester.show_all_tests(true);

    {
      const bool result = tester.compare(*probedInverse, *invLo, Teuchos::ptrFromRef(out));
      if (!result) {
        out << "Apply: FAILURE" << std::endl;
        success = false;
      } else
        out << "Apply: SUCCESS" << std::endl;
    }
  }

  {
    RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tOp =
        Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(lo, true);
    RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > tMat =
        Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT> >(
            tOp->getConstTpetraOperator(), true);
    Teuchos::RCP<const Tpetra::CrsGraph<LO, GO, NT> > theGraph = tMat->getCrsGraph();

    Teuchos::ParameterList pl;
    pl.set("Inverse Type", "Ifpack2");
    pl.set("Probing Graph", theGraph);

    Teuchos::RCP<Teko::ProbingPreconditionerFactory> probeFact =
        rcp(new Teko::ProbingPreconditionerFactory);
    probeFact->initializeFromParameterList(pl);

    RCP<Teko::InverseFactory> invFact =
        Teuchos::rcp(new Teko::PreconditionerInverseFactory(probeFact, Teuchos::null));

    Teko::LinearOp probedInverse = Teko::buildInverse(*invFact, lo);
    Teko::LinearOp invLo         = Teko::buildInverse(*directSolveFactory, lo);

    Thyra::LinearOpTester<double> tester;
    tester.dump_all(true);
    tester.show_all_tests(true);

    {
      const bool result = tester.compare(*probedInverse, *invLo, Teuchos::ptrFromRef(out));
      if (!result) {
        out << "Apply: FAILURE" << std::endl;
        success = false;
      } else
        out << "Apply: SUCCESS" << std::endl;
    }
  }
}

TEUCHOS_UNIT_TEST(tProbingFactory, invlib_constr) {
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  Teko::LinearOp lo = buildSystem(Comm, 10);

  Teuchos::ParameterList subList;
  subList.set("Type", "Probing Preconditioner");
  subList.set("Inverse Type", "Ifpack2");
  subList.set("Probing Graph Operator", lo);

  Teuchos::ParameterList pl;
  pl.set("Prober", subList);

  Teuchos::RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(pl);
  Teuchos::RCP<Teko::InverseFactory> proberFactory      = invLib->getInverseFactory("Prober");
  Teuchos::RCP<Teko::InverseFactory> directSolveFactory = invLib->getInverseFactory("Ifpack2");

  {
    Teko::LinearOp probedInverse = Teko::buildInverse(*proberFactory, lo);
    Teko::LinearOp invLo         = Teko::buildInverse(*directSolveFactory, lo);

    Thyra::LinearOpTester<double> tester;
    tester.dump_all(true);
    tester.show_all_tests(true);

    {
      const bool result = tester.compare(*probedInverse, *invLo, Teuchos::ptrFromRef(out));
      if (!result) {
        out << "Apply: FAILURE" << std::endl;
        success = false;
      } else
        out << "Apply: SUCCESS" << std::endl;
    }
  }
}
