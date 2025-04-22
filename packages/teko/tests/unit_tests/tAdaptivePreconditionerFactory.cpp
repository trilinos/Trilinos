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

#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

// Teko-Package includes
#include "Teko_ConfigDefs.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_AdaptivePreconditionerFactory.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_PreconditionerLinearOp.hpp"

#include "Thyra_TpetraLinearOp.hpp"

// Test-rig

typedef Teko::ST ST;
typedef Teko::LO LO;
typedef Teko::GO GO;
typedef Teko::NT NT;

using Teuchos::ParameterList;
using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

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

TEUCHOS_UNIT_TEST(tAdaptivePreconditionerFactory, relative_residual_test) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  const double targetReduction = 1e-4;

  Teuchos::ParameterList pl;
  Teuchos::ParameterList& adaptivePl = pl.sublist("adapt");
  adaptivePl.set("Type", "Adaptive");
  adaptivePl.set("Inverse Type 1", "Ifpack2");
  adaptivePl.set("Inverse Type 2", "Amesos2");

  RCP<Teko::InverseLibrary> invLibrary = Teko::InverseLibrary::buildFromParameterList(pl);
  RCP<Teko::InverseFactory> invFact    = invLibrary->getInverseFactory("adapt");
  Teko::LinearOp A                     = buildSystem(Comm, 1000);
  Teko::ModifiableLinearOp invA        = Teko::buildInverse(*invFact, A);

  Teko::MultiVector b = Thyra::createMember(A->domain());
  Teko::MultiVector x = Thyra::createMember(A->range());
  Thyra::randomize(-1.0, 1.0, b.ptr());

  Teko::MultiVector residual = Teko::deepcopy(b);

  Teko::applyOp(invA, b, x);
  Teko::applyOp(A, x, residual, -1.0, 1.0);

  const auto resNorm    = Teko::norm_2(residual, 0);
  const auto rhsNorm    = Teko::norm_2(b, 0);
  const auto relResNorm = resNorm / rhsNorm;

  success = relResNorm <= targetReduction;
}
