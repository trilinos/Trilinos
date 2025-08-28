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

#include "Teko_DiagnosticPreconditionerFactory.hpp"
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

auto extract_diagnostic_preconditioner_factory(Teuchos::RCP<Teko::InverseFactory>& inverse) {
  auto factory = Teuchos::rcp_dynamic_cast<Teko::PreconditionerInverseFactory>(inverse, true);
  auto diagPrecFactory = Teuchos::rcp_dynamic_cast<Teko::DiagnosticPreconditionerFactory>(
      factory->getPrecFactory(), true);
  return diagPrecFactory;
}

auto extract_adaptive_preconditioner_factory(Teuchos::RCP<Teko::InverseFactory>& inverse) {
  auto factory = Teuchos::rcp_dynamic_cast<Teko::PreconditionerInverseFactory>(inverse, true);
  auto adaptPrecFactory = Teuchos::rcp_dynamic_cast<Teko::AdaptivePreconditionerFactory>(
      factory->getPrecFactory(), true);
  return adaptPrecFactory;
}

auto extract_build_counts(Teuchos::RCP<Teko::InverseFactory>& inverse,
                          RCP<Teko::InverseLibrary>& invLibrary) {
  auto adapt      = extract_adaptive_preconditioner_factory(inverse);
  auto jacobi_inv = adapt->get_inverses()[0];
  auto direct_inv = adapt->get_inverses()[1];

  auto jacobi_diag = extract_diagnostic_preconditioner_factory(jacobi_inv);
  auto direct_diag = extract_diagnostic_preconditioner_factory(direct_inv);
  return std::make_pair(jacobi_diag->numInitialBuilds(), direct_diag->numInitialBuilds());
}

TEUCHOS_UNIT_TEST(tAdaptivePreconditionerFactory, switch_to_new_solver) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  const double targetReduction = 1e-8;

  Teuchos::ParameterList pl;
  Teuchos::ParameterList& adaptivePl = pl.sublist("adapt");
  adaptivePl.set("Type", "Adaptive");
  adaptivePl.set("Inverse Type 1", "diag_jacobi");
  adaptivePl.set("Inverse Type 2", "diag_direct");
  adaptivePl.set("Target Residual Reduction", targetReduction);

  Teuchos::ParameterList& diag_ifpack2 = pl.sublist("diag_jacobi");
  diag_ifpack2.set("Type", "Diagnostic Inverse");
  diag_ifpack2.set("Inverse Factory", "Jacobi");
  diag_ifpack2.set("Descriptive Label", "Ifpack2_preconditioner");

  Teuchos::ParameterList& jacobi_pl = pl.sublist("Jacobi");
  jacobi_pl.set("Type", "Ifpack2");
  jacobi_pl.set("Prec Type", "Relaxation");
  jacobi_pl.sublist("Ifpack2 Settings").set("relaxation: type", "Jacobi");

  Teuchos::ParameterList& diag_amesos2 = pl.sublist("diag_direct");
  diag_amesos2.set("Type", "Diagnostic Inverse");
  diag_amesos2.set("Inverse Factory", "Amesos2");
  diag_amesos2.set("Descriptive Label", "Amesos2_preconditioner");

  RCP<Teko::InverseLibrary> invLibrary = Teko::InverseLibrary::buildFromParameterList(pl);
  RCP<Teko::InverseFactory> invFact    = invLibrary->getInverseFactory("adapt");
  Teko::LinearOp A                     = buildSystem(Comm, 1000);
  Teko::ModifiableLinearOp invA        = Teko::buildInverse(*invFact, A);

  Teko::MultiVector b = Thyra::createMember(A->domain());
  Teko::MultiVector x = Thyra::createMember(A->range());
  Thyra::randomize(-1.0, 1.0, b.ptr());

  Teko::MultiVector residual = Teko::deepcopy(b);

  success = true;
  {
    auto [jacobi_build_count, direct_build_count] = extract_build_counts(invFact, invLibrary);
    success &= jacobi_build_count == 1;
    success &= direct_build_count == 0;
  }

  Teko::applyOp(invA, b, x);
  Teko::applyOp(A, x, residual, -1.0, 1.0);

  // relative residual from a single sweep of Jacobi isn't sufficient to reach 1e-8 relative
  // residual criterion hence, adaptive solver will switch to the direct solver
  {
    auto [jacobi_build_count, direct_build_count] = extract_build_counts(invFact, invLibrary);
    success &= jacobi_build_count == 1;
    success &= direct_build_count == 1;
  }

  const auto resNorm    = Teko::norm_2(residual, 0);
  const auto rhsNorm    = Teko::norm_2(b, 0);
  const auto relResNorm = resNorm / rhsNorm;

  success &= relResNorm <= targetReduction;
}

TEUCHOS_UNIT_TEST(tAdaptivePreconditionerFactory, use_initial_solver) {
  // build global (or serial communicator)
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  const double targetReduction = 0.9;

  Teuchos::ParameterList pl;
  Teuchos::ParameterList& adaptivePl = pl.sublist("adapt");
  adaptivePl.set("Type", "Adaptive");
  adaptivePl.set("Inverse Type 1", "diag_jacobi");
  adaptivePl.set("Inverse Type 2", "diag_direct");
  adaptivePl.set("Target Residual Reduction", targetReduction);

  Teuchos::ParameterList& diag_ifpack2 = pl.sublist("diag_jacobi");
  diag_ifpack2.set("Type", "Diagnostic Inverse");
  diag_ifpack2.set("Inverse Factory", "Jacobi");
  diag_ifpack2.set("Descriptive Label", "Ifpack2_preconditioner");

  Teuchos::ParameterList& jacobi_pl = pl.sublist("Jacobi");
  jacobi_pl.set("Type", "Ifpack2");
  jacobi_pl.set("Prec Type", "Relaxation");
  jacobi_pl.sublist("Ifpack2 Settings").set("relaxation: type", "Jacobi");

  Teuchos::ParameterList& diag_amesos2 = pl.sublist("diag_direct");
  diag_amesos2.set("Type", "Diagnostic Inverse");
  diag_amesos2.set("Inverse Factory", "Amesos2");
  diag_amesos2.set("Descriptive Label", "Amesos2_preconditioner");

  RCP<Teko::InverseLibrary> invLibrary = Teko::InverseLibrary::buildFromParameterList(pl);
  RCP<Teko::InverseFactory> invFact    = invLibrary->getInverseFactory("adapt");
  Teko::LinearOp A                     = buildSystem(Comm, 1000);
  Teko::ModifiableLinearOp invA        = Teko::buildInverse(*invFact, A);

  Teko::MultiVector b = Thyra::createMember(A->domain());
  Teko::MultiVector x = Thyra::createMember(A->range());
  Thyra::randomize(-1.0, 1.0, b.ptr());

  Teko::MultiVector residual = Teko::deepcopy(b);

  success = true;
  {
    auto [jacobi_build_count, direct_build_count] = extract_build_counts(invFact, invLibrary);
    success &= jacobi_build_count == 1;
    success &= direct_build_count == 0;
  }

  Teko::applyOp(invA, b, x);
  Teko::applyOp(A, x, residual, -1.0, 1.0);

  // In this instance, single sweep Jacobi is sufficient to reduce the relative residual.
  {
    auto [jacobi_build_count, direct_build_count] = extract_build_counts(invFact, invLibrary);
    success &= jacobi_build_count == 1;
    success &= direct_build_count == 0;
  }

  const auto resNorm    = Teko::norm_2(residual, 0);
  const auto rhsNorm    = Teko::norm_2(b, 0);
  const auto relResNorm = resNorm / rhsNorm;

  success &= relResNorm <= targetReduction;
}