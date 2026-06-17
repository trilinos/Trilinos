// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Teuchos includes
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

// Tpetra includes
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_TpetraBlockPreconditioner.hpp"
#include "Teko_ConfigDefs.hpp"

#include "ExamplePreconditionerFactory.cpp"

#include <iostream>
#include <vector>

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

// utility function to construct Tpetra operators
RCP<Tpetra::CrsMatrix<Teko::ST, Teko::LO, Teko::GO, Teko::NT> > build2x2(
    double a, double b, double c, double d, const RCP<const Teuchos::Comm<int> >& comm) {
  using ST = Teko::ST;
  using LO = Teko::LO;
  using GO = Teko::GO;
  using NT = Teko::NT;

  using map_t = Tpetra::Map<LO, GO, NT>;
  using crs_t = Tpetra::CrsMatrix<ST, LO, GO, NT>;

  RCP<const map_t> map = rcp(new map_t(2, 0, comm));
  RCP<crs_t> matrix    = rcp(new crs_t(map, 2));

  Teuchos::Array<GO> indices(2);
  Teuchos::Array<ST> values(2);

  indices[0] = 0;
  indices[1] = 1;

  // build first row
  if (map->isNodeGlobalElement(0)) {
    values[0] = a;
    values[1] = b;
    matrix->insertGlobalValues(0, indices(), values());
  }

  // build second row
  if (map->isNodeGlobalElement(1)) {
    values[0] = c;
    values[1] = d;
    matrix->insertGlobalValues(1, indices(), values());
  }

  matrix->fillComplete();

  return matrix;
}

int main(int argc, char* argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

  using ST = Teko::ST;
  using LO = Teko::LO;
  using GO = Teko::GO;
  using NT = Teko::NT;

  using vec_t = Tpetra::Vector<ST, LO, GO, NT>;
  using mv_t  = Tpetra::MultiVector<ST, LO, GO, NT>;

  auto comm = Tpetra::getDefaultComm();

  // build the sub blocks
  auto mat_00         = build2x2(1, 2, 2, 1, comm);
  Teko::LinearOp A_00 = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat_00->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat_00->getDomainMap()), mat_00);

  auto mat_01         = build2x2(0, -1, 3, 4, comm);
  Teko::LinearOp A_01 = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat_01->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat_01->getDomainMap()), mat_01);

  auto mat_10         = build2x2(1, 6, -3, 2, comm);
  Teko::LinearOp A_10 = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat_10->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat_10->getDomainMap()), mat_10);

  auto mat_11         = build2x2(2, 1, 1, 2, comm);
  Teko::LinearOp A_11 = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat_11->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(mat_11->getDomainMap()), mat_11);

  // build the Tpetra operator wrapper
  Teuchos::RCP<Teko::TpetraHelpers::TpetraOperatorWrapper> A = Teuchos::rcp(
      new Teko::TpetraHelpers::TpetraOperatorWrapper(Teko::block2x2(A_00, A_01, A_10, A_11)));

  // build the Tpetra vectors
  RCP<vec_t> b = rcp(new vec_t(A->getRangeMap()));
  RCP<vec_t> x = rcp(new vec_t(A->getDomainMap()));
  x->putScalar(0.0);

  // build the RHS vector
  b->replaceGlobalValue(0, 1.0);
  b->replaceGlobalValue(1, 2.0);
  b->replaceGlobalValue(2, 3.0);
  b->replaceGlobalValue(3, 4.0);

  // Build the preconditioner
  /////////////////////////////////////////////////////////

  RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();

  RCP<const Teko::InverseFactory> inverse = invLib->getInverseFactory("Amesos2");

  RCP<Teko::BlockPreconditionerFactory> precFact =
      rcp(new ExamplePreconditionerFactory(inverse, 0.9));

  Teko::TpetraHelpers::TpetraBlockPreconditioner prec(precFact);
  prec.buildPreconditioner(A);

  // apply the preconditioner
  RCP<mv_t> B = b;
  RCP<mv_t> X = x;

  prec.apply(*B, *X);

  x->describe(*(Teuchos::VerboseObjectBase::getDefaultOStream()),
              Teuchos::EVerbosityLevel::VERB_EXTREME);

  return 0;
}
