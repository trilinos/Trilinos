// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_DiagonalPreconditionerFactory.hpp"
#include "Teko_DiagonalPreconditionerOp.hpp"

#include "Teko_TpetraHelpers.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "TpetraExt_PointToBlockDiagPermute_decl.hpp"

using Teuchos::rcp;
using Teuchos::RCP;

namespace Teko {

DiagonalPrecondState::DiagonalPrecondState() {}

/*****************************************************/

DiagonalPreconditionerFactory::DiagonalPreconditionerFactory() {}

RCP<PreconditionerState> DiagonalPreconditionerFactory::buildPreconditionerState() const {
  return Teuchos::rcp(new DiagonalPrecondState());
}

LinearOp DiagonalPreconditionerFactory::buildPreconditionerOperator(
    LinearOp& lo, PreconditionerState& state) const {
  if (diagonalType_ == BlkDiag) {
    DiagonalPrecondState& MyState = Teuchos::dyn_cast<DiagonalPrecondState>(state);

    if (TpetraHelpers::isTpetraLinearOp(lo)) {
      RCP<const Thyra::TpetraLinearOp<ST, LO, GO, NT> > tlo =
          Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<ST, LO, GO, NT> >(lo, true);
      RCP<const Tpetra::Operator<ST, LO, GO, NT> > top = tlo->getConstTpetraOperator();
      RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > MAT =
          Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT> >(top, true);

      RCP<Tpetra::Ext::PointToBlockDiagPermute<ST, LO, GO, NT> > BDP;
      if (MyState.BDP_ == Teuchos::null) {
        BDP = Teuchos::rcp(new Tpetra::Ext::PointToBlockDiagPermute<ST, LO, GO, NT>(*MAT));
        BDP->setParameters(List_);
        BDP->compute();
        MyState.BDP_ = BDP;
      } else {
        BDP = MyState.BDP_;
      }

      RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > Hcrs = BDP->createCrsMatrix();
      return Thyra::tpetraLinearOp<ST, LO, GO, NT>(
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(Hcrs->getRangeMap()),
          Thyra::tpetraVectorSpace<ST, LO, GO, NT>(Hcrs->getDomainMap()), Hcrs);
    }

    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                               "DiagonalPreconditionerFactory::buildPreconditionerOperator: "
                               "BlkDiag requested, but operator is not Tpetra-supported.");
  }

  return getInvDiagonalOp(lo, diagonalType_);
}

void DiagonalPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList& pl) {
  List_ = pl;

  diagonalType_ = BlkDiag;
  if (pl.isParameter("Diagonal Type")) {
    diagonalType_ = getDiagonalType(pl.get<std::string>("Diagonal Type"));
    TEUCHOS_TEST_FOR_EXCEPT(diagonalType_ == NotDiag);
  }

  if (diagonalType_ == BlkDiag) {
    // Reset default to invert mode if the user hasn't specified something else
    Teuchos::ParameterList& SubList = List_.sublist("blockdiagmatrix: list");
    SubList.set("apply mode", SubList.get("apply mode", "invert"));
  }
}

}  // end namespace Teko
