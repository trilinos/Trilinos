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
    TEUCHOS_TEST_FOR_EXCEPTION(TpetraHelpers::isTpetraLinearOp(lo), std::runtime_error,
                               "BlkDiag not implemented for Tpetra operators");
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
