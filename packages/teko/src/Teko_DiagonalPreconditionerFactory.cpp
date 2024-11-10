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
#ifdef TEKO_HAVE_EPETRA
#include "Thyra_get_Epetra_Operator.hpp"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_PointToBlockDiagPermute.h"
#endif

#include "Teko_TpetraHelpers.hpp"
#ifdef TEKO_HAVE_EPETRA
#include "Thyra_EpetraLinearOp.hpp"
#endif
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
#ifdef TEKO_HAVE_EPETRA
    // Sanity check the state
    DiagonalPrecondState& MyState = Teuchos::dyn_cast<DiagonalPrecondState>(state);

    // Get the underlying Epetra_CrsMatrix, if we have one
    Teuchos::RCP<const Epetra_Operator> eo = Thyra::get_Epetra_Operator(*lo);
    TEUCHOS_ASSERT(eo != Teuchos::null);
    Teuchos::RCP<const Epetra_CrsMatrix> MAT =
        Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(eo);
    TEUCHOS_ASSERT(MAT != Teuchos::null);

    // Create a new EpetraExt_PointToBlockDiagPermute for the state object, if we don't have one
    Teuchos::RCP<EpetraExt_PointToBlockDiagPermute> BDP;
    if (MyState.BDP_ == Teuchos::null) {
      BDP = Teuchos::rcp(new EpetraExt_PointToBlockDiagPermute(*MAT));
      BDP->SetParameters(List_);
      BDP->Compute();
      MyState.BDP_ = BDP;
    }

    RCP<Epetra_FECrsMatrix> Hcrs = rcp(MyState.BDP_->CreateFECrsMatrix());
    return Thyra::epetraLinearOp(Hcrs);

    // Build the LinearOp object  (NTS: swapping the range and domain)
    // LinearOp MyOp = Teuchos::rcp(new
    // DiagonalPreconditionerOp(MyState.BDP_,lo->domain(),lo->range()));
#endif
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
