// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_DiagonalPreconditionerOp.hpp"

#include "TpetraExt_PointToBlockDiagPermute_decl.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Tpetra_MultiVector.hpp"

using Teuchos::RCP;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

using Thyra::MultiVectorBase;

namespace Teko {

DiagonalPreconditionerOp::DiagonalPreconditionerOp(
    Teuchos::RCP<Tpetra::Ext::PointToBlockDiagPermute<ST, LO, GO, NT>> BDP, const VectorSpace range,
    const VectorSpace domain)
    : BDP_(BDP), range_(range), domain_(domain) {}

void DiagonalPreconditionerOp::implicitApply(const MultiVector& x, MultiVector& y,
                                             const double alpha, const double beta) const {
  TEUCHOS_TEST_FOR_EXCEPTION(BDP_ == Teuchos::null, std::runtime_error,
                             "DiagonalPreconditionerOp::implicitApply: null BDP_");

  RCP<const Tpetra::MultiVector<ST, LO, GO, NT>> x_ =
      Thyra::TpetraOperatorVectorExtraction<ST, LO, GO, NT>::getConstTpetraMultiVector(x);
  RCP<Tpetra::MultiVector<ST, LO, GO, NT>> y_ =
      Thyra::TpetraOperatorVectorExtraction<ST, LO, GO, NT>::getTpetraMultiVector(y);

  TEUCHOS_ASSERT(x_ != Teuchos::null);
  TEUCHOS_ASSERT(y_ != Teuchos::null);

  if (beta == 0.0) {
    BDP_->applyInverse(*x_, *y_);
    scale(alpha, y);
  } else {
    MultiVector y0 = deepcopy(y);
    BDP_->applyInverse(*x_, *y_);
    update(alpha, y, beta, y0);
  }
}

void DiagonalPreconditionerOp::describe(Teuchos::FancyOStream& out_arg,
                                        const Teuchos::EVerbosityLevel verbLevel) const {
  switch (verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW: out_arg << this->description() << std::endl; break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME: {
      if (BDP_ != Teuchos::null) {
        RCP<Tpetra::CrsMatrix<ST, LO, GO, NT>> H = BDP_->createCrsMatrix();
        if (H != Teuchos::null) H->describe(out_arg, verbLevel);
      }
      break;
    }
    default: TEUCHOS_TEST_FOR_EXCEPT(true);
  }
}

}  // end namespace Teko
