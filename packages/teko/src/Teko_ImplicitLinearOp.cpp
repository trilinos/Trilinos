// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_ImplicitLinearOp.hpp"

namespace Teko {

using Teuchos::RCP;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

using Thyra::MultiVectorBase;

bool ImplicitLinearOp::opSupportedImpl(const Thyra::EOpTransp M_trans) const {
  return (M_trans == Thyra::NOTRANS);
}

void ImplicitLinearOp::applyImpl(const Thyra::EOpTransp M_trans,
                                 const Thyra::MultiVectorBase<double>& x,
                                 const Teuchos::Ptr<Thyra::MultiVectorBase<double> >& y,
                                 const double alpha, const double beta) const {
  TEUCHOS_TEST_FOR_EXCEPTION(M_trans != Thyra::NOTRANS, std::runtime_error,
                             "Linear operators of inherited type Teko::ImplicitLinearOp "
                             "cannot handle conjugation (yet!)");

  MultiVector srcX  = rcp_const_cast<MultiVectorBase<double> >(rcpFromRef(x));
  MultiVector destY = rcp_dynamic_cast<MultiVectorBase<double> >(rcpFromRef(*y));

  // call apply
  implicitApply(srcX, destY, alpha, beta);
}

}  // end namespace Teko
