// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_BlockImplicitLinearOp.hpp"

namespace Teko {

using Teuchos::RCP;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

using Thyra::ProductMultiVectorBase;

void BlockImplicitLinearOp::implicitApply(const Thyra::EOpTransp M_trans,
                                          const BlockedMultiVector& x, BlockedMultiVector& y,
                                          const double alpha, const double beta) const {
  TEUCHOS_TEST_FOR_EXCEPTION(M_trans != Thyra::NOTRANS, std::runtime_error,
                             "Linear operators of inherited type BlockImplicitLinearOp "
                             "cannot handle conjugation (yet!)");

  // call apply
  implicitApply(x, y, alpha, beta);
}

bool BlockImplicitLinearOp::opSupportedImpl(const Thyra::EOpTransp M_trans) const {
  return (M_trans == Thyra::NOTRANS);
}

void BlockImplicitLinearOp::applyImpl(const Thyra::EOpTransp M_trans,
                                      const Thyra::MultiVectorBase<double>& x,
                                      const Teuchos::Ptr<Thyra::MultiVectorBase<double> >& y,
                                      const double alpha, const double beta) const {
  // cast source vector
  RCP<const ProductMultiVectorBase<double> > src =
      rcp_dynamic_cast<const ProductMultiVectorBase<double> >(rcpFromRef(x));
  BlockedMultiVector srcX = rcp_const_cast<ProductMultiVectorBase<double> >(src);

  // cast destination vector
  BlockedMultiVector destY = rcp_dynamic_cast<ProductMultiVectorBase<double> >(rcpFromPtr(y));

  // call apply
  implicitApply(M_trans, srcX, destY, alpha, beta);
}

}  // end namespace Teko
