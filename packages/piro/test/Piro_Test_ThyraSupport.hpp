// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_TEST_THYRASUPPORT_HPP
#define PIRO_TEST_THYRASUPPORT_HPP

#include "Thyra_VectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DetachedVectorView.hpp"

#include "Thyra_VectorStdOps.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

namespace Piro {

namespace Test {

Teuchos::Array<double> arrayFromVector(const Thyra::VectorBase<double> &v)
{
  const Thyra::ConstDetachedVectorView<double> view(v, Thyra::Range1D(), true);
  return Teuchos::Array<double>(view.values(), view.values() + view.subDim());
}

Teuchos::RCP<Thyra::VectorBase<double> > vectorFromLinOp(const Thyra::LinearOpBase<double> &op, int col)
{
  const Teuchos::RCP<Thyra::VectorBase<double> > result = Thyra::createMember(op.range());
  const Teuchos::RCP<Thyra::VectorBase<double> > v = Thyra::createMember(op.domain());
  Thyra::put_scalar(0.0, v.ptr());
  Thyra::set_ele(col, 1.0, v.ptr());
  Thyra::apply(op, Thyra::NOTRANS, *v, result.ptr(), 1.0, 0.0);
  return result;
}

Teuchos::Array<double> arrayFromLinOp(const Thyra::LinearOpBase<double> &op, int col)
{
  return arrayFromVector(*vectorFromLinOp(op, col));
}

} // namespace Test

} // namespace Piro

#endif /* PIRO_TEST_THYRASUPPORT_HPP */
