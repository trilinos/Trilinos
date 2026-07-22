// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFUALT_LINEAR_OP_SOURCE_HPP
#define THYRA_DEFUALT_LINEAR_OP_SOURCE_HPP


#include "Thyra_DefaultLinearOpSource_decl.hpp"
#include "Thyra_LinearOpBase.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template <class Scalar>
DefaultLinearOpSource<Scalar>::DefaultLinearOpSource()
{}


template <class Scalar>
DefaultLinearOpSource<Scalar>::DefaultLinearOpSource(
  const Teuchos::RCP<LinearOpBase<Scalar> >    &op
  )
{
  op_.initialize(op);
}


template <class Scalar>
DefaultLinearOpSource<Scalar>::DefaultLinearOpSource(
  const Teuchos::RCP<const LinearOpBase<Scalar> >    &op
  )
{
  op_.initialize(op);
}


template <class Scalar>
void DefaultLinearOpSource<Scalar>::initialize(
  const Teuchos::RCP<LinearOpBase<Scalar> >    &op
  )
{
  op_.initialize(op);
}


template <class Scalar>
void DefaultLinearOpSource<Scalar>::initialize(
  const Teuchos::RCP<const LinearOpBase<Scalar> >    &op
  )
{
  op_.initialize(op);
}


template <class Scalar>
void DefaultLinearOpSource<Scalar>::uninitialize()
{
  op_.uninitialize();
}


// Overridden from LinearOpSourceBase


template <class Scalar>
bool DefaultLinearOpSource<Scalar>::isOpConst() const
{
  return op_.isConst();
}


template <class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
DefaultLinearOpSource<Scalar>::getNonconstOp()
{
  return op_.getNonconstObj();
}


template <class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultLinearOpSource<Scalar>::getOp() const
{
  return op_.getConstObj();
}


} // namespace Thyra


#endif // THYRA_DEFUALT_LINEAR_OP_SOURCE_HPP
