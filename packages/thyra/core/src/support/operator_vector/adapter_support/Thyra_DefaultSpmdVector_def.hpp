// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_SPMD_VECTOR_DEF_HPP
#define THYRA_DEFAULT_SPMD_VECTOR_DEF_HPP


#include "Thyra_DefaultSpmdVector_decl.hpp"
#include "Thyra_SpmdVectorDefaultBase.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
DefaultSpmdVector<Scalar>::DefaultSpmdVector()
  :stride_(0)
{}


template<class Scalar>
DefaultSpmdVector<Scalar>::DefaultSpmdVector(
  const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdSpace_in,
  const ArrayRCP<Scalar> &localValues,
  const Ordinal stride
  )
{
  initialize(spmdSpace_in, localValues, stride);
}


template<class Scalar>
void DefaultSpmdVector<Scalar>::initialize(
  const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdSpace_in
  ,const ArrayRCP<Scalar> &localValues
  ,const Ordinal stride
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(is_null(spmdSpace_in));
  TEUCHOS_TEST_FOR_EXCEPT(spmdSpace_in->localSubDim() > 0 && localValues.get()==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(stride==0);
#endif
  spmdSpace_ = spmdSpace_in;
  localValues_ = localValues;
  stride_ = stride;
  this->updateSpmdSpace();
}


template<class Scalar>
void DefaultSpmdVector<Scalar>::uninitialize(
  RCP<const SpmdVectorSpaceBase<Scalar> > *spmdSpace_in
  ,ArrayRCP<Scalar> *localValues
  ,Ordinal *stride
  )
{
  if(spmdSpace_in) *spmdSpace_in = spmdSpace_;
  if(localValues) *localValues = localValues_;
  if(stride) *stride = stride_;

  spmdSpace_ = Teuchos::null;
  localValues_ = Teuchos::null;
  stride_ = 0;

  this->updateSpmdSpace();
}


// Overridden from SpmdMultiVectorBase


template<class Scalar>
RCP<const SpmdVectorSpaceBase<Scalar> >
DefaultSpmdVector<Scalar>::spmdSpaceImpl() const
{
  return spmdSpace_;
}


// Overridden from SpmdVectorBase


template<class Scalar>
void DefaultSpmdVector<Scalar>::getNonconstLocalVectorDataImpl(
  const Ptr<ArrayRCP<Scalar> > &localValues )
{
  *localValues = localValues_;
}


template<class Scalar>
void DefaultSpmdVector<Scalar>::getLocalVectorDataImpl(
  const Ptr<ArrayRCP<const Scalar> > &localValues ) const
{
  *localValues = localValues_;
}


} // end namespace Thyra


#endif // THYRA_DEFAULT_SPMD_VECTOR_DEF_HPP
