// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_LINEAR_OP_SCALAR_PROD_DEF_HPP
#define THYRA_LINEAR_OP_SCALAR_PROD_DEF_HPP

#include "Thyra_LinearOpScalarProd_decl.hpp"
#include "Thyra_ScalarProdBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"


namespace Thyra {


// Constructors, initializers, accessors


template<class Scalar>
LinearOpScalarProd<Scalar>::LinearOpScalarProd()
{}


template<class Scalar>
LinearOpScalarProd<Scalar>::LinearOpScalarProd(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &op_in )
{
  this->initialize(op_in);
}


template<class Scalar>
void LinearOpScalarProd<Scalar>::initialize(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &op_in
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(is_null(op_in));
#endif
  op_ = op_in;
}


template<class Scalar>
void LinearOpScalarProd<Scalar>::uninitialize(
  const Ptr<RCP<const LinearOpBase<Scalar> > > &op_out
  )
{
  if (!is_null(op_out)) *op_out = op_;
  op_ = Teuchos::null;
}


// Overridden from ScalarProdBase


template<class Scalar>
bool LinearOpScalarProd<Scalar>::isEuclideanImpl() const
{
  return false;
}


template<class Scalar>
void LinearOpScalarProd<Scalar>::scalarProdsImpl(
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
  const ArrayView<Scalar> &scalarProds_out
  ) const
{
  Teuchos::RCP<MultiVectorBase<Scalar> >
    T = createMembers(Y.range() ,Y.domain()->dim());
  Thyra::apply(*op_, NOTRANS,Y, T.ptr());
  dots(X, *T, scalarProds_out);
}


template<class Scalar>
RCP<const LinearOpBase<Scalar> >
LinearOpScalarProd<Scalar>::getLinearOpImpl() const
{
  return op_;
}


} // end namespace Thyra


#endif  // THYRA_LINEAR_OP_SCALAR_PROD_DEF_HPP
