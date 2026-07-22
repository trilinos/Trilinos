// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_LINEAR_OP_BASE_DEF_HPP
#define THYRA_LINEAR_OP_BASE_DEF_HPP

#include "Thyra_LinearOpBase_decl.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"


namespace Thyra {


// Public interface functions


template<class Scalar>
RCP<const LinearOpBase<Scalar> > 
LinearOpBase<Scalar>::clone() const
{
  return Teuchos::null;
}


}	// end namespace Thyra


// ToDo: You can move this back to the decl file after you have refactored
// apply(...) to not use raw pointers.  Otherwise the Y.ptr() call needs to have
// the definition of MultiVectorBase.


template<class Scalar>
void Thyra::apply(
  const LinearOpBase<Scalar> &M,
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  )
{
  M.apply(M_trans, X, Y, alpha, beta);
}


//
// Explicit instantiation macro
//

#define THYRA_LINEAR_OP_BASE_INSTANT(SCALAR) \
  \
  template class LinearOpBase<SCALAR >; \
  \
  template void apply(  \
    const LinearOpBase<SCALAR > &M,  \
    const EOpTransp M_trans,  \
    const MultiVectorBase<SCALAR > &X,  \
    const Ptr<MultiVectorBase<SCALAR > > &Y,  \
    const SCALAR alpha,  \
    const SCALAR beta  \
    );


#endif // THYRA_LINEAR_OP_BASE_DEF_HPP
