// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_PRECONDITIONER_FACTORY_BASE_DEF_HPP
#define THYRA_PRECONDITIONER_FACTORY_BASE_DEF_HPP

#include "Thyra_PreconditionerFactoryBase_decl.hpp"
#include "Teuchos_ScalarTraits.hpp"


namespace Thyra {


template<class Scalar>
bool PreconditionerFactoryBase<Scalar>::applySupportsConj(EConj conj) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return ( ST::isComplex ? ( conj==NONCONJ_ELE ) : true );
}


template<class Scalar>
bool
PreconditionerFactoryBase<Scalar>::applyTransposeSupportsConj(EConj /* conj */) const
{
  return false;
}


} // namespace Thyra


#endif // THYRA_PRECONDITIONER_FACTORY_BASE_DEF_HPP
