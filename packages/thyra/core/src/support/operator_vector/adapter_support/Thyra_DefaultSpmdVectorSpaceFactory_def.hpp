// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_SPMD_VECTOR_SPACE_FACTORY_HPP
#define THYRA_DEFAULT_SPMD_VECTOR_SPACE_FACTORY_HPP

#include "Thyra_DefaultSpmdVectorSpaceFactory_decl.hpp"
#include "Thyra_VectorSpaceFactoryBase.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"

namespace Thyra {

template<class Scalar>
DefaultSpmdVectorSpaceFactory<Scalar>::DefaultSpmdVectorSpaceFactory(
  const Teuchos::RCP<const Teuchos::Comm<Ordinal> > &comm
  )
  :comm_(comm)
{}

template<class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DefaultSpmdVectorSpaceFactory<Scalar>::createVecSpc(int dim) const
{
  return locallyReplicatedDefaultSpmdVectorSpace<Scalar>(comm_, dim);
}

} // end namespace Thyra

#endif // THYRA_DEFAULT_SPMD_VECTOR_SPACE_FACTORY_HPP
