// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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
  return defaultSpmdVectorSpace<Scalar>(comm_,dim,dim);
}

} // end namespace Thyra

#endif // THYRA_DEFAULT_SPMD_VECTOR_SPACE_FACTORY_HPP
