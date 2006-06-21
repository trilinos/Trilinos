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

#ifndef THYRA_MPI_VECTOR_STD_HPP
#define THYRA_MPI_VECTOR_STD_HPP

#include "Thyra_DefaultMPIVectorDecl.hpp"
#include "Thyra_MPIVectorBase.hpp"
#include "Thyra_MPIVectorSpaceDefaultBase.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template<class Scalar>
DefaultMPIVector<Scalar>::DefaultMPIVector()
  :stride_(0)
{}

template<class Scalar>
DefaultMPIVector<Scalar>::DefaultMPIVector(
  const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    &mpiSpace
  ,const Teuchos::RefCountPtr<Scalar>                                  &localValues
  ,const Index                                                         stride
  )
{
  initialize(mpiSpace,localValues,stride);
}

template<class Scalar>
void DefaultMPIVector<Scalar>::initialize(
  const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    &mpiSpace
  ,const Teuchos::RefCountPtr<Scalar>                                  &localValues
  ,const Index                                                         stride
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(mpiSpace.get()==NULL);
  TEST_FOR_EXCEPT(mpiSpace->localSubDim() > 0 && localValues.get()==NULL);
  TEST_FOR_EXCEPT(stride==0);
#endif
  mpiSpace_      = mpiSpace;
  localValues_   = localValues;
  stride_        = stride;
  this->updateMpiSpace();
}

template<class Scalar>
void DefaultMPIVector<Scalar>::uninitialize(
  Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    *mpiSpace
  ,Teuchos::RefCountPtr<Scalar>                                  *localValues
  ,Index                                                         *stride
  )
{
  if(mpiSpace)      *mpiSpace       = mpiSpace_;
  if(localValues)   *localValues    = localValues_;
  if(stride)        *stride         = stride_;

  mpiSpace_       = Teuchos::null;
  localValues_    = Teuchos::null;
  stride_         = 0;

  this->updateMpiSpace();
}

// Overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultMPIVector<Scalar>::description() const
{
  return (std::string("DefaultMPIVector<") + Teuchos::ScalarTraits<Scalar>::name() + std::string(">"));
}

// Overridden from MPIVectorBase

template<class Scalar>
Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >
DefaultMPIVector<Scalar>::mpiSpace() const
{
  return mpiSpace_;
}

template<class Scalar>
void DefaultMPIVector<Scalar>::getLocalData( Scalar** localValues, Index* stride )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( localValues==NULL );
  TEST_FOR_EXCEPT( stride==NULL );
#endif
  *localValues = &*localValues_;
  *stride = stride_;
}

template<class Scalar>
void DefaultMPIVector<Scalar>::commitLocalData( Scalar* localValues )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( localValues!=&*localValues_ );
#endif
  // Nothing to commit!
}

template<class Scalar>
void DefaultMPIVector<Scalar>::getLocalData( const Scalar** localValues, Index* stride ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( localValues==NULL );
  TEST_FOR_EXCEPT( stride==NULL );
#endif
  *localValues = &*localValues_;
  *stride = stride_;
}

template<class Scalar>
void DefaultMPIVector<Scalar>::freeLocalData( const Scalar* localValues ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( localValues!=&*localValues_ );
#endif
  // Nothing to free!
}

} // end namespace Thyra

#endif // THYRA_MPI_VECTOR_STD_HPP
