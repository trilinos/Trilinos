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

#ifndef THYRA_MPI_VECTOR_SPACE_BASE_HPP
#define THYRA_MPI_VECTOR_SPACE_BASE_HPP

#include "Thyra_MPIVectorSpaceDefaultBaseDecl.hpp"
#include "Thyra_ScalarProdVectorSpaceBase.hpp"
#include "Thyra_DefaultMPIVectorSpaceFactory.hpp"
#ifdef RTOp_USE_MPI
#  include "Thyra_MPIVectorSpaceUtilities.hpp"
#  include "Teuchos_RawMPITraits.hpp"
#endif

namespace Thyra {

/** \brief . */
template<class Scalar>
MPIVectorSpaceDefaultBase<Scalar>::MPIVectorSpaceDefaultBase()
  :mapCode_(-1),defaultLocalOffset_(-1),defaultGlobalDim_(-1),localSubDim_(-1)
{}

// Virtual methods with default implementations

template<class Scalar>
Index MPIVectorSpaceDefaultBase<Scalar>::localOffset() const
{
  return defaultLocalOffset_;
}

template<class Scalar>
Index MPIVectorSpaceDefaultBase<Scalar>::mapCode() const
{
  return mapCode_;
}

// Overridden from VectorSpaceBase

template<class Scalar>
Index MPIVectorSpaceDefaultBase<Scalar>::dim() const
{
  return defaultGlobalDim_;
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceFactoryBase<Scalar> >
MPIVectorSpaceDefaultBase<Scalar>::smallVecSpcFcty() const
{
  return smallVecSpcFcty_;
}

template<class Scalar>
bool MPIVectorSpaceDefaultBase<Scalar>::isCompatible( const VectorSpaceBase<Scalar>& vecSpc ) const
{
  if( this->hasInCoreView() && vecSpc.hasInCoreView() )
    return this->dim() == vecSpc.dim();
  const MPIVectorSpaceDefaultBase<Scalar>
    *mpiVecSpc = dynamic_cast<const MPIVectorSpaceDefaultBase<Scalar>*>(&vecSpc);
  if(mpiVecSpc)
    return mapCode() == mpiVecSpc->mapCode();
  return false;
}

// protected

template<class Scalar>
void MPIVectorSpaceDefaultBase<Scalar>::updateState( const Index globalDim )
{
  localSubDim_ = this->localSubDim(); 
  const MPI_Comm mpiComm  = MPI_COMM_NULL;
  if( localSubDim_ >= 0 ) {
#ifdef RTOp_USE_MPI
    typedef Teuchos::RawMPITraits<Index> IRMT;
    MPI_Comm mpiComm = this->mpiComm();
    int numProc = 1;
    int procRank = 0;
    if( mpiComm != MPI_COMM_NULL ) {
      MPI_Comm_size( mpiComm, &numProc );
      MPI_Comm_rank( mpiComm, &procRank );
    }
    if( numProc > 1 && (localSubDim_ < globalDim || globalDim < 0) ) {
      mapCode_ = MPIVectorSpaceUtilities::computeMapCode(mpiComm,localSubDim_);
      defaultLocalOffset_ = MPIVectorSpaceUtilities::computeLocalOffset(mpiComm,localSubDim_);
      if( globalDim < 1 ) {
        defaultGlobalDim_ = MPIVectorSpaceUtilities::computeGlobalDim(mpiComm,localSubDim_);
      }
      else {
        defaultGlobalDim_ = globalDim;
        // ToDo: Perform global reduction to check that this is correct in debug build
      }
    }
    else {
#endif // RTOp_USE_MPI
      // This is a serial or a locally-replicated parallel
      // vector space.
      mapCode_ = localSubDim_;
      defaultLocalOffset_ = 0;
      defaultGlobalDim_ = localSubDim_;
#ifdef RTOp_USE_MPI
    }
#endif
  }
  else {
    mapCode_  = -1;     // Uninitialized!
    defaultLocalOffset_ = -1;
    defaultGlobalDim_ = -1;
  }
  smallVecSpcFcty_ = Teuchos::rcp(new DefaultMPIVectorSpaceFactory<Scalar>(mpiComm));
}
  
} // end namespace Thyra

#endif // THYRA_MPI_VECTOR_SPACE_BASE_HPP
