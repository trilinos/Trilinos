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

#include "Thyra_MPIVectorSpaceBaseDecl.hpp"
#include "Thyra_ScalarProdVectorSpaceBase.hpp"
#include "Thyra_MPIVectorSpaceFactoryStd.hpp"
#ifdef RTOp_USE_MPI
#  include "Teuchos_RawMPITraits.hpp"
#endif

namespace Thyra {

/** \brief . */
template<class Scalar>
MPIVectorSpaceBase<Scalar>::MPIVectorSpaceBase()
  :mapCode_(-1),isInCore_(false),defaultLocalOffset_(-1),defaultGlobalDim_(-1)
{}

// Virtual methods with default implementations

template<class Scalar>
Index MPIVectorSpaceBase<Scalar>::localOffset() const
{
  return defaultLocalOffset_;
}

template<class Scalar>
Index MPIVectorSpaceBase<Scalar>::mapCode() const
{
  return mapCode_;
}

// Overridden from VectorSpaceBase

template<class Scalar>
Index MPIVectorSpaceBase<Scalar>::dim() const
{
  return defaultGlobalDim_;
}

template<class Scalar>
bool MPIVectorSpaceBase<Scalar>::isInCore() const
{
  return isInCore_;
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceFactoryBase<Scalar> >
MPIVectorSpaceBase<Scalar>::smallVecSpcFcty() const
{
  return smallVecSpcFcty_;
}

template<class Scalar>
bool MPIVectorSpaceBase<Scalar>::isCompatible( const VectorSpaceBase<Scalar>& vecSpc ) const
{
  if( isInCore() && vecSpc.isInCore() )
    return this->dim() == vecSpc.dim();
  const MPIVectorSpaceBase<Scalar>
    *mpiVecSpc = dynamic_cast<const MPIVectorSpaceBase<Scalar>*>(&vecSpc);
  if(mpiVecSpc)
    return mapCode() == mpiVecSpc->mapCode();
  return false;
}

// protected

template<class Scalar>
void MPIVectorSpaceBase<Scalar>::updateState( const Index globalDim )
{
  const Index localSubDim = this->localSubDim(); 
  const MPI_Comm mpiComm  = MPI_COMM_NULL;
  if( localSubDim >= 0 ) {
#ifdef RTOp_USE_MPI
    typedef Teuchos::RawMPITraits<Index> IRMT;
    MPI_Comm mpiComm = this->mpiComm();
    int numProc = 1;
    int procRank = 0;
    if( mpiComm != MPI_COMM_NULL ) {
      MPI_Comm_size( mpiComm, &numProc );
      MPI_Comm_rank( mpiComm, &procRank );
    }
    if( numProc > 1 && (localSubDim < globalDim || globalDim < 0) ) {
      //
      // Here we will make a map code out of just the local
      // sub-dimension on each processor.  If each processor
      // has the same number of local elements, then the maps
      // will be the same and this is all you need for
      // RTOp compatibility unless the operations are not
      // coordinate invariant.  I will work on this issue
      // if it becomes a problem.
      //
      Index localCode = localSubDim % (procRank+1) + localSubDim;
      MPI_Allreduce(
        &localCode                              // sendbuf
        ,&mapCode_                              // recvbuf
        ,IRMT::adjustCount(1)                   // count
        ,IRMT::type()                           // datatype
        ,IRMT::sumOp()                          // op
        ,mpiComm                                // comm
        );
      // Set the default localOffset automatically
      Index localOffset = localSubDim;
      MPI_Scan(
        &localOffset                            // sendbuf
        ,&defaultLocalOffset_                   // recvbuf
        ,IRMT::adjustCount(1)                   // count
        ,IRMT::type()                           // datatype
        ,IRMT::sumOp()                          // op
        ,mpiComm                                // comm
        );
      defaultLocalOffset_ -= localSubDim;
      //int procRank; MPI_Comm_rank( mpiComm, &procRank );
      //std::cout << "\nMPIVectorSpaceBase<Scalar>::updateState(): procRank = " << procRank << ", defaultLocalOffset = " << defaultLocalOffset_ << std::endl;
      // 
      isInCore_ = false;  // This is not an inCore vector
      // Determine the global dimension
      if( globalDim < 0 ) {
        MPI_Allreduce(
          (void*)&localSubDim                     // sendbuf
          ,&defaultGlobalDim_                     // recvbuf
          ,IRMT::adjustCount(1)                   // count
          ,IRMT::type()                           // datatype
          ,IRMT::sumOp()                          // op
          ,mpiComm                                // comm
          );
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
      mapCode_ = localSubDim;
      isInCore_ = true;
      defaultLocalOffset_ = 0;
      defaultGlobalDim_ = localSubDim;
#ifdef RTOp_USE_MPI
    }
#endif
  }
  else {
    mapCode_  = -1;     // Uninitialized!
    isInCore_ = false;
    defaultLocalOffset_ = -1;
    defaultGlobalDim_ = -1;
  }
  smallVecSpcFcty_ = Teuchos::rcp(new MPIVectorSpaceFactoryStd<Scalar>(mpiComm));
}
  
} // end namespace Thyra

#endif // THYRA_MPI_VECTOR_SPACE_BASE_HPP
