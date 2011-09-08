// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_SPMD_VECTOR_SPACE_BASE_DEF_HPP
#define THYRA_SPMD_VECTOR_SPACE_BASE_DEF_HPP

#include "Thyra_SpmdVectorSpaceDefaultBase_decl.hpp"
#include "Thyra_ScalarProdVectorSpaceBase.hpp"
#include "Thyra_DefaultSpmdVectorSpaceFactory.hpp"
#include "Thyra_SpmdVectorSpaceUtilities.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"


namespace Thyra {


template<class Scalar>
SpmdVectorSpaceDefaultBase<Scalar>::SpmdVectorSpaceDefaultBase()
  :mapCode_(-1), defaultLocalOffset_(-1), defaultGlobalDim_(-1), localSubDim_(-1)
{}


// Virtual methods with default implementations


template<class Scalar>
Ordinal SpmdVectorSpaceDefaultBase<Scalar>::localOffset() const
{
  return defaultLocalOffset_;
}


template<class Scalar>
Ordinal SpmdVectorSpaceDefaultBase<Scalar>::mapCode() const
{
  return mapCode_;
}


template<class Scalar>
std::string SpmdVectorSpaceDefaultBase<Scalar>::description() const
{
  using Teuchos::RCP; using Teuchos::Comm; using Teuchos::null;
  std::ostringstream ostr;
  ostr << Teuchos::typeName(*this) << "{";
  ostr << "globalDim="<<this->dim();
  ostr << ",localSubDim="<<this->localSubDim();
  ostr << ",localOffset="<<this->localOffset();
  ostr << ",comm=";
  RCP<const Comm<Ordinal> > comm;
  if ( (comm=this->getComm())!=null ) {
    ostr << comm->description();
  }
  else {
    ostr << "NULL";
  }
  ostr << "}";
  return ostr.str();
}


// Overridden from VectorSpaceBase


template<class Scalar>
Ordinal SpmdVectorSpaceDefaultBase<Scalar>::dim() const
{
  return defaultGlobalDim_;
}


template<class Scalar>
Teuchos::RCP< const VectorSpaceFactoryBase<Scalar> >
SpmdVectorSpaceDefaultBase<Scalar>::smallVecSpcFcty() const
{
  return smallVecSpcFcty_;
}


template<class Scalar>
bool SpmdVectorSpaceDefaultBase<Scalar>::isCompatible(
  const VectorSpaceBase<Scalar>& vecSpc
  ) const
{

  using Teuchos::ptrFromRef;
  using Teuchos::ptr_dynamic_cast;
  
  // Check for exact match of vector space
  const Ptr<const SpmdVectorSpaceBase<Scalar> >
    spmdVecSpc = ptr_dynamic_cast<const SpmdVectorSpaceBase<Scalar> >(ptrFromRef(vecSpc));
  if (nonnull(spmdVecSpc)) {
    return mapCode() == spmdVecSpc->mapCode();
  }

  // Check for in-core views
  if( this->hasInCoreView() && vecSpc.hasInCoreView() && this->dim() == vecSpc.dim() )
    return true;
  // 2009/05/11: rabartl: ToDo: Remove this!
  
  // Check for product vector interface
  const Ptr<const ProductVectorSpaceBase<Scalar> > pvsb =
    ptr_dynamic_cast<const ProductVectorSpaceBase<Scalar> >(ptrFromRef(vecSpc));
  
  if (nonnull(pvsb)) {
    if (pvsb->numBlocks() == 1 ) {
      return pvsb->getBlock(0)->isCompatible(*this);
    }
    else {
      return false;
    }
  }
  
  // If we get here, we are not compatible!
  return false;
  
}


// protected


template<class Scalar>
void SpmdVectorSpaceDefaultBase<Scalar>::updateState( const Ordinal globalDim )
{
  localSubDim_ = this->localSubDim(); 
  const Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = this->getComm();
  if( localSubDim_ >= 0 ) {
    int numProc = 1;
    int procRank = 0;
    if( comm.get() ) {
      numProc = comm->getSize();
      procRank = comm->getRank();
    }
    if( numProc > 1 && (localSubDim_ < globalDim || globalDim < 0) ) {
      mapCode_ = SpmdVectorSpaceUtilities::computeMapCode(*comm,localSubDim_);
      defaultLocalOffset_
        = SpmdVectorSpaceUtilities::computeLocalOffset(*comm,localSubDim_);
      if( globalDim < 1 ) {
        defaultGlobalDim_
          = SpmdVectorSpaceUtilities::computeGlobalDim(*comm,localSubDim_);
      }
      else {
        defaultGlobalDim_ = globalDim;
        // ToDo: Perform global reduction to check that this is correct in
        // debug build
      }
    }
    else {
      // This is a serial or a locally-replicated parallel
      // vector space.
      mapCode_ = localSubDim_;
      defaultLocalOffset_ = 0;
      defaultGlobalDim_ = localSubDim_;
    }
  }
  else {
    mapCode_  = -1;     // Uninitialized!
    defaultLocalOffset_ = -1;
    defaultGlobalDim_ = -1;
  }
  smallVecSpcFcty_ = defaultSpmdVectorSpaceFactory<Scalar>(comm);
}

 
} // end namespace Thyra


#endif // THYRA_SPMD_VECTOR_SPACE_BASE_DEF_HPP
