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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
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
  :mapCode_(-1), defaultLocalOffset_(-1), defaultGlobalDim_(-1),
   localSubDim_(-1), isLocallyReplicated_(false)
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
bool SpmdVectorSpaceDefaultBase<Scalar>::isLocallyReplicated() const
{
  return isLocallyReplicated_;
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
void SpmdVectorSpaceDefaultBase<Scalar>::updateState(const Ordinal globalDim,
  const bool isLocallyReplicated_in)
{
  namespace SVSU = SpmdVectorSpaceUtilities;

  //
  // A) Get the comm, comm info, and local subdim
  //

  localSubDim_ = this->localSubDim(); 

  const Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = this->getComm();

  int numProcs = 1;
  if (nonnull(comm)) {
    numProcs = comm->getSize();
  }

  //
  // B) Determine the type of vector space
  //

  Ordinal sumLocalSubDims = localSubDim_;

  bool isSerialOrLocallyReplicated = false;
  bool isEmpty = false;
  bool isDistributed = false;

  if (isLocallyReplicated_in) {
    // Avoid communication when we know this is locally replicated
    isSerialOrLocallyReplicated = true;
    if (sumLocalSubDims == 0) {
      isEmpty = true;
    }
    TEUCHOS_ASSERT_EQUALITY(localSubDim_, globalDim);
  }
  else {
    if (numProcs > 1) {
      sumLocalSubDims = SVSU::computeGlobalDim(*comm, localSubDim_);
    }
    if (sumLocalSubDims == 0) {
      // This is an uninitialized space (zero on every process)
      isEmpty = true;
    }
    else if (
      numProcs == 1
      ||
      (
        sumLocalSubDims / numProcs == globalDim
        &&
        sumLocalSubDims % numProcs == 0
        )
      )
    {
      // This is a serial or a locally-replicated parallel
      // vector space.
      isSerialOrLocallyReplicated = true;
      //TEUCHOS_TEST_FOR_EXCEPTION(numProcs > 1, std::logic_error,
      //  "Error, you are creating a locally replicated vector space implicitly which"
      //  " is very inefficient.  Please pass in isLocallyReplicated=true to avoid"
      //  " unnecessary global communication!");
      // ToDo: Add runtime option to assert that an implicit VS is not being
      // created which is a performance problem.
    }
    else {
      // This is a regular distributed vector space
      isDistributed = true;
    }
  }

  //
  // C) Set the state of the vector space for the given type
  //

  if (isEmpty) {
    mapCode_  = 0;
    defaultLocalOffset_ = 0;
    defaultGlobalDim_ = 0;
  }
  else if (isSerialOrLocallyReplicated) {
    isLocallyReplicated_ = true;
    mapCode_ = localSubDim_;
    defaultLocalOffset_ = 0;
    defaultGlobalDim_ = localSubDim_;
  }
  else {
    TEUCHOS_ASSERT(isDistributed);
    defaultGlobalDim_ = sumLocalSubDims;
    mapCode_ = SVSU::computeMapCode(*comm, localSubDim_);
    defaultLocalOffset_ = SVSU::computeLocalOffset(*comm, localSubDim_);
  }

  smallVecSpcFcty_ = defaultSpmdVectorSpaceFactory<Scalar>(comm);

}

 
} // end namespace Thyra


#endif // THYRA_SPMD_VECTOR_SPACE_BASE_DEF_HPP
