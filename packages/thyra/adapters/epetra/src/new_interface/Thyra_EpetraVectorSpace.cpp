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

#include "Thyra_EpetraVectorSpace.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraVector.hpp"
#include "Thyra_EpetraMultiVector.hpp"

#include "Thyra_EuclideanScalarProd.hpp"

#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

namespace Thyra {

RCP<EpetraVectorSpace>
EpetraVectorSpace::create()
{
  const RCP<EpetraVectorSpace> vs(new EpetraVectorSpace);
  vs->weakSelfPtr_ = vs.create_weak();
  return vs;
}

void EpetraVectorSpace::initialize(const RCP<const Epetra_Map> &epetraMap)
{
  comm_ = convertEpetraToThyraComm(epetraMap->Comm());
  epetraMap_ = epetraMap;
  const bool isLongLong = epetraMap->GlobalIndicesLongLong();
  if (isLongLong) {
    this->updateState(epetraMap->NumGlobalElements64(),
                      !epetraMap->DistributedGlobal());
  } else {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    this->updateState(epetraMap->NumGlobalElements(),
                      !epetraMap->DistributedGlobal());
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "Error! The indices in the Epetra_Map are not long long, "
                               "and yet 32bits indices are not enabled.\n");
#endif
  }
  this->setScalarProd(Teuchos::rcp( new EuclideanScalarProd<double>() ));
}


// Overridden from VectorSpace


RCP<VectorBase<double> >
EpetraVectorSpace::createMember() const
{
  return epetraVector(
    weakSelfPtr_.create_strong().getConst(),
    Teuchos::rcp(
      new Epetra_Vector(*epetraMap_, false)
      )
    );
}


RCP< MultiVectorBase<double> >
EpetraVectorSpace::createMembers(int numMembers) const
{
  return epetraMultiVector(
    weakSelfPtr_.create_strong().getConst(),
    epetraVectorSpace(Teuchos::rcp( new Epetra_LocalMap(numMembers, 0, epetraMap_->Comm()) )),
    Teuchos::rcp( new Epetra_MultiVector(*epetraMap_, numMembers, false) )
    );
}


bool EpetraVectorSpace::hasInCoreView(
  const Range1D& rng_in, const EViewType viewType, const EStrideType strideType
  ) const
{
  const Range1D rng = full_range(rng_in,0,this->dim()-1);
  const Ordinal l_localOffset = this->localOffset();

  const Ordinal myLocalSubDim = epetraMap_.is_null () ?
    static_cast<Ordinal> (0) : epetraMap_->NumMyElements();

  return ( l_localOffset<=rng.lbound() && rng.ubound()<l_localOffset+myLocalSubDim );
}


RCP< const VectorSpaceBase<double> >
EpetraVectorSpace::clone() const
{
  return epetraVectorSpace(epetraMap_);
}

RCP<const Epetra_Map>
EpetraVectorSpace::getEpetraMap() const
{
  return epetraMap_;
}

// Overridden from SpmdVectorSpaceDefaultBase


RCP<const Teuchos::Comm<Ordinal> >
EpetraVectorSpace::getComm() const
{
  return comm_;
}


Ordinal EpetraVectorSpace::localSubDim() const
{
  return epetraMap_.is_null () ? static_cast<Ordinal> (0) :
    static_cast<Ordinal> (epetraMap_->NumMyElements());
}

// private


EpetraVectorSpace::EpetraVectorSpace()
{
  // The base classes should automatically default initialize to a safe
  // uninitialized state.
}

} // end namespace Thyra
