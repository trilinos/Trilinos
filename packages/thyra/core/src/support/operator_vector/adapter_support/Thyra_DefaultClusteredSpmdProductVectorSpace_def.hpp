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

#ifndef THYRA_DEFAULT_CLUSTERED_SPMD_PRODUCT_VECTOR_SPACE_HPP
#define THYRA_DEFAULT_CLUSTERED_SPMD_PRODUCT_VECTOR_SPACE_HPP

#include "Thyra_DefaultClusteredSpmdProductVectorSpace_decl.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_DefaultClusteredSpmdProductVector.hpp"
#include "Thyra_VectorSpaceDefaultBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_SpmdVectorSpaceUtilities.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_CommHelpers.hpp"

namespace Thyra {

// Constructors/Intializers/Accessors

template<class Scalar>
DefaultClusteredSpmdProductVectorSpace<Scalar>::DefaultClusteredSpmdProductVectorSpace()
  :isEuclidean_(false),globalDim_(0),clusterSubDim_(-1),clusterOffset_(-1)
{}

template<class Scalar>
DefaultClusteredSpmdProductVectorSpace<Scalar>::DefaultClusteredSpmdProductVectorSpace(
  const Teuchos::RCP<const Teuchos::Comm<Ordinal> >          &intraClusterComm_in
  ,const int                                                       clusterRootRank_in
  ,const Teuchos::RCP<const Teuchos::Comm<Ordinal> >         &interClusterComm_in
  ,const int                                                       numBlocks_in
  ,const Teuchos::RCP<const VectorSpaceBase<Scalar> >      vecSpaces[]
  )
{
  initialize(intraClusterComm_in,clusterRootRank_in,interClusterComm_in,numBlocks_in,vecSpaces);
}

template<class Scalar>
void DefaultClusteredSpmdProductVectorSpace<Scalar>::initialize(
  const Teuchos::RCP<const Teuchos::Comm<Ordinal> >          &intraClusterComm_in
  ,const int                                                       clusterRootRank_in
  ,const Teuchos::RCP<const Teuchos::Comm<Ordinal> >         &interClusterComm_in
  ,const int                                                       numBlocks_in
  ,const Teuchos::RCP<const VectorSpaceBase<Scalar> >      vecSpaces[]
  )
{
  // Set state
  intraClusterComm_ = intraClusterComm_in.assert_not_null();
  clusterRootRank_ = clusterRootRank_in;
  interClusterComm_ = interClusterComm_in; // This can be NULL!
  vecSpaces_.resize(numBlocks_in);
  isEuclidean_ = true;
  Ordinal l_clusterSubDim = 0;
  for( int k = 0; k < numBlocks_in; ++k ) {
    l_clusterSubDim += vecSpaces[k]->dim();
    if(!vecSpaces[k]->isEuclidean())
      isEuclidean_ = false;
    vecSpaces_[k] = vecSpaces[k];
  }
  // We must compute the offsets between clusters and the global dimension by
  // only involving the root process in each cluster.
  if(interClusterComm_.get()) {
    clusterOffset_ = SpmdVectorSpaceUtilities::computeLocalOffset(
      *interClusterComm_,l_clusterSubDim
      );
    globalDim_ = SpmdVectorSpaceUtilities::computeGlobalDim(
      *interClusterComm_,l_clusterSubDim
    );
  }
  // Here must then broadcast the values to all processes within each cluster.
  {
    const Ordinal num = 2;
    Ordinal buff[num] = { clusterOffset_, globalDim_ };
    Teuchos::broadcast<Ordinal>(*intraClusterComm_, clusterRootRank_, num, &buff[0]);
    clusterOffset_ = buff[0];
    globalDim_     = buff[1];

  }
  //
  clusterSubDim_ = l_clusterSubDim;
  // ToDo: Do a global communication across all clusters to see if all vector
  // spaces are all Euclidean.  It is unlikely to be the case where all of the
  // clusters do not have the same vector spaces so I do not think this will
  // even come up.  But just in case, we should keep this in mind!
}

// Overridden form Teuchos::Describable

template<class Scalar>
std::string DefaultClusteredSpmdProductVectorSpace<Scalar>::description() const
{
  std::ostringstream oss;
  oss << "DefaultClusteredSpmdProductVectorSpace{";
  oss << "numBlocks="<<vecSpaces_.size();
  oss << ",globalDim="<<globalDim_;
  oss << ",clusterOffset="<<clusterOffset_;
  oss << "}";
  return oss.str();
}

// Public overridden from VectorSpaceBase

template<class Scalar>
Ordinal DefaultClusteredSpmdProductVectorSpace<Scalar>::dim() const
{
  return globalDim_;
}

template<class Scalar>
bool DefaultClusteredSpmdProductVectorSpace<Scalar>::isCompatible(
  const VectorSpaceBase<Scalar>& vecSpc
  ) const
{
  if( &vecSpc == this )
    return true;
  // For now, I will just do the dynamic cast but in the future, we could get
  // more sophisticated.
  TEST_FOR_EXCEPT(true);
  return false;
}

template<class Scalar>
Teuchos::RCP< const VectorSpaceFactoryBase<Scalar> >
DefaultClusteredSpmdProductVectorSpace<Scalar>::smallVecSpcFcty() const
{
  if(!vecSpaces_.size())
    return Teuchos::null;
  return vecSpaces_[0]->smallVecSpcFcty();
}

template<class Scalar>
Scalar DefaultClusteredSpmdProductVectorSpace<Scalar>::scalarProd(
  const VectorBase<Scalar>& x, const VectorBase<Scalar>& y
  ) const
{
  Scalar scalarProds_out[1];
  this->scalarProds(x, y, &scalarProds_out[0]);
  return scalarProds_out[0];
}

template<class Scalar>
void DefaultClusteredSpmdProductVectorSpace<Scalar>::scalarProdsImpl(
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
  const ArrayView<Scalar> &scalarProds_out ) const
{
  TEST_FOR_EXCEPTION(
    !isEuclidean_, std::logic_error
    ,"Error, have not implemented support for none Euclidean scalar products yet!"
    );
  return dots(X, Y, scalarProds_out);
  // ToDo:
  // * Create DefaultClusteredSpmdProductMultiVector subclass
  // * Cast X and Y this type
  // * Accumulate the scalar products across all of the blocks in this cluster
  // * Accumulate the full scalar products across all of the clusters
  //   using *interClusterComm
  // * Broadcast the final scalar products to all of the processes in
  //   a cluster using *intraClusterComm
}

template<class Scalar>
bool DefaultClusteredSpmdProductVectorSpace<Scalar>::isEuclidean() const
{
  return isEuclidean_;
}

template<class Scalar>
bool DefaultClusteredSpmdProductVectorSpace<Scalar>::hasInCoreView(
  const Range1D& rng, const EViewType viewType, const EStrideType strideType
  ) const
{
  return false; // ToDo: Figure this out for real!
}

template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultClusteredSpmdProductVectorSpace<Scalar>::clone() const
{
  return Teuchos::rcp(new DefaultClusteredSpmdProductVectorSpace<Scalar>(*this));
}

// Protected overridden from ProductVectorSpaceBase

template<class Scalar>
int DefaultClusteredSpmdProductVectorSpace<Scalar>::numBlocks() const
{
  return vecSpaces_.size();
}

template<class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DefaultClusteredSpmdProductVectorSpace<Scalar>::getBlock(const int k) const
{
  using Teuchos::implicit_cast;
  TEST_FOR_EXCEPT( !( 0 <= k && k < implicit_cast<int>(vecSpaces_.size()) ) );
  return vecSpaces_[k];
} 

// Protected overridden from VectorSpaceBase

template<class Scalar>
Teuchos::RCP<VectorBase<Scalar> >
DefaultClusteredSpmdProductVectorSpace<Scalar>::createMember() const
{
  using Teuchos::rcp;
  return rcp(new DefaultClusteredSpmdProductVector<Scalar>(rcp(this,false),NULL));
}

template<class Scalar>
Teuchos::RCP<MultiVectorBase<Scalar> >
DefaultClusteredSpmdProductVectorSpace<Scalar>::createMembers(int numMembers) const
{
  return VectorSpaceDefaultBase<Scalar>::createMembers(numMembers);
  // ToDo: Provide an optimized multi-vector implementation when needed!
}

} // end namespace Thyra

#endif // THYRA_DEFAULT_CLUSTERED_SPMD_PRODUCT_VECTOR_SPACE_HPP
