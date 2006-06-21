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

#ifndef THYRA_PRODUCT_VECTOR_SPACE_STD_HPP
#define THYRA_PRODUCT_VECTOR_SPACE_STD_HPP

#include "Thyra_DefaultProductVectorSpaceDecl.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template<class Scalar>
DefaultProductVectorSpace<Scalar>::DefaultProductVectorSpace(
  const int                                                       numBlocks
  ,const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >     vecSpaces[]
  )
{
  initialize(numBlocks,vecSpaces);
}

template<class Scalar>
void DefaultProductVectorSpace<Scalar>::initialize(
  const int                                                       numBlocks
  ,const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >     vecSpaces[]
  )
{
  //
  // Check preconditions and compute cached quantities
  //
  TEST_FOR_EXCEPT( numBlocks < 0 );
  TEST_FOR_EXCEPT( vecSpaces == NULL );
  bool  hasInCoreView = true;
  for( int k = 0; k < numBlocks; ++k ) {
    TEST_FOR_EXCEPTION(
      vecSpaces[k].get() == NULL, std::invalid_argument
      ,"Error, the smart pointer vecSpaces["<<k<<"] can not be NULL!"
      );
    if( !vecSpaces[k]->hasInCoreView() ) hasInCoreView = false;
  }
  //
  // Setup private data members (should not throw an exception from here)
  //
  numBlocks_ = numBlocks;
  vecSpaces_ = Teuchos::rcp(new vecSpaces_t(numBlocks));
  std::copy( vecSpaces, vecSpaces+numBlocks, &(*vecSpaces_)[0] );
  vecSpacesOffsets_ = Teuchos::rcp(new vecSpacesOffsets_t(numBlocks+1));
  (*vecSpacesOffsets_)[0] = 0;
  dim_ = 0;
  for( int k = 1; k <= numBlocks; ++k ) {
    const Index dim_km1 = vecSpaces[k-1]->dim();
    (*vecSpacesOffsets_)[k] = (*vecSpacesOffsets_)[k-1] + dim_km1;
    dim_ += dim_km1;
  }
  isInCore_ = hasInCoreView;
}

template<class Scalar>
void DefaultProductVectorSpace<Scalar>::uninitialize(
  int                                                             *numBlocks
  ,Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >           vecSpaces[]
  )
{
  vecSpaces_        = Teuchos::null;
  vecSpacesOffsets_ = Teuchos::null;
  dim_              = 0;
  isInCore_         = false;
}

template<class Scalar>
void DefaultProductVectorSpace<Scalar>::getVecSpcPoss(
  Index i, int* kth_vector_space, Index* kth_global_offset
  ) const
{
  // Validate the preconditions
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    !(0 <= i && i < this->dim()), std::out_of_range
    ,"VectorSpaceBlocked::get_vector_space_position(...): Error, i = "
    << i << " is not in range [0,"<<(this->dim()-1)<<"]"
    );
#endif
  *kth_vector_space  = 0;
  *kth_global_offset = 0;
  while( *kth_vector_space < numBlocks_ ) {
    const Index off_kp1 = (*vecSpacesOffsets_)[*kth_vector_space+1];
    if( off_kp1 > i ) {
      *kth_global_offset = (*vecSpacesOffsets_)[*kth_vector_space];
      break;
    }
    ++(*kth_vector_space);
  }
  TEST_FOR_EXCEPT( !(*kth_vector_space < numBlocks_) );
}

// Overridden from DefaultProductVectorSpace

template<class Scalar>
int DefaultProductVectorSpace<Scalar>::numBlocks() const
{
  return numBlocks_;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultProductVectorSpace<Scalar>::getBlock(const int k) const
{
  TEST_FOR_EXCEPT( k < 0 || numBlocks_ < k );
  return (*vecSpaces_)[k];
}

// Overridden from VectorSpaceBase

template<class Scalar>
Index DefaultProductVectorSpace<Scalar>::dim() const
{
  return dim_;
}

template<class Scalar>
bool DefaultProductVectorSpace<Scalar>::isCompatible( const VectorSpaceBase<Scalar>& vecSpc ) const
{
  // Check for in-core
  if( this->hasInCoreView(Range1D(),VIEW_TYPE_DETACHED,STRIDE_TYPE_NONUNIT) && vecSpc.hasInCoreView() && ( this->dim() == vecSpc.dim() ) )
    return true;
  // Check for product vector interface
  const ProductVectorSpaceBase<Scalar> *pvsb = dynamic_cast<const ProductVectorSpaceBase<Scalar>*>(&vecSpc);
  if( !pvsb )
    return false;
  // Validate that constituent vector spaces are compatible
  const int numBlocks = this->numBlocks(); 
  if( numBlocks != pvsb->numBlocks() )
    return false;
  for( int i = 0; i < numBlocks; ++i ) {
    if( !this->getBlock(i)->isCompatible(*pvsb->getBlock(i)) )
      return false;
  }
  return true;
}

template<class Scalar>
Teuchos::RefCountPtr< VectorBase<Scalar> >
DefaultProductVectorSpace<Scalar>::createMember() const
{
  return Teuchos::rcp(
    new DefaultProductVector<Scalar>(Teuchos::rcp(this,false))
    );
}

template<class Scalar>
Scalar DefaultProductVectorSpace<Scalar>::scalarProd(
  const VectorBase<Scalar>   &x_in
  ,const VectorBase<Scalar>  &y_in
  ) const
{
  const int numBlocks = this->numBlocks(); 
  const ProductVectorBase<Scalar>
    &x = Teuchos::dyn_cast<const ProductVectorBase<Scalar> >(x_in),
    &y = Teuchos::dyn_cast<const ProductVectorBase<Scalar> >(y_in);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( numBlocks!=x.productSpace()->numBlocks() || numBlocks!=y.productSpace()->numBlocks() );
#endif
  Scalar scalarProd = Teuchos::ScalarTraits<Scalar>::zero();
  for( int k = 0; k < numBlocks; ++k )
    scalarProd += (*vecSpaces_)[k]->scalarProd(*x.getBlock(k),*y.getBlock(k));
  return scalarProd;
}

template<class Scalar>
void DefaultProductVectorSpace<Scalar>::scalarProds(
  const MultiVectorBase<Scalar>    &X_in
  ,const MultiVectorBase<Scalar>   &Y_in
  ,Scalar                          scalar_prods[]
  ) const
{
  using Teuchos::Workspace;
  const VectorSpaceBase<Scalar> &domain = *X_in.domain();
  const Index m = domain.dim();
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(scalar_prods==NULL);
  TEST_FOR_EXCEPT( !domain.isCompatible(*Y_in.domain()) );
#endif
  if(m==1) {
    scalar_prods[0] = this->scalarProd(*X_in.col(0),*Y_in.col(0));
    return;
    // ToDo: Remove this if(...) block once we have a DefaultProductMultiVector implementation!
  }
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  const int numBlocks = this->numBlocks(); 
  const ProductMultiVectorBase<Scalar>
    &X = Teuchos::dyn_cast<const ProductMultiVectorBase<Scalar> >(X_in),
    &Y = Teuchos::dyn_cast<const ProductMultiVectorBase<Scalar> >(Y_in);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( numBlocks!=X.productSpace()->numBlocks() || numBlocks!=Y.productSpace()->numBlocks() );
#endif
  Workspace<Scalar> _scalar_prods(wss,m,false);
  std::fill_n( scalar_prods, m, Teuchos::ScalarTraits<Scalar>::zero() );
  for( int k = 0; k < numBlocks; ++k ) {
    (*vecSpaces_)[k]->scalarProds(*X.getMultiVectorBlock(k),*Y.getMultiVectorBlock(k),&_scalar_prods[0]);
    for( int j = 0; j < m; ++j ) scalar_prods[j] += _scalar_prods[j];
  }
}

template<class Scalar>
bool DefaultProductVectorSpace<Scalar>::hasInCoreView(const Range1D& rng_in, const EViewType viewType, const EStrideType strideType) const
{
  const Range1D rng = full_range(rng_in,0,dim_-1);
  // First see if rng fits in a single constituent vector
  int    kth_vector_space  = -1;
  Index  kth_global_offset = 0;
  this->getVecSpcPoss(rng.lbound(),&kth_vector_space,&kth_global_offset);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
  if( rng.lbound() + rng.size() <= kth_global_offset + (*vecSpaces_)[kth_vector_space]->dim() ) {
    return (*vecSpaces_)[kth_vector_space]->hasInCoreView(rng_in-kth_global_offset,viewType,strideType);
  }
  // If we get here, rng does not fit in a single constituent vector which
  // also means that numBlocks_ > 1 must also be true!
  //
  // Next, if the client is asking for a direct view then we have to return
  // false since this range spans more than one constituent vector.
  if( viewType == VIEW_TYPE_DIRECT )
    return false;
  // If we get here then hasDirectView==false and therefore we are allowed to
  // create a copy.  Therefore, if all of the constituent vectors are "in
  // core" then we can return true.
  if(isInCore_)
    return true;
  // Finally, loop through all of the constituent vectors spaned by rng and
  // see if they are each in core.
  //
  // Todo: Implement this if you have to!
  //
  // We must give up and return false
  return false;
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceFactoryBase<Scalar> >
DefaultProductVectorSpace<Scalar>::smallVecSpcFcty() const
{
  if( dim_ ) return (*vecSpaces_)[0]->smallVecSpcFcty(); // They should all be compatible!
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr< MultiVectorBase<Scalar> >
DefaultProductVectorSpace<Scalar>::createMembers(int numMembers) const
{
  return VectorSpaceDefaultBase<Scalar>::createMembers(numMembers); // ToDo: Specialize for ProductMultiVector when needed!
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultProductVectorSpace<Scalar>::clone() const
{
  // Warning! If the client uninitialized this object then changes the
  // constituent vector spaces then we are in trouble!  The client is warned
  // in documentation!
  Teuchos::RefCountPtr<DefaultProductVectorSpace<Scalar> >
    pvs = Teuchos::rcp(new DefaultProductVectorSpace<Scalar>());
  pvs->numBlocks_          = numBlocks_;
  pvs->vecSpaces_          = vecSpaces_;
  pvs->vecSpacesOffsets_   = vecSpacesOffsets_;
  pvs->dim_                = dim_;
  pvs->isInCore_           = isInCore_;
  return pvs;
}

} // namespace Thyra

#endif // THYRA_PRODUCT_VECTOR_SPACE_STD_HPP
