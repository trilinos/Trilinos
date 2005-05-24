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

#include "Thyra_ProductVectorSpaceDecl.hpp"
#include "Thyra_ProductVector.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template<class Scalar>
ProductVectorSpace<Scalar>::ProductVectorSpace(
  const int                                                       numBlocks
  ,const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >     vecSpaces[]
  )
{
  initialize(numBlocks,vecSpaces);
}

template<class Scalar>
void ProductVectorSpace<Scalar>::initialize(
  const int                                                       numBlocks
  ,const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >     vecSpaces[]
  )
{
  //
  // Check preconditions and compute cached quantities
  //
  TEST_FOR_EXCEPT( numBlocks < 0 );
  TEST_FOR_EXCEPT( vecSpaces == NULL );
  bool  isInCore = true;
  for( int k = 0; k < numBlocks; ++k ) {
    TEST_FOR_EXCEPTION(
      vecSpaces[k].get() == NULL, std::invalid_argument
      ,"Error, the smart pointer vecSpaces["<<k<<"] can not be NULL!"
      );
    if( !vecSpaces[k]->isInCore() ) isInCore = false;
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
  isInCore_ = isInCore;
}

template<class Scalar>
void ProductVectorSpace<Scalar>::uninitialize(
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
void ProductVectorSpace<Scalar>::getVecSpcPoss(
  Index i, int* kth_vector_space, Index* kth_global_offset
  ) const
{
  // Validate the preconditions
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(
    i < 1 || this->dim() < i, std::out_of_range
    ,"VectorSpaceBlocked::get_vector_space_position(...): Error, i = "
    << i << " is not in range [1,"<<this->dim()<<"]"
    );
#endif
  *kth_vector_space  = 0;
  *kth_global_offset = 0;
  while( *kth_vector_space < numBlocks_ ) {
    const Index off_kp1 = (*vecSpacesOffsets_)[*kth_vector_space+1];
    if( off_kp1 + 1 > i ) {
      *kth_global_offset = (*vecSpacesOffsets_)[*kth_vector_space];
      break;
    }
    ++(*kth_vector_space);
  }
  TEST_FOR_EXCEPT( !(*kth_vector_space < numBlocks_) );
}

// Overridden from ProductVectorSpace

template<class Scalar>
int ProductVectorSpace<Scalar>::numBlocks() const
{
  return numBlocks_;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
ProductVectorSpace<Scalar>::getBlock(const int k) const
{
  TEST_FOR_EXCEPT( k < 0 || numBlocks_ < k );
  return (*vecSpaces_)[k];
}

// Overridden from VectorSpaceBase

template<class Scalar>
Index ProductVectorSpace<Scalar>::dim() const
{
  return dim_;
}

template<class Scalar>
bool ProductVectorSpace<Scalar>::isCompatible( const VectorSpaceBase<Scalar>& vecSpc ) const
{
  // Check for in-core
  if( this->isInCore() && vecSpc.isInCore() && ( this->dim() == vecSpc.dim() ) )
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
ProductVectorSpace<Scalar>::createMember() const
{
  return Teuchos::rcp(new ProductVector<Scalar>(Teuchos::rcp(this,false),NULL));
}

template<class Scalar>
Scalar ProductVectorSpace<Scalar>::scalarProd(
  const VectorBase<Scalar>   &x_in
  ,const VectorBase<Scalar>  &y_in
  ) const
{
  const int numBlocks = this->numBlocks(); 
  const ProductVectorBase<Scalar>
    &x = Teuchos::dyn_cast<const ProductVectorBase<Scalar> >(x_in),
    &y = Teuchos::dyn_cast<const ProductVectorBase<Scalar> >(y_in);
  TEST_FOR_EXCEPT( numBlocks!=x.productSpace()->numBlocks() || numBlocks!=y.productSpace()->numBlocks() );
  Scalar scalarProd = Teuchos::ScalarTraits<Scalar>::zero();
  for( int k = 0; k < numBlocks; ++k )
    scalarProd += (*vecSpaces_)[k]->scalarProd(*x.getBlock(k),*y.getBlock(k));
  return scalarProd;
}

template<class Scalar>
void ProductVectorSpace<Scalar>::scalarProds(
  const MultiVectorBase<Scalar>    &X_in
  ,const MultiVectorBase<Scalar>   &Y_in
  ,Scalar                      scalar_prods[]
  ) const
{
   using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  const int numBlocks = this->numBlocks(); 
  const ProductMultiVectorBase<Scalar>
    &X = Teuchos::dyn_cast<const ProductMultiVectorBase<Scalar> >(X_in),
    &Y = Teuchos::dyn_cast<const ProductMultiVectorBase<Scalar> >(Y_in);
  TEST_FOR_EXCEPT( numBlocks!=X.productSpace()->numBlocks() || numBlocks!=Y.productSpace()->numBlocks() );
  const VectorSpaceBase<Scalar> &domain = *X.domain();
  TEST_FOR_EXCEPT( !domain.isCompatible(*Y.domain()) );
  const Index m = domain.dim();
  Workspace<Scalar> _scalar_prods(wss,m,false);
  std::fill_n( scalar_prods, m, Teuchos::ScalarTraits<Scalar>::zero() );
  for( int k = 0; k < numBlocks; ++k ) {
    (*vecSpaces_)[k]->scalarProds(*X.getBlock(k),*Y.getBlock(k),&_scalar_prods[0]);
    for( int j = 0; j < m; ++j ) scalar_prods[j] += _scalar_prods[j];
  }
}

template<class Scalar>
bool ProductVectorSpace<Scalar>::isInCore() const
{
  return isInCore_;
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceFactoryBase<Scalar> >
ProductVectorSpace<Scalar>::smallVecSpcFcty() const
{
  if( dim_ ) return (*vecSpaces_)[0]->smallVecSpcFcty(); // They should all be compatible!
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr< MultiVectorBase<Scalar> >
ProductVectorSpace<Scalar>::createMembers(int numMembers) const
{
  return VectorSpaceBase<Scalar>::createMembers(numMembers); // ToDo: Specialize for ProductMultiVector when needed!
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
ProductVectorSpace<Scalar>::clone() const
{
  // Warning! If the client uninitialized this object then changes the
  // constituent vector spaces then we are in trouble!  The client is warned
  // in documentation!
  Teuchos::RefCountPtr<ProductVectorSpace<Scalar> >
    pvs = Teuchos::rcp(new ProductVectorSpace<Scalar>());
  pvs->numBlocks_          = numBlocks_;
  pvs->vecSpaces_          = vecSpaces_;
  pvs->vecSpacesOffsets_   = vecSpacesOffsets_;
  pvs->dim_                = dim_;
  pvs->isInCore_           = isInCore_;
  return pvs;
}

} // namespace Thyra

#endif // THYRA_PRODUCT_VECTOR_SPACE_STD_HPP
