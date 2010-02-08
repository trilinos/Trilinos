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

#ifndef THYRA_DEFAULT_MULTI_VECTOR_PRODUCT_VECTOR_SPACE_DEF_HPP
#define THYRA_DEFAULT_MULTI_VECTOR_PRODUCT_VECTOR_SPACE_DEF_HPP


#include "Thyra_DefaultMultiVectorProductVectorSpace_decl.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
DefaultMultiVectorProductVectorSpace<Scalar>::DefaultMultiVectorProductVectorSpace()
  : numColumns_(-1)
{}


template<class Scalar>
void DefaultMultiVectorProductVectorSpace<Scalar>::initialize(
  const Teuchos::RCP<const VectorSpaceBase<Scalar> > &space,
  const int numColumns
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(space));
  TEST_FOR_EXCEPT(numColumns <= 0);
#endif
  space_ = space;
  numColumns_ = numColumns;
  defaultProdVecSpc_ = productVectorSpace(space,numColumns);
}


template<class Scalar>
void DefaultMultiVectorProductVectorSpace<Scalar>::uninitialize(
  Teuchos::RCP<const VectorSpaceBase<Scalar> > *space,
  int *numColumns
  )
{
  TEST_FOR_EXCEPT("ToDo: Implement when needed!");
}
  
  
// Overridden from DefaultMultiVectorProductVectorSpace


template<class Scalar>
int DefaultMultiVectorProductVectorSpace<Scalar>::numBlocks() const
{
  return numColumns_;
}


template<class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DefaultMultiVectorProductVectorSpace<Scalar>::getBlock(const int k) const
{
  TEST_FOR_EXCEPT( k < 0 || numColumns_ < k );
  return space_;
}


// Overridden from VectorSpaceBase


template<class Scalar>
Ordinal DefaultMultiVectorProductVectorSpace<Scalar>::dim() const
{
  if (nonnull(space_))
    return numColumns_ * space_->dim();
  return -1;
}


template<class Scalar>
bool DefaultMultiVectorProductVectorSpace<Scalar>::isCompatible(
  const VectorSpaceBase<Scalar>& vecSpc
  ) const
{
  const DefaultMultiVectorProductVectorSpace<Scalar> *multiVecProdVecSpc
    = dynamic_cast<const DefaultMultiVectorProductVectorSpace<Scalar>*>(&vecSpc);
  if ( multiVecProdVecSpc != 0 ) {
    return (
      ( numColumns_ == multiVecProdVecSpc->numColumns_ )
      && 
      ( space_->isCompatible(*multiVecProdVecSpc->space_) )
      );
  }
  return false;
}


template<class Scalar>
Teuchos::RCP< VectorBase<Scalar> >
DefaultMultiVectorProductVectorSpace<Scalar>::createMember() const
{
  return multiVectorProductVector<Scalar>(
    Teuchos::rcp(new DefaultMultiVectorProductVectorSpace<Scalar>(*this))
    );
}


template<class Scalar>
Scalar DefaultMultiVectorProductVectorSpace<Scalar>::scalarProd(
  const VectorBase<Scalar> &x_in,
  const VectorBase<Scalar> &y_in
  ) const
{
  return defaultProdVecSpc_->scalarProd(x_in,y_in);
  // 2007/05/23: rabartl: ToDo: Implement this in a more efficient way using
  // the block scalar product using a single global reduction!
}

template<class Scalar>
void DefaultMultiVectorProductVectorSpace<Scalar>::scalarProdsImpl(
  const MultiVectorBase<Scalar> &X_in,
  const MultiVectorBase<Scalar> &Y_in,
  const ArrayView<Scalar> &scalarProds_out
  ) const
{
  defaultProdVecSpc_->scalarProds(X_in, Y_in, scalarProds_out);
  // 2007/05/23: rabartl: ToDo: Implement this in a more efficient way once
  // you have a specialized multi-vector implementation.
}

template<class Scalar>
bool DefaultMultiVectorProductVectorSpace<Scalar>::hasInCoreView(
  const Range1D& rng_in, const EViewType viewType, const EStrideType strideType
  ) const
{
  return defaultProdVecSpc_->hasInCoreView(rng_in,viewType,strideType);
}

template<class Scalar>
Teuchos::RCP< const VectorSpaceFactoryBase<Scalar> >
DefaultMultiVectorProductVectorSpace<Scalar>::smallVecSpcFcty() const
{
  if (!is_null(space_))
    return space_->smallVecSpcFcty();
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RCP< MultiVectorBase<Scalar> >
DefaultMultiVectorProductVectorSpace<Scalar>::createMembers(int numMembers) const
{
  return VectorSpaceDefaultBase<Scalar>::createMembers(numMembers);
  // 2007/05/23: rabartl: ToDo: Return MultiVectorProductMultiVector object
  // once MultiVectorProductMultiVector is created when needed!
}


template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultMultiVectorProductVectorSpace<Scalar>::clone() const
{
  // Warning! If the client uninitialized this object then changes the
  // constituent vector spaces then we are in trouble!  The client is warned
  // in documentation!
  Teuchos::RCP<DefaultMultiVectorProductVectorSpace<Scalar> >
    mvpvs = Teuchos::rcp(new DefaultMultiVectorProductVectorSpace<Scalar>());
  mvpvs->numColumns_ = numColumns_;
  mvpvs->space_ = space_;
  mvpvs->defaultProdVecSpc_ = defaultProdVecSpc_;
  return mvpvs;
}


// Overridden from Teuchos::Describable

                                                
template<class Scalar>
std::string DefaultMultiVectorProductVectorSpace<Scalar>::description() const
{
  std::ostringstream oss;
  oss
    << Teuchos::Describable::description() << "{"
    << "dim="<<this->dim()
    << ", numBlocks="<<numColumns_
    << "}";
  return oss.str();
}


template<class Scalar>
void DefaultMultiVectorProductVectorSpace<Scalar>::describe(
  Teuchos::FancyOStream &out_arg,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  using Teuchos::RCP;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  RCP<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      *out
        << this->description() << std::endl;
      if (nonnull(space_)) {
        OSTab tab2(out);
        *out
          <<  "Constituent vector space 'space' is the same for all spaces V[0],V[1],,,V[numBlocks-1]:\n";
        tab.incrTab();
        *out << "space = " << Teuchos::describe(*space_,verbLevel);
      }
      break;
    }
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


} // namespace Thyra


#endif // THYRA_DEFAULT_MULTI_VECTOR_PRODUCT_VECTOR_SPACE_DEF_HPP
