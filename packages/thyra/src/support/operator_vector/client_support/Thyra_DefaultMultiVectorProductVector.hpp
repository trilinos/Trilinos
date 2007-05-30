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

#ifndef THYRA_MULTI_VECTOR_PRODUCT_VECTOR_HPP
#define THYRA_MULTI_VECTOR_PRODUCT_VECTOR_HPP


#include "Thyra_DefaultMultiVectorProductVectorDecl.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template <class Scalar>
DefaultMultiVectorProductVector<Scalar>::DefaultMultiVectorProductVector()
{
  uninitialize();
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::initialize(
  const Teuchos::RefCountPtr<const DefaultMultiVectorProductVectorSpace<Scalar> > &productSpace,
  const Teuchos::RefCountPtr<MultiVectorBase<Scalar> > &multiVec
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(productSpace));
  TEST_FOR_EXCEPT(is_null(multiVec));
  THYRA_ASSERT_VEC_SPACES(
    "DefaultMultiVectorProductVector<Scalar>::initialize(productSpace,multiVec)",
    *multiVec->range(), *productSpace->getBlock(0)
    );
  TEST_FOR_EXCEPT( multiVec->domain()->dim() != productSpace->numBlocks() );
#endif

  numBlocks_ = productSpace->numBlocks();

  productSpace_ = productSpace;

  multiVec_ = multiVec;

}


template <class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getNonconstMultiVector()
{
  return multiVec_.getNonconstObj();
}


template <class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getMultiVector() const
{
  return multiVec_.getConstObj();
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::uninitialize()
{
  numBlocks_ = 0;
  productSpace_ = Teuchos::null;
  multiVec_.uninitialize();
}


// Overridden from Teuchos::Describable

                                                
template<class Scalar>
std::string DefaultMultiVectorProductVector<Scalar>::description() const
{
  std::ostringstream oss;
  oss
    << Teuchos::Describable::description()
    << "{"
    << "dim="<<this->space()->dim()
    << ",numColumns = "<<numBlocks_
    << "}";
  return oss.str();
}

template<class Scalar>
void DefaultMultiVectorProductVector<Scalar>::describe(
  Teuchos::FancyOStream &out_arg,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  using Teuchos::RefCountPtr;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::describe;
  RefCountPtr<FancyOStream> out = rcp(&out_arg,false);
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
        << Teuchos::Describable::description() << "{"
        << "dim=" << this->space()->dim()
        << "}\n";
      OSTab tab(out);
      *out <<  "multiVec = " << describe(*multiVec_.getConstObj(),verbLevel);
      break;
    }
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


// Overridden from ProductVectorBase


template <class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getNonconstVectorBlock(const int k)
{
  TEST_FOR_EXCEPT( k < 0 || numBlocks_-1 < k);
  return multiVec_.getNonconstObj()->col(k);
}


template <class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getVectorBlock(const int k) const
{
  TEST_FOR_EXCEPT( k < 0 || numBlocks_-1 < k);
  return multiVec_.getConstObj()->col(k);
}


// Overridden from ProductMultiVectorBase


template <class Scalar>
Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::productSpace() const
{
  return productSpace_;
}


template <class Scalar>
bool DefaultMultiVectorProductVector<Scalar>::blockIsConst(const int k) const
{
  TEST_FOR_EXCEPT( k < 0 || numBlocks_-1 < k);
  return multiVec_.isConst();
}


template <class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getNonconstMultiVectorBlock(const int k)
{
  return getNonconstVectorBlock(k);
}


template <class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getMultiVectorBlock(const int k) const
{
  return getVectorBlock(k);
}


// Overridden from VectorBase


template <class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::space() const
{
  return productSpace_;
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::applyOp(
  const RTOpPack::RTOpT<Scalar>    &op
  ,const int                       num_vecs
  ,const VectorBase<Scalar>*const  vecs[]
  ,const int                       num_targ_vecs
  ,VectorBase<Scalar>*const        targ_vecs[]
  ,RTOpPack::ReductTarget          *reduct_obj
  ,const Index                     first_ele_offset_in
  ,const Index                     sub_dim_in
  ,const Index                     global_offset_in
  ) const
{
  this->getDefaultProductVector()->applyOp(
    op, num_vecs, vecs, num_targ_vecs, targ_vecs,
    reduct_obj,
    first_ele_offset_in, sub_dim_in, global_offset_in
    );
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::acquireDetachedView(
  const Range1D& rng_in, RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
  this->getDefaultProductVector()->acquireDetachedView(rng_in,sub_vec);
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::releaseDetachedView(
  RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
  this->getDefaultProductVector()->releaseDetachedView(sub_vec);
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::acquireDetachedView(
  const Range1D& rng_in, RTOpPack::SubVectorView<Scalar>* sub_vec
  )
{
  TEST_FOR_EXCEPT("ToDo: Implement DefaultMultiVectorProductVector<Scalar>::acquireDetachedView(...)!");
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::commitDetachedView(
  RTOpPack::SubVectorView<Scalar>* sub_vec
  )
{
  TEST_FOR_EXCEPT("ToDo: Implement DefaultMultiVectorProductVector<Scalar>::commitDetachedView(...)!");
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::setSubVector(
  const RTOpPack::SparseSubVectorT<Scalar>& sub_vec
  )
{
  TEST_FOR_EXCEPT("ToDo: Implement DefaultMultiVectorProductVector<Scalar>::setSubVector(...)!");
}


// private


template <class Scalar>
Teuchos::RefCountPtr<const DefaultProductVector<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getDefaultProductVector() const
{

  // This function exists since in general we can not create views of a column
  // vectors and expect the changes to be mirrored in the mulit-vector
  // automatically.  Later, we might be able to change this once we have a
  // Thyra::MultiVectorBase::hasDirectColumnVectorView() function and it
  // returns true.  Until then, this is the safe way to do this ...

  using Teuchos::Array; using Teuchos::RefCountPtr;

  Array<RefCountPtr<const VectorBase<Scalar> > > vecArray;
  for ( int k = 0; k < numBlocks_; ++k) {
    vecArray.push_back(multiVec_.getConstObj()->col(k));
  }

  return Thyra::defaultProductVector<Scalar>(
    productSpace_->getDefaultProductVectorSpace(),
    &vecArray[0]
    );

}


} // namespace Thyra


#endif // THYRA_MULTI_VECTOR_PRODUCT_VECTOR_HPP
