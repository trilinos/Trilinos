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

#ifndef THYRA_VECTOR_SPACE_DEFAULT_BASE_HPP
#define THYRA_VECTOR_SPACE_DEFAULT_BASE_HPP

#include "Thyra_VectorSpaceDefaultBaseDecl.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_SerialVectorSpaceFactoryStd.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_MultiVectorCols.hpp"

namespace Thyra {

// Helper classes

template<class Scalar>
class CopyVectorViewBack {
public:
  CopyVectorViewBack( const VectorBase<Scalar> *v, const RTOpPack::MutableSubVectorT<Scalar>  &raw_v )
    :v_(v), raw_v_(raw_v)
    {}
  ~CopyVectorViewBack()
    {
      RTOpPack::SubVectorT<Scalar> sv;
      v_->getSubVector(Range1D(),&sv);
      RTOpPack::assign_entries( &raw_v_, sv );
      v_->freeSubVector(&sv);
    }
private:
  const VectorBase<Scalar>                   *v_;
  const RTOpPack::MutableSubVectorT<Scalar>  raw_v_;
};

template<class Scalar>
class CopyMultiVectorViewBack {
public:
  CopyMultiVectorViewBack( const MultiVectorBase<Scalar> *mv, const RTOpPack::MutableSubMultiVectorT<Scalar>  &raw_mv )
    :mv_(mv), raw_mv_(raw_mv)
    {}
  ~CopyMultiVectorViewBack()
    {
      RTOpPack::SubMultiVectorT<Scalar> smv;
      mv_->getSubMultiVector(Range1D(),Range1D(),&smv);
      RTOpPack::assign_entries( &raw_mv_, smv );
      mv_->freeSubMultiVector(&smv);
    }
private:
  const MultiVectorBase<Scalar>                       *mv_;
  const RTOpPack::MutableSubMultiVectorT<Scalar>  raw_mv_;
};

// Overridden from VectorSpaceBase

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceFactoryBase<Scalar> >
VectorSpaceDefaultBase<Scalar>::smallVecSpcFcty() const
{
  return Teuchos::rcp(new SerialVectorSpaceFactoryStd<Scalar>());
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> > 
VectorSpaceDefaultBase<Scalar>::createMembers(int numMembers) const
{
  return Teuchos::rcp(new MultiVectorCols<Scalar> (Teuchos::rcp(this,false),this->smallVecSpcFcty()->createVecSpc(numMembers)));
}

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
VectorSpaceDefaultBase<Scalar>::createMemberView( const RTOpPack::MutableSubVectorT<Scalar> &raw_v ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( raw_v.subDim() != this->dim() );
#endif
  // Create a vector
  Teuchos::RefCountPtr<VectorBase<Scalar> > v = this->createMember();
  // Copy initial values in raw_v into vector
  RTOpPack::MutableSubVectorT<Scalar> sv;
  v->getSubVector(Range1D(),&sv);
  RTOpPack::assign_entries( &sv, raw_v );
  v->commitSubVector(&sv);
  // Setup smart pointer to vector to copy view back out just before vector is destroyed
  Teuchos::set_extra_data(
    Teuchos::rcp(new CopyVectorViewBack<Scalar>(&*v,raw_v))
    ,"CopyVectorViewBack"
    ,&v
    ,true
    ,Teuchos::PRE_DESTROY
    );
  return v;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
VectorSpaceDefaultBase<Scalar>::createMemberView( const RTOpPack::SubVectorT<Scalar> &raw_v ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( raw_v.subDim() != this->dim() );
#endif
  // Create a vector
  Teuchos::RefCountPtr<VectorBase<Scalar> > v = this->createMember();
  // Copy initial values in raw_v into vector
  RTOpPack::MutableSubVectorT<Scalar> sv;
  v->getSubVector(Range1D(),&sv);
  RTOpPack::assign_entries( &sv, raw_v );
  v->commitSubVector(&sv);
  return v;
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
VectorSpaceDefaultBase<Scalar>::createMembersView( const RTOpPack::MutableSubMultiVectorT<Scalar> &raw_mv ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( raw_mv.subDim() != this->dim() );
#endif
  // Create a multi-vector
  Teuchos::RefCountPtr< MultiVectorBase<Scalar> > mv = this->createMembers(raw_mv.numSubCols());
  // Copy initial values in raw_mv into multi-vector
  RTOpPack::MutableSubMultiVectorT<Scalar> smv;
  mv->getSubMultiVector(Range1D(),Range1D(),&smv);
  RTOpPack::assign_entries( &smv, raw_mv );
  mv->commitSubMultiVector(&smv);
  // Setup smart pointer to multi-vector to copy view back out just before multi-vector is destroyed
  Teuchos::set_extra_data(
    Teuchos::rcp(new CopyMultiVectorViewBack<Scalar>(&*mv,raw_mv))
    ,"CopyMultiVectorViewBack"
    ,&mv
    ,true
    ,Teuchos::PRE_DESTROY
    );
  return mv;
}

template<class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
VectorSpaceDefaultBase<Scalar>::createMembersView( const RTOpPack::SubMultiVectorT<Scalar> &raw_mv ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( raw_mv.subDim() != this->dim() );
#endif
  // Create a multi-vector
  Teuchos::RefCountPtr< MultiVectorBase<Scalar> > mv = this->createMembers(raw_mv.numSubCols());
  // Copy values in raw_mv into multi-vector
  RTOpPack::MutableSubMultiVectorT<Scalar> smv;
  mv->getSubMultiVector(Range1D(),Range1D(),&smv);
  RTOpPack::assign_entries( &smv, raw_mv );
  mv->commitSubMultiVector(&smv);
  return mv;
}

} // end namespace Thyra

#endif // THYRA_VECTOR_SPACE_DEFAULT_BASE_HPP
