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

#ifndef THYRA_VECTOR_SPACE_DEFAULT_BASE_DEF_HPP
#define THYRA_VECTOR_SPACE_DEFAULT_BASE_DEF_HPP

#include "Thyra_VectorSpaceDefaultBase_decl.hpp"
#include "Thyra_VectorSpaceFactoryBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DefaultColumnwiseMultiVector.hpp"


namespace Thyra {


// Helper classes


template<class Scalar>
class CopyVectorViewBack {
public:
  CopyVectorViewBack( const VectorBase<Scalar> *v, const RTOpPack::SubVectorView<Scalar>  &raw_v )
    :v_(v), raw_v_(raw_v)
    {}
  ~CopyVectorViewBack()
    {
      RTOpPack::ConstSubVectorView<Scalar> sv;
      v_->acquireDetachedView(Range1D(),&sv);
      RTOpPack::assign_entries<Scalar>( Teuchos::outArg(raw_v_), sv );
      v_->releaseDetachedView(&sv);
    }
private:
  const VectorBase<Scalar>                   *v_;
  const RTOpPack::SubVectorView<Scalar>  raw_v_;
};


template<class Scalar>
class CopyMultiVectorViewBack {
public:
  CopyMultiVectorViewBack( const MultiVectorBase<Scalar> *mv, const RTOpPack::SubMultiVectorView<Scalar>  &raw_mv )
    :mv_(mv), raw_mv_(raw_mv)
    {}
  ~CopyMultiVectorViewBack()
    {
      RTOpPack::ConstSubMultiVectorView<Scalar> smv;
      mv_->acquireDetachedView(Range1D(),Range1D(),&smv);
      RTOpPack::assign_entries<Scalar>( Teuchos::outArg(raw_mv_), smv );
      mv_->releaseDetachedView(&smv);
    }
private:
  const MultiVectorBase<Scalar>                       *mv_;
  const RTOpPack::SubMultiVectorView<Scalar>  raw_mv_;
};


// Overridden protected functions from VectorSpaceBase


template<class Scalar>
Teuchos::RCP<MultiVectorBase<Scalar> > 
VectorSpaceDefaultBase<Scalar>::createMembers(int numMembers) const
{
  return Teuchos::rcp(
    new DefaultColumnwiseMultiVector<Scalar>(
      Teuchos::rcp(this,false),
      this->smallVecSpcFcty()->createVecSpc(numMembers)));
  // ToDo: Use the "self object-reference" idiom ot fix this!
}


template<class Scalar>
Teuchos::RCP<VectorBase<Scalar> >
VectorSpaceDefaultBase<Scalar>::createMemberView(
  const RTOpPack::SubVectorView<Scalar> &raw_v ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( raw_v.subDim() != this->dim() );
#endif
  // Create a vector
  RCP<VectorBase<Scalar> > v = this->createMember();
  // Copy initial values in raw_v into vector
  RTOpPack::SubVectorView<Scalar> sv;
  v->acquireDetachedView(Range1D(),&sv);
  RTOpPack::assign_entries<Scalar>(
    Teuchos::ptr_implicit_cast<const RTOpPack::SubVectorView<Scalar> >(Teuchos::outArg(sv)),
    RTOpPack::ConstSubVectorView<Scalar>(raw_v) );
  v->commitDetachedView(&sv);
  // Setup smart pointer to vector to copy view back out just before vector is destroyed
  Teuchos::set_extra_data(
    Teuchos::rcp(new CopyVectorViewBack<Scalar>(&*v,raw_v)),
    "CopyVectorViewBack",
    Teuchos::outArg(v),
    Teuchos::PRE_DESTROY
    );
  return v;
}


template<class Scalar>
Teuchos::RCP<const VectorBase<Scalar> >
VectorSpaceDefaultBase<Scalar>::createMemberView( const RTOpPack::ConstSubVectorView<Scalar> &raw_v ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( raw_v.subDim() != this->dim() );
#endif
  // Create a vector
  RCP<VectorBase<Scalar> > v = this->createMember();
  // Copy initial values in raw_v into vector
  RTOpPack::SubVectorView<Scalar> sv;
  v->acquireDetachedView(Range1D(),&sv);
  RTOpPack::assign_entries<Scalar>(
    Teuchos::ptr_implicit_cast<const RTOpPack::SubVectorView<Scalar> >(Teuchos::outArg(sv)),
    raw_v );
  v->commitDetachedView(&sv);
  return v;
}


template<class Scalar>
Teuchos::RCP<MultiVectorBase<Scalar> >
VectorSpaceDefaultBase<Scalar>::createMembersView(
  const RTOpPack::SubMultiVectorView<Scalar> &raw_mv ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( raw_mv.subDim() != this->dim() );
#endif
  // Create a multi-vector
  RCP< MultiVectorBase<Scalar> > mv = this->createMembers(raw_mv.numSubCols());
  // Copy initial values in raw_mv into multi-vector
  RTOpPack::SubMultiVectorView<Scalar> smv;
  mv->acquireDetachedView(Range1D(),Range1D(),&smv);
  RTOpPack::assign_entries<Scalar>(
    Ptr<const RTOpPack::SubMultiVectorView<Scalar> >(Teuchos::outArg(smv)),
    raw_mv
    );
  mv->commitDetachedView(&smv);
  // Setup smart pointer to multi-vector to copy view back out just before multi-vector is destroyed
  Teuchos::set_extra_data(
    Teuchos::rcp(new CopyMultiVectorViewBack<Scalar>(&*mv,raw_mv)),
    "CopyMultiVectorViewBack",
    Teuchos::outArg(mv),
    Teuchos::PRE_DESTROY
    );
  return mv;
}


template<class Scalar>
Teuchos::RCP<const MultiVectorBase<Scalar> >
VectorSpaceDefaultBase<Scalar>::createMembersView(
  const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( raw_mv.subDim() != this->dim() );
#endif
  // Create a multi-vector
  RCP< MultiVectorBase<Scalar> > mv =
    this->createMembers(raw_mv.numSubCols());
  // Copy values in raw_mv into multi-vector
  RTOpPack::SubMultiVectorView<Scalar> smv;
  mv->acquireDetachedView(Range1D(),Range1D(),&smv);
  RTOpPack::assign_entries<Scalar>(
    Ptr<const RTOpPack::SubMultiVectorView<Scalar> >(Teuchos::outArg(smv)),
    raw_mv );
  mv->commitDetachedView(&smv);
  return mv;
}


} // end namespace Thyra


#endif // THYRA_VECTOR_SPACE_DEFAULT_BASE_DEF_HPP
