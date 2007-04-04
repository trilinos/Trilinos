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

#ifndef THYRA_VECTOR_SPACE_BASE_HPP
#define THYRA_VECTOR_SPACE_BASE_HPP

#include "Thyra_VectorSpaceBaseDecl.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

namespace Thyra {

// Virtual functions with default implementations

template<class Scalar>
bool VectorSpaceBase<Scalar>::isEuclidean() const
{
  return false;
}

template<class Scalar>
bool VectorSpaceBase<Scalar>::hasInCoreView(const Range1D& rng, const EViewType viewType, const EStrideType strideType) const
{
  return false;
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
VectorSpaceBase<Scalar>::clone() const
{
  return Teuchos::null;
}

} // end namespace Thyra

// Nonmember functions

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >
Thyra::makeHaveOwnership( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs )
{
  if( vs.has_ownership() ) return vs;
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > _vs = vs->clone();
  TEST_FOR_EXCEPTION(
    _vs.get() == NULL, std::logic_error
    ,"Thyra::makeHaveOwnership(vs): Error, the concrete VectorSpaceBase object identified as \'"
    << vs->description() << "\' does not support the clone() function!"
    );
  return _vs;
}

template<class Scalar>
Teuchos::RefCountPtr< Thyra::VectorBase<Scalar> >
Thyra::createMember(
  const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs,
  const std::string &label
  )
{
  Teuchos::RefCountPtr<VectorBase<Scalar> > v = vs->createMember();
#ifdef THYRA_INITIALIZE_VECS_MULTIVECS_WITH_NANS
  assign(&*v,Teuchos::ScalarTraits<Scalar>::nan());
#endif  
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase", &v );
  if(label.length()) v->setObjectLabel(label);
  return v;
}
  
template<class Scalar>
Teuchos::RefCountPtr< Thyra::VectorBase<Scalar> >
Thyra::createMember(
  const VectorSpaceBase<Scalar> &vs, const std::string &label
  )
{
  return createMember(Teuchos::rcp(&vs,false),label);
}

template<class Scalar>
Teuchos::RefCountPtr< Thyra::MultiVectorBase<Scalar> >
Thyra::createMembers(
  const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs,
  int numMembers,  const std::string &label
  )
{
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
    mv = vs->createMembers(numMembers);
#ifdef THYRA_INITIALIZE_VECS_MULTIVECS_WITH_NANS
  assign(&*mv,Teuchos::ScalarTraits<Scalar>::nan());
#endif  
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase", &mv );
  if(label.length()) mv->setObjectLabel(label);
  return mv;
}

template<class Scalar>
Teuchos::RefCountPtr< Thyra::MultiVectorBase<Scalar> >
Thyra::createMembers(
  const VectorSpaceBase<Scalar> &vs, int numMembers,
  const std::string &label
  )
{
  return createMembers(Teuchos::rcp(&vs,false),numMembers,label);
}

template<class Scalar>
Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
Thyra::createMemberView(
  const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs,
  const RTOpPack::SubVectorView<Scalar> &raw_v,
  const std::string &label
  )
{
  Teuchos::RefCountPtr<VectorBase<Scalar> >
    v = vs->createMemberView(raw_v);
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase", &v );
  if(label.length()) v->setObjectLabel(label);
  return v;
}

template<class Scalar>
Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
Thyra::createMemberView(
  const VectorSpaceBase<Scalar> &vs,
  const RTOpPack::SubVectorView<Scalar> &raw_v,
  const std::string &label
  )
{
  return createMemberView(Teuchos::rcp(&vs,false),raw_v,label);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> >
Thyra::createMemberView(
  const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs,
  const RTOpPack::ConstSubVectorView<Scalar> &raw_v,
  const std::string &label
  )
{
  Teuchos::RefCountPtr<const VectorBase<Scalar> >
    v = vs->createMemberView(raw_v);
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase", &v );
  if(label.length())
    Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(v)->setObjectLabel(label);
  return v;
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> >
Thyra::createMemberView(
  const VectorSpaceBase<Scalar> &vs,
  const RTOpPack::ConstSubVectorView<Scalar> &raw_v,
  const std::string &label
  )
{
  return createMemberView(Teuchos::rcp(&vs,false),raw_v,label);
}

template<class Scalar>
Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >
Thyra::createMembersView(
  const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs,
  const RTOpPack::SubMultiVectorView<Scalar> &raw_mv,
  const std::string &label
  )
{
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
    mv = vs->createMembersView(raw_mv);
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase", &mv );
  if(label.length()) mv->setObjectLabel(label);
  return mv;
}

template<class Scalar>
Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >
Thyra::createMembersView(
  const VectorSpaceBase<Scalar> &vs,
  const RTOpPack::SubMultiVectorView<Scalar> &raw_mv,
  const std::string &label
  )
{
  return createMembersView(Teuchos::rcp(&vs,false),raw_mv,label);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::MultiVectorBase<Scalar> >
Thyra::createMembersView(
  const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs,
  const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv,
  const std::string &label
  )
{
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
    mv = vs->createMembersView(raw_mv);
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase", &mv );
  if(label.length())
    Teuchos::rcp_const_cast<MultiVectorBase<Scalar> >(mv)->setObjectLabel(label);
  return mv;
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::MultiVectorBase<Scalar> >
Thyra::createMembersView( const VectorSpaceBase<Scalar> &vs,
  const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv,
  const std::string &label
  )
{
  return createMembersView(Teuchos::rcp(&vs,false),raw_mv,label);
}

#endif // THYRA_VECTOR_SPACE_BASE_HPP
