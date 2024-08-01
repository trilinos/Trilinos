// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_VECTOR_SPACE_BASE_DEF_HPP
#define THYRA_VECTOR_SPACE_BASE_DEF_HPP

#include "Thyra_VectorSpaceBase_decl.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Teuchos_Tuple.hpp"


#ifdef TEUCHOS_DEBUG
#  define THYRA_INITIALIZE_VECS_MULTIVECS_WITH_NANS
#endif


namespace Thyra {


//
// VectorSpaceBase
//



// Virtual functions with default implementations


template<class Scalar>
bool VectorSpaceBase<Scalar>::isEuclidean() const
{
  return false;
}


template<class Scalar>
bool VectorSpaceBase<Scalar>::hasInCoreView(const Range1D& /* rng */,
  const EViewType /* viewType */, const EStrideType /* strideType */) const
{
  return false;
}


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
VectorSpaceBase<Scalar>::clone() const
{
  return Teuchos::null;
}


} // end namespace Thyra


//
// Nonmember functions
//


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Thyra::makeHaveOwnership( const RCP<const VectorSpaceBase<Scalar> > &vs_in )
{
  if (vs_in.has_ownership())
    return vs_in;
  const RCP<const VectorSpaceBase<Scalar> > vs = vs_in->clone();
  TEUCHOS_TEST_FOR_EXCEPTION(
    is_null(vs), std::logic_error
    ,"Thyra::makeHaveOwnership(vs): Error, the concrete VectorSpaceBase object identified as \'"
    << vs->description() << "\' does not support the clone() function!"
    );
  return vs;
}


template<class Scalar>
Teuchos::RCP< Thyra::VectorBase<Scalar> >
Thyra::createMember(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  const std::string &label
  )
{
  RCP<VectorBase<Scalar> > v = vs->createMember();
#ifdef THYRA_INITIALIZE_VECS_MULTIVECS_WITH_NANS
  if (vs->dim()) {
    v->assign(ScalarTraits<Scalar>::nan());
  }
#endif  
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase",
    Teuchos::outArg(v) );
  if (label.length()) v->setObjectLabel(label);
  return v;
}
  

template<class Scalar>
Teuchos::RCP< Thyra::VectorBase<Scalar> >
Thyra::createMember(
  const VectorSpaceBase<Scalar> &vs, const std::string &label
  )
{
  return createMember(Teuchos::rcpFromRef(vs), label);
}


template<class Scalar>
Teuchos::RCP< Thyra::MultiVectorBase<Scalar> >
Thyra::createMembers(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  int numMembers,  const std::string &label
  )
{
  RCP<MultiVectorBase<Scalar> >
    mv = vs->createMembers(numMembers);
#ifdef THYRA_INITIALIZE_VECS_MULTIVECS_WITH_NANS
  if (vs->dim()) {
    mv->assign(ScalarTraits<Scalar>::nan());
  }
#endif  
  Teuchos::set_extra_data(makeHaveOwnership(vs), "VectorSpaceBase",
    Teuchos::outArg(mv));
  if(label.length()) mv->setObjectLabel(label);
  return mv;
}


template<class Scalar>
Teuchos::RCP< Thyra::MultiVectorBase<Scalar> >
Thyra::createMembers(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  const RCP<const VectorSpaceBase<Scalar> > &domain,
  const std::string &label
  )
{
  return createMembers(vs, domain->dim(), label);
}


template<class Scalar>
Teuchos::RCP< Thyra::MultiVectorBase<Scalar> >
Thyra::createMembers(
  const VectorSpaceBase<Scalar> &vs, int numMembers,
  const std::string &label
  )
{
  return createMembers(Teuchos::rcp(&vs,false), numMembers, label);
}


template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
Thyra::createMemberView(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  const RTOpPack::SubVectorView<Scalar> &raw_v,
  const std::string &label
  )
{
  RCP<VectorBase<Scalar> >
    v = vs->createMemberView(raw_v);
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase",
    Teuchos::outArg(v) );
  if (label.length()) v->setObjectLabel(label);
  return v;
}


template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
Thyra::createMemberView(
  const VectorSpaceBase<Scalar> &vs,
  const RTOpPack::SubVectorView<Scalar> &raw_v,
  const std::string &label
  )
{
  return createMemberView(Teuchos::rcp(&vs,false),raw_v,label);
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
Thyra::createMemberView(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  const RTOpPack::ConstSubVectorView<Scalar> &raw_v,
  const std::string &label
  )
{
  RCP<const VectorBase<Scalar> >
    v = vs->createMemberView(raw_v);
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase",
    Teuchos::outArg(v) );
  if (label.length())
    Teuchos::rcp_const_cast<VectorBase<Scalar> >(v)->setObjectLabel(label);
  return v;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
Thyra::createMemberView(
  const VectorSpaceBase<Scalar> &vs,
  const RTOpPack::ConstSubVectorView<Scalar> &raw_v,
  const std::string &label
  )
{
  return createMemberView(Teuchos::rcp(&vs,false),raw_v,label);
}


template<class Scalar>
Teuchos::RCP<Thyra::MultiVectorBase<Scalar> >
Thyra::createMembersView(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  const RTOpPack::SubMultiVectorView<Scalar> &raw_mv,
  const std::string &label
  )
{
  RCP<MultiVectorBase<Scalar> >
    mv = vs->createMembersView(raw_mv);
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase",
    Teuchos::outArg(mv) );
  if (label.length()) mv->setObjectLabel(label);
  return mv;
}


template<class Scalar>
Teuchos::RCP<Thyra::MultiVectorBase<Scalar> >
Thyra::createMembersView(
  const VectorSpaceBase<Scalar> &vs,
  const RTOpPack::SubMultiVectorView<Scalar> &raw_mv,
  const std::string &label
  )
{
  return createMembersView(Teuchos::rcp(&vs,false),raw_mv,label);
}


template<class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> >
Thyra::createMembersView(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv,
  const std::string &label
  )
{
  RCP<const MultiVectorBase<Scalar> >
    mv = vs->createMembersView(raw_mv);
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase",
    Teuchos::outArg(mv) );
  if (label.length())
    Teuchos::rcp_const_cast<MultiVectorBase<Scalar> >(mv)->setObjectLabel(label);
  return mv;
}


template<class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> >
Thyra::createMembersView( const VectorSpaceBase<Scalar> &vs,
  const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv,
  const std::string &label
  )
{
  return createMembersView(Teuchos::rcp(&vs,false),raw_mv,label);
}



//
// Explicit instantiation macro
//
// Must be expanded from within the Thyra namespace!
//


#define THYRA_VECTOR_SPACE_BASE_INSTANT(SCALAR) \
  \
  template class VectorSpaceBase<SCALAR >; \
   \
  template RCP< VectorBase<SCALAR > > \
  createMember( \
    const RCP<const VectorSpaceBase<SCALAR > > &vs, \
    const std::string &label \
    ); \
   \
  template RCP< VectorBase<SCALAR > > \
  createMember( \
    const VectorSpaceBase<SCALAR > &vs, const std::string &label \
    ); \
   \
  template RCP< MultiVectorBase<SCALAR > > \
  createMembers( \
    const RCP<const VectorSpaceBase<SCALAR > > &vs, \
    int numMembers,  const std::string &label \
    ); \
   \
  template RCP< Thyra::MultiVectorBase<SCALAR > > \
  createMembers( \
    const RCP<const VectorSpaceBase<SCALAR > > &vs, \
    const RCP<const VectorSpaceBase<SCALAR > > &domain, \
    const std::string &label \
    ); \
  \
  template RCP< MultiVectorBase<SCALAR > > \
  createMembers( \
    const VectorSpaceBase<SCALAR > &vs, int numMembers, \
    const std::string &label \
    ); \
   \
  template RCP<VectorBase<SCALAR > > \
  createMemberView( \
    const RCP<const VectorSpaceBase<SCALAR > > &vs, \
    const RTOpPack::SubVectorView<SCALAR > &raw_v, \
    const std::string &label \
    ); \
   \
  template RCP<VectorBase<SCALAR > > \
  createMemberView( \
    const VectorSpaceBase<SCALAR > &vs, \
    const RTOpPack::SubVectorView<SCALAR > &raw_v, \
    const std::string &label \
    ); \
   \
  template RCP<const VectorBase<SCALAR > > \
  createMemberView( \
    const RCP<const VectorSpaceBase<SCALAR > > &vs, \
    const RTOpPack::ConstSubVectorView<SCALAR > &raw_v, \
    const std::string &label \
    ); \
   \
  template RCP<const VectorBase<SCALAR > > \
  createMemberView( \
    const VectorSpaceBase<SCALAR > &vs, \
    const RTOpPack::ConstSubVectorView<SCALAR > &raw_v, \
    const std::string &label \
    ); \
   \
  template RCP<MultiVectorBase<SCALAR > > \
  createMembersView( \
    const RCP<const VectorSpaceBase<SCALAR > > &vs, \
    const RTOpPack::SubMultiVectorView<SCALAR > &raw_mv, \
    const std::string &label \
    ); \
   \
  template RCP<MultiVectorBase<SCALAR > > \
  createMembersView( \
    const VectorSpaceBase<SCALAR > &vs, \
    const RTOpPack::SubMultiVectorView<SCALAR > &raw_mv, \
    const std::string &label \
    ); \
   \
  template RCP<const MultiVectorBase<SCALAR > > \
  createMembersView( \
    const RCP<const VectorSpaceBase<SCALAR > > &vs, \
    const RTOpPack::ConstSubMultiVectorView<SCALAR > &raw_mv, \
    const std::string &label \
    ); \
   \
  template RCP<const MultiVectorBase<SCALAR > > \
  createMembersView( const VectorSpaceBase<SCALAR > &vs, \
    const RTOpPack::ConstSubMultiVectorView<SCALAR > &raw_mv, \
    const std::string &label \
    );


#endif // THYRA_VECTOR_SPACE_BASE_DEF_HPP
