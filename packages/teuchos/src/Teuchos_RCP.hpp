// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

#ifndef TEUCHOS_RCP_HPP
#define TEUCHOS_RCP_HPP


/*! \file Teuchos_RCP.hpp
    \brief Reference-counted pointer class and non-member templated function implementations.
*/

/** \example example/RCP/cxx_main.cpp
    This is an example of how to use the <tt>Teuchos::RCP</tt> class.
*/

/** \example test/RCP/cxx_main.cpp
    This is a more detailed testing program that uses all of the <tt>Teuchos::RCP</tt> class.
*/

#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Exceptions.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_map.hpp"
#include "Teuchos_TypeNameTraits.hpp"


namespace Teuchos {


// Constructors/destructors/initializers


template<class T>
inline
RCP<T>::RCP( ENull )
  : ptr_(NULL)
{}


template<class T>
inline
RCP<T>::RCP( T* p, ENull null_arg )
  : ptr_(p)
#ifndef TEUCHOS_DEBUG
  , node_(RCP_createNewRCPNodeRawPtrNonowned(p))
#endif // TEUCHOS_DEBUG
{
#ifdef TEUCHOS_DEBUG
  if (p) {
    node_ = RCPNodeHandle(
      RCP_createNewRCPNodeRawPtrNonowned(p),
      p, typeName(*p), concreteTypeName(*p),
      false
      );
  }
#endif // TEUCHOS_DEBUG
}


template<class T>
inline
RCP<T>::RCP( T* p, bool has_ownership_in )
  : ptr_(p)
#ifndef TEUCHOS_DEBUG
  , node_(RCP_createNewRCPNodeRawPtr(p, has_ownership_in))
#endif // TEUCHOS_DEBUG
{
#ifdef TEUCHOS_DEBUG
  if (p) {
    node_ = RCPNodeHandle(
      RCP_createNewRCPNodeRawPtr(p, has_ownership_in),
      p, typeName(*p), concreteTypeName(*p),
      has_ownership_in
      );
  }
#endif // TEUCHOS_DEBUG
}


template<class T>
REFCOUNTPTR_INLINE
template<class Dealloc_T>
RCP<T>::RCP( T* p, Dealloc_T dealloc, bool has_ownership_in )
  : ptr_(p)
#ifndef TEUCHOS_DEBUG
  , node_(RCP_createNewDeallocRCPNodeRawPtr(p, dealloc, has_ownership_in))
#endif // TEUCHOS_DEBUG
{
#ifdef TEUCHOS_DEBUG
  if (p) {
    node_ = RCPNodeHandle(
      RCP_createNewDeallocRCPNodeRawPtr(p, dealloc, has_ownership_in),
      p, typeName(*p), concreteTypeName(*p),
      has_ownership_in
      );
  }
#endif // TEUCHOS_DEBUG
}


template<class T>
REFCOUNTPTR_INLINE
RCP<T>::RCP(const RCP<T>& r_ptr)
  : ptr_(r_ptr.ptr_), node_(r_ptr.node_)
{}


template<class T>
REFCOUNTPTR_INLINE
template<class T2>
RCP<T>::RCP(const RCP<T2>& r_ptr)
  : ptr_(r_ptr.get()), // will not compile if T is not base class of T2
    node_(r_ptr.access_private_node())
{}


template<class T>
REFCOUNTPTR_INLINE
RCP<T>::~RCP()
{}


template<class T>
REFCOUNTPTR_INLINE
RCP<T>& RCP<T>::operator=(const RCP<T>& r_ptr)
{
  if( this == &r_ptr )
    return *this; // Assignment to self
  node_ = r_ptr.access_private_node(); // May throw in debug mode!
  ptr_ = r_ptr.ptr_;
  return *this;
  // NOTE: It is critical that the assignment of ptr_ come *after* the
  // assignment of node_ since node_ might throw an exception!
}


template<class T>
REFCOUNTPTR_INLINE
template<class T2>
void RCP<T>::reset(T2* p, bool has_ownership_in)
{
  *this = rcp(p, has_ownership_in);
}


// Object query and access functions


template<class T>
inline
bool RCP<T>::is_null() const
{
  return ptr_ == 0;
}


template<class T>
inline
T* RCP<T>::operator->() const
{
  debug_assert_not_null();
  debug_assert_valid_ptr();
  return ptr_;
}


template<class T>
inline
T& RCP<T>::operator*() const
{
  debug_assert_not_null();
  debug_assert_valid_ptr();
  return *ptr_;
}

template<class T>
inline
T* RCP<T>::get() const
{
  debug_assert_valid_ptr();
  return ptr_;
}


template<class T>
inline
T* RCP<T>::getRawPtr() const
{
  return this->get();
}


template<class T>
inline
Ptr<T> RCP<T>::ptr() const
{
#ifdef TEUCHOS_DEBUG
  return Ptr<T>(this->create_weak());
#else
  return Ptr<T>(getRawPtr());
#endif
}


// Reference counting


template<class T>
inline
ERCPStrength RCP<T>::strength() const
{
  return node_.strength();
}


template<class T>
inline
bool RCP<T>::is_valid_ptr() const
{
  if (ptr_)
    return node_.is_valid_ptr();
  return true;
}


template<class T>
inline
int RCP<T>::strong_count() const
{
  return node_.strong_count();
}


template<class T>
inline
int RCP<T>::weak_count() const
{
  return node_.weak_count();
}


template<class T>
inline
int RCP<T>::total_count() const
{
  return node_.total_count();
}


template<class T>
REFCOUNTPTR_INLINE
void RCP<T>::set_has_ownership()
{
  node_.has_ownership(true);
}


template<class T>
REFCOUNTPTR_INLINE
bool RCP<T>::has_ownership() const
{
  return node_.has_ownership();
}


template<class T>
REFCOUNTPTR_INLINE
Ptr<T> RCP<T>::release()
{
  debug_assert_valid_ptr();
  node_.has_ownership(false);
  return Ptr<T>(ptr_);
}


template<class T>
inline
RCP<T> RCP<T>::create_weak() const
{
  debug_assert_valid_ptr();
  return RCP<T>(ptr_, node_.create_weak());
}


template<class T>
inline
RCP<T> RCP<T>::create_strong() const
{
  debug_assert_valid_ptr();
  return RCP<T>(ptr_, node_.create_strong());
}


template<class T>
REFCOUNTPTR_INLINE
template <class T2>
bool RCP<T>::shares_resource(const RCP<T2>& r_ptr) const
{
  return node_.same_node(r_ptr.access_private_node());
  // Note: above, r_ptr is *not* the same class type as *this so we can not
  // access its node_ member directly!  This is an interesting detail to the
  // C++ protected/private protection mechanism!
}


// Assertions


template<class T>
inline
const RCP<T>& RCP<T>::assert_not_null() const
{
  if (!ptr_)
    throw_null_ptr_error(typeName(*this));
  return *this;
}


template<class T>
inline
const RCP<T>& RCP<T>::assert_valid_ptr() const
{
  if (ptr_)
    node_.assert_valid_ptr(*this);
  return *this;
}


// Deprecated


template<class T>
REFCOUNTPTR_INLINE
int RCP<T>::count() const
{
  return node_.count();
}


// very bad public functions


template<class T>
inline
RCPNode* RCP_createNewRCPNodeRawPtrNonowned( T* p )
{
  return new RCPNodeTmpl<T,DeallocNull<T> >(
    p, DeallocNull<T>(), false
    );
}


template<class T>
inline
RCPNode* RCP_createNewRCPNodeRawPtr( T* p, bool has_ownership_in )
{
  return new RCPNodeTmpl<T,DeallocDelete<T> >(
    p, DeallocDelete<T>(), has_ownership_in
    );
}


template<class T, class Dealloc_T>
inline
RCPNode* RCP_createNewDeallocRCPNodeRawPtr(
  T* p, Dealloc_T dealloc, bool has_ownership_in
  )
{
  return new RCPNodeTmpl<T,Dealloc_T>(p, dealloc, has_ownership_in);
}


template<class T>
inline
RCP<T>::RCP( T* p, const RCPNodeHandle& node)
  : ptr_(p), node_(node)
{}


template<class T>
inline
T* RCP<T>::access_private_ptr() const
{  return ptr_; }


template<class T>
inline
RCPNodeHandle& RCP<T>::nonconst_access_private_node()
{  return node_; }


template<class T>
inline
const RCPNodeHandle& RCP<T>::access_private_node() const
{  return node_; }


}  // end namespace Teuchos


// /////////////////////////////////////////////////////////////////////////////////
// Inline non-member functions for RCP


template<class T>
inline
Teuchos::RCP<T>
Teuchos::rcp( T* p, bool owns_mem )
{
  return RCP<T>(p, owns_mem);
}


template<class T, class Dealloc_T>
inline
Teuchos::RCP<T>
Teuchos::rcp( T* p, Dealloc_T dealloc, bool owns_mem )
{
  return RCP<T>(p, dealloc, owns_mem);
}


template<class T>
Teuchos::RCP<T>
Teuchos::rcpFromRef( T& r )
{
  return RCP<T>(&r, null);
}


template<class T, class Embedded>
Teuchos::RCP<T>
Teuchos::rcpWithEmbeddedObjPreDestroy(
  T* p, const Embedded &embedded, bool owns_mem
  )
{
  return rcp(
    p, embeddedObjDeallocDelete<T>(embedded,PRE_DESTROY), owns_mem
    );
}


template<class T, class Embedded>
Teuchos::RCP<T>
Teuchos::rcpWithEmbeddedObjPostDestroy(
  T* p, const Embedded &embedded, bool owns_mem
  )
{
  return rcp( p, embeddedObjDeallocDelete<T>(embedded,POST_DESTROY), owns_mem );
}


template<class T, class Embedded>
Teuchos::RCP<T>
Teuchos::rcpWithEmbeddedObj( T* p, const Embedded &embedded, bool owns_mem )
{
  return rcpWithEmbeddedObjPostDestroy<T,Embedded>(p,embedded,owns_mem);
}


template<class T>
REFCOUNTPTR_INLINE
bool Teuchos::is_null( const RCP<T> &p )
{
  return p.is_null();
}


template<class T>
REFCOUNTPTR_INLINE
bool Teuchos::nonnull( const RCP<T> &p )
{
  return !p.is_null();
}


template<class T>
REFCOUNTPTR_INLINE
bool Teuchos::operator==( const RCP<T> &p, ENull )
{
  return p.get() == NULL;
}


template<class T>
REFCOUNTPTR_INLINE
bool Teuchos::operator!=( const RCP<T> &p, ENull )
{
  return p.get() != NULL;
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator==( const RCP<T1> &p1, const RCP<T2> &p2 )
{
  return p1.access_private_node().same_node(p2.access_private_node());
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator!=( const RCP<T1> &p1, const RCP<T2> &p2 )
{
  return !p1.access_private_node().same_node(p2.access_private_node());
}


template<class T2, class T1>
REFCOUNTPTR_INLINE
Teuchos::RCP<T2>
Teuchos::rcp_implicit_cast(const RCP<T1>& p1)
{
  // Make the compiler check if the conversion is legal
  T2 *check = p1.get();
  return RCP<T2>(check, p1.access_private_node());
}


template<class T2, class T1>
REFCOUNTPTR_INLINE
Teuchos::RCP<T2>
Teuchos::rcp_static_cast(const RCP<T1>& p1)
{
  // Make the compiler check if the conversion is legal
  T2 *check = static_cast<T2*>(p1.get());
  return RCP<T2>(check, p1.access_private_node());
}


template<class T2, class T1>
REFCOUNTPTR_INLINE
Teuchos::RCP<T2>
Teuchos::rcp_const_cast(const RCP<T1>& p1)
{
  // Make the compiler check if the conversion is legal
  T2 *check = const_cast<T2*>(p1.get());
  return RCP<T2>(check, p1.access_private_node());
}


template<class T2, class T1>
REFCOUNTPTR_INLINE
Teuchos::RCP<T2>
Teuchos::rcp_dynamic_cast(const RCP<T1>& p1, bool throw_on_fail)
{
  if (!is_null(p1)) {
    T2 *p = NULL;
    if (throw_on_fail) {
      p = &dyn_cast<T2>(*p1);
    }
    else {
      // Make the compiler check if the conversion is legal
      p = dynamic_cast<T2*>(p1.get());
    }
    if (p) {
      return RCP<T2>(p, p1.access_private_node());
    }
  }
  return null;
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
void Teuchos::set_extra_data( const T1 &extra_data, const std::string& name,
  const Ptr<RCP<T2> > &p, EPrePostDestruction destroy_when, bool force_unique )
{
  p->assert_not_null();
  p->nonconst_access_private_node().set_extra_data(
    any(extra_data), name, destroy_when,
    force_unique );
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
const T1& Teuchos::get_extra_data( const RCP<T2>& p, const std::string& name )
{
  p.assert_not_null();
  return any_cast<T1>(
    p.access_private_node().get_extra_data(
      TypeNameTraits<T1>::name(), name
      )
    );
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
T1& Teuchos::get_nonconst_extra_data( RCP<T2>& p, const std::string& name )
{
  p.assert_not_null();
  return any_cast<T1>(
    p.nonconst_access_private_node().get_extra_data(
      TypeNameTraits<T1>::name(), name
      )
    );
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
Teuchos::Ptr<const T1>
Teuchos::get_optional_extra_data( const RCP<T2>& p, const std::string& name )
{
  p.assert_not_null();
  const any *extra_data = p.access_private_node().get_optional_extra_data(
    TypeNameTraits<T1>::name(), name);
  if (extra_data)
    return Ptr<const T1>(&any_cast<T1>(*extra_data));
  return null;
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
Teuchos::Ptr<T1>
Teuchos::get_optional_nonconst_extra_data( RCP<T2>& p, const std::string& name )
{
  p.assert_not_null();
  any *extra_data = p.nonconst_access_private_node().get_optional_extra_data(
    TypeNameTraits<T1>::name(), name);
  if (extra_data)
    return Ptr<T1>(&any_cast<T1>(*extra_data));
  return null;
}


template<class Dealloc_T, class T>
inline
const Dealloc_T& Teuchos::get_dealloc( const RCP<T>& p )
{
  return get_nonconst_dealloc<Dealloc_T>(const_cast<RCP<T>&>(p));
}


template<class Dealloc_T, class T>
REFCOUNTPTR_INLINE
Dealloc_T& Teuchos::get_nonconst_dealloc( const RCP<T>& p )
{
  typedef RCPNodeTmpl<typename Dealloc_T::ptr_t,Dealloc_T>  requested_type;
  p.assert_not_null();
  RCPNodeTmpl<typename Dealloc_T::ptr_t,Dealloc_T>
    *dnode = dynamic_cast<RCPNodeTmpl<typename Dealloc_T::ptr_t,Dealloc_T>*>(
      p.access_private_node().node_ptr());
  TEST_FOR_EXCEPTION(
    dnode==NULL, NullReferenceError
    ,"get_dealloc<" << TypeNameTraits<Dealloc_T>::name()
    << "," << TypeNameTraits<T>::name() << ">(p): "
    << "Error, requested type \'" << TypeNameTraits<requested_type>::name()
    << "\' does not match actual type of the node \'"
    << typeName(*p.access_private_node().node_ptr()) << "!"
    );
  return dnode->get_nonconst_dealloc();
}


template<class Dealloc_T, class T>
REFCOUNTPTR_INLINE
Teuchos::Ptr<Dealloc_T>
Teuchos::get_optional_nonconst_dealloc( const RCP<T>& p )
{
  p.assert_not_null();
  typedef RCPNodeTmpl<typename Dealloc_T::ptr_t,Dealloc_T> RCPNT;
  RCPNT *dnode = dynamic_cast<RCPNT*>(p.access_private_node().node_ptr());
  if(dnode)
    return ptr(&dnode->get_nonconst_dealloc());
  return null;
}


template<class Dealloc_T, class T>
inline
Teuchos::Ptr<const Dealloc_T>
Teuchos::get_optional_dealloc( const RCP<T>& p )
{
  return get_optional_nonconst_dealloc<Dealloc_T>(const_cast<RCP<T>&>(p));
}


template<class TOrig, class Embedded, class T>
const Embedded& Teuchos::getEmbeddedObj( const RCP<T>& p )
{
  typedef EmbeddedObjDealloc<TOrig,Embedded,DeallocDelete<TOrig> > Dealloc_t;
  return get_dealloc<Dealloc_t>(p).getObj();
}


template<class TOrig, class Embedded, class T>
Embedded& Teuchos::getNonconstEmbeddedObj( const RCP<T>& p )
{
  typedef EmbeddedObjDealloc<TOrig,Embedded,DeallocDelete<TOrig> > Dealloc_t;
  return get_nonconst_dealloc<Dealloc_t>(p).getNonconstObj();
}


template<class TOrig, class Embedded, class T>
Teuchos::Ptr<const Embedded>
Teuchos::getOptionalEmbeddedObj( const RCP<T>& p )
{
  typedef EmbeddedObjDealloc<TOrig,Embedded,DeallocDelete<TOrig> > Dealloc_t;
  const Ptr<const Dealloc_t> dealloc = get_optional_dealloc<Dealloc_t>(p);
  if (!is_null(dealloc)) {
    return ptr(&dealloc->getObj());
  }
  return null;
}


template<class TOrig, class Embedded, class T>
Teuchos::Ptr<Embedded>
Teuchos::getOptionalNonconstEmbeddedObj( const RCP<T>& p )
{
  typedef EmbeddedObjDealloc<TOrig,Embedded,DeallocDelete<TOrig> > Dealloc_t;
  const Ptr<Dealloc_t> dealloc = get_optional_nonconst_dealloc<Dealloc_t>(p);
  if (!is_null(dealloc)) {
    return ptr(&dealloc->getNonconstObj());
  }
  return null;
}


template<class T>
std::ostream& Teuchos::operator<<( std::ostream& out, const RCP<T>& p )
{
  out
    << typeName(p) << "{"
    << "ptr="<<(const void*)(p.get()) // I can't find any alternative to this C cast :-(
    <<",node="<<p.access_private_node()
    <<",count="<<p.count()
    <<"}";
  return out;
}


#endif // TEUCHOS_RCP_HPP
