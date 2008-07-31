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

#ifndef TEUCHOS_ARRAY_RCP_HPP
#define TEUCHOS_ARRAY_RCP_HPP


#include "Teuchos_ArrayRCPDecl.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_map.hpp"


namespace Teuchos {


// Constructors/Initializers


template<class T>
inline
ArrayRCP<T>::ArrayRCP( ENull )
  : ptr_(NULL)
  , node_(NULL)
  , lowerOffset_(0)
  , upperOffset_(-1)
{}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>::ArrayRCP(const ArrayRCP<T>& r_ptr)
  :ptr_(r_ptr.ptr_), node_(r_ptr.node_),
   lowerOffset_(r_ptr.lowerOffset_),
   upperOffset_(r_ptr.upperOffset_)
{
  if(node_) node_->incr_count();
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>::~ArrayRCP()
{
  if(node_ && node_->deincr_count() == 0 ) {
#ifdef TEUCHOS_DEBUG
    local_printActiveRCPNodes.foo(); // Make sure this object is used!
    remove_RCPNode(node_);
#endif
    delete node_;
  }
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>& ArrayRCP<T>::operator=(const ArrayRCP<T>& r_ptr)
{
  if( this == &r_ptr )
    return *this; // Assignment to self
  if( node_ && !node_->deincr_count() ) {
#ifdef TEUCHOS_DEBUG
    remove_RCPNode(node_);
#endif
    delete node_;
  }
  ptr_   = r_ptr.ptr_;
  node_  = r_ptr.node_;
  lowerOffset_ = r_ptr.lowerOffset_;
  upperOffset_ = r_ptr.upperOffset_;
  if(node_) node_->incr_count();
  return *this;
}


// Object/Pointer Access Functions


template<class T>
inline
T* ArrayRCP<T>::operator->() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assert_in_range(0,1);
#endif
  return ptr_;
}


template<class T>
inline
T& ArrayRCP<T>::operator*() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assert_in_range(0,1);
#endif
  return *ptr_;
}


template<class T>
inline
T* ArrayRCP<T>::get() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  if(ptr_) {
    assert_in_range(0,1);
  }
#endif
  return ptr_;
}


template<class T>
inline
T* ArrayRCP<T>::getRawPtr() const
{
  return this->get();
}


template<class T>
REFCOUNTPTR_INLINE
T& ArrayRCP<T>::operator[](Ordinal offset) const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assert_in_range(offset,1);
#endif
  return ptr_[offset];
}


// Pointer Arithmetic Functions


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>& ArrayRCP<T>::operator++()
{
  if(ptr_) {
    ++ptr_;
    --lowerOffset_;
    --upperOffset_;
  }
  return *this;
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T> ArrayRCP<T>::operator++(int)
{
  ArrayRCP<T> r_ptr = *this;
  ++(*this);
  return r_ptr;
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>& ArrayRCP<T>::operator--()
{
  if(ptr_) {
    --ptr_;
    ++lowerOffset_;
    ++upperOffset_;
  }
  return *this;
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T> ArrayRCP<T>::operator--(int)
{
  ArrayRCP<T> r_ptr = *this;
  --(*this);
  return r_ptr;
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>& ArrayRCP<T>::operator+=(Ordinal offset)
{
  if(ptr_) {
    ptr_ += offset;
    lowerOffset_ -= offset;
    upperOffset_ -= offset;
  }
  return *this;
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>& ArrayRCP<T>::operator-=(Ordinal offset)
{
  if(ptr_) {
    ptr_ -= offset;
    lowerOffset_ += offset;
    upperOffset_ += offset;
  }
  return *this;
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T> ArrayRCP<T>::operator+(Ordinal offset) const
{
  ArrayRCP<T> r_ptr = *this;
  r_ptr+=(offset);
  return r_ptr;
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T> ArrayRCP<T>::operator-(Ordinal offset) const
{
  ArrayRCP<T> r_ptr = *this;
  r_ptr-=offset;
  return r_ptr;
}


// ArrayView views


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<const T> ArrayRCP<T>::getConst() const
{
  const T *cptr = ptr_; // Will not compile if not legal!
  return ArrayRCP<const T>(cptr,lowerOffset_,upperOffset_,node_);
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>
ArrayRCP<T>::persistingView( Ordinal lowerOffset_in, Ordinal size_in ) const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assert_in_range(lowerOffset_in,size_in);
#endif
  ArrayRCP<T> ptr = *this;
  ptr.ptr_ = ptr.ptr_ + lowerOffset_in;
  ptr.lowerOffset_ = 0;
  ptr.upperOffset_ = size_in-1;
  return ptr;
}


// General query functions


template<class T>
REFCOUNTPTR_INLINE
int ArrayRCP<T>::count() const {
  if(node_)
    return node_->count();
  return 0;
}


template<class T>
REFCOUNTPTR_INLINE
template <class T2>
bool ArrayRCP<T>::shares_resource(const ArrayRCP<T2>& r_ptr) const
{
  return node_ == r_ptr.access_node();
  // Note: above, r_ptr is *not* the same class type as *this so we can not
  // access its node_ member directly!  This is an interesting detail to the
  // C++ protected/private protection mechanism!
}


template<class T>
REFCOUNTPTR_INLINE
typename ArrayRCP<T>::Ordinal
ArrayRCP<T>::lowerOffset() const
{
  return lowerOffset_;
}


template<class T>
REFCOUNTPTR_INLINE
typename ArrayRCP<T>::Ordinal
ArrayRCP<T>::upperOffset() const
{
  return upperOffset_;
}


template<class T>
REFCOUNTPTR_INLINE
typename ArrayRCP<T>::Ordinal
ArrayRCP<T>::size() const
{
  return upperOffset_-lowerOffset_+1;
}


// Standard Container-Like Functions


template<class T>
REFCOUNTPTR_INLINE
typename ArrayRCP<T>::const_iterator ArrayRCP<T>::begin() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return *this;
#else
  return ptr_;
#endif
}


template<class T>
REFCOUNTPTR_INLINE
typename ArrayRCP<T>::const_iterator
ArrayRCP<T>::end() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return *this + (upperOffset_ + 1);
#else
  return ptr_ + (upperOffset_ + 1);
#endif
}


// Views 


template<class T> inline
ArrayView<T> ArrayRCP<T>::view( Ordinal lowerOffset_in, Ordinal size_in ) const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assert_in_range(lowerOffset_in,size_in);
#endif
  return arrayView(ptr_ + lowerOffset_in, size_in );
}


template<class T> inline
ArrayView<T> ArrayRCP<T>::operator()( Ordinal lowerOffset_in, Ordinal size_in ) const
{
  return view(lowerOffset_in,size_in);
}


template<class T> inline
ArrayView<T> ArrayRCP<T>::operator()() const
{
  if (!size())
    return null;
  return arrayView(ptr_ + lowerOffset_, size() );
}


template<class T> inline
ArrayRCP<T>::operator ArrayView<T>() const
{
  return this->operator()();
}


template<class T> inline
ArrayRCP<T>::operator ArrayRCP<const T>() const
{
  return ArrayRCP<const T>(ptr_,lowerOffset_,upperOffset_,node_);
}


// Ownership


template<class T>
REFCOUNTPTR_INLINE
T* ArrayRCP<T>::release()
{
  if(node_)
    node_->has_ownership(false);
  return ptr_;
}


template<class T>
REFCOUNTPTR_INLINE
void ArrayRCP<T>::set_has_ownership()
{
  if(node_)
    node_->has_ownership(true);
}


template<class T>
REFCOUNTPTR_INLINE
bool ArrayRCP<T>::has_ownership() const
{
  if(node_)
    return node_->has_ownership();
  return false;
}


// Assertion Functions.


template<class T>
inline
const ArrayRCP<T>&
ArrayRCP<T>::assert_not_null() const
{
  if(!ptr_)
    throw_null_ptr_error(TypeNameTraits<T>::name());
  return *this;
}


template<class T>
inline
const ArrayRCP<T>&
ArrayRCP<T>::assert_in_range( Ordinal lowerOffset_in, Ordinal size_in ) const
{
  assert_not_null();
  TEST_FOR_EXCEPTION(
    !( lowerOffset_ <= lowerOffset_in && lowerOffset_in+size_in-1 <= upperOffset_ ),
    Teuchos::RangeError,
    "Teuchos::ArrayRCP<"<<TypeNameTraits<T>::name()<<">::assert_in_range:"
    " Error, [lowerOffset,lowerOffset+size-1] = ["
    <<lowerOffset_in<<","<<(lowerOffset_in+size_in-1)<<"] does not lie in the"
    " range ["<<lowerOffset_<<","<<upperOffset_<<"]!"
    );
  return *this;
}


// very bad public functions


template<class T>
inline
ArrayRCP<T>::ArrayRCP(
  T* p, Ordinal lowerOffset_in, Ordinal upperOffset_in, bool has_ownership_in
  )
  : ptr_(p)
  , node_(
    p
    ? new RCPNodeTmpl<T,DeallocArrayDelete<T> >(
      p,DeallocArrayDelete<T>(),has_ownership_in
      )
    : NULL
    )
  ,lowerOffset_(lowerOffset_in)
  ,upperOffset_(upperOffset_in)
{
#ifdef TEUCHOS_DEBUG
  if(node_ && isTracingActiveRCPNodes()) {
    std::ostringstream os;
    os << "{T=\'"<<TypeNameTraits<T>::name()<<"\',Concrete T=\'"
       <<typeName(*p)<<"\',p="<<p<<",has_ownership="<<has_ownership_in<<"}";
    add_new_RCPNode(node_,os.str());
  }
#endif
}


template<class T>
REFCOUNTPTR_INLINE
template<class Dealloc_T>
ArrayRCP<T>::ArrayRCP(
  T* p, Ordinal lowerOffset_in, Ordinal upperOffset_in,
  Dealloc_T dealloc, bool has_ownership_in
  )
  : ptr_(p)
  , node_( p
    ? new RCPNodeTmpl<T,Dealloc_T>(p,dealloc,has_ownership_in)
    : NULL )
  ,lowerOffset_(lowerOffset_in)
  ,upperOffset_(upperOffset_in)
{
#ifdef TEUCHOS_DEBUG
  if(node_ && isTracingActiveRCPNodes()) {
    std::ostringstream os;
    os << "{T=\'"<<TypeNameTraits<T>::name()<<"\',Concrete T=\'"
       <<typeName(*p)<<"\',p="<<p<<",has_ownership="<<has_ownership_in<<"}";
    add_new_RCPNode(node_,os.str());
  }
#endif
}


template<class T>
inline
ArrayRCP<T>::ArrayRCP(
  T* p, Ordinal lowerOffset_in, Ordinal upperOffset_in, RCPNode* node
  )
  :ptr_(p),
   node_(node),
   lowerOffset_(lowerOffset_in),
   upperOffset_(upperOffset_in)
{
  if(node_)
    node_->incr_count();
}


template<class T>
inline
T*& ArrayRCP<T>::access_ptr()
{ 
  return ptr_;
}


template<class T>
inline
T* ArrayRCP<T>::access_ptr() const
{
  return ptr_;
}


template<class T>
inline
RCPNode*& ArrayRCP<T>::access_node()
{
  return node_;
}


template<class T>
inline
RCPNode* ArrayRCP<T>::access_node() const
{
  return node_;
}


}  // end namespace Teuchos


// ///////////////////////////////////////////
// Non-member functions for ArrayRCP


namespace Teuchos {
namespace Utilities {
template<class T1, class T2>
inline void assert_shares_resource(
  const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2
  )
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_FOR_EXCEPTION(
    !p1.shares_resource(p2), IncompatibleIteratorsError,
    "Error, these iterators are *not* pointing to the same valid memory!"
    );
#endif
}
} // namespace Utilities
} // namespace Teuchos


template<class T>
inline
Teuchos::ArrayRCP<T>
Teuchos::arcp(
T* p, typename ArrayRCP<T>::Ordinal lowerOffset
  ,typename ArrayRCP<T>::Ordinal size_in
  ,bool owns_mem
  )
{
  return ArrayRCP<T>(p,lowerOffset,lowerOffset+size_in-1,owns_mem);
}


template<class T, class Dealloc_T>
inline
Teuchos::ArrayRCP<T>
Teuchos::arcp(
T* p, typename ArrayRCP<T>::Ordinal lowerOffset
  ,typename ArrayRCP<T>::Ordinal size_in
  ,Dealloc_T dealloc, bool owns_mem
  )
{
  return ArrayRCP<T>(p,lowerOffset,lowerOffset+size_in-1,dealloc,owns_mem);
}


template<class T>
inline
Teuchos::ArrayRCP<T>
Teuchos::arcp( typename ArrayRCP<T>::Ordinal size )
{
  return ArrayRCP<T>(new T[size],0,size-1,true);
}


template<class T>
inline
Teuchos::ArrayRCP<T>
Teuchos::arcpClone( const ArrayView<const T> &v )
{
  const ArrayRCP<T> new_arcp = arcp<T>(v.size());
  std::copy( v.begin(), v.end(), new_arcp.begin() );
  return new_arcp;
}


template<class T, class Embedded>
Teuchos::ArrayRCP<T>
Teuchos::arcpWithEmbeddedObjPreDestroy(
  T* p,
  typename ArrayRCP<T>::Ordinal lowerOffset,
  typename ArrayRCP<T>::Ordinal size,
  const Embedded &embedded,
  bool owns_mem
  )
{
  return arcp(
    p, lowerOffset, size,
    embeddedObjDeallocDelete<T>(embedded,PRE_DESTROY),
    owns_mem
    );
}


template<class T, class Embedded>
Teuchos::ArrayRCP<T>
Teuchos::arcpWithEmbeddedObjPostDestroy(
  T* p,
  typename ArrayRCP<T>::Ordinal lowerOffset,
  typename ArrayRCP<T>::Ordinal size,
  const Embedded &embedded,
  bool owns_mem
  )
{
  return arcp(
    p, lowerOffset, size,
    embeddedObjDeallocDelete<T>(embedded,POST_DESTROY),
    owns_mem
    );
}


template<class T, class Embedded>
Teuchos::ArrayRCP<T>
Teuchos::arcpWithEmbeddedObj(
  T* p,
  typename ArrayRCP<T>::Ordinal lowerOffset,
  typename ArrayRCP<T>::Ordinal size,
  const Embedded &embedded,
  bool owns_mem
  )
{
  return arcpWithEmbeddedObjPostDestroy<T,Embedded>(
    p, lowerOffset, size, embedded, owns_mem );
}


template<class T>
REFCOUNTPTR_INLINE
Teuchos::ArrayRCP<T>
Teuchos::arcp( const RCP<std::vector<T> > &v )
{
  if ( is_null(v) || !v->size() )
    return null;
  return arcpWithEmbeddedObjPostDestroy<T,RCP<std::vector<T> > >(
    &(*v)[0], 0, v->size(),
    v, false
    );
}


template<class T>
REFCOUNTPTR_INLINE
Teuchos::ArrayRCP<const T>
Teuchos::arcp( const RCP<const std::vector<T> > &v )
{
  if ( is_null(v) || !v->size() )
    return null;
  return arcpWithEmbeddedObjPostDestroy<const T,RCP<const std::vector<T> > >(
    &(*v)[0], 0, v->size(),
    v, false
    );
}


template<class T>
REFCOUNTPTR_INLINE
Teuchos::RCP<std::vector<T> >
Teuchos::get_std_vector( const ArrayRCP<T> &ptr )
{
  return getEmbeddedObj<T,RCP<std::vector<T> > >(ptr);
}


template<class T>
REFCOUNTPTR_INLINE
Teuchos::RCP<const std::vector<T> >
Teuchos::get_std_vector( const ArrayRCP<const T> &ptr )
{
  return getEmbeddedObj<const T,RCP<const std::vector<T> > >(ptr);
}


template<class T>
REFCOUNTPTR_INLINE
bool Teuchos::is_null( const ArrayRCP<T> &p )
{
  return p.access_ptr() == NULL;
}


template<class T>
REFCOUNTPTR_INLINE
bool Teuchos::operator==( const ArrayRCP<T> &p, ENull )
{
  return p.access_ptr() == NULL;
}


template<class T>
REFCOUNTPTR_INLINE
bool Teuchos::operator!=( const ArrayRCP<T> &p, ENull )
{
  return p.access_ptr() != NULL;
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator==( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  return p1.access_ptr() == p2.access_ptr();
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator!=( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  return p1.access_ptr() != p2.access_ptr();
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator<( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  return p1.access_ptr() < p2.access_ptr();
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator<=( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  Utilities::assert_shares_resource(p1,p2);
  return p1.access_ptr() <= p2.access_ptr();
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator>( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  Utilities::assert_shares_resource(p1,p2);
  return p1.access_ptr() > p2.access_ptr();
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator>=( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  Utilities::assert_shares_resource(p1,p2);
  return p1.access_ptr() >= p2.access_ptr();
}


template<class T>
typename Teuchos::ArrayRCP<T>::difference_type
Teuchos::operator-( const ArrayRCP<T> &p1, const ArrayRCP<T> &p2 )
{
  Utilities::assert_shares_resource(p1,p2);
  return p1.access_ptr() - p2.access_ptr();
}


template<class T2, class T1>
REFCOUNTPTR_INLINE
Teuchos::ArrayRCP<T2>
Teuchos::arcp_reinterpret_cast(const ArrayRCP<T1>& p1)
{
  typedef typename ArrayRCP<T1>::Ordinal Ordinal;
  const int sizeOfT2ToT1 = sizeof(T2)/sizeof(T1);
  Ordinal lowerOffset2 = p1.lowerOffset() / sizeOfT2ToT1;
  Ordinal upperOffset2 = (p1.upperOffset()+1) / sizeOfT2ToT1 -1;
  T2 *ptr2 = reinterpret_cast<T2*>(p1.get());
  return ArrayRCP<T2>(
    ptr2, lowerOffset2, upperOffset2,
    p1.access_node()
    );
  // Note: Above is just fine even if p1.get()==NULL!
}


template<class T2, class T1>
REFCOUNTPTR_INLINE
Teuchos::ArrayRCP<T2>
Teuchos::arcp_const_cast(const ArrayRCP<T1>& p1)
{
  typedef typename ArrayRCP<T1>::Ordinal Ordinal;
  T2 *ptr2 = const_cast<T2*>(p1.get());
  return ArrayRCP<T2>(
    ptr2, p1.lowerOffset(), p1.upperOffset(),
    p1.access_node()
    );
  // Note: Above is just fine even if p1.get()==NULL!
}


template<class T2, class T1>
REFCOUNTPTR_INLINE
Teuchos::ArrayRCP<T2>
Teuchos::arcp_implicit_cast(const ArrayRCP<T1>& p1)
{
  typedef typename ArrayRCP<T1>::Ordinal Ordinal;
  T2 * raw_ptr2 = p1.get();
  return ArrayRCP<T2>(
    raw_ptr2,p1.lowerOffset(),p1.upperOffset()
    ,p1.access_node()
    );
  // Note: Above is just fine even if p1.get()==NULL!
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
void Teuchos::set_extra_data(
  const T1 &extra_data, const std::string& name, Teuchos::ArrayRCP<T2> *p
  ,EPrePostDestruction destroy_when, bool force_unique
  )
{
  p->assert_not_null();
  p->access_node()->set_extra_data( any(extra_data), name, destroy_when, force_unique );
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
T1& Teuchos::get_extra_data( ArrayRCP<T2>& p, const std::string& name )
{
  p.assert_not_null();
  return any_cast<T1>(p.access_node()->get_extra_data(TypeNameTraits<T1>::name(),name));
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
const T1& Teuchos::get_extra_data( const ArrayRCP<T2>& p, const std::string& name )
{
  p.assert_not_null();
  return any_cast<T1>(p.access_node()->get_extra_data(TypeNameTraits<T1>::name(),name));
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
T1* Teuchos::get_optional_extra_data( ArrayRCP<T2>& p, const std::string& name )
{
  p.assert_not_null();
  any *extra_data = p.access_node()->get_optional_extra_data(TypeNameTraits<T1>::name(),name);
  if( extra_data ) return &any_cast<T1>(*extra_data);
  return NULL;
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
const T1* Teuchos::get_optional_extra_data( const ArrayRCP<T2>& p, const std::string& name )
{
  p.assert_not_null();
  any *extra_data = p.access_node()->get_optional_extra_data(TypeNameTraits<T1>::name(),name);
  if( extra_data ) return &any_cast<T1>(*extra_data);
  return NULL;
}


template<class Dealloc_T, class T>
REFCOUNTPTR_INLINE
const Dealloc_T&
Teuchos::get_dealloc( const ArrayRCP<T>& p )
{
  return get_nonconst_dealloc<Dealloc_T>(p);
}


template<class Dealloc_T, class T>
inline
Dealloc_T& 
Teuchos::get_nonconst_dealloc( const Teuchos::ArrayRCP<T>& p )
{
  typedef RCPNodeTmpl<typename Dealloc_T::ptr_t,Dealloc_T>  requested_type;
  p.assert_not_null();
  RCPNodeTmpl<typename Dealloc_T::ptr_t,Dealloc_T>
    *dnode = dynamic_cast<RCPNodeTmpl<typename Dealloc_T::ptr_t,Dealloc_T>*>(
      p.access_node());
  TEST_FOR_EXCEPTION(
    dnode==NULL, NullReferenceError
    ,"get_dealloc<" << TypeNameTraits<Dealloc_T>::name()
    << "," << TypeNameTraits<T>::name() << ">(p): "
    << "Error, requested type \'" << TypeNameTraits<requested_type>::name()
    << "\' does not match actual type of the node \'" << typeName(*p.access_node()) << "!"
    );
  return dnode->get_nonconst_dealloc();
}


template<class Dealloc_T, class T>
REFCOUNTPTR_INLINE
const Dealloc_T*
Teuchos::get_optional_dealloc( const ArrayRCP<T>& p )
{
  return get_optional_dealloc<Dealloc_T>(p);
}


template<class Dealloc_T, class T>
inline
Dealloc_T*
Teuchos::get_optional_nonconst_dealloc( const Teuchos::ArrayRCP<T>& p )
{
  p.assert_not_null();
  typedef RCPNodeTmpl<typename Dealloc_T::ptr_t,Dealloc_T>
    RCPNT;
  RCPNT *dnode = dynamic_cast<RCPNT*>(p.access_node());
  if(dnode)
    return &dnode->get_nonconst_dealloc();
  return 0;
}


template<class TOrig, class Embedded, class T>
const Embedded& Teuchos::getEmbeddedObj( const ArrayRCP<T>& p )
{
  typedef EmbeddedObjDealloc<TOrig,Embedded,DeallocDelete<TOrig> > Dealloc_t;
  return get_dealloc<Dealloc_t>(p).getObj();
}


template<class TOrig, class Embedded, class T>
Embedded& Teuchos::getNonconstEmbeddedObj( const ArrayRCP<T>& p )
{
  typedef EmbeddedObjDealloc<TOrig,Embedded,DeallocDelete<TOrig> > Dealloc_t;
  return get_nonconst_dealloc<Dealloc_t>(p).getNonconstObj();
}


template<class T>
std::ostream& Teuchos::operator<<( std::ostream& out, const ArrayRCP<T>& p )
{
  out
    << TypeNameTraits<ArrayRCP<T> >::name() << "{"
    << "ptr="<<(const void*)(p.get()) // I can't find any alternative to this C cast :-(
    <<",lowerOffset="<<p.lowerOffset()
    <<",upperOffset="<<p.upperOffset()
    <<",size="<<p.size()
    <<",node="<<p.access_node()
    <<",count="<<p.count()
    <<"}";
  return out;
}


#endif // TEUCHOS_ARRAY_RCP_HPP
