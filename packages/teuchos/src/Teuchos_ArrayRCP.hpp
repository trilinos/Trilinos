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
#include "Teuchos_as.hpp"


namespace Teuchos {


// Helper functions


template<class T>
inline
RCPNode* ArrayRCP_createNewRCPNodeRawPtr( T* p, bool has_ownership_in )
{
  return new RCPNodeTmpl<T,DeallocArrayDelete<T> >(
    p, DeallocArrayDelete<T>(), has_ownership_in
    );
}


template<class T, class Dealloc_T>
inline
RCPNode* ArrayRCP_createNewDeallocRCPNodeRawPtr(
  T* p, Dealloc_T dealloc, bool has_ownership_in
  )
{
    return new RCPNodeTmpl<T,Dealloc_T>(p, dealloc, has_ownership_in);
}


// Constructors/Destructors/Initializers 


template<class T>
inline
ArrayRCP<T>::ArrayRCP( ENull )
  : ptr_(NULL), lowerOffset_(0), upperOffset_(-1)
{}


template<class T>
inline
ArrayRCP<T>::ArrayRCP(Ordinal n, const T& val)
  : ptr_(0), lowerOffset_(0), upperOffset_(-1)
{
  *this = arcp<T>(n);
  std::fill_n(begin(), n, val);
}


template<class T>
inline
ArrayRCP<T>::ArrayRCP(
  T* p, Ordinal lowerOffset_in, Ordinal upperOffset_in, bool has_ownership_in
  )
  : ptr_(p),
#ifndef TEUCHOS_DEBUG
    node_(ArrayRCP_createNewRCPNodeRawPtr(p, has_ownership_in)),
#endif // TEUCHOS_DEBUG
    lowerOffset_(lowerOffset_in),
    upperOffset_(upperOffset_in)
{
#ifdef TEUCHOS_DEBUG
  if (p) {
    node_ = RCPNodeHandle(
      ArrayRCP_createNewRCPNodeRawPtr(p, has_ownership_in),
      p, typeName(*p), concreteTypeName(*p),
      has_ownership_in
      );
  }
#endif // TEUCHOS_DEBUG
}


template<class T>
REFCOUNTPTR_INLINE
template<class Dealloc_T>
ArrayRCP<T>::ArrayRCP(
  T* p, Ordinal lowerOffset_in, Ordinal upperOffset_in,
  Dealloc_T dealloc, bool has_ownership_in
  )
  : ptr_(p),
#ifndef TEUCHOS_DEBUG
    node_(ArrayRCP_createNewDeallocRCPNodeRawPtr(p, dealloc, has_ownership_in)),
#endif // TEUCHOS_DEBUG
    lowerOffset_(lowerOffset_in),
    upperOffset_(upperOffset_in)
{
#ifdef TEUCHOS_DEBUG
  if (p) {
    node_ = RCPNodeHandle(
      ArrayRCP_createNewDeallocRCPNodeRawPtr(p, dealloc, has_ownership_in),
      p, typeName(*p), concreteTypeName(*p),
      has_ownership_in
      );
  }
#endif // TEUCHOS_DEBUG
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>::ArrayRCP(const ArrayRCP<T>& r_ptr)
  :ptr_(r_ptr.ptr_),
   node_(r_ptr.node_),
   lowerOffset_(r_ptr.lowerOffset_),
   upperOffset_(r_ptr.upperOffset_)
{}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>::~ArrayRCP()
{}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>& ArrayRCP<T>::operator=(const ArrayRCP<T>& r_ptr)
{
  if( this == &r_ptr )
    return *this; // Assignment to self
  node_ = r_ptr.access_private_node(); // May throw in debug mode!
  ptr_ = r_ptr.ptr_;
  lowerOffset_ = r_ptr.lowerOffset_;
  upperOffset_ = r_ptr.upperOffset_;
  return *this;
  // NOTE: It is critical that the assignment of ptr_ come *after* the
  // assignment of node_ since node_ might throw an exception!
}


// Object/Pointer Access Functions


template<class T>
inline
bool ArrayRCP<T>::is_null() const
{
  return ptr_ == 0;
}


template<class T>
inline
T* ArrayRCP<T>::operator->() const
{
  debug_assert_valid_ptr();
  debug_assert_in_range(0,1);
  return ptr_;
}


template<class T>
inline
T& ArrayRCP<T>::operator*() const
{
  debug_assert_valid_ptr();
  debug_assert_in_range(0,1);
  return *ptr_;
}


template<class T>
inline
T* ArrayRCP<T>::get() const
{
  if(ptr_) {
    debug_assert_valid_ptr();
    debug_assert_in_range(0,1);
  }
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
  debug_assert_valid_ptr();
  debug_assert_in_range(offset,1);
  return ptr_[offset];
}


// Pointer Arithmetic Functions


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>& ArrayRCP<T>::operator++()
{
  if(ptr_) {
    debug_assert_valid_ptr();
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
  debug_assert_valid_ptr();
  ArrayRCP<T> r_ptr = *this;
  ++(*this);
  return r_ptr;
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>& ArrayRCP<T>::operator--()
{
  if(ptr_) {
    debug_assert_valid_ptr();
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
  debug_assert_valid_ptr();
  ArrayRCP<T> r_ptr = *this;
  --(*this);
  return r_ptr;
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>& ArrayRCP<T>::operator+=(Ordinal offset)
{
  if(ptr_) {
    debug_assert_valid_ptr();
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
    debug_assert_valid_ptr();
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


// Standard Container-Like Functions


template<class T>
REFCOUNTPTR_INLINE
typename ArrayRCP<T>::const_iterator ArrayRCP<T>::begin() const
{
  debug_assert_valid_ptr();
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
  debug_assert_valid_ptr();
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return *this + (upperOffset_ + 1);
#else
  return ptr_ + (upperOffset_ + 1);
#endif
}


// ArrayRCP Views 


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<const T> ArrayRCP<T>::getConst() const
{
  if (ptr_) {
    debug_assert_valid_ptr();
    const T *cptr = ptr_; // Will not compile if not legal!
    return ArrayRCP<const T>(cptr, lowerOffset_, upperOffset_, node_);
  }
  return null;
}


template<class T>
REFCOUNTPTR_INLINE
ArrayRCP<T>
ArrayRCP<T>::persistingView( Ordinal lowerOffset_in, Ordinal size_in ) const
{
  debug_assert_valid_ptr();
  debug_assert_in_range(lowerOffset_in,size_in);
  ArrayRCP<T> ptr = *this;
  ptr.ptr_ = ptr.ptr_ + lowerOffset_in;
  ptr.lowerOffset_ = 0;
  ptr.upperOffset_ = size_in-1;
  return ptr;
}


// Size and extent query functions


template<class T>
REFCOUNTPTR_INLINE
typename ArrayRCP<T>::Ordinal
ArrayRCP<T>::lowerOffset() const
{
  debug_assert_valid_ptr();
  return lowerOffset_;
}


template<class T>
REFCOUNTPTR_INLINE
typename ArrayRCP<T>::Ordinal
ArrayRCP<T>::upperOffset() const
{
  debug_assert_valid_ptr();
  return upperOffset_;
}


template<class T>
REFCOUNTPTR_INLINE
typename ArrayRCP<T>::Ordinal
ArrayRCP<T>::size() const
{
  debug_assert_valid_ptr();
  return upperOffset_-lowerOffset_+1;
}


// ArrayView views 


template<class T> inline
ArrayView<T> ArrayRCP<T>::view( Ordinal lowerOffset_in, Ordinal size_in ) const
{
  debug_assert_valid_ptr();
  debug_assert_in_range(lowerOffset_in,size_in);
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return ArrayView<T>(persistingView(lowerOffset_in, size_in).create_weak());
#else
  return arrayView(ptr_ + lowerOffset_in, size_in);
#endif
  // ToDo: Implement checks for dangling references!
}


template<class T> inline
ArrayView<T> ArrayRCP<T>::operator()( Ordinal lowerOffset_in, Ordinal size_in ) const
{
  return view(lowerOffset_in,size_in);
}


template<class T> inline
ArrayView<T> ArrayRCP<T>::operator()() const
{
  if (size())
    return view(lowerOffset_, upperOffset_-lowerOffset_+1);
  return null;
}


// Implicit conversions


template<class T> inline
ArrayRCP<T>::operator ArrayView<T>() const
{
  return this->operator()();
}


template<class T> inline
ArrayRCP<T>::operator ArrayRCP<const T>() const
{
  if (size())
    return ArrayRCP<const T>(ptr_, lowerOffset_, upperOffset_, node_);
  return null;
}


// std::vector like functions


template<class T>
inline
void ArrayRCP<T>::assign(Ordinal n, const T &val)
{
  *this = arcp<T>(n);
  std::fill_n(this->begin(), n, val);
}


template<class T>
template<class Iter>
inline
void ArrayRCP<T>::assign(Iter first, Iter last)
{
  const Ordinal new_n = std::distance(first, last);
  if (new_n != size())
    *this = arcp<T>(new_n);
  std::copy( first, last, begin() );
}


template<class T>
inline
void ArrayRCP<T>::resize(const Ordinal n, const T &val)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(lowerOffset(), 0);
#endif
  if (n == 0) {
    clear();
    return;
  }
  const Ordinal orig_n = size();
  if (n != orig_n) {
    ArrayRCP<T> tmp = *this;
    *this = arcp<T>(n);
    const Ordinal small_n = std::min(n, orig_n);
    for (Ordinal i = 0; i < small_n; ++i)
      (*this)[i] = tmp[i];
    for (Ordinal i = orig_n; i < n; ++i)
      (*this)[i] = val;
    upperOffset_ = n-1;
  }
}


template<class T>
inline
void ArrayRCP<T>::clear()
{
  *this = null;
}


// Misc functions

template<class T>
inline
void ArrayRCP<T>::deepCopy(const ArrayView<const T>& av)
{
  if (av.size() == 0) {
    *this = null;
    return;
  }
  assign(av.begin(), av.end());
}


// Reference counting


template<class T>
inline
ERCPStrength ArrayRCP<T>::strength() const
{
  return node_.strength();
}


template<class T>
inline
bool ArrayRCP<T>::is_valid_ptr() const
{
  if (ptr_)
    return node_.is_valid_ptr();
  return true;
}


template<class T>
inline
int ArrayRCP<T>::strong_count() const
{
  return node_.strong_count();
}


template<class T>
inline
int ArrayRCP<T>::weak_count() const
{
  return node_.weak_count();
}


template<class T>
inline
int ArrayRCP<T>::total_count() const
{
  return node_.total_count();
}


template<class T>
REFCOUNTPTR_INLINE
void ArrayRCP<T>::set_has_ownership()
{
  node_.has_ownership(true);
}


template<class T>
REFCOUNTPTR_INLINE
bool ArrayRCP<T>::has_ownership() const
{
  return node_.has_ownership();
}


template<class T>
REFCOUNTPTR_INLINE
T* ArrayRCP<T>::release()
{
  debug_assert_valid_ptr();
  node_.has_ownership(false);
  return ptr_;
}


template<class T>
inline
ArrayRCP<T> ArrayRCP<T>::create_weak() const
{
  debug_assert_valid_ptr();
  return ArrayRCP<T>(ptr_, lowerOffset_, upperOffset_, node_.create_weak());
}


template<class T>
inline
ArrayRCP<T> ArrayRCP<T>::create_strong() const
{
  debug_assert_valid_ptr();
  return ArrayRCP<T>(ptr_, lowerOffset_, upperOffset_, node_.create_strong());
}


template<class T>
REFCOUNTPTR_INLINE
template <class T2>
bool ArrayRCP<T>::shares_resource(const ArrayRCP<T2>& r_ptr) const
{
  return node_.same_node(r_ptr.access_private_node());
  // Note: above, r_ptr is *not* the same class type as *this so we can not
  // access its node_ member directly!  This is an interesting detail to the
  // C++ protected/private protection mechanism!
}


// Assertion Functions


template<class T>
inline
const ArrayRCP<T>&
ArrayRCP<T>::assert_not_null() const
{
  if(!ptr_)
    throw_null_ptr_error(typeName(*this));
  return *this;
}


template<class T>
inline
const ArrayRCP<T>& ArrayRCP<T>::assert_valid_ptr() const
{
  if (ptr_)
    node_.assert_valid_ptr(*this);
  return *this;
}


template<class T>
inline
const ArrayRCP<T>&
ArrayRCP<T>::assert_in_range( Ordinal lowerOffset_in, Ordinal size_in ) const
{
  assert_not_null();
  TEST_FOR_EXCEPTION(
    !(
      (lowerOffset_ <= lowerOffset_in && lowerOffset_in+size_in-1 <= upperOffset_)
      &&
      size_in >= 0
      ),
    Teuchos::RangeError,
    typeName(*this)<<"::assert_in_range:"
    " Error, [lowerOffset,lowerOffset+size-1] = ["
    <<lowerOffset_in<<","<<(lowerOffset_in+size_in-1)<<"] does not lie in the"
    " range ["<<lowerOffset_<<","<<upperOffset_<<"]!"
    );
  return *this;
}


// Deprecated


template<class T>
REFCOUNTPTR_INLINE
int ArrayRCP<T>::count() const {
  return node_.count();
}


// very bad public functions


template<class T>
inline
ArrayRCP<T>::ArrayRCP(
  T* p, Ordinal lowerOffset_in, Ordinal upperOffset_in,
  const RCPNodeHandle& node
  )
  :ptr_(p),
   node_(node),
   lowerOffset_(lowerOffset_in),
   upperOffset_(upperOffset_in)
{}


template<class T>
inline
T* ArrayRCP<T>::access_private_ptr() const
{
  return ptr_;
}


template<class T>
inline
RCPNodeHandle& ArrayRCP<T>::nonconst_access_private_node()
{
  return node_;
}


template<class T>
inline
const RCPNodeHandle& ArrayRCP<T>::access_private_node() const
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
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_INEQUALITY( size, >, 0 );
#endif
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
  return p.is_null();
}


template<class T>
REFCOUNTPTR_INLINE
bool Teuchos::nonnull( const ArrayRCP<T> &p )
{
  return !p.is_null();
}


template<class T>
REFCOUNTPTR_INLINE
bool Teuchos::operator==( const ArrayRCP<T> &p, ENull )
{
  return p.is_null();
}


template<class T>
REFCOUNTPTR_INLINE
bool Teuchos::operator!=( const ArrayRCP<T> &p, ENull )
{
  return !p.is_null();
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator==( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  return p1.access_private_ptr() == p2.access_private_ptr();
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator!=( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  return p1.access_private_ptr() != p2.access_private_ptr();
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator<( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  return p1.access_private_ptr() < p2.access_private_ptr();
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator<=( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  Utilities::assert_shares_resource(p1,p2);
  return p1.access_private_ptr() <= p2.access_private_ptr();
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator>( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  Utilities::assert_shares_resource(p1,p2);
  return p1.access_private_ptr() > p2.access_private_ptr();
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator>=( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  Utilities::assert_shares_resource(p1,p2);
  return p1.access_private_ptr() >= p2.access_private_ptr();
}


template<class T>
typename Teuchos::ArrayRCP<T>::difference_type
Teuchos::operator-( const ArrayRCP<T> &p1, const ArrayRCP<T> &p2 )
{
  Utilities::assert_shares_resource(p1,p2);
  return p1.access_private_ptr() - p2.access_private_ptr();
}


template<class T2, class T1>
REFCOUNTPTR_INLINE
Teuchos::ArrayRCP<T2>
Teuchos::arcp_reinterpret_cast(const ArrayRCP<T1>& p1)
{
  typedef typename ArrayRCP<T1>::Ordinal Ordinal;
  const int sizeOfT1 = sizeof(T1);
  const int sizeOfT2 = sizeof(T2);
  Ordinal lowerOffset2 = (p1.lowerOffset()*sizeOfT1) / sizeOfT2;
  Ordinal upperOffset2 = ((p1.upperOffset()+1)*sizeOfT1) / sizeOfT2 - 1;
  T2 *ptr2 = reinterpret_cast<T2*>(p1.get());
  return ArrayRCP<T2>(
    ptr2, lowerOffset2, upperOffset2,
    p1.access_private_node()
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
    p1.access_private_node()
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
    raw_ptr2,p1.lowerOffset(),p1.upperOffset(),
    p1.access_private_node()
    );
  // Note: Above is just fine even if p1.get()==NULL!
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
void Teuchos::set_extra_data(
  const T1 &extra_data, const std::string& name,
  const Ptr<ArrayRCP<T2> > &p, EPrePostDestruction destroy_when,
  bool force_unique
  )
{
  p->assert_not_null();
  p->nonconst_access_private_node().set_extra_data( any(extra_data), name, destroy_when,
    force_unique );
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
T1& Teuchos::get_extra_data( ArrayRCP<T2>& p, const std::string& name )
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
const T1& Teuchos::get_extra_data( const ArrayRCP<T2>& p, const std::string& name )
{
  p.assert_not_null();
  return any_cast<T1>(
    p.access_private_node().get_extra_data(
      TypeNameTraits<T1>::name() ,name
      )
    );
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
T1* Teuchos::get_optional_extra_data( ArrayRCP<T2>& p, const std::string& name )
{
  p.assert_not_null();
  any *extra_data = p.nonconst_access_private_node().get_optional_extra_data(
    TypeNameTraits<T1>::name(), name);
  if( extra_data ) return &any_cast<T1>(*extra_data);
  return NULL;
}


template<class T1, class T2>
REFCOUNTPTR_INLINE
const T1* Teuchos::get_optional_extra_data( const ArrayRCP<T2>& p, const std::string& name )
{
  p.assert_not_null();
  any *extra_data = p.access_private_node().get_optional_extra_data(
    TypeNameTraits<T1>::name(), name);
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
  RCPNT *dnode = dynamic_cast<RCPNT*>(p.access_private_node().node_ptr());
  if (dnode)
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
    <<",node="<<p.access_private_node()
    <<",count="<<p.count()
    <<"}";
  return out;
}


#endif // TEUCHOS_ARRAY_RCP_HPP
