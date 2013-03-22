// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_ARRAY_RCP_HPP
#define TEUCHOS_ARRAY_RCP_HPP


#include "Teuchos_ArrayRCPDecl.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_as.hpp"


namespace Teuchos {


// Helper code (not for general clients)


template<class T> inline
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


template<class T2, class T1>
class ArcpReinterpretCastEmbeddedObj
{
public:
  typedef T2 ptr_t;
  ArcpReinterpretCastEmbeddedObj()
    : arcp_pod_(null)
    {}
  ArcpReinterpretCastEmbeddedObj(const ArrayRCP<T1> &arcp_pod)
    : arcp_pod_(arcpCloneNode(arcp_pod)) // Unique reference count!
    {}
  // NOTE: The default copy constructor is allowed and does the right thing
  ~ArcpReinterpretCastEmbeddedObj()
    { freeMemory(); }
  ArcpReinterpretCastEmbeddedObj&
  operator=(const ArcpReinterpretCastEmbeddedObj& arceo)
    {
      assert(is_null(arceo.arcp_pod_)); // Can only be a catestrophic programming error!
      freeMemory();
      return *this;
    }
private:
  ArrayRCP<T1> arcp_pod_;
  void freeMemory()
    {
      typedef typename ArrayRCP<T2>::iterator itr_t;
      if (arcp_pod_.strong_count() == 1) {
        ArrayRCP<T2> arcp2 = arcp_reinterpret_cast<T2>(arcp_pod_);
        for (itr_t itr = arcp2.begin(); itr != arcp2.end(); ++itr) {
          itr->~T2();
        }
        arcp_pod_ = null;
      }
    }
};


// Constructors/Destructors/Initializers

template<class T> inline
ArrayRCP<T>::ArrayRCP( ENull )
  : ptr_(NULL), lowerOffset_(0), upperOffset_(-1)
{}


template<class T> inline
ArrayRCP<T>::ArrayRCP(size_type n, const T& val)
  : ptr_(0), lowerOffset_(0), upperOffset_(-1)
{
  *this = arcp<T>(n);
  std::fill_n(begin(), n, val);
}


template<class T> inline
ArrayRCP<T>::ArrayRCP(
  T* p, size_type lowerOffset_in, size_type size_in,
  bool has_ownership_in, const ERCPNodeLookup rcpNodeLookup
  )
  : ptr_(p),
#ifndef TEUCHOS_DEBUG
    node_(ArrayRCP_createNewRCPNodeRawPtr(p, has_ownership_in)),
#endif // TEUCHOS_DEBUG
    lowerOffset_(lowerOffset_in),
    upperOffset_(size_in + lowerOffset_in - 1)
{
#ifdef TEUCHOS_DEBUG
  if (p) {
    RCPNode* existing_RCPNode = 0;
    if (!has_ownership_in && rcpNodeLookup==RCP_ENABLE_NODE_LOOKUP) {
      existing_RCPNode = RCPNodeTracer::getExistingRCPNode(p);
    }
    if (existing_RCPNode) {
      // Will not call add_new_RCPNode(...)
      node_ = RCPNodeHandle(existing_RCPNode, RCP_WEAK, false);
    }
    else {
      // Will call add_new_RCPNode(...)
      RCPNodeThrowDeleter nodeDeleter(ArrayRCP_createNewRCPNodeRawPtr(p, has_ownership_in));
      node_ = RCPNodeHandle(
        nodeDeleter.get(),
        p, typeName(*p), concreteTypeName(*p),
        has_ownership_in
        );
      nodeDeleter.release();
    }
  }
#else // NOT TEUCHOS_DEBUG
  (void) rcpNodeLookup; // Silence "unused variable" compiler warning.
#endif // TEUCHOS_DEBUG
}


template<class T>
template<class Dealloc_T>
inline
ArrayRCP<T>::ArrayRCP(
  T* p, size_type lowerOffset_in, size_type size_in,
  Dealloc_T dealloc, bool has_ownership_in
  )
  : ptr_(p),
#ifndef TEUCHOS_DEBUG
    node_(ArrayRCP_createNewDeallocRCPNodeRawPtr(p, dealloc, has_ownership_in)),
#endif // TEUCHOS_DEBUG
    lowerOffset_(lowerOffset_in),
    upperOffset_(size_in + lowerOffset_in - 1)
{
#ifdef TEUCHOS_DEBUG
  if (p) {
    node_ = RCPNodeHandle(
      ArrayRCP_createNewDeallocRCPNodeRawPtr(p, dealloc, has_ownership_in),
      p, typeName(*p), concreteTypeName(*p),
      has_ownership_in
      //, RCP_STRONG, false
      );
  }
#endif // TEUCHOS_DEBUG
}


template<class T> inline
ArrayRCP<T>::ArrayRCP(const ArrayRCP<T>& r_ptr)
  :ptr_(r_ptr.ptr_),
   node_(r_ptr.node_),
   lowerOffset_(r_ptr.lowerOffset_),
   upperOffset_(r_ptr.upperOffset_)
{}


template<class T> inline
ArrayRCP<T>::~ArrayRCP()
{}


template<class T> inline
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


template<class T> inline
bool ArrayRCP<T>::is_null() const
{
  return ptr_ == 0;
}


template<class T> inline
T* ArrayRCP<T>::operator->() const
{
  debug_assert_valid_ptr();
  debug_assert_in_range(0,1);
  return ptr_;
}


template<class T> inline
T& ArrayRCP<T>::operator*() const
{
  debug_assert_valid_ptr();
  debug_assert_in_range(0,1);
  return *ptr_;
}


template<class T> inline
T* ArrayRCP<T>::get() const
{
  if(ptr_) {
    debug_assert_valid_ptr();
    debug_assert_in_range(0,1);
  }
  return ptr_;
}


template<class T> inline
T* ArrayRCP<T>::getRawPtr() const
{
  return this->get();
}


template<class T> inline
T& ArrayRCP<T>::operator[](size_type offset) const
{
  debug_assert_valid_ptr();
  debug_assert_in_range(offset,1);
  return ptr_[offset];
}


// Pointer Arithmetic Functions


template<class T> inline
ArrayRCP<T>& ArrayRCP<T>::operator++()
{
  debug_assert_valid_ptr();
  ++ptr_;
  --lowerOffset_;
  --upperOffset_;
  return *this;
}


template<class T> inline
ArrayRCP<T> ArrayRCP<T>::operator++(int)
{
  debug_assert_valid_ptr();
  ArrayRCP<T> r_ptr = *this;
  ++(*this);
  return r_ptr;
}


template<class T> inline
ArrayRCP<T>& ArrayRCP<T>::operator--()
{
  debug_assert_valid_ptr();
  --ptr_;
  ++lowerOffset_;
  ++upperOffset_;
  return *this;
}


template<class T> inline
ArrayRCP<T> ArrayRCP<T>::operator--(int)
{
  debug_assert_valid_ptr();
  ArrayRCP<T> r_ptr = *this;
  --(*this);
  return r_ptr;
}


template<class T> inline
ArrayRCP<T>& ArrayRCP<T>::operator+=(size_type offset)
{
  debug_assert_valid_ptr();
  ptr_ += offset;
  lowerOffset_ -= offset;
  upperOffset_ -= offset;
  return *this;
}


template<class T> inline
ArrayRCP<T>& ArrayRCP<T>::operator-=(size_type offset)
{
  debug_assert_valid_ptr();
  ptr_ -= offset;
  lowerOffset_ += offset;
  upperOffset_ += offset;
  return *this;
}


template<class T> inline
ArrayRCP<T> ArrayRCP<T>::operator+(size_type offset) const
{
  ArrayRCP<T> r_ptr = *this;
  r_ptr+=(offset);
  return r_ptr;
}


template<class T> inline
ArrayRCP<T> ArrayRCP<T>::operator-(size_type offset) const
{
  ArrayRCP<T> r_ptr = *this;
  r_ptr-=offset;
  return r_ptr;
}


// Standard Container-Like Functions


template<class T> inline
typename ArrayRCP<T>::iterator ArrayRCP<T>::begin() const
{
  debug_assert_valid_ptr();
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return *this;
#else
  return ptr_;
#endif
}


template<class T> inline
typename ArrayRCP<T>::iterator ArrayRCP<T>::end() const
{
  debug_assert_valid_ptr();
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return *this + (upperOffset_ + 1);
#else
  return ptr_ + (upperOffset_ + 1);
#endif
}


// ArrayRCP Views


template<class T> inline
ArrayRCP<const T> ArrayRCP<T>::getConst() const
{
  if (ptr_) {
    debug_assert_valid_ptr();
    const T *cptr = ptr_; // Will not compile if not legal!
    return ArrayRCP<const T>(cptr, lowerOffset_, size(), node_);
  }
  return null;
}


template<class T> inline
ArrayRCP<T>
ArrayRCP<T>::persistingView( size_type lowerOffset_in, size_type size_in ) const
{
  if (size_in == 0) {
    return null;
  }
  debug_assert_valid_ptr();
  debug_assert_in_range(lowerOffset_in, size_in);
  ArrayRCP<T> ptr = *this;
  ptr.ptr_ = ptr.ptr_ + lowerOffset_in;
  ptr.lowerOffset_ = 0;
  ptr.upperOffset_ = size_in - 1;
  return ptr;
}


// Size and extent query functions


template<class T> inline
typename ArrayRCP<T>::size_type
ArrayRCP<T>::lowerOffset() const
{
  debug_assert_valid_ptr();
  return lowerOffset_;
}


template<class T> inline
typename ArrayRCP<T>::size_type
ArrayRCP<T>::upperOffset() const
{
  debug_assert_valid_ptr();
  return upperOffset_;
}


template<class T> inline
typename ArrayRCP<T>::size_type
ArrayRCP<T>::size() const
{
  debug_assert_valid_ptr();
  return upperOffset_ - lowerOffset_ + 1;
}


// ArrayView views


template<class T> inline
ArrayView<T> ArrayRCP<T>::view( size_type lowerOffset_in, size_type size_in ) const
{
  if (size_in == 0) {
    return null;
  }
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
ArrayView<T> ArrayRCP<T>::operator()( size_type lowerOffset_in, size_type size_in ) const
{
  return view(lowerOffset_in, size_in);
}


template<class T> inline
ArrayView<T> ArrayRCP<T>::operator()() const
{
  if (size()) {
    return view(lowerOffset_, size());
  }
  return null;
}


// Implicit conversions


template<class T> inline
ArrayRCP<T>::operator ArrayRCP<const T>() const
{
  if (size()) {
    return ArrayRCP<const T>(ptr_, lowerOffset_, size(), node_);
  }
  return null;
}


// std::vector like functions


template<class T> inline
void ArrayRCP<T>::assign(size_type n, const T &val)
{
  *this = arcp<T>(n);
  std::fill_n(this->begin(), n, val);
}


template<class T>
template<class Iter>
inline
void ArrayRCP<T>::assign(Iter first, Iter last)
{
  const size_type new_n = std::distance(first, last);
  if (new_n != size())
    *this = arcp<T>(new_n);
  std::copy( first, last, begin() );
}


template<class T> inline
void ArrayRCP<T>::resize(const size_type n, const T &val)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(lowerOffset(), 0);
#endif
  if (n == 0) {
    clear();
    return;
  }
  const size_type orig_n = size();
  if (n != orig_n) {
    ArrayRCP<T> tmp = *this;
    *this = arcp<T>(n);
    const size_type small_n = std::min(n, orig_n);
    for (size_type i = 0; i < small_n; ++i)
      (*this)[i] = tmp[i];
    for (size_type i = orig_n; i < n; ++i)
      (*this)[i] = val;
    upperOffset_ = n-1;
  }
}


template<class T> inline
void ArrayRCP<T>::clear()
{
  *this = null;
}


// Misc functions

template<class T> inline
void ArrayRCP<T>::deepCopy(const ArrayView<const T>& av)
{
  if (av.size() == 0) {
    *this = null;
    return;
  }
  assign(av.begin(), av.end());
}


// Reference counting


template<class T> inline
ERCPStrength ArrayRCP<T>::strength() const
{
  return node_.strength();
}


template<class T> inline
bool ArrayRCP<T>::is_valid_ptr() const
{
  if (ptr_)
    return node_.is_valid_ptr();
  return true;
}


template<class T> inline
int ArrayRCP<T>::strong_count() const
{
  return node_.strong_count();
}


template<class T> inline
int ArrayRCP<T>::weak_count() const
{
  return node_.weak_count();
}


template<class T> inline
int ArrayRCP<T>::total_count() const
{
  return node_.total_count();
}


template<class T> inline
void ArrayRCP<T>::set_has_ownership()
{
  node_.has_ownership(true);
}


template<class T> inline
bool ArrayRCP<T>::has_ownership() const
{
  return node_.has_ownership();
}


template<class T> inline
T* ArrayRCP<T>::release()
{
  debug_assert_valid_ptr();
  node_.has_ownership(false);
  return ptr_;
}


template<class T> inline
ArrayRCP<T> ArrayRCP<T>::create_weak() const
{
  debug_assert_valid_ptr();
  return ArrayRCP<T>(ptr_, lowerOffset_, size(), node_.create_weak());
}


template<class T> inline
ArrayRCP<T> ArrayRCP<T>::create_strong() const
{
  debug_assert_valid_ptr();
  return ArrayRCP<T>(ptr_, lowerOffset_, size(), node_.create_strong());
}


template<class T>
template <class T2>
inline
bool ArrayRCP<T>::shares_resource(const ArrayRCP<T2>& r_ptr) const
{
  return node_.same_node(r_ptr.access_private_node());
  // Note: above, r_ptr is *not* the same class type as *this so we can not
  // access its node_ member directly!  This is an interesting detail to the
  // C++ protected/private protection mechanism!
}


// Assertion Functions


template<class T> inline
const ArrayRCP<T>&
ArrayRCP<T>::assert_not_null() const
{
  if(!ptr_)
    throw_null_ptr_error(typeName(*this));
  return *this;
}


template<class T> inline
const ArrayRCP<T>& ArrayRCP<T>::assert_valid_ptr() const
{
  if (ptr_)
    node_.assert_valid_ptr(*this);
  return *this;
}


template<class T> inline
const ArrayRCP<T>&
ArrayRCP<T>::assert_in_range( size_type lowerOffset_in, size_type size_in ) const
{
  assert_not_null();
  TEUCHOS_TEST_FOR_EXCEPTION(
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


template<class T> inline
int ArrayRCP<T>::count() const {
  return node_.count();
}


// very bad public functions


template<class T> inline
ArrayRCP<T>::ArrayRCP(
  T* p, size_type lowerOffset_in, size_type size_in,
  const RCPNodeHandle& node
  )
  :ptr_(p),
   node_(node),
   lowerOffset_(lowerOffset_in),
   upperOffset_(size_in + lowerOffset_in - 1)
{}


template<class T> inline
T* ArrayRCP<T>::access_private_ptr() const
{
  return ptr_;
}


template<class T> inline
RCPNodeHandle& ArrayRCP<T>::nonconst_access_private_node()
{
  return node_;
}


template<class T> inline
const RCPNodeHandle& ArrayRCP<T>::access_private_node() const
{
  return node_;
}

// Array<void> and Array<const void> specializations


ArrayRCP<void>::ArrayRCP()
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
}


ArrayRCP<const void>::ArrayRCP()
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
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
  TEUCHOS_TEST_FOR_EXCEPTION(
    !p1.shares_resource(p2), IncompatibleIteratorsError,
    "Error, these iterators are *not* pointing to the same valid memory!"
    );
#endif
}
} // namespace Utilities
} // namespace Teuchos


template<class T> inline
Teuchos::ArrayRCP<T>
Teuchos::arcp(
T* p, typename ArrayRCP<T>::size_type lowerOffset
  ,typename ArrayRCP<T>::size_type size_in
  ,bool owns_mem
  )
{
  return ArrayRCP<T>(p, lowerOffset, size_in, owns_mem);
}


template<class T, class Dealloc_T>
inline
Teuchos::ArrayRCP<T>
Teuchos::arcp(
T* p, typename ArrayRCP<T>::size_type lowerOffset
  ,typename ArrayRCP<T>::size_type size_in
  ,Dealloc_T dealloc, bool owns_mem
  )
{
  return ArrayRCP<T>(p, lowerOffset, size_in, dealloc, owns_mem);
}


template<class T> inline
Teuchos::ArrayRCP<T>
Teuchos::arcp( typename ArrayRCP<T>::size_type size )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_INEQUALITY( size, >=, 0 );
#endif
  if (size == 0) {
    return null;
  }
  return ArrayRCP<T>(new T[size], 0, size, true);
}


template<class T> inline
Teuchos::ArrayRCP<T>
Teuchos::arcpCloneNode(const ArrayRCP<T> &a)
{
  if (is_null(a)) {
    return null;
  }
  return arcpWithEmbeddedObj(a.getRawPtr(), a.lowerOffset(), a.size(),
    a, false);
}


template<class T> inline
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
  typename ArrayRCP<T>::size_type lowerOffset,
  typename ArrayRCP<T>::size_type size,
  const Embedded &embedded,
  bool owns_mem
  )
{
  return arcp(
    p, lowerOffset, size,
    embeddedObjDeallocArrayDelete<T>(embedded, PRE_DESTROY),
    owns_mem
    );
}


template<class T, class Embedded>
Teuchos::ArrayRCP<T>
Teuchos::arcpWithEmbeddedObjPostDestroy(
  T* p,
  typename ArrayRCP<T>::size_type lowerOffset,
  typename ArrayRCP<T>::size_type size,
  const Embedded &embedded,
  bool owns_mem
  )
{
  return arcp(
    p, lowerOffset, size,
    embeddedObjDeallocArrayDelete<T>(embedded, POST_DESTROY),
    owns_mem
    );
}


template<class T, class Embedded>
Teuchos::ArrayRCP<T>
Teuchos::arcpWithEmbeddedObj(
  T* p,
  typename ArrayRCP<T>::size_type lowerOffset,
  typename ArrayRCP<T>::size_type size,
  const Embedded &embedded,
  bool owns_mem
  )
{
  return arcpWithEmbeddedObjPostDestroy<T,Embedded>(
    p, lowerOffset, size, embedded, owns_mem );
}


template<class T> inline
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


template<class T> inline
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


template<class T> inline
Teuchos::ArrayRCP<T>
Teuchos::arcpFromArrayView(const ArrayView<T> &av)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return av.access_private_arcp();
#else
  return arcp(av.getRawPtr(), 0, av.size(), false);
#endif
}


template<class T> inline
Teuchos::RCP<std::vector<T> >
Teuchos::get_std_vector( const ArrayRCP<T> &ptr )
{
  return getEmbeddedObj<T, RCP<std::vector<T> > >(ptr);
}


template<class T> inline
Teuchos::RCP<const std::vector<T> >
Teuchos::get_std_vector( const ArrayRCP<const T> &ptr )
{
  return getEmbeddedObj<const T, RCP<const std::vector<T> > >(ptr);
}


template<class T> inline
bool Teuchos::is_null( const ArrayRCP<T> &p )
{
  return p.is_null();
}


template<class T> inline
bool Teuchos::nonnull( const ArrayRCP<T> &p )
{
  return !p.is_null();
}


template<class T> inline
bool Teuchos::operator==( const ArrayRCP<T> &p, ENull )
{
  return p.is_null();
}


template<class T> inline
bool Teuchos::operator!=( const ArrayRCP<T> &p, ENull )
{
  return !p.is_null();
}


template<class T1, class T2>
inline
bool Teuchos::operator==( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  return p1.access_private_ptr() == p2.access_private_ptr();
}


template<class T1, class T2>
inline
bool Teuchos::operator!=( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  return p1.access_private_ptr() != p2.access_private_ptr();
}


template<class T1, class T2>
inline
bool Teuchos::operator<( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  return p1.access_private_ptr() < p2.access_private_ptr();
}


template<class T1, class T2>
inline
bool Teuchos::operator<=( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  Utilities::assert_shares_resource(p1,p2);
  return p1.access_private_ptr() <= p2.access_private_ptr();
}


template<class T1, class T2>
inline
bool Teuchos::operator>( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 )
{
  Utilities::assert_shares_resource(p1,p2);
  return p1.access_private_ptr() > p2.access_private_ptr();
}


template<class T1, class T2>
inline
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
inline
Teuchos::ArrayRCP<T2>
Teuchos::arcp_reinterpret_cast(const ArrayRCP<T1>& p1)
{
  typedef typename ArrayRCP<T1>::size_type size_type;
  const int sizeOfT1 = sizeof(T1);
  const int sizeOfT2 = sizeof(T2);
  size_type lowerOffset2 = (p1.lowerOffset()*sizeOfT1) / sizeOfT2;
  size_type upperOffset2 = ((p1.upperOffset()+1)*sizeOfT1) / sizeOfT2 - 1;
  T2 *ptr2 = reinterpret_cast<T2*>(p1.get());
  return ArrayRCP<T2>(
    ptr2, lowerOffset2, upperOffset2 - lowerOffset2 + 1,
    p1.access_private_node()
    );
  // Note: Above is just fine even if p1.get()==NULL!
}


template<class T2, class T1>
Teuchos::ArrayRCP<T2>
Teuchos::arcp_reinterpret_cast_nonpod(const ArrayRCP<T1>& p1, const T2& val)
{
  typedef typename ArrayRCP<T2>::iterator itr_t;
  ArrayRCP<T2> arcp2 = arcp_reinterpret_cast<T2>(p1);
  for (itr_t itr = arcp2.begin(); itr != arcp2.end(); ++itr) {
    new (&*itr) T2(val);
  }
  return arcpWithEmbeddedObj(
    arcp2.getRawPtr(), 0, arcp2.size(),
    ArcpReinterpretCastEmbeddedObj<T2, T1>(p1),
    false);
  // Above, the ownership of the memory is totally owned by the embedded
  // object and the default deallocator policy object does not do anything.
  // This is just fine.
}


template<class T2, class T1>
inline
Teuchos::ArrayRCP<T2>
Teuchos::arcp_const_cast(const ArrayRCP<T1>& p1)
{
  typedef typename ArrayRCP<T1>::size_type size_type;
  T2 *ptr2 = const_cast<T2*>(p1.get());
  return ArrayRCP<T2>(
    ptr2, p1.lowerOffset(), p1.size(),
    p1.access_private_node()
    );
  // Note: Above is just fine even if p1.get()==NULL!
}


template<class T2, class T1>
inline
Teuchos::ArrayRCP<T2>
Teuchos::arcp_implicit_cast(const ArrayRCP<T1>& p1)
{
  typedef typename ArrayRCP<T1>::size_type size_type;
  T2 * raw_ptr2 = p1.get();
  return ArrayRCP<T2>(
    raw_ptr2, p1.lowerOffset(), p1.size(),
    p1.access_private_node()
    );
  // Note: Above is just fine even if p1.get()==NULL!
}


template<class T1, class T2>
inline
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
inline
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
inline
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
inline
T1* Teuchos::get_optional_extra_data( ArrayRCP<T2>& p, const std::string& name )
{
  p.assert_not_null();
  any *extra_data = p.nonconst_access_private_node().get_optional_extra_data(
    TypeNameTraits<T1>::name(), name);
  if( extra_data ) return &any_cast<T1>(*extra_data);
  return NULL;
}


template<class T1, class T2>
inline
const T1* Teuchos::get_optional_extra_data( const ArrayRCP<T2>& p, const std::string& name )
{
  p.assert_not_null();
  any *extra_data = p.access_private_node().get_optional_extra_data(
    TypeNameTraits<T1>::name(), name);
  if( extra_data ) return &any_cast<T1>(*extra_data);
  return NULL;
}


template<class Dealloc_T, class T>
inline
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
  TEUCHOS_TEST_FOR_EXCEPTION(
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
inline
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
  typedef EmbeddedObjDealloc<TOrig,Embedded,DeallocArrayDelete<TOrig> > Dealloc_t;
  return get_dealloc<Dealloc_t>(p).getObj();
}


template<class TOrig, class Embedded, class T>
Embedded& Teuchos::getNonconstEmbeddedObj( const ArrayRCP<T>& p )
{
  typedef EmbeddedObjDealloc<TOrig,Embedded,DeallocArrayDelete<TOrig> > Dealloc_t;
  return get_nonconst_dealloc<Dealloc_t>(p).getNonconstObj();
}


template<class T>
std::ostream& Teuchos::operator<<( std::ostream& out, const ArrayRCP<T>& p )
{
  out
    << TypeNameTraits<ArrayRCP<T> >::name() << "{"
    << "ptr="<<(const void*)(p.access_private_ptr())
    <<",lowerOffset="<<p.lowerOffset()
    <<",upperOffset="<<p.upperOffset()
    <<",size="<<p.size()
    <<",node="<<p.access_private_node()
    <<",strong_count="<<p.strong_count()
    <<",weak_count="<<p.weak_count()
    <<"}";
  return out;
  // NOTES:
  // * I can't find any alternative to this C cast (problems with char data)
  // * Don't range check the pointer since this code does not dereference it.
  //   This is needed to allow printing the end() or past end() for debugging.
}


#endif // TEUCHOS_ARRAY_RCP_HPP
