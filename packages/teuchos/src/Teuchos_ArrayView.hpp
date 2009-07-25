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

#ifndef TEUCHOS_ARRAY_VIEW_HPP
#define TEUCHOS_ARRAY_VIEW_HPP


#include "Teuchos_ArrayViewDecl.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_as.hpp"


namespace Teuchos {


// Constructors/Destructors


template<class T> inline
ArrayView<T>::ArrayView( ENull )
  :ptr_(0), size_(0)
{
  setUpIterators();
}


template<class T> inline
ArrayView<T>::ArrayView( T* p, Ordinal size_in )
  :ptr_(p), size_(size_in)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_FOR_EXCEPT( p != 0 && size_in <= 0 );
  TEST_FOR_EXCEPT( p == 0 && size_in != 0 );
  setUpIterators();
#endif
}


template<class T> inline
ArrayView<T>::ArrayView(const ArrayView<T>& array)
  :ptr_(array.ptr_), size_(array.size_)
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  ,arcp_(array.arcp_)
#endif
{}


template<class T> inline
ArrayView<T>::ArrayView(
  std::vector<typename ConstTypeTraits<T>::NonConstType>& vec
  )
  : ptr_( vec.empty() ? 0 : &vec[0] ), size_(vec.size())
{
  setUpIterators();
}


template<class T> inline
ArrayView<T>::ArrayView(
  const std::vector<typename ConstTypeTraits<T>::NonConstType>& vec
  )
  : ptr_( vec.empty() ? 0 : &vec[0] ), size_(vec.size())
{
  setUpIterators();
}


template<class T> inline
ArrayView<T>& ArrayView<T>::operator=(const ArrayView<T>& array)
{
  ptr_ = array.ptr_;
  size_ = array.size_;
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  arcp_ = array.arcp_;
#endif
  return *this;
}


template<class T> inline
ArrayView<T>::~ArrayView()
{}


// General query functions 


template<class T>
inline
bool ArrayView<T>::is_null() const
{
  return ptr_ == 0;
}


template<class T> inline
typename ArrayView<T>::Ordinal ArrayView<T>::size() const
{
  debug_assert_valid_ptr();
  return size_;
}


template<typename T>
std::string ArrayView<T>::toString() const
{

  using Teuchos::as;
  std::ostringstream ss;

  debug_assert_valid_ptr();

  ss << "{";

  for (int i=0; i < as<int>(size()); ++i)
  {
    ss << operator[](i);
    if (i < size()-1) ss << ", ";
  }
  ss << "}";

  return ss.str();

}


// Element Access Functions


template<class T> inline
T* ArrayView<T>::getRawPtr() const
{
  debug_assert_valid_ptr();
  return ptr_;
}


template<class T> inline
T& ArrayView<T>::operator[](Ordinal i) const
{
  debug_assert_valid_ptr();
  debug_assert_in_range(i,1);
  return ptr_[i];
}


template<class T> inline
T& ArrayView<T>::front() const
{
  debug_assert_not_null();
  debug_assert_valid_ptr();
  return *ptr_;
}


template<class T> inline
T& ArrayView<T>::back() const
{
  debug_assert_not_null();
  debug_assert_valid_ptr();
  return *(ptr_+size_-1);
}


// Views 


template<class T> inline
ArrayView<T> ArrayView<T>::view(Ordinal offset, Ordinal size_in) const
{
  debug_assert_valid_ptr();
  debug_assert_in_range(offset, size_in);
  return ArrayView<T>(ptr_+offset,size_in);
  // WARNING: The above code had better be correct since we are using raw
  // pointer arithmetic!
}


template<class T> inline
ArrayView<T> ArrayView<T>::operator()( Ordinal offset, Ordinal size_in ) const
{
  return view(offset, size_in);
}


template<class T> inline
const ArrayView<T>& ArrayView<T>::operator()() const
{
  debug_assert_valid_ptr();
  return *this;
}


template<class T> inline
ArrayView<const T> ArrayView<T>::getConst() const
{
  debug_assert_valid_ptr();
  return ArrayView<const T>(ptr_,size_);
}


template<class T> inline
ArrayView<T>::operator ArrayView<const T>() const
{
  return getConst();
}


// Assignment


template<class T>
void ArrayView<T>::assign(const ArrayView<const T>& array) const
{
  debug_assert_valid_ptr();
  debug_assert_not_null();
  if (this->getRawPtr()==array.getRawPtr() && this->size()==array.size())
    return; // Assignment to self
  debug_assert_in_range(0,array.size());
  std::copy( array.begin(), array.end(), this->begin() );
  // Note: Above, in debug mode, the iterators are range checked!  In
  // optimized mode, these are raw pointers which should run very fast!
}


// Standard Container-Like Functions 


template<class T>
typename ArrayView<T>::iterator ArrayView<T>::begin() const
{
  debug_assert_valid_ptr();
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return arcp_.create_weak();
#else
  return ptr_;
#endif
}


template<class T>
typename ArrayView<T>::iterator ArrayView<T>::end() const
{
  debug_assert_valid_ptr();
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return arcp_.create_weak() + size_;
#else
  return ptr_ + size_;
#endif
}


// Assertion Functions. 


template<class T>
const ArrayView<T>& ArrayView<T>::assert_not_null() const
{
  if(!ptr_)
    throw_null_ptr_error(typeName(*this));
  return *this;
}


template<class T>
const ArrayView<T>&
ArrayView<T>::assert_in_range( Ordinal offset, Ordinal size_in ) const
{
  
  assert_not_null();
  TEST_FOR_EXCEPTION(
    !(
      ( 0 <= offset && offset+size_in <= this->size() )
      &&
      size_in >= 0
      ),
    Teuchos::RangeError,
    typeName(*this)<<"::assert_in_range():"
    " Error, [offset,offset+size) = ["<<offset<<","<<(offset+size_in)<<")"
    " does not lie in the range [0,"<<this->size()<<")!"
    );
  return*this;
}


#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

template<class T>
ArrayView<T>::ArrayView( const ArrayRCP<T> &arcp )
  : ptr_(arcp.getRawPtr()), size_(arcp.size()), arcp_(arcp)
{}


template<class T>
ArrayView<T>::ArrayView( T* p, Ordinal size_in, const ArrayRCP<T> &arcp )
  : ptr_(p), size_(size_in), arcp_(arcp)
{}


#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


// private


template<class T>
void ArrayView<T>::setUpIterators()
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  if (ptr_)
    arcp_ = arcp(ptr_, 0, size_, false);
  // 2008/09/16: rabartl: Above, this will catch improper usage of the array
  // for out-of-bounds errors but it will not catch dangling references.
#endif
}


} // namespace Teuchos


//
// Nonmember helper functions
//


template<class T> inline
Teuchos::ArrayView<T>
Teuchos::arrayView( T* p, typename ArrayView<T>::Ordinal size )
{
  return ArrayView<T>(p,size);
}


template<class T> inline
Teuchos::ArrayView<T> Teuchos::arrayViewFromVector( std::vector<T>& vec )
{
  return ArrayView<T>(vec);
}


template<class T> inline
Teuchos::ArrayView<const T> Teuchos::arrayViewFromVector( const std::vector<T>& vec )
{
  return ArrayView<const T>(vec);
}


#ifndef __sun

template<class T> inline
std::vector<T> Teuchos::createVector( const ArrayView<T> &av )
{
  std::vector<T> v(av.begin(), av.end());
  return v;
}

#endif // __sun


template<class T> inline
std::vector<T> Teuchos::createVector( const ArrayView<const T> &av )
{
  std::vector<T> v(av.begin(), av.end());
  return v;
}


template<class T> inline
bool Teuchos::is_null( const ArrayView<T> &av )
{
  return av.is_null();
}


template<class T>
std::ostream& Teuchos::operator<<( std::ostream& out, const ArrayView<T>& p )
{
  return out << p.toString();
}


template<class T2, class T1>
REFCOUNTPTR_INLINE
Teuchos::ArrayView<T2>
Teuchos::av_reinterpret_cast(const ArrayView<T1>& p1)
{
  typedef typename ArrayView<T1>::Ordinal Ordinal;
  const int sizeOfT1 = sizeof(T1);
  const int sizeOfT2 = sizeof(T2);
  Ordinal size2 = (p1.size()*sizeOfT1) / sizeOfT2;
  T2 *ptr2 = reinterpret_cast<T2*>(p1.getRawPtr());
  return ArrayView<T2>(
    ptr2, size2
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    ,arcp_reinterpret_cast<T2>(p1.access_private_arcp())
#endif
    );
  // Note: Above is just fine even if p1.get()==NULL!
}


#endif	// TEUCHOS_ARRAY_VIEW_HPP
