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
ArrayView<T>::ArrayView( ENull /*null_arg*/ )
  :ptr_(0), size_(0)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  setUpIterators();
#endif
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
  ,arcp_(array.arcp_) // It is okay to share the same iterator implementation!
#endif
{}


template<class T> inline
ArrayView<T>::ArrayView(
  std::vector<typename ConstTypeTraits<T>::NonConstType>& vec
  )
  : ptr_( vec.empty() ? 0 : &vec[0] ), size_(vec.size())
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  setUpIterators();
#endif
}


template<class T> inline
ArrayView<T>::ArrayView(
  const std::vector<typename ConstTypeTraits<T>::NonConstType>& vec
  )
  : ptr_( vec.empty() ? 0 : &vec[0] ), size_(vec.size())
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  setUpIterators();
#endif
}


template<class T> inline
ArrayView<T>::~ArrayView()
{}


// General query functions 


template<class T> inline
typename ArrayView<T>::Ordinal ArrayView<T>::size() const
{
  return size_;
}


template<typename T>
std::string ArrayView<T>::toString() const
{
  using Teuchos::as;
  std::ostringstream ss;
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
  return ptr_;
}


template<class T> inline
T& ArrayView<T>::operator[](Ordinal i) const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assert_in_range(i,1);
#endif
  return ptr_[i];
}


template<class T> inline
T& ArrayView<T>::front() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assert_not_null();
#endif
  return *ptr_;
}


template<class T> inline
T& ArrayView<T>::back() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assert_not_null();
#endif
  return *(ptr_+size_-1);
}


// Views 


template<class T> inline
ArrayView<T> ArrayView<T>::view(Ordinal offset, Ordinal size_in) const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assert_in_range(offset,size_in);
#endif
  return ArrayView<T>(ptr_+offset,size_in);
  // WARNING: The above code had better be correct since we are using raw
  // pointer arithmetic!
}


template<class T> inline
ArrayView<T> ArrayView<T>::operator()( Ordinal offset, Ordinal size_in ) const
{
  return view(offset,size_in);
}


template<class T> inline
const ArrayView<T>& ArrayView<T>::operator()() const
{
  return *this;
}


template<class T> inline
ArrayView<const T> ArrayView<T>::getConst() const
{
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
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assert_not_null();
#endif
  if (this->getRawPtr()==array.getRawPtr() && this->size()==array.size())
    return; // Assignment to self
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assert_in_range(0,array.size());
#endif
  std::copy( array.begin(), array.end(), this->begin() );
  // Note: Above, in debug mode, the iterators are range checked!  In
  // optimized mode, these are raw pointers which should run very fast!
}

template<class T>
ArrayView<T>& ArrayView<T>::operator=(const ArrayView<T>& array) 
{
  ptr_ = array.ptr_;
  size_ = array.size_;
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  arcp_ = array.arcp_;
#endif
  return *this;
}

// Standard Container-Like Functions 


template<class T>
typename ArrayView<T>::iterator ArrayView<T>::begin() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return arcp_;
#else
  return ptr_;
#endif
}


template<class T>
typename ArrayView<T>::iterator ArrayView<T>::end() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return arcp_ + size_;
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
    !( 0 <= offset && offset+size_in <= this->size() ), Teuchos::RangeError,
    typeName(*this)<<"::assert_in_range():"
    " Error, [offset,offset+size) = ["<<offset<<","<<(offset+size_in)<<")"
    " does not lie in the range [0,"<<this->size()<<")!"
    );
  return*this;
}


// private


#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

template<class T>
void ArrayView<T>::setUpIterators()
{
  if (ptr_)
    arcp_ = arcp(ptr_,0,size_,false);
}

#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


template<class T> inline
void* ArrayView<T>::operator new(size_t)
{
  // Should never be called!
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_FOR_EXCEPT(true);
#endif
  return 0;
}


#ifndef TEUCHOS_PRIVIATE_DELETE_NOT_SUPPORTED
template<class T> inline
void ArrayView<T>::operator delete(void*)
{
  // Should never be called!
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_FOR_EXCEPT(true);
#endif
}
#endif // TEUCHOS_PRIVIATE_DELETE_NOT_SUPPORTED


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
std::vector<T> Teuchos::createVector( const ArrayView<T> &ptr )
{
  std::vector<T> v(ptr.begin(),ptr.end());
  return v;
}

#endif // __sun


template<class T> inline
std::vector<T> Teuchos::createVector( const ArrayView<const T> &ptr )
{
  std::vector<T> v(ptr.begin(),ptr.end());
  return v;
}


template<class T>
std::ostream& Teuchos::operator<<( std::ostream& out, const ArrayView<T>& p )
{
  return out << p.toString();
}


#endif	// TEUCHOS_ARRAY_VIEW_HPP
