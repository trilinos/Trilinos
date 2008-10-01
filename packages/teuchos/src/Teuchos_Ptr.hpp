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


#ifndef TEUCHOS_PTR_HPP
#define TEUCHOS_PTR_HPP


#include "Teuchos_PtrDecl.hpp"
#include "Teuchos_RCP.hpp"


namespace Teuchos {


namespace PtrPrivateUtilityPack {
void throw_null( const std::string &type_name );
} // namespace PtrPrivateUtilityPack


template<class T> inline
Ptr<T>::Ptr( ENull null_in )
  : ptr_(0)
{}


template<class T> inline
Ptr<T>::Ptr( T *ptr )
  : ptr_(ptr)
{}


template<class T> inline
Ptr<T>::Ptr(const Ptr<T>& ptr)
  :ptr_(ptr.ptr_)
{}


template<class T>
template<class T2> inline
Ptr<T>::Ptr(const Ptr<T2>& ptr)
  :ptr_(ptr.get())
{}


template<class T> inline
Ptr<T>& Ptr<T>::operator=(const Ptr<T>& ptr)
{
  ptr_ = ptr.get();
  return *this;
}


template<class T> inline
T* Ptr<T>::operator->() const
{
  debug_assert_not_null();
  debug_assert_valid_ptr();
  return ptr_;
}


template<class T> inline
T& Ptr<T>::operator*() const
{
  debug_assert_not_null();
  debug_assert_valid_ptr();
  return *ptr_;
}


template<class T> inline
T* Ptr<T>::get() const
{
  debug_assert_valid_ptr();
  return ptr_;
}


template<class T> inline
T* Ptr<T>::getRawPtr() const
{
  return get();
}


template<class T> inline
const Ptr<T>& Ptr<T>::assert_not_null() const
{
  if(!ptr_)
    PtrPrivateUtilityPack::throw_null(TypeNameTraits<T>::name());
  return *this;
}


template<class T> inline
void Ptr<T>::debug_assert_valid_ptr() const
{
#ifdef TEUCHOS_DEBUG
  rcp_.access_private_node().assert_valid_ptr(*this);
#endif
}


#ifdef TEUCHOS_DEBUG


template<class T> inline
Ptr<T>::Ptr( const RCP<T> &p )
  : ptr_(p.getRawPtr()), rcp_(p)
{}


#endif // TEUCHOS_DEBUG


} // namespace Teuchos


template<class T>
std::ostream& Teuchos::operator<<( std::ostream& out, const Ptr<T>& p )
{
  out
    << TypeNameTraits<RCP<T> >::name() << "{"
    << "ptr="<<(const void*)(p.get()) // I can't find any alternative to this C cast :-(
    <<"}";
  return out;
}


#endif // TEUCHOS_PTR_HPP
