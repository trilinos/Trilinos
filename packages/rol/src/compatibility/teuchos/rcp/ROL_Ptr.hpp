// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once

#include <memory>
#include <type_traits>

/* \file ROL_Ptr.hpp
 * \brief Provides unified interface Teuchos::RCP for legacy support.
 *        ROL will be build with this implementation by default  
 *        unless ROL_ENABLE_STD_SHARED_PTR:BOOL=ON
 */

#include <cstddef>
#include <utility>

#include "Teuchos_RCP.hpp"

namespace ROL {

template<class T> using Ptr = Teuchos::RCP<T>;

static const Teuchos::ENull nullPtr = Teuchos::null;

}

namespace ROL {

// Special case handling until C++14 (auto type deduction) is required
// because template type deduction does not work for initializer_list here. 

template<class T, class U>
inline
Ptr<T> makePtr( std::initializer_list<U>&& u ) {
  return Teuchos::rcp( new T(std::forward<std::initializer_list<U>>(u)) );
}


template<class T, class... Args>
inline 
Ptr<T> makePtr( Args&&... args ) {
  return Teuchos::rcp( new T(std::forward<Args>(args)...) );
}

template<class T>
inline
Ptr<T> makePtrFromRef( T& obj ) {
  return Teuchos::rcpFromRef(obj);
}

template<class T, class U> 
inline
Ptr<T> staticPtrCast( const Ptr<U>& r ) noexcept {
  return Teuchos::rcp_static_cast<T>(r);
}

template< class T, class U > 
inline
Ptr<T> constPtrCast( const Ptr<U>& r ) noexcept {
  return Teuchos::rcp_const_cast<T>(r);
}

template< class T, class U > 
inline
Ptr<T> dynamicPtrCast( const Ptr<U>& r ) noexcept {
  return Teuchos::rcp_dynamic_cast<T>(r);
}

template<class T>
inline
const T* getRawPtr( const Ptr<const T>& x ) {
  return x.get();
}

template<class T>
inline
T* getRawPtr( const Ptr<T>& x ) {
  return x.get();
}

template<class T>
inline
int getCount( const Ptr<T>& x ) {
  return x.strong_count();
}

template<class T>
inline
bool is_nullPtr( const Ptr<T>& x ) {
  return x.is_null();
}

template<typename T>
inline 
Ptr<T> toPtr( const Ptr<T>& ptr ) { 
  return ptr;
}

template<typename T>
inline 
Ptr<const T> toPtr( const Ptr<const T>& ptr ) { 
  return ptr;
}

template<class T>
struct IsSharedPtr : public std::false_type {};

template<class T>
struct IsSharedPtr<std::shared_ptr<T>> : public std::true_type {};

} // namespace ROL

