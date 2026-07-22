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
#include <cstddef>
#include <utility>

/* \file  ROL_Ptr.hpp
 * \brief Wraps the C++11 std::shared_ptr
 *        ROL will be build with this implementation if CMake is
 *        configured with ROL::Ptr='shared_ptr'
 */

namespace ROL {

template<class T> using Ptr = std::shared_ptr<T>;

static std::nullptr_t nullPtr = nullptr;

template<class T, class... Args>
inline
Ptr<T> makePtr( Args&&... args ) {
  return std::make_shared<T>(args...);
}

template<class T, class U>
inline
Ptr<T> makePtr( std::initializer_list<U>&& u )  {
  return std::make_shared<T>(std::forward<std::initializer_list<U>>(u));
}

template<class T>
inline
Ptr<T> makePtrFromRef( T& obj ) {
  return std::shared_ptr<T>(&obj,[](void*){});
}

template<class T>
inline
Ptr<const T> makePtrFromRef( const T& obj ) {
  return std::shared_ptr<const T>(&obj,[](const void*){});
}

template< class T, class U > 
inline
Ptr<T> staticPtrCast( const Ptr<U>& r ) noexcept {
  return std::static_pointer_cast<T>(r);
}

template< class T, class U > 
inline
Ptr<T> constPtrCast( const Ptr<U>& r ) noexcept {
  return std::const_pointer_cast<T>(r);
}

template< class T, class U > 
inline
Ptr<T> dynamicPtrCast( const Ptr<U>& r ) noexcept {
  return std::dynamic_pointer_cast<T>(r);
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
  return x.use_count();
}

template<class T>
inline
bool is_nullPtr( const Ptr<T>& x ) {
  return x == nullPtr;
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

