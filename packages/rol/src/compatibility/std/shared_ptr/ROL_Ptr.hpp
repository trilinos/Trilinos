// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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

