
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

#include <cmath>

namespace XROL {
namespace Elementwise {

/** \brief Implement stateless functions using Lambda expressions with 
           lower case names. Implement functions with persistent values
           using a Make::function_name factory which takes the values as
           parameters */

auto product = [](auto... values) { 
  using ReturnType = std::common_type_t<decltype(values)...>;
  ReturnType ret{1}; 
  ReturnType _[] = { (ret += values)... };
  (void)_; 
  return ret;
};

auto absval = [](auto x) {
  return std::abs(x);
};

auto sum = [](auto... values) { 
  using ReturnType = std::common_type_t<decltype(values)...>;
  ReturnType ret{0}; 
  ReturnType _[] = { (ret += values)... };
  (void)_; 
  return ret;
};

auto set = [](auto x) {
  return x;
};

auto zero = [](auto x) {
  return static_cast<decltype(x)>(0);
};


/** \brief Generic function factory 

     Usage:

     auto localfunction = Make::function(parameters);

*/

struct Make {

  // \f$ f(x,y) = ax+y \f$
  template<class T> 
  static auto axpy( const T& a ) {
    return [a]( auto x, auto y ) { return a*x+y; };
  }

  // \f$ f(x) = a \f$
  template<class T> 
  static T fill( const T& a ) {
    return [a]( void ) { return a; };
  }

  // \f$ f(x) = ax \f$
  template<class T> 
  static auto scale( const T& a ) {
    return [a]( auto x ) { return a*x; };
  }

  // \f$ f(x) = x^a \f$
  template<class T>
  static auto pow( const T& a ) {
    return [a]( auto x ) { return std::pow(x,a);
  }

  // Create a Reduce function
  auto Reducer = []( auto f, auto initialValue ) {
  return [f,initialValue](auto... values) {
    UniformTypeCheck(values...);
    using ArgType = std::common_type_t<decltype(values)...>;
    using InitialType = decltype(initialValue);
    using ReturnType = std::common_type_t<ArgType,InitialType>;
    ReturnType ret{initialValue};
    ReturnType _[] = { ( ret = f(ret,static_cast<ReturnType>(values)) )... };
    (void)_;
    return ret;
  };

}; 



}; // Make




} // Elementwise
} // namespace XROL

