
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
           using templated functor classes with upper case names. */


auto sum = [](auto... values) { 
  using ReturnType = std::common_type_t<decltype(values)...>;
  ReturnType ret{0}; 
  ReturnType _[] = { (ret += values)... };
  (void)_; // To silence possible warning
  return ret;
};

template<class T>
class Axpy {
private:
  T a_;
public:
  Axpy( T a ) : a_(a) {}
  auto operator() ( auto y, auto x ) const {
    return a_*x+y;
  }
};

template<class T>
class Scale {
private:
  T a_;
public:
  Scale( T a ) : a_(a) {}
  auto operator() ( auto y ) const {
    return a_*y;
  }
};

auto set = [](auto x) {
  return x;
};

auto zero = [](auto x) {
  return static_cast<decltype(x)>(0);
};

template<class T>
class Fill {
private:
  T a_;
public:
  Fill( T a ) : a_(a) {}
  auto operator() ( auto x ) const {
    return a_;
  }
};

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

template<class T>
class Power {
private:
  T p_;
public:
  Power( T p ) : p_(p) {}
  auto operator() ( auto x ) const {
    return std::pow(x,p_);
  }
};


} // Elementwise
} // namespace XROL

