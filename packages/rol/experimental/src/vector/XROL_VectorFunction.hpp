
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

#include "XROL.hpp"

namespace XROL {

template<class V>
// std::unique_ptr<V> 
auto clone( const V& x ) { 
  static_assert("clone() must be specialized for type " +
   typeid(x).name() );
  return nullptr; 
}

template<class V> 
// std::unique_ptr<V> 
auto basis( const V& x, ElementTraits<V>::IndexType i ) {
  static_assert("basis() must be specialized for type " +
   typeid(x).name() );
  return nullptr
}

template<class V>
// ElementTraits<V>::IndexType 
auto dimension( const V& x ) {
  static_assert("dimension() be specialized for type " +
   typeid(x).name() );
  return 0;
}

template<class V>
void plus( V& x, const V& y ) {
  VectorFunctionTraits<V>::transform(x,[](auto xe, auto ye){xe+=ye;},y);
}

template<class V>
void scale( V& x, const ElementTraits<V>::ElementType alpha ) {
  VectorFunctionTraits::transform(x,[alpha](auto xe){xs*=alpha;},y);
}

template<class V>
// ElementTraits<V>::ElementType
auto dot( const V& x, const V &y ) {
  using ElementType = ElementTraits<V>::ElementType;
  Sum<ElementType> sum;
  return VectorFunctionTraits<V>::transform_and_reduce( 
    [](auto xe,auto ye){return xe*ye;}, sum, x, y); 
}

template<class V> 
// ElementTraits<V>::MagnitudeType
auto norm( const V& x ) {
  using ElementType   = ElementTraits<V>::ElementType;
  using MagnitudeType = ElementTraits<V>::MagnitudeType;
  Sum<ElementType> sum;
  auto norm2 = VFT::transform_and_reduce( [](auto xe){return xe*xe;}, sum, x);  
  return std::sqrt(norm2);
}

template<class V>
void axpy( V& x, const ElementTraits<V>::ElementType



} // namespace XROL

