m 
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

#include "XROL_VectorTraits.hpp"
#include <vector>
#include <memory>

/** @ingroup la_group
    \file XROL::StdVectorTraits
    \brief Specializes the XROL::VectorTraits interface for std::vector 
*/


namespace XROL {

template<class Element> 
struct ElementTraits<std::vector<Element>> {

  using IndexType     = typename std::vector<Element>::size_type;
  using ElementType   = Element;
  using MagnitudeType = MagnitudeTypeTraits<Element>::type;

};



template<class Element>
struct VectorCreationTraits<std::vector<Element>> {

  using V = std::vector<Element>;

  using IndexType     = typename ElementTraits<V>::IndexType;
  using ElementType   = typename ElementTraits<V>::ElementType;
  using MagnitudeType = typename ElementTraits<V>::MagnitudeType;

  static auto clone( const V& v ) {
    return std::make_shared<V>( V(v.size()) ); 
  }

  static auto basis( const V& v, IndexType i ) { 
    auto b = std::make_shared<V>( V(v.size() );
    for( auto &e : *b ) e = ElementType(0);
    (*b)[i] = ElementType(1);
    return b;
  }
}; // VectorCreationTraits

 
template<class Element>
struct VectorFunctionTraits<std::vector<Element>> {

  using V = std::vector<Element>;
  using VST = VectorSpaceTraits<V>;

  using IndexType     = typename ElementTraits<V>::IndexType;
  using ElementType   = typename ElementTraits<V>::ElementType;
  using MagnitudeType = typename ElementTraits<V>::MagnitudeType;

  template<class T> Reduce = Elementwise::Reduce<T>;

  template<class Function, class ...Vs>
  static void transform( V& x, const Function &f, const Vs&... vs ) {
    VST::CheckDimensions(vs...);
    // TODO: Add x to check
    for(IndexType i=0; i<x.size(); ++i)
      x[i] = Elementwise::call(f,std::make_tuple(vs[i]...));
  } 

  template<class R>
  static ElementType reduce( const Reduce<R>& r, const V& x ) {
    ElementType result = r(); // Initial value
    for( auto e : x ) result = r(e,result);
    return result;
  }

  template<class Function, class R, class ...Vs>
  static ElementType reduce( const Function &f, const Reduce<R> &r, const Vs&... vs ) {
    VST::CheckDimensions(vs...);
    ElementType result = r(); 
    IndexType dim = VST::dimension(std::get<0>(std::make_tuple(vs...)));
    for(IndexType i=0; i<dim; ++i) {
      result = r(result, Elementwise::call(f,std::make_tuple(vs[i]...)));
    }
    return result; 
  }

}; // VectorFunctionTraits

 

template<class Element>
struct VectorOptionalTraits {

   
  using V = std::vector<Element>;

  // TODO: Check dimensions and make check const-invariant
  using VST = VectorSpaceTraits<V>;

  using IndexType     = typename ElementTraits<V>::IndexType;
  using ElementType   = typename ElementTraits<V>::ElementType;
  using MagnitudeType = typename ElementTraits<V>::MagnitudeType;

  static void plus( V& x, const V& y ) { 
    std::transform( x.begin(), x.end(), y.cbegin(), x.begin(), 
      [](auto xe, auto ye){return xe+ye} );
  }
  
  static void set( V& x, const V& y ) { 
    std::transform( y.cbegin(), y.cend(), x.begin(),
      [](auto ye){return ye} );
  }
 
  static void scale( V& x, const ElementType alpha ) {
    for( auto &e : x ) e * =alpha;
  }

  static void zero( V &x ) {
    for( auto &e : x ) x = ElementType(0);
  }

  static void axpy( V &x, const ElementType alpha, const V &y) {
    std::transform( x.begin(), x.end(), y.cbegin(), x.begin(), 
      [alpha]( auto xe, auto ye ) { return alpha*xe+ye; } );
  }

  static ElementType dot( const V& x, const V& y ) {
    return std::inner_product(x.cbegin(), x.cend(), y.begin(), 0); 
  }

};


template<class Element>
struct VectorPrintTraits<std::vector<Element>> {

  using V = std::vector<Element>;

  static void print( const V& x, std::ostream &os ) {
    os << "\n";
    for( auto e: x ) os << e << " ";
    os << std::endl;
  }
};



} // namespace XROL
