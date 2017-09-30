 
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

#include "cxxstd.hpp"
#include "XROL_VectorTraits.hpp"

/** @ingroup la_group
    \file XROL::StdVector
    \brief Overload functions for std::vector
*/

namespace XROL {

// Traits specialization

template<class T>
struct VectorIndex<std::vector<T>> {
  using type = typename std::vector<T>::size_type;
};

template<class T>
struct VectorElement<std::vector<T>> {
  using type = T;
};

template<class T>
struct VectorMagnitude<std::vector<T>> {
  using type = typename Magnitude<T>::type;
};

template<class T>
struct VectorDual<std::vector<T>> {
  using type = std::vector<T>;
};

template<class T>
std::unique_ptr<std::vector<T>>
clone( const std::vector<T>& v ) {
  return std::move(std::make_unique<std::vector<T>>( v.size() )); 
}

template<class T>
index_t<std::vector<T>>
dimension( const std::vector<T>& x ) {
  return x.size();
}

template<class T>  
void set( std::vector<T>& x, const std::vector<T>& y ) { 
  std::transform( y.cbegin(), y.cend(), x.begin(),
  [](auto ye){ return ye; } );
}

template<class T>
void dual( dual_t<std::vector<T>>& xdual, 
            const std::vector<T>& xprim ) {
  set(xdual,xprim);
}

template<class T>
void plus( std::vector<T>& x, const std::vector<T>& y ) { 
  std::transform( x.begin(), x.end(), y.cbegin(), x.begin(), 
  [](auto xe, auto ye){ return xe+ye; } );
}

template<class T> 
void scale( std::vector<T>& x, 
            const element_t<std::vector<T>> alpha ) {
  for( auto &e : x ) e *= alpha;
}

template<class T> 
void fill( std::vector<T>& x,
           const element_t<std::vector<T>> alpha ) {
  for( auto &e : x ) e = alpha;
}

template<class T>
void axpy( std::vector<T> &x, 
           const element_t<std::vector<T>> alpha, 
           const std::vector<T> &y ) {
  std::transform( x.begin(), x.end(), y.cbegin(), x.begin(), 
  [alpha]( auto xe, auto ye ) { return alpha*ye+xe; } );
}

template<class T>
std::unique_ptr<std::vector<T>>
basis( const std::vector<T>& v, 
       index_t<std::vector<T>> i ) { 
  auto b = std::make_unique<std::vector<T>>( v.size() ); 
  fill(*b,0);
  (*b)[i] = T(1.0);
  return std::move(b);
}

template<class T>
element_t<std::vector<T>>
dot( const std::vector<T>& x, const std::vector<T>& y ) {
//  return std::inner_product(x.cbegin(), x.cend(), y.begin(), 0); 
  using V = std::vector<T>;
  element_t<V> result = 0;
  for( index_t<V> i=0; i<x.size(); ++i ) result += x[i]*y[i];
  return result;
}

template<class T>
magnitude_t<std::vector<T>>
norm( const std::vector<T>& x ) {
//  return std::inner_product(x.cbegin(), x.cend(), x.begin(), 0); 
  using V = std::vector<T>;
  magnitude_t<V> sum2 = 0;
  for( auto e: x ) sum2 += e*e;
  return std::sqrt(sum2);
}


template<class T>
void print( std::vector<T> &x, std::ostream &os ) {
  for( auto e : x ) 
    os << e << " ";
  os << std::endl;
}


template<class R, class T>
auto reduce( const R& r, const std::vector<T>& x ) {
  auto result = r(); 
  for( auto e : x ) result = r(e,result);
  return result;
}


template<class Generator, class Distribution, class T>
void randomize( Generator& g, Distribution& d, std::vector<T>& x ) {
  for( auto &e : x ) e = d(g);
}


// Forward declaration of evaluate
namespace Elementwise {
template<class F, class Tuple>
decltype(auto) evaluate( const F& f, Tuple && t );
}



template<class T, class Function, class ...Vs>
void eval_function( std::vector<T>& x, const Function &f, const Vs&... vs ) {
  check_dimensions(x,vs...);
  for(auto i=0; i<x.size(); ++i)
    x[i] = Elementwise::evaluate(f,std::make_tuple(vs[i]...));
} 


template<class R, class F, class T, class ...Vs>
auto eval_function_and_reduce( const R& r, const F& f, 
  const std::vector<T>& x, const Vs&... vs ) {

  check_dimensions(x,vs...);
  auto result = r(); 
  auto dim = dimension(x);
  for(auto i=0; i<dim; ++i) {
    result = r(result, Elementwise::evaluate(f,std::make_tuple(x,vs[i]...)));
  }
  return result; 
}


} // namespace XROL
