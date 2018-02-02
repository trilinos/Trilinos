 
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
template<class T> using StdVector = std::vector<T>;

template<class T>
struct VectorIndex<StdVector<T>> {
  using type = typename StdVector<T>::size_type;
};

template<class T>
struct VectorElement<StdVector<T>> {
  using type = T;
};

template<class T>
struct VectorMagnitude<StdVector<T>> {
  using type = typename Magnitude<T>::type;
};

template<class T>
struct VectorDual<StdVector<T>> {
  using type = StdVector<T>;
};

template<class T>
struct implements_elementwise<StdVector<T>> : std::true_type {};

template<class T>
struct implements_core<StdVector<T>> : std::true_type {};

// Functions



template<class T>
std::unique_ptr<StdVector<T>>
clone( const StdVector<T>& v ) {
  return std::move(std::make_unique<StdVector<T>>( v.size() )); 
}

template<class T>
index_t<StdVector<T>>
dimension( const StdVector<T>& x ) {
  return x.size();
}

template<class T>  
void set( StdVector<T>& x, const StdVector<T>& y ) { 
  std::transform( y.cbegin(), y.cend(), x.begin(),
  [](auto ye){ return ye; } );
}

template<class T>
void dual( dual_t<StdVector<T>>& xdual, 
            const StdVector<T>& xprim ) {
  set(xdual,xprim);
}

template<class T>
void plus( StdVector<T>& x, const StdVector<T>& y ) { 
  std::transform( x.begin(), x.end(), y.cbegin(), x.begin(), 
  [](auto xe, auto ye){ return xe+ye; } );
}

template<class T> 
void scale( StdVector<T>& x, 
            const element_t<StdVector<T>> alpha ) {
  for( auto &e : x ) e *= alpha;
}

template<class T> 
void fill( StdVector<T>& x,
           const element_t<StdVector<T>> alpha ) {
  for( auto &e : x ) e = alpha;
}

template<class T>
void axpy( StdVector<T> &x, 
           const element_t<StdVector<T>> alpha, 
           const StdVector<T> &y ) {
  std::transform( x.begin(), x.end(), y.cbegin(), x.begin(), 
  [alpha]( auto xe, auto ye ) { return alpha*ye+xe; } );
}

template<class T>
void basis( StdVector<T>& b, 
            index_t<StdVector<T>> i ) { 
  fill(b,0);
  b[i] = T(1.0);
}

template<class T>
element_t<StdVector<T>>
dot( const StdVector<T>& x, const StdVector<T>& y ) {
  using V = StdVector<T>;
  element_t<V> result = 0;
  for( index_t<V> i=0; i<x.size(); ++i ) result += x[i]*y[i];
  return result;
}

template<class T>
magnitude_t<StdVector<T>>
norm( const StdVector<T>& x ) {
  using V = StdVector<T>;
  magnitude_t<V> sum2 = 0;
  for( auto e: x ) sum2 += e*e;
  return std::sqrt(sum2);
} 


template<class T>
void print( StdVector<T>& x, std::ostream& os ) {
  for( auto e : x ) 
    os << e << " ";
  os << std::endl;
}


template<class R, class T>
auto reduce( const R& r, const StdVector<T>& x ) {
  auto result = r(); 
  for( auto e : x ) result = r(e,result);
  return result;
}


template<class Generator, class Distribution, class T>
void randomize( Generator& g, Distribution& d, StdVector<T>& x ) {
  for( auto &e : x ) e = d(g);
}


// Forward declaration of evaluate
namespace Elementwise {
template<class F, class Tuple>
decltype(auto) evaluate( const F& f, Tuple && t );
}



template<class T, class Function, class ...Vs>
void eval_function( StdVector<T>& x, const Function &f, const Vs&... vs ) {
  check_dimensions(x,vs...);
  for(auto i=0; i<x.size(); ++i)
    x[i] = Elementwise::evaluate(f,std::make_tuple(vs[i]...));
} 


template<class R, class F, class T, class ...Vs>
auto eval_function_and_reduce( const R& r, const F& f, 
  const StdVector<T>& x, const Vs&... vs ) {

  check_dimensions(x,vs...);
  auto result = r(); 
  auto dim = dimension(x);
  for(auto i=0; i<dim; ++i) {
    result = r(result, Elementwise::evaluate(f,std::make_tuple(x,vs[i]...)));
  }
  return result; 
}


} // namespace XROL
