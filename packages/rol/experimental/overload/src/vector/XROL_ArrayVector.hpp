 
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
    \file XROL::ArrayVector
    \brief Overload functions for std::array
*/

namespace XROL {

template<class T, std::size_t N> using ArrayVector = std::array<T,N>;

// Traits specialization

template<class T,std::size_t N>
struct VectorIndex<ArrayVector<T,N>> {
  using type = std::size_t;
};

template<class T, std::size_t N>
struct VectorElement<ArrayVector<T,N>> {
  using type = T;
};

template<class T, std::size_t N>
struct VectorMagnitude<ArrayVector<T,N>> {
  using type = typename Magnitude<T>::type;
};

template<class T, std::size_t N>
struct VectorDual<ArrayVector<T,N>> {
  using type = ArrayVector<T,N>;
};

template<class T, std::size_t N>
struct implements_elementwise<ArrayVector<T,N>> : std::true_type {};

template<class T, std::size_t N>
struct implements_core<ArrayVector<T,N>> : std::false_type {};

// Functions



template<class T, std::size_t N>
std::unique_ptr<ArrayVector<T,N>>
clone( const ArrayVector<T,N>& v ) {
  return std::move(std::make_unique<ArrayVector<T,N>>); 
}

template<class T, std::size_t N>
index_t<ArrayVector<T,N>>
dimension( const ArrayVector<T,N>& x ) {
  return N;
}

template<class T, std::size_t N>  
void set( ArrayVector<T,N>& x, const ArrayVector<T,N>& y ) { 
  std::transform( y.cbegin(), y.cend(), x.begin(),
  [](auto ye){ return ye; } );
}

template<class T, std::size_t N>
void dual( dual_t<ArrayVector<T,N>>& xdual, 
            const ArrayVector<T,N>& xprim ) {
  set(xdual,xprim);
}

template<class T, std::size_t N>
void plus( ArrayVector<T,N>& x, const ArrayVector<T,N>& y ) { 
  std::transform( x.begin(), x.end(), y.cbegin(), x.begin(), 
  [](auto xe, auto ye){ return xe+ye; } );
}

template<class T, std::size_t N> 
void scale( ArrayVector<T,N>& x, 
            const element_t<ArrayVector<T,N>> alpha ) {
  for( auto &e : x ) e *= alpha;
}

template<class T, std::size_t N> 
void fill( ArrayVector<T,N>& x,
           const element_t<ArrayVector<T,N>> alpha ) {
  for( auto &e : x ) e = alpha;
}

template<class T, std::size_t N>
void axpy( ArrayVector<T,N> &x, 
           const element_t<ArrayVector<T,N>> alpha, 
           const ArrayVector<T,N> &y ) {
  std::transform( x.begin(), x.end(), y.cbegin(), x.begin(), 
  [alpha]( auto xe, auto ye ) { return alpha*ye+xe; } );
}

template<class T, std::size_t N>
void basis( ArrayVector<T,N>& b, 
            index_t<ArrayVector<T,N>> i ) { 
  fill(b,0);
  b[i] = T(1.0);
}

template<class T, std::size_t N>
element_t<ArrayVector<T,N>>
dot( const ArrayVector<T,N>& x, const ArrayVector<T,N>& y ) {
  using V = ArrayVector<T,N>;
  element_t<V> result = 0;
  for( index_t<V> i=0; i<N; ++i ) result += x[i]*y[i];
  return result;
}

template<class T, std::size_t N>
magnitude_t<ArrayVector<T,N>>
norm( const ArrayVector<T,N>& x ) {
  using V = ArrayVector<T,N>;
  magnitude_t<V> sum2 = 0;
  for( auto e: x ) sum2 += e*e;
  return std::sqrt(sum2);
} 


template<class T, std::size_t N>
void print( ArrayVector<T,N>& x, std::ostream& os ) {
  for( auto e : x ) 
    os << e << " ";
  os << std::endl;
}


template<class R, class T, std::size_t N>
auto reduce( const R& r, const ArrayVector<T,N>& x ) {
  auto result = r(); 
  for( auto e : x ) result = r(e,result);
  return result;
}


template<class Generator, class Distribution, class T, std::size_t N>
void randomize( Generator& g, Distribution& d, ArrayVector<T,N>& x ) {
  for( auto &e : x ) e = d(g);
}


// Forward declaration of evaluate
namespace Elementwise {
template<class F, class Tuple>
decltype(auto) evaluate( const F& f, Tuple && t );
}



template<class T, std::size_t N, class Function, class ...Vs>
void eval_function( ArrayVector<T,N>& x, const Function &f, const Vs&... vs ) {
  check_dimensions(x,vs...);
  for(auto i=0; i<N; ++i)
    x[i] = Elementwise::evaluate(f,std::make_tuple(vs[i]...));
} 


template<class R, class F, class T, std::size_t N, class ...Vs>
auto eval_function_and_reduce( const R& r, const F& f, 
  const ArrayVector<T,N>& x, const Vs&... vs ) {

  check_dimensions(x,vs...);
  auto result = r(); 
  auto dim = dimension(x);
  for(auto i=0; i<dim; ++i) {
    result = r(result, Elementwise::evaluate(f,std::make_tuple(x,vs[i]...)));
  }
  return result; 
}


} // namespace XROL
