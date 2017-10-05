
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

namespace XROL {

// Forward Declaration
template<class...Vs> void check_dimension(Vs...);


namespace Elementwise {

namespace details {
 
using namespace std;

template <class F, class Tuple, bool Done, size_t Total, size_t... N>
struct call_impl {
  static decltype(auto) call(const F& f, Tuple && t) {
    constexpr auto narg = sizeof...(N);
    return call_impl<F, Tuple, Total == 1+narg, Total, N..., narg>::call(f, forward<Tuple>(t));
  }
};

// Specialize terminating case
template <class F, class Tuple, size_t Total, size_t... N>
struct call_impl<F, Tuple, true, Total, N...> {
  static decltype(auto) call(const F& f, Tuple && t) {
    return f(get<N>(forward<Tuple>(t))...);
  }
};

} // namespace details 

/** \fn         evaluate 
    \brief      Evaluate a variadic function using an input tuple and return the result
    @param[in]  f  A function taking a nonzero number of arguments
    @param[in]  t  A tuple of values 
    \return     The value \f$f(t_1,t_2,\hdots)\f$
    Note: This can be replaced by std::apply in C++17 */
template<class F, class Tuple>
decltype(auto) evaluate( const F& f, Tuple && t ) {
  using namespace std;
  constexpr auto tsize = tuple_size<decay_t<Tuple>>::value;
  return details::call_impl<F,Tuple,0 == tsize,tsize>::call(f,forward<Tuple>(t));
}

/** \fn         evaluate_unary_composition 
    \brief      Apply a unary function to every element a tuple and 
                a variadic function to each of the resulting values
    @param[in]  f  A function taking a nonzero number of arguments
    @param[in]  u  A function taking a single argument
    @param[in]  t  A tuple of values 
    \return     The value \f$ f(u(t_1),u(t_2),...) \f$
 
*/
template<class F, class U, class Tuple>
decltype(auto) evaluate_unary_composition( const F&f, const U& u, Tuple && t ) {
  using namespace std;
  constexpr auto tsize = tuple_size<decay_t<Tuple>>::value;
  return details::call_impl<F,Tuple,0 == tsize,tsize>::call(f,forward<Tuple>(u(t)));
}

/** \fn         evaluate_all
    \brief      Apply a unary modify-in-place function to every element of all of 
                the supplied vectors
    @param[in]      f   A function taking a non-const scalar argument 
    @param[in,out]  vs  A set of vectors to be modified by f
*/
template<class F, class... Vs>
void evaluate_all( const F& f, Vs&... vs ) {
  using expand_type = int[];
  expand_type{ (f(vs),0)...};
} 

/** CRTP implementation of common reduce types 

    Alternately, we could redefine a reduce operation purely in terms
    of elementwise lambdas with a captured initial value

*/
template<class T, template<class> class ReduceType>
struct Reduce {

  // Initial value
  T operator()(void) const {
    return static_cast<const ReduceType<T>&>(*this).initialValue_;
  }

  T operator()( const T& x, const T& y ) const {
    return static_cast<const ReduceType<T>&>(*this)(x,y);
  }
}; // Reduce

template<class T> 
struct Minimum : Reduce<T,Minimum> {
  static constexpr T initialValue_ = std::numeric_limits<T>::max();
  T operator() ( const T& x, const T& y ) const { return x<y?x:y; }
};
template<class T> constexpr T Minimum<T>::initialValue_;

template<class T> 
struct Maximum : Reduce<T,Maximum> {
  static constexpr T initialValue_ = std::numeric_limits<T>::lowest();
  T operator() ( const T& x, const T& y ) const { return x>y?x:y; }
};
template<class T> constexpr T Maximum<T>::initialValue_;

template<class T>
struct Sum : Reduce<T,Sum> {
  static constexpr T initialValue_ = T(0);
  T operator()( const T& x, const T& y )  const { return x+y; }
};
template<class T> constexpr T Sum<T>::initialValue_;

template<class T>
struct Product : Reduce<T,Product> {
  static constexpr T initialValue_ = T(1);
  T operator()( const T& x, const T& y ) const { return x*y; }
};
template<class T> constexpr T Product<T>::initialValue_;


} // namespace Elementwise


template<class V>
auto make_min( const V& x ) {
  return Elementwise::Minimum<element_t<V>>(); 
}

template<class V>
auto make_max( const V& x ) {
  return Elementwise::Maximum<element_t<V>>(); 
}

template<class V>
auto make_sum( const V& x ) {
  return Elementwise::Sum<element_t<V>>();
}

template<class V>
auto make_product( const V& x ) {
  return Elementwise::Sum<element_t<V>>();
}


} // namespace XROL

