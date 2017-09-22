
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

namespace FunctionDetails {
/**  \brief Converts tuple into variadic arguments for generic function calling
   */
  using std::size_t;
  template <class F,class Tuple, bool Done, size_t Total, size_t... N>
  struct call_impl {
    static auto call(F f, Tuple && t) {
        constexpr auto narg = sizeof...(N);
        return call_impl<F, Tuple, Total == 1+narg, Total, N..., narg>::call(f, std::forward<Tuple>(t));
    }
  };
  // Specialize terminating case
  template <class F, class Tuple, size_t Total, size_t... N>
  struct call_impl<F, Tuple, true, Total, N...> {
    static auto call(const F& f, Tuple && t) {
        return f(std::get<N>(std::forward<Tuple>(t))...);
    }
  };
}

namespace Elementwise {

/* This will be replaced by std::apply in C++17 */
template <class F, class Tuple>
auto evaluate = []( const F& f, Tuple && t) {
  // TODO: Type check the tuple
  constexpr auto tsize = std::tuple_size<std::decay_t<Tuple>>::value;    
  return detail::call_impl<F, Tuple, 0 == tsize, tsize>::call(f, std::forward<Tuple>(t));
};


// Extract template type of template template class
template<class> struct ValueTraits;
template<class T, template<class> class C> struct ValueTraits<C<T>> { using type = T; };


/* Implement reduction operations using the Curiously Recurring Template Pattern
*/ 

template<class ReduceType>
struct Reduce {
  using T = typename ValueTraits<ReduceType>::type;

  // Initial value
  T operator()(void) const {
    return static_cast<const ReduceType&>(*this).initialValue_;
  }
  
  T operator( const T& x, const T& y ) const {
    return static_cast<const ReduceType&>(*this)(x,y);
  }
};

template<class T> 
struct Minimum : Reducer<Minimum<T>> {
  static constexpr T initialValue_ = std::numeric_limits<T>::max();
  T operator() ( const T& x, const T& y ) const { return x<y?x:y; }
};
template<class T> constexpr T Minimum<T>::initialValue_;

template<class T> 
struct Maximum : Reducer<Minimum<T>> {
  static constexpr T initialValue_ = std::numeric_limits<T>::max();
  T operator() ( const T& x, const T& y ) const { return x<y?x:y; }
};
template<class T> constexpr T Minimum<T>::initialValue_;

template<class T>
struct Sum : Reducer<Sum<T>> {
  static constexpr T initialValue_ = T(0);
  T operator()( const T& x, const T& y )  const { return x+y; }
};
template class<T> constexpr T Sum<T>::initialValue_;

template<class T>
struct Product : Reducer<Product<T>> {
  static constexpr T initialValue_ = T(1);
  T operator()( const T& x, const T& y ) const { return x*y; }
};
template<class T> constexpr T Product<T>::initialValue_;

} // namespace Elementwise
} // namespace XROL

