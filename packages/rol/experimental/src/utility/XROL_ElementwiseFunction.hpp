
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
#include <algorithm>
#include <limits>
#include <tuple>

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


template<class V> 
struct Elementwise { 

  using ElementType   = typename ElementTraits<V>::ElementType;
  using MagnitudeType = typename ElementTraits<V>::MagnitudeType;

  template<class ...Es>
  static ElementType product( Es... es ) {
    UniformTypeCheck(es...);
    ElementType ret{1.0};
    ElementType _[] = { (ret *= es)... };
    (void)_;
    return ret;
  }
 
  template<class ...Es>
  static ElementType sum( Es... es ) {
    UniformTypeCheck(es...);
    ElementType ret{0.0};
    ElementType _[] = { (ret *= es)... };
    (void)_;
    return ret;
  }
 
  template<class ...Es>
  static ElementType min( Es... es ) {
    UniformTypeCheck(es...);
    ElementType ret{std::numeric_limits<ElementType>::max()};
    ElementType _[] = { (ret = std::min(ret,es))... };
    (void)_;
    return ret;
  }
 
  template<class ...Es>
  static ElementType max( Es... es ) {
    UniformTypeCheck(es...);
    ElementType ret{std::numeric_limits<ElementType>::lowest()};
    ElementType _[] = { (ret = std::max(ret,es))... };
    (void)_;
    return ret;
  }

  template <class F, class Tuple>
  static ElementType call( const F& f, Tuple && t) {
    // TODO: Type check the tuple
    constexpr auto tsize = std::tuple_size<std::decay_t<Tuple>>::value;    
    return detail::call_impl<F, Tuple, 0 == tsize, tsize>::call(f, std::forward<Tuple>(t));
  }
}; 


} // namespace XROL

