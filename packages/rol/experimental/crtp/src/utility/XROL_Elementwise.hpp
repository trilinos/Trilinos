
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

#include "XROL_Core.hpp"

#include <utility>
#include <tuple>

namespace XROL {

namespace details {
using namespace std;

template<class Tuple>
using Indices = make_index_sequence<tuple_size<decay_t<Tuple>>::value>;

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

template<class F, class Tuple>
decltype(auto) evaluate( F&& f, Tuple && t ) {
  constexpr auto tsize = tuple_size<decay_t<Tuple>>::value;
  return call_impl<F,Tuple,0 == tsize,tsize>::call(forward<F>(f),forward<Tuple>(t));
}


template<class F, class Tuple, size_t... I>
decltype(auto) unary_apply_impl(F&& f, Tuple&& t, index_sequence<I...>) {
  return forward<F>(f)(get<I>(forward<Tuple>(t))...);
}

template<class F, class Tuple>
decltype(auto) unary_apply(F&& f, Tuple&& t) {
  return unary_apply_impl(forward<F>(f), forward<Tuple>(t), Indices<Tuple>{});
}

} // namespace details 

// Common elementwise functions
auto sum = [](auto... args) {
  using ReturnType = std::common_type_t<decltype(args)...>;
  ReturnType ret{0};
  ReturnType  _[] = { (ret += args )... };
  (void)_;
  return ret;
};

auto product = [](auto... args) {
  using ReturnType = std::common_type_t<decltype(args)...>;
  ReturnType ret{1};
  ReturnType  _[] = { (ret *= args )... };
  (void)_;
  return ret;
};

auto divide = [] ( auto num, auto den ) { return num/den; };




} // namespace XROL

