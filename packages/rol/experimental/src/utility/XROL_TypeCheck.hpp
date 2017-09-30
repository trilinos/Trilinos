
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

namespace XROL {

namespace details {

using namespace std;

/* Generate compile time error if there are any different template types */
template<class T, class... Ts>
struct CheckType {
  using u_type = typename CheckType<Ts...>::type;
  static_assert(is_same<T,u_type>::value, "Incompatible types!");
  using type = T;
};  

/* Terminating case */
template<class T, class U>
struct CheckType<U,T> {
  static_assert(is_same<T,U>::value, "Incompatible types!");
  using type = T;
};

/* Generate compile time error if there are any different template types
    Ignores differences in const, volatile, raw pointer, or reference modifier
 */
template<class T, class... Ts>
struct CheckDecayType {
  using u_type = decay_t<typename CheckDecayType<Ts...>::type>;
  static_assert(is_same<decay_t<T>,u_type>::value,"Incompatible decay types!");
};

/* Terminating case */
template<class T, class U>
struct CheckDecayType<U,T> {
  static_assert(is_same<decay_t<T>,decay_t<U>>::value, "Incompatible decay types!");
  using u_type = decay_t<U>;
};

} // namespace details


template<class... Ts>
using check_type_t = typename details::CheckType<Ts...>::type;

template<class... Ts>
using check_decay_type_t = typename details::CheckDecayType<Ts...>::type;

template<class... Ts>
constexpr void check_types( Ts... ts ) {
  check_type_t<Ts...> _;
  (void)_;
}

template<class... Ts>
constexpr void check_decay_types( Ts... ts ) {
  check_decay_type_t<Ts...> _;
  (void)_;
}

// Forward declaration
struct IncompatibleDimensions;

template<class ...Vs>
void check_dimensions( Vs... vecs ) {

  check_decay_type(vecs...);
  auto first = std::get<0>(std::make_tuple(vecs...)).dimension();
  auto cond = [first](auto i){ return i.dimension() == first; };
  bool match{true};
  bool _1[] = { ( match &= cond(vecs) )... }; (void)_1;

  if(!match) {
    std::string dims;
    bool _2[] = { ( dims += std::to_string( vecs.dimension() ) + " ")...}; (void)_2;
    throw IncompatibleDimensions(dims);
  }
}




} // namespace XROL

