
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

namespace XROL {

namespace tuple_functions {

using namespace std;

template<class Tuple>
using Indices = make_index_sequence<tuple_size<decay_t<Tuple>>::value>;

/*-------------------------------------------------------------*/
/*                                                             */
/*   unary_apply(Function,Tuple)                               */ 
/*                                                             */
/*   Apply a function f to every element e_i in a tuple and    */
/*   return a new tuple containing the elements f_i = f(e_i)   */
/*                                                             */
/*   returns new tuple                                         */
/*   ( f_1, f_2, ... ) <- ( f(e_1), f(e_2), ... )              */
/*                                                             */
/*                                                             */
/*-------------------------------------------------------------*/

template<class Function, class Tuple, size_t... Is>
constexpr auto unary_apply_impl( const Function& f, 
                                 Tuple&& x, 
                                 index_sequence<Is...>) {

  return make_tuple( f(get<Is>(x) )... );
}

template<class Function, class Tuple>
constexpr auto unary_apply( Function&& f, 
                            Tuple&& x ) {

  return unary_apply_impl(forward<Function>(f),
                          forward<Tuple>(x),
                          Indices<Tuple>{});
}


/*---------------------------------------------------------------*/
/*                                                               */
/*   unary_eval(Function,Tuple)                                  */ 
/*                                                               */
/*   Apply a unary function f to every element e_i in a tuple    */
/*   and does one of two things depending on the nature of f:    */
/*                                                               */
/*   1) If f takes arbitrary arguments, modifies them in-place   */
/*   and does not return a value, this function does the same    */
/*                                                               */
/*   modifies tuple in place                                     */
/*   ( e_1, e_2, ... ) <- ( f(e_1), f(e_2), ... )                */
/*                                                               */
/*   2) If f reduces its arguments to a single value, so will    */
/*   this function.                                              */                                
/*                                                               */
/*   returns a common_type quantity                              */
/*   r = reduce(r,e_i) for all i                                 */
/*                                                               */
/*---------------------------------------------------------------*/


template<class Function, class Tuple, size_t... Is>
decltype(auto) unary_eval_impl( Function&& f, 
                                Tuple&& t, 
                                index_sequence<Is...> ) {

  return forward<Function>(f)(get<Is>(forward<Tuple>(t))...);
}

template<class Function, class Tuple>
decltype(auto) unary_eval( Function&& f, 
                           Tuple&& t ) {

  return unary_eval_impl(forward<Function>(f), 
                         forward<Tuple>(t), 
                         Indices<Tuple>{});
}

} // namespace tuple_functions


// Aliases for brevity

template<class Function, class Tuple>
constexpr decltype(auto) unary_apply = tuple_functions::unary_apply<Function,Tuple>;

template<class Function, class Tuple>
constexpr decltype(auto) unary_eval = tuple_functions::unary_eval<Function,Tuple>;




} // namespace XROL

