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

#ifndef PDEOPT_TEMPLATE_TOOLS_HPP
#define PDEOPT_TEMPLATE_TOOLS_HPP

#include "Sacado.hpp"
#include <type_traits>
#include <utility>


/*! \file  template_tools.hpp
    \brief Provides generic programming components with emphasis on 
           automatic differentiation
*/

namespace TemplateTools {

template<class...> struct voider { using type=void; };
template<class...Ts> using void_t = typename voider<Ts...>::type;

template<template<class...>class Z, class, class...Ts>
struct has_method: public std::false_type {};

template<template<class...>class Z, class...Ts>
struct has_method<Z, void_t<Z<Ts...>>, Ts...> : public std::true_type {};


// Detect existence of .size()
template <class Container>
using c_size_t = decltype(std::declval<Container>().size());

template <class Container> 
using has_size = has_method<c_size_t,Container>;

// Detect subscript operator []
template <class Container, class Index>
using c_subscript_t = decltype(std::declval<Container>()[std::declval<Index>()]);

template <class Container, class Index>
using has_subscript = has_method<c_subscript_t, Container, Index>;

// Call this method if the container has a .size()
//template<class Container> 
//typename std::enable_if<has_size<Container>::value,decltype(std::declval<Container>().size())>::type 
//size( const Container &c ) {
//  return c.size();  
//}

// Container type deduction
template<class Container>
struct ContainerTypes {
    typename std::add_rvalue_reference<Container>::type arg = std::declval<Container>();
    
  static constexpr bool hasSize      = has_size<Container>::value;
   
  using IndexType = typename std::conditional<hasSize,
                      decltype(std::declval<Container>().size()),void>::type;

  static constexpr bool hasSubscript = has_subscript<Container,IndexType>::value;

  using ElementType = typename std::conditional<hasSubscript,
                        decltype(std::declval<Container>()[std::declval<IndexType>()]),void>::type; 

}; // ContainerInterface

template<class Container>
typename std::enable_if<has_size<Container>::value,decltype(std::declval<Container>().size())>::type 
size( const Container &c ) {
  return c.size();  
}

template<class Container>
typename std::enable_if<!has_size<Container>::value,decltype(std::declval<Container>().size())>::type 
size( const Container &c ) {
  return 0;  
}


} // namespace TemplateTools

namespace AD {

/*  ResultType<T...> deduces the type that would be needed to add 
 *  the constituent types together. For example if DFad = Sacado::Fad::DFad<Real>
 *  then ReturnType<Real,DFad> would be DFad since adding the two together
 *  would yield a DFad.
 */
template<typename T,typename ...Args>
using ResultType = typename std::common_type<T,Args...>::type;


/* Determine at compile time, whether the ScalarT has both of the required
 * properties:
 *
 *  1) ScalarT has a dx() member 
 *  2) dx() takes an integer argument
 *
 */

// Detect dx() methjod
//template <class Container, class Index>
//using dx_t = decltype(std::declval<Container>().dx(std::declval<Index>()));

//template <class Container, class Index>
//using has_dx = has_method<dx_t, Container, Index>;




template <class ScalarT> struct has_dx
{
  std::add_rvalue_reference<int>::type arg = std::declval<int>();

  // checks whether there is a method .dx(int)
  template <typename C> static constexpr decltype(std::declval<C>().dx(arg), bool()) test(int)
  { return true;  }
  template <typename C> static constexpr bool test(...)
  { return false; }
  static constexpr bool value = test<ScalarT>(int());
};

/* Determine whether ScalarT has vale() method */
template <class ScalarT> struct has_val
{
  template <typename C> static constexpr decltype(std::declval<C>().val(), bool()) test(int)
  { return true;  }
  template <typename C> static constexpr bool test(...)
  { return false; }
  static constexpr bool value = test<ScalarT>(int());
};




/* Automatic differentiation function for arbitrary types
 *
 * If the given object has a .dx(int) function, call it, otherwise
 * return 0
 */
template<class ScalarT> 
typename std::enable_if<has_dx<ScalarT>::value,decltype(std::declval<ScalarT>().dx(0))>::type 
derivative( const ScalarT &f, int i ) {
  return f.dx(i);  // Does not check validity of i value!
}


// If the method does not exist, return integer 0
template<class ScalarT> 
typename std::enable_if<!has_dx<ScalarT>::value,int>::type 
derivative( const ScalarT &f, int i ) {
  return int(0);  
}

template<class ScalarT> 
typename std::enable_if<has_val<ScalarT>::value,decltype(std::declval<ScalarT>().val())>::type 
value( ScalarT &f ) {
  return f.val();  
}

// Return a copy if the argument is const 
template<class ScalarT> 
typename std::enable_if<has_val<ScalarT>::value,decltype(std::declval<ScalarT>().val())>::type 
value( const ScalarT &f ) {
  return ScalarT(f.val());  
}


template<class ScalarT> 
typename std::enable_if<!has_val<ScalarT>::value,const ScalarT>::type 
value( const ScalarT &f ) {
  return f;  
}


} // namespace AD




#endif // PDEOPT_TEMPLATE_TOOLS_HPP
