// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PDEOPT_TEMPLATE_TOOLS_HPP
#define PDEOPT_TEMPLATE_TOOLS_HPP

#include "Sacado.hpp"
#include <limits>
#include <type_traits>
#include <tuple>
#include <utility>


/*! \file  template_tools.hpp
    \brief Provides generic programming components with emphasis on 
           automatic differentiation
*/

namespace TemplateTools {

// Is there a Variadic solution?
#define DV(A)                std::declval<A>()
#define METHOD_0(A,B)        decltype( DV(A).B())
#define METHOD_1(A,B,C)      decltype( DV(A).B(C))


// Need to create this template parameter pack until C++17 is adopted
template<class...> struct voider { using type=void; };
template<class...Ts> using void_t = typename voider<Ts...>::type;
template<template<class...>class Z, class, class...Ts>
struct has_method_: public std::false_type {};

template<template<class...>class Z, class...Ts>
struct has_method_<Z, void_t<Z<Ts...>>, Ts...> : public std::true_type {};

template<template<class...>class Z, class...Ts>
using has_method=has_method_<Z, void, Ts...>;

template<class> struct element{};
template<template<class...>class Container, class T, class...Rest> 
  struct element<Container<T, Rest...>>{ using type=T; };

template<class Container> using element_t=typename element<Container>::type;

template<bool HasMethod, class TypeIf, class TypeElse> using CondType = 
    typename std::conditional<HasMethod, TypeIf, TypeElse>::type;


// Detect existence of .size()
template <class Container> using c_size_t = METHOD_0(Container,size); 
template <class Container> using has_size = has_method<c_size_t,Container>;

template<class Container>
using index_t = CondType<has_size<Container>::value,METHOD_0(Container,size),void>;

template<class Container>
typename std::enable_if<has_size<Container>::value,index_t<Container>>::type 
inline size( const Container &c ) { return c.size(); }


template<class Container>
typename std::enable_if<!has_size<Container>::value,index_t<Container>>::type 
inline size( const Container &c ) { return 0; }

// Detect subscript operator []
//template <class Container>
//using c_subscript_t = decltype(DV(Container)[c_size_t<Container>::type]);

//template <class Container>
//using has_subscript = has_method<c_subscript_t, Container>;

template <class T, class Index>
using subscript_t = decltype(std::declval<T>()[std::declval<Index>()]);

template <class T, class Index=decltype(std::declval<T>().size())>
using has_subscript = has_method<subscript_t, T, Index>;


template<class Container,class Index=index_t<Container>>
typename std::enable_if<has_subscript<Container,Index>::value,element_t<Container>>::type
inline get( const Container &c, const Index i ) {
  return c[i];
}

template<class Container,class Index=index_t<Container>>
typename std::enable_if<!has_subscript<Container>::value,void>::type
inline get( const Container &c, const Index i ) {
}


// Returns the element x_j such that x_j >= x_i for all i in [0,size-1]
template<class Container> 
typename std::enable_if<has_subscript<Container>::value,element_t<Container>>::type
inline max( const Container &c ) {
  auto result = std::numeric_limits<element_t<Container>>::lowest();
  for( index_t<Container> i=0; i<size(c); ++i ) {
    auto val = get(c,i);
    result = val > result ? val : result;
  } 
  return result;
}

template<class Container> 
typename std::enable_if<!has_subscript<Container>::value,void>::type
inline max( const Container &c ) {
}


// Returns the element x_j such that x_j <= x_i for all i in [0,size-1]
template<class Container> 
typename std::enable_if<has_subscript<Container>::value,element_t<Container>>::type
inline min( const Container &c ) {
  auto result = std::numeric_limits<element_t<Container>>::max();
  for( index_t<Container> i=0; i<size(c); ++i ) {
    auto val = get(c,i);
    result = val < result ? val : result;
  } 
  return result;
}

template<class Container> 
typename std::enable_if<!has_subscript<Container>::value,void>::type
inline min( const Container &c ) {
}



} // namespace TemplateTools

namespace AD {

  using namespace TemplateTools;

/*  ResultType<T...> deduces the type that would be needed to add 
 *  the constituent types together. For example if DFad = Sacado::Fad::DFad<Real>
 *  then ReturnType<Real,DFad> would be DFad since adding the two together
 *  would yield a DFad.
 */
template<typename T,typename ...Args>
using ResultType = typename std::common_type<T,Args...>::type;

// Detect existence of .val()
template <class ScalarT>
using val_t = METHOD_0(ScalarT,val); 

template <class ScalarT> 
using has_val = has_method<val_t,ScalarT>;

template<class ScalarT> 
typename std::enable_if<has_val<ScalarT>::value,METHOD_0(ScalarT,val)>::type 
inline value( ScalarT &f ) { return f.val(); }

// Return a copy if the argument is const 
template<class ScalarT> 
typename std::enable_if<has_val<ScalarT>::value,METHOD_0(ScalarT,val)>::type 
inline value( const ScalarT &f ) {
  return ScalarT(f.val());  
}

template<class ScalarT> 
typename std::enable_if<!has_val<ScalarT>::value,const ScalarT>::type 
inline value( const ScalarT &f ) {
  return f;  
}

template <class ScalarT> struct has_dx
{
  std::add_rvalue_reference<int>::type arg = DV(int);

  // checks whether there is a method .dx(int)
  template <typename C> static constexpr decltype(DV(C).dx(arg), bool()) test(int)
  { return true;  }
  template <typename C> static constexpr bool test(...)
  { return false; }
  static constexpr bool value = test<ScalarT>(int());
}; 

/* // TODO: Convert this the above approach
template <class ScalarT>
using dx_t = METHOD_1(ScalarT,dx,int); 

template <class ScalarT> 
using has_dx = has_method<dx_t,ScalarT, int>;
*/

/* Automatic differentiation function for arbitrary types
 *
 * If the given object has a .dx(int) function, call it, otherwise
 * return 0
 */
template<class ScalarT> 
typename std::enable_if<has_dx<ScalarT>::value,METHOD_1(ScalarT,dx,0)>::type 
derivative( const ScalarT &f, int i ) {
  return f.dx(i);  // Does not check validity of i value!
}

// If the method does not exist, return integer 0
template<class ScalarT> 
typename std::enable_if<!has_dx<ScalarT>::value,int>::type 
derivative( const ScalarT &f, int i ) {
  return int(0);  
}

template<class Param, template<class> class Array>
struct ScalarEvaluator {

  using Real  = element_t<Param>;

  template<template<class,class> class ScalarFunction>
  static Real eval( const ScalarFunction<Param, Array<Real>> &F, const Param &param, 
                    const Array<Real> &x, const Array<int> &deriv ) {
    return 0;
  }


};


} // namespace AD




#endif // PDEOPT_TEMPLATE_TOOLS_HPP
