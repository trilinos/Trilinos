#pragma once
#ifndef ROL2_TYPEDEDUCTION_HPP
#define ROL2_TYPEDEDUCTION_HPP

namespace ROL2 {
namespace detail {

template<template<typename...> class Z, typename, typename...Ts>
struct has_member : std::false_type {};

template<template<typename...> class Z, typename...Ts>
struct has_member<Z, void_t<Z<Ts...>>, Ts...> : std::true_type {};

} // namespace detail

template<template<typename...> class Z, typename...Ts>
using has_member = detail::has_member<Z,void,Ts...>;

template<template<typename...> class Z, typename...Ts>
constexpr bool has_member_v = has_member<Z,Ts...>::value;

template<typename> struct function_traits;

template<typename ret_type, typename...arg_types>
struct function_traits<ret_type(arg_types...)> {

   static constexpr std::size_t arity = sizeof...(arg_types);

   template<std::size_t N>
   static constexpr auto argument( std::integral_constant<std::size_t,N> ) ->
                                   std::tuple_element_t<N,std::tuple<arg_types...>>;
};

template<typename ret_type, typename...arg_types>
struct function_traits<ret_type(*)(arg_types...)> :
       function_traits<ret_type(arg_types...)> {};

// member function pointer
template<class C, class R, class... Args>
struct function_traits<R(C::*)(Args...)> : function_traits<R(C&,Args...)> {};

// const member function pointer
template<class C, class R, class... Args>
struct function_traits<R(C::*)(Args...) const> : function_traits<R(C&,Args...)> {};

// member object pointer
template<class C, class R>
struct function_traits<R(C::*)> : function_traits<R(C&)> {};

template<std::size_t N,typename F>
using argument_t =
decltype( function_traits<F>::argument(
          std::declval<std::integral_constant<std::size_t,N>>() ) );

} // namespace ROL2

#endif //ROL2_TYPEDEDUCTION_HPP

