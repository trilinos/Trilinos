#pragma once
#ifndef ROL2_TYPEDEDUCTION_HPP
#define ROL2_TYPEDEDUCTION_HPP

namespace ROL2 {

#if __cplusplus >= 201703L
using std::void_t;
#else
template<typename...> struct make_void { using type = void; }
template<typename...Ts> using void_t = typename make_void<Ts...>:type;
#endif

namespace detail {

template<template<typename...> class Z, typename, typename...Ts>
struct has_member : std::false_type {};

template<template<...> class Z, typename...Ts>
struct has_member<Z, std::void_t<Z<Ts...>>, Ts...> : std::true_type {};

} // namespace detail

template<template<typename...> class Z, typename...Ts>
using has_member = details::has_member<Z,void,Ts...>;

template<template<typename...> class Z, typename...Ts>
constexpr bool has_member_v = has_member<Z,Ts...>::value;

} // namespace ROL2

#endif //ROL2_TYPEDEDUCTION_HPP

