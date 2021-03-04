#pragma once
#ifndef ROL2_TYPETRAITS_HPP
#define ROL2_TYPETRAITS_HPP

/**
 *    \file ROL2_TypeTraits.hpp
 *    \brief Alias or implement C++17 
 */

namespace ROL2 { 

#if __cplusplus >= 201703L // Have C++17

using std::conjunction;
using std::disjunction;
using std::void_t;

#else  // Implement C++17 type traits using C++14 

// Implement std::conjunction
template<class...> struct conjunction : std::true_type { };
template<class B1> struct conjunction<B1> : B1 { };
template<class B1, class... Bn>
struct conjunction<B1, Bn...> 
    : std::conditional_t<bool(B1::value), conjunction<Bn...>, B1> {};

// Implement std::disjunction
template<class...> struct disjunction : std::false_type { };
template<class B1> struct disjunction<B1> : B1 { };
template<class B1, class... Bn>
struct disjunction<B1, Bn...> 
    : std::conditional_t<bool(B1::value), B1, disjunction<Bn...>>  { };


namespace detail {
template<typename...> struct make_void { using type = void; };
} // namespace detail

template<typename...Ts> using void_t = typename detail::make_void<Ts...>::type;

#endif // __cplusplus >= 201703L


template<class... B>
constexpr bool conjunction_v = conjunction<B...>::value;

template<class... B>
constexpr bool disjunction_v = disjunction<B...>::value;

using std::enable_if_t;

using std::is_convertible;

template<typename From, typename To>
constexpr bool is_convertible_v = std::is_convertible<From,To>::value;

template<typename T>
constexpr bool is_enum_v = std::is_enum<T>::value;

using std::underlying_type_t;

} // namespace ROL2


#endif // ROL2_TYPETRAITS_HPP

