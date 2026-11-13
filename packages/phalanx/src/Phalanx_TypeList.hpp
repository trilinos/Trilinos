// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_TYPE_LIST_HPP
#define PHALANX_TYPE_LIST_HPP

#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "Phalanx_config.hpp"

// Only needed for interoperability with Sacado::mpl during the
// deprecation period (see the end of this file).
#include "Sacado_mpl_vector.hpp"

/** \file Phalanx_TypeList.hpp
    \brief Compile-time type list and algorithms used by Phalanx
    (e.g. for Traits::EvalTypes and the scalar types of each
    evaluation type).

    The query algorithms (size_v, at_t, contains_v, find_v, found_v,
    find_if_v, for_each_type) accept either a PHX::mpl::type_list or a
    Sacado::mpl::vector. Passing a Sacado::mpl::vector is deprecated.
    The algorithms that build new lists (push_back_t, remove_t, ...)
    only accept PHX::mpl::type_list.

    For backward compatibility, Sacado::mpl::size, find, at and
    for_each also accept a PHX::mpl::type_list, so existing code that
    applies them to a Traits::EvalTypes keeps compiling.
*/

namespace PHX::mpl {

// ---- Core type list --------------------------------------------------
//! Tag used to dispatch Sacado::mpl algorithms on a type_list.
struct type_list_tag {};

template <typename... Ts>
struct type_list {
  static constexpr std::size_t size = sizeof...(Ts);
  using tag = type_list_tag;
};

// ---- Concept: "is a type_list" ---------------------------------------
template <typename T>
struct is_type_list : std::false_type {};

template <typename... Ts>
struct is_type_list<type_list<Ts...>> : std::true_type {};

template <typename T>
concept TypeList = is_type_list<T>::value;

// ---- Conversion to type_list -------------------------------------------
//! Converts a supported type sequence to a PHX::mpl::type_list.
template <typename L>
struct to_type_list;

template <typename... Ts>
struct to_type_list<type_list<Ts...>> {
  using type = type_list<Ts...>;
};

// Named so the deprecation warning is self-explanatory even when the
// compiler drops the message.
template <typename... Ts>
struct PHALANX_DEPRECATED_MSG("Use PHX::mpl::type_list instead of Sacado::mpl::vector for Phalanx type sequences (e.g. Traits::EvalTypes and eval_scalar_types).")
Sacado_mpl_vector_is_deprecated_use_PHX_mpl_type_list {
  using type = type_list<Ts...>;
};

template <typename... Ts>
struct to_type_list<Sacado::mpl::vector<Ts...>> {
  using type = typename Sacado_mpl_vector_is_deprecated_use_PHX_mpl_type_list<Ts...>::type;
};

template <typename L>
using to_type_list_t = typename to_type_list<L>::type;

namespace detail {
  template <typename T>
  struct is_sacado_vector : std::false_type {};

  template <typename... Ts>
  struct is_sacado_vector<Sacado::mpl::vector<Ts...>> : std::true_type {};
}

//! A type list, or a type sequence that can be converted to one.
template <typename L>
concept TypeSequence = TypeList<L> || detail::is_sacado_vector<L>::value;

// ---- size --------------------------------------------------------------
template <TypeSequence L>
inline constexpr std::size_t size_v = to_type_list_t<L>::size;

// ---- at<I> (replacement for Sacado::mpl::at) ----------------------------
namespace detail {
  template <std::size_t I, typename... Ts>
  struct at_impl;

  template <std::size_t I, typename T, typename... Ts>
  struct at_impl<I, T, Ts...> : at_impl<I - 1, Ts...> {};

  template <typename T, typename... Ts>
  struct at_impl<0, T, Ts...> { using type = T; };

  template <std::size_t I, typename L>
  struct at;

  template <std::size_t I, typename... Ts>
  struct at<I, type_list<Ts...>> {
    using type = typename at_impl<I, Ts...>::type;
  };
}

template <TypeSequence L, std::size_t I>
  requires (I < size_v<L>)
using at_t = typename detail::at<I, to_type_list_t<L>>::type;

// ---- contains --------------------------------------------------------
namespace detail {
  template <typename T, typename L>
  struct contains;

  template <typename T, typename... Ts>
  struct contains<T, type_list<Ts...>>
    : std::bool_constant<(std::is_same_v<T, Ts> || ...)> {};
}

template <typename T, TypeSequence L>
inline constexpr bool contains_v = detail::contains<T, to_type_list_t<L>>::value;

// ---- find (replacement for Sacado::mpl::find) --------------------------
// Returns the index of the first occurrence of T in the list, or
// size_v<L> (see npos) if T isn't present.
namespace detail {
  template <typename T, typename... Ts>
  constexpr std::size_t find_index() {
    constexpr bool matches[] = { std::is_same_v<T, Ts>..., false }; // +false avoids zero-size array
    for (std::size_t i = 0; i < sizeof...(Ts); ++i) {
      if (matches[i]) return i;
    }
    return sizeof...(Ts);
  }

  template <typename T, typename L>
  struct find;

  template <typename T, typename... Ts>
  struct find<T, type_list<Ts...>> {
    static constexpr std::size_t value = find_index<T, Ts...>();
  };
}

template <typename T, TypeSequence L>
inline constexpr std::size_t find_v = detail::find<T, to_type_list_t<L>>::value;

//! Sentinel returned by find_v and find_if_v when no type matches.
template <TypeSequence L>
inline constexpr std::size_t npos = size_v<L>;

template <typename T, TypeSequence L>
inline constexpr bool found_v = find_v<T, L> != npos<L>;

/** \brief Index of T in L. Unlike find_v, a T that is not in L is a
    compile-time error. Use this when the index is used to access
    per-type storage. */
template <typename T, TypeSequence L>
  requires found_v<T, L>
inline constexpr std::size_t index_v = find_v<T, L>;

// ---- find_if (replacement for Sacado::mpl::find_if) ---------------------
// Pred is expected to expose Pred<T>::value, e.g. std::is_floating_point.
namespace detail {
  template <template <typename> class Pred, typename... Ts>
  constexpr std::size_t find_if_index() {
    constexpr bool matches[] = { Pred<Ts>::value..., false };
    for (std::size_t i = 0; i < sizeof...(Ts); ++i) {
      if (matches[i]) return i;
    }
    return sizeof...(Ts);
  }

  template <template <typename> class Pred, typename L>
  struct find_if;

  template <template <typename> class Pred, typename... Ts>
  struct find_if<Pred, type_list<Ts...>> {
    static constexpr std::size_t value = find_if_index<Pred, Ts...>();
  };
}

template <template <typename> class Pred, TypeSequence L>
inline constexpr std::size_t find_if_v = detail::find_if<Pred, to_type_list_t<L>>::value;

// ---- push_back (replacement for Sacado::mpl::push_back) ----------------
template <TypeList L, typename New>
struct push_back;

template <typename New, typename... Ts>
struct push_back<type_list<Ts...>, New> {
  using type = type_list<Ts..., New>;
};

template <TypeList L, typename New>
using push_back_t = typename push_back<L, New>::type;

// ---- push_front --------------------------------------------------------
template <TypeList L, typename New>
struct push_front;

template <typename New, typename... Ts>
struct push_front<type_list<Ts...>, New> {
  using type = type_list<New, Ts...>;
};

template <TypeList L, typename New>
using push_front_t = typename push_front<L, New>::type;

// ---- pop_front ---------------------------------------------------------
template <TypeList L>
  requires (L::size > 0)
struct pop_front;

template <typename T, typename... Ts>
struct pop_front<type_list<T, Ts...>> {
  using type = type_list<Ts...>;
};

template <TypeList L>
using pop_front_t = typename pop_front<L>::type;

// ---- pop_back ----------------------------------------------------------
template <TypeList L>
  requires (L::size > 0)
struct pop_back;

template <typename T>
struct pop_back<type_list<T>> {
  using type = type_list<>;
};

template <typename T, typename... Ts>
  requires (sizeof...(Ts) > 0)
struct pop_back<type_list<T, Ts...>> {
  using type = typename push_front<typename pop_back<type_list<Ts...>>::type, T>::type;
};

template <TypeList L>
using pop_back_t = typename pop_back<L>::type;

// ---- concat ------------------------------------------------------------
template <TypeList... Ls>
struct concat;

template <>
struct concat<> { using type = type_list<>; };

template <TypeList L>
struct concat<L> { using type = L; };

template <typename... As, typename... Bs, TypeList... Rest>
struct concat<type_list<As...>, type_list<Bs...>, Rest...> {
  using type = typename concat<type_list<As..., Bs...>, Rest...>::type;
};

template <TypeList... Ls>
using concat_t = typename concat<Ls...>::type;

// ---- remove (replacement for Sacado::mpl::remove/erase) ----------------
// Removes every occurrence of T from the list.
template <typename T, TypeList L>
struct remove;

template <typename T, typename... Ts>
struct remove<T, type_list<Ts...>> {
  using type = concat_t<
    std::conditional_t<std::is_same_v<T, Ts>, type_list<>, type_list<Ts>>...
  >;
};

template <typename T, TypeList L>
using remove_t = typename remove<T, L>::type;

// ---- unique ------------------------------------------------------------
// Deduplicates the list, keeping the first occurrence of each type in
// its original relative order.
template <TypeList L>
struct unique;

template <>
struct unique<type_list<>> { using type = type_list<>; };

template <typename T, typename... Ts>
struct unique<type_list<T, Ts...>> {
private:
  using rest             = typename unique<type_list<Ts...>>::type;
  using rest_without_dup = remove_t<T, rest>;
public:
  using type = typename push_front<rest_without_dup, T>::type;
};

template <TypeList L>
using unique_t = typename unique<L>::type;

// ---- for_each_type (replacement for Sacado::mpl::for_each) -------------
/** \brief Calls `f.template operator()<T>()` for each type T in L, in order.

    Typically used with a generic lambda that has an explicit template
    parameter:
    \code
    PHX::mpl::for_each_type<Traits::EvalTypes>([&]<typename EvalT>() { ... });
    \endcode
    Unlike Sacado::mpl::for_each, no object of type T is constructed.
*/
template <TypeSequence L, typename F>
constexpr void for_each_type(F&& f) {
  [&]<typename... Ts>(type_list<Ts...>) {
    (f.template operator()<Ts>(), ...);
  }(to_type_list_t<L>{});
}

// ---- for_each_index (replacement for for_each over Sacado::mpl::range_c)
//! Calls `f.template operator()<I>()` for each I in [0,N), in order.
template <std::size_t N, typename F>
constexpr void for_each_index(F&& f) {
  [&]<std::size_t... Is>(std::index_sequence<Is...>) {
    (f.template operator()<Is>(), ...);
  }(std::make_index_sequence<N>{});
}

} // namespace PHX::mpl

// ---- Sacado::mpl interoperability ----------------------------------------
// Lets Sacado::mpl::size, begin/end, find, at and for_each operate on a
// PHX::mpl::type_list, so downstream code that applies Sacado::mpl to a
// Traits::EvalTypes keeps compiling after the switch to type_list.
namespace Sacado::mpl {

  template <>
  struct size_impl<PHX::mpl::type_list_tag> {
    template <class L>
    struct apply {
      static const int value = static_cast<int>(L::size);
    };
  };

  template <>
  struct begin_impl<PHX::mpl::type_list_tag> {
    template <class L>
    struct apply {
      typedef vector_iterator<L,0> type;
    };
  };

  template <>
  struct end_impl<PHX::mpl::type_list_tag> {
    template <class L>
    struct apply {
      typedef vector_iterator<L,static_cast<int>(L::size)> type;
    };
  };

  template <int Pos>
  struct at_impl<PHX::mpl::type_list_tag, Pos> {
    template <class L>
    struct apply {
      typedef PHX::mpl::at_t<L,static_cast<std::size_t>(Pos)> type;
    };
  };

}

#endif
