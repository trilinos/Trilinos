//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef _KOKKOSKERNELS_PREDICATES_HPP
#define _KOKKOSKERNELS_PREDICATES_HPP

#include "Kokkos_ArithTraits.hpp"

/*! \file KokkosKernels_Predicates.hpp
 * Define predicates for KokkosKernels search functions
 */

namespace KokkosKernels {

/**
 * @brief Struct template for a greater-than predicate
 * @tparam T Type to be compared
 */
template <typename T>
struct GT {
  using value_type = T;
  static_assert(!Kokkos::ArithTraits<T>::is_complex, "Please define custom predicates for ordering complex types");

  /**
   * @brief Return true if a is greater than b
   * @param a First value to be compared
   * @param b Second value to be compared
   */
  KOKKOS_INLINE_FUNCTION constexpr bool operator()(const value_type &a, const value_type &b) const noexcept {
    return a > b;
  }
};

/*! \brief "Greater-than-or-equal" predicate, a >= b
    \tparam T the type to compare
*/
template <typename T>
struct GTE {
  using value_type = T;
  static_assert(!Kokkos::ArithTraits<T>::is_complex, "Please define custom predicates for ordering complex types");

  /// \brief return a >= b
  KOKKOS_INLINE_FUNCTION constexpr bool operator()(const value_type &a, const value_type &b) const noexcept {
    return a >= b;
  }
};

/*! \brief "Less-than" predicate, a < b
    \tparam T the type to compare
*/
template <typename T>
struct LT {
  using value_type = T;
  static_assert(!Kokkos::ArithTraits<T>::is_complex, "Please define custom predicates for ordering complex types");

  /// \brief return a < b
  KOKKOS_INLINE_FUNCTION constexpr bool operator()(const value_type &a, const value_type &b) const noexcept {
    return a < b;
  }
};

/*! \brief "Less-than-or-equal" predicate, a <= b
    \tparam T the type to compare
*/
template <typename T>
struct LTE {
  using value_type = T;
  static_assert(!Kokkos::ArithTraits<T>::is_complex, "Please define custom predicates for ordering complex types");

  /// \brief return a <= b
  KOKKOS_INLINE_FUNCTION constexpr bool operator()(const value_type &a, const value_type &b) const noexcept {
    return a <= b;
  }
};

/*! \brief "Equal" predicate, a == b
    \tparam T the type to compare
*/
template <typename T>
struct Equal {
  using value_type = T;

  /// \brief return a == b
  KOKKOS_INLINE_FUNCTION constexpr bool operator()(const value_type &a, const value_type &b) const { return a == b; }
};

/**
 * @brief Struct template for inverting a predicate
 * @tparam Pred Predicate type to be inverted
 */
template <typename Pred>
struct Neg {
  using value_type = typename Pred::value_type;

  /**
   * @brief Constructor
   * @param pred Predicate object to be inverted
   */
  KOKKOS_INLINE_FUNCTION
  constexpr Neg(const Pred &pred) : pred_(pred) {}

  /**
   * @brief Return the boolean inverse of the underlying predicate
   * @param a First value to be compared by the predicate
   * @param b Second value to be compared by the predicate
   * @return Boolean inverse of the result of the predicate applied to a and b
   */
  KOKKOS_INLINE_FUNCTION constexpr bool operator()(const value_type &a, const value_type &b) const {
    return !pred_(a, b);
  }

 private:
  Pred pred_;  //< Underlying predicate object
};

/*! \brief Reflect a predicate, pred(b, a)
    \tparam Pred the type of the predicate to reflect
*/
template <typename Pred>
struct Refl {
  using value_type = typename Pred::value_type;

  KOKKOS_INLINE_FUNCTION
  constexpr Refl(const Pred &pred) : pred_(pred) {}

  /// \brief return the underlying binary predicate with reversed arguments
  KOKKOS_INLINE_FUNCTION constexpr bool operator()(const value_type &a, const value_type &b) const {
    return pred_(b, a);
  }

 private:
  Pred pred_;
};

}  // namespace KokkosKernels

#endif  // _KOKKOSKERNELS_PREDICATES_HPP