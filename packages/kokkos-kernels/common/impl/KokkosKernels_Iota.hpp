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

#ifndef _KOKKOSKERNELS_IOTA_HPP
#define _KOKKOSKERNELS_IOTA_HPP

#include <type_traits>

#include <Kokkos_Core.hpp>

#include "KokkosKernels_Error.hpp"

/*! \file KokkosKernels_Iota.hpp
 * Define an Iota struct that implements a small subset of Kokkos::View and
 * related utilities.
 */

namespace KokkosKernels {
namespace Impl {

/*! \class Iota
  \brief A class that mimics a small subset of Kokkos::View

  \tparam T the type returned by operator()
  \tparam SizeType a custom offset type

  \typedef size_type SizeType
  \typedef value_type T
  \typedef non_const_value_type non-const T
  \typedef device_type void
  \typedef data_type const value_type *
  \enum rank always 1

  Iota::operator() returns offset + i
  Meant to be used in place of a Kokkos::View where entry i holds i + offset.
  Unlike a Kokkos::View, Iota is not materialized in memory.

  Constructing with a size less than 0 yeilds a 0-size Iota
*/
template <typename T, typename SizeType = size_t>
class Iota {
 public:
  using size_type            = SizeType;
  using value_type           = T;
  using non_const_value_type = std::remove_const_t<value_type>;
  using device_type          = void;
  using data_type            = const value_type *;

  /*! \brief construct an Iota where iota(i) -> offset + i

      \param[in] size the number of entries
      \param[in] offset the offset of the first entry

      Constructing with size < 0 yeilds a 0-size Iota
  */
  KOKKOS_INLINE_FUNCTION
  constexpr Iota(const size_type &size, const value_type offset) : size_(size), offset_(offset) {
    if constexpr (std::is_signed_v<size_type>) {
      if (size_ < size_type(0)) {
        size_ = 0;
      }
    }
  }

  /*! \brief construct an Iota where iota(i) ->  i

      \param[in] size the number of entries
  */
  KOKKOS_INLINE_FUNCTION
  explicit constexpr Iota(const size_type &size) : Iota(size, 0) {}

  /*! \brief construct a zero-sized iota
   */
  KOKKOS_INLINE_FUNCTION
  constexpr Iota() : size_(0), offset_(0) {}

  /*! \brief Construct Iota subview

    Like the Kokkos::View 1D subview constructor:
    \verbatim
    Kokkos::View a(10); // size = 10
    Kokkos::View b(a, Kokkos::pair{3,7}); // entries 3,4,5,6 of a

    Iota a(10);
    Iota b(a, Kokkos::pair{3,7}); // entries // 3,4,5,6 of a
    \endverbatim

    Creating a subview outside of the base Iota yeilds undefined behavior
  */
  template <typename P1, typename P2>
  KOKKOS_INLINE_FUNCTION constexpr Iota(const Iota &base, const Kokkos::pair<P1, P2> &range)
      : Iota(range.second - range.first, base.offset_ + range.first) {}

  /*! \brief Construct Iota subview

     i >= size() or i < 0 yields undefined behavior.
  */
  KOKKOS_INLINE_FUNCTION
  constexpr T operator()(size_type i) const noexcept { return value_type(i + offset_); };

  /// \brief return the size of the iota
  KOKKOS_INLINE_FUNCTION
  constexpr size_t size() const noexcept { return size_; }

  /// \brief Iotas are always like a rank-1 Kokkos::View
  enum { rank = 1 };

 private:
  size_type size_;
  value_type offset_;
};

/// \class is_iota
/// \brief is_iota<T>::value is true if T is a Iota<...>, false otherwise
template <typename>
struct is_iota : public std::false_type {};
template <typename... P>
struct is_iota<Iota<P...>> : public std::true_type {};
template <typename... P>
struct is_iota<const Iota<P...>> : public std::true_type {};
template <typename... P>
inline constexpr bool is_iota_v = is_iota<P...>::value;

}  // namespace Impl
}  // namespace KokkosKernels

#endif  // _KOKKOSKERNELS_IOTA_HPP
