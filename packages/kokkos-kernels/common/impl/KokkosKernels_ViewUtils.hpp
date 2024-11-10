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

#ifndef KOKKOSKERNELS_VIEWUTILS_HPP
#define KOKKOSKERNELS_VIEWUTILS_HPP
#include "Kokkos_Core.hpp"

namespace KokkosKernels::Impl {

/*! \brief Yields a type that is View with Kokkos::Unmanaged added to the memory
 * traits
 */
template <typename View>
class with_unmanaged {
  using data_type    = typename View::data_type;
  using layout_type  = typename View::array_layout;
  using memory_space = typename View::memory_space;

  using orig_traits                    = typename View::memory_traits;
  static constexpr unsigned new_traits = orig_traits::impl_value | Kokkos::Unmanaged;

 public:
  using type = Kokkos::View<data_type, layout_type, memory_space, Kokkos::MemoryTraits<new_traits> >;
};

/*! \brief A type that is View with Kokkos::Unmanaged added to the memory traits

    \tparam View the type to add Kokkos::Unmanaged to
 */
template <typename View>
using with_unmanaged_t = typename with_unmanaged<View>::type;

/*! \brief Returns an unmanaged version of v

    \tparam View the type of the input view v
 */
template <typename View>
auto make_unmanaged(const View &v) {
  return typename with_unmanaged<View>::type(v);
}

}  // namespace KokkosKernels::Impl

#endif
