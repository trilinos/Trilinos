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

#include "Kokkos_Core.hpp"
#ifndef _KOKKOSKERNELS_VECTORUTILS_HPP
#define _KOKKOSKERNELS_VECTORUTILS_HPP

namespace KokkosKernels {

namespace Impl {

template <typename out_array_t, typename in_array_t, typename scalar_1, typename scalar_2>
struct A_times_X_plus_B {
  out_array_t out_view;
  in_array_t in_view;
  const scalar_1 a;
  const scalar_2 b;
  A_times_X_plus_B(out_array_t out_view_, in_array_t in_view_, scalar_1 a_, scalar_2 b_)
      : out_view(out_view_), in_view(in_view_), a(a_), b(b_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii) const { out_view(ii) = in_view(ii) * a + b; }
};

template <typename out_array_type, typename in_array_type>
struct ModularView {
  typedef typename in_array_type::value_type vt;
  out_array_type out_view;
  in_array_type in_view;
  const int modular_constant;
  ModularView(out_array_type out_view_, in_array_type in_view_, int mod_factor_)
      : out_view(out_view_), in_view(in_view_), modular_constant(mod_factor_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii) const { out_view(ii) = in_view(ii) % modular_constant; }
};

template <typename from_vector, typename to_vector>
struct CopyVectorFunctor {
  from_vector from;
  to_vector to;

  CopyVectorFunctor(from_vector &from_, to_vector to_) : from(from_), to(to_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i) const { to[i] = from[i]; }
};

/**
 * \brief given a input view in_arr, sets the corresponding index of out_arr by
 * out_arr(ii) = in_arr(ii) * a + b;
 * \param num_elements: number of elements in input and output arrays.
 * \param out_arr: output arr, can be same as input array.
 * \param in_arr: input arr.
 * \param a: scalar for multiplication
 * \param b: scalar for addition
 */
template <typename out_array_t, typename in_array_t, typename scalar_1, typename scalar_2, typename MyExecSpace>
inline void kk_a_times_x_plus_b(typename in_array_t::value_type num_elements, out_array_t out_arr, in_array_t in_arr,
                                scalar_1 a, scalar_2 b) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_for("KokkosKernels::Common::ATimesXPlusB", my_exec_space(0, num_elements),
                       A_times_X_plus_B<out_array_t, in_array_t, scalar_1, scalar_2>(out_arr, in_arr, a, b));
}

/**
 * \brief calculates the modular of each entry input array and writes it to
 * corresponding vector. \param num_elements: number of elements in input and
 * output arrays. \param out_arr: output arr, can be same as input array. \param
 * in_arr: input arr. \param mod_factor_: for what value the modular will be
 * applied.
 */
template <typename out_array_type, typename in_array_type, typename MyExecSpace>
inline void kk_modular_view(typename in_array_type::value_type num_elements, out_array_type out_arr,
                            in_array_type in_arr, int mod_factor_) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_for("KokkosKernels::Common::ModularView", my_exec_space(0, num_elements),
                       ModularView<out_array_type, in_array_type>(out_arr, in_arr, mod_factor_));
}

template <typename from_vector, typename to_vector, typename MyExecSpace>
void kk_copy_vector(size_t num_elements, from_vector from, to_vector to) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_for("KokkosKernels::Common::CopyVector", my_exec_space(0, num_elements),
                       CopyVectorFunctor<from_vector, to_vector>(from, to));
}
}  // namespace Impl
}  // namespace KokkosKernels

#endif
