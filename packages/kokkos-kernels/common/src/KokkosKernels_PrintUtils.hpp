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

#ifndef _KOKKOSKERNELS_PRINTUTILS_HPP
#define _KOKKOSKERNELS_PRINTUTILS_HPP
#include "Kokkos_Core.hpp"
#include <ostream>

namespace KokkosKernels {

namespace Impl {

template <typename in_lno_view_t, typename out_lno_view_t>
struct Histogram {
  in_lno_view_t inview;
  out_lno_view_t outview;
  Histogram(in_lno_view_t inview_, out_lno_view_t outview_) : inview(inview_), outview(outview_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t& ii) const {
    typedef typename std::remove_reference<decltype(outview(0))>::type atomic_incr_type;
    Kokkos::atomic_fetch_add(&(outview(inview(ii))), atomic_incr_type(1));
  }
};

/**
 * \brief given an integer input view, it fills the histogram that
 * represents how many of each value exists.
 * \param in_elements: number of the elements in input view.
 * \param in_view: the input view. Has to be integer-like.
 * \param histogram: the output histogram. User is responsible from initializing
 * them with 0, and size must be big enough to hold all values in input view.
 */
template <typename in_lno_view_t, typename out_lno_view_t, typename MyExecSpace>
inline void kk_get_histogram(typename in_lno_view_t::size_type in_elements, in_lno_view_t in_view,
                             out_lno_view_t histogram /*must be initialized with 0s*/) {
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_for("KokkosKernels::Common::GetHistogram", my_exec_space(0, in_elements),
                       Histogram<in_lno_view_t, out_lno_view_t>(in_view, histogram));
  MyExecSpace().fence();
}

/**
 * \brief Prints the given 1D view.
 * \param os: Stream to print to. To print to stdout use std::cout, stderr,
 * std::cerr, or a file use an ofstream object. \param view: input view to
 * print. \param print_all: whether to print all elements or not. If it is
 * false, print print_size/2 first and last elements. \param sep: Element
 * separator. Default is a single space: " " \param print_size: Total elements
 * to print if print_all is false print_size/2 first and last elements are
 * pritned. This parameter is not used if print_all is set to true.
 */
template <typename idx_array_type>
inline std::enable_if_t<idx_array_type::rank <= 1> kk_print_1Dview(std::ostream& os, idx_array_type view,
                                                                   bool print_all = false, const char* sep = " ",
                                                                   size_t print_size = 40) {
  typedef typename idx_array_type::HostMirror host_type;
  typedef typename idx_array_type::size_type idx;
  host_type host_view = Kokkos::create_mirror_view(view);
  Kokkos::deep_copy(host_view, view);
  const auto print_range = [&](idx begin, idx end) {
    for (idx i = begin; i < end; ++i) os << host_view.access(i) << sep;
  };
  idx nr = host_view.extent(0);
  if (print_all || nr <= print_size) {
    print_range(0, nr);
  } else {
    idx n = print_size / 2;
    print_range(0, n);
    os << "... ... ..." << sep;
    print_range(nr - n, nr);
  }
  os << std::endl;
}

/**
 * \brief Multi-vector variant (see rank-1 for param description). Prints Nx1
 * rank-2 vectors same like rank-1 vectors and prints multi-vector dimensions.
 */
template <typename idx_array_type>
inline std::enable_if_t<idx_array_type::rank >= 2> kk_print_1Dview(std::ostream& os, idx_array_type view,
                                                                   bool print_all = false, const char* sep = " ",
                                                                   size_t print_size = 40) {
  if (idx_array_type::rank == 2 && view.extent(1) == 1) {
    kk_print_1Dview(os, subview(view, Kokkos::ALL, 0), print_all, sep, print_size);
    return;
  }
  os << "[" << view.extent(0);
  // ::rank is a Kokkos::...::integral_constant, not appropriate for `i`
  for (int i = 1; i < int(idx_array_type::rank); ++i) {
    os << "x" << view.extent(i);
  }
  os << " multi-vector]" << std::endl;
}

/**
 * \brief Prints the given 1D view.
 * \param view: input view to print.
 * \param print_all: whether to print all elements or not. If it is false,
 * only first and last 20 elements are printed.
 *
 * This interface is provided for backwards compatiblity.
 */
template <typename idx_array_type>
inline void kk_print_1Dview(idx_array_type view, bool print_all = false, size_t print_size = 40) {
  kk_print_1Dview(std::cout, view, print_all, " ", print_size);
}

}  // namespace Impl
}  // namespace KokkosKernels

#endif
