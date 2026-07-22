// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSGRAPH_RCB_HPP
#define KOKKOSGRAPH_RCB_HPP

#include "KokkosGraph_RCB_impl.hpp"

namespace KokkosGraph {
namespace Experimental {

// The recursive coordinate bisection (RCB) algorithm partitions the graph
// according to the coordinates of the mesh points. This function returns
// a vector containing sizes of partitions, the coordinate list (organized in RCB
// order), a permutation array describing the mapping from the original order
// to RCB order, and a reverse permutation array describing the mapping from
// the RCB order to the original order

template <typename coors_view_type, typename perm_view_type>
std::vector<typename perm_view_type::value_type> recursive_coordinate_bisection(coors_view_type &coordinates,
                                                                                perm_view_type &perm,
                                                                                perm_view_type &reverse_perm,
                                                                                const int &n_levels) {
  static_assert(Kokkos::is_view_v<coors_view_type>,
                "KokkosGraph::Experimental::recursive_coordinate_bisection: coors_view_type must be a Kokkos::View.");
  static_assert(Kokkos::is_view_v<perm_view_type>,
                "KokkosGraph::Experimental::recursive_coordinate_bisection: perm_view_type must be a Kokkos::View.");

  static_assert(static_cast<int>(coors_view_type::rank()) == 2,
                "Kokkos::Experimental::recursive_coordinate_bisection: coors_view_type must have rank 2.");
  static_assert(static_cast<int>(perm_view_type::rank()) == 1,
                "Kokkos::Experimental::recursive_coordinate_bisection: perm_view_type must have rank 1.");

  if (n_levels < 2) {
    std::ostringstream os;
    os << "KokkosGraph::Experimental::recursive_coordinate_bisection only works with more than 1 level of bisection "
          "(i.e., 2 partitions).";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  if (static_cast<int>(coordinates.extent(1)) > 3) {
    std::ostringstream os;
    os << "KokkosGraph::Experimental::recursive_coordinate_bisection currently only supports 1-D, 2-D, or 3-D "
          "coordinates.";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  return KokkosGraph::Impl::rcb<coors_view_type, perm_view_type>(coordinates, perm, reverse_perm, n_levels);
}

}  // namespace Experimental
}  // namespace KokkosGraph

#endif
