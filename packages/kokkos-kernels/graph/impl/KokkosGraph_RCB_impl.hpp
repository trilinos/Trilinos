// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSGRAPH_RCB_IMPL_HPP
#define KOKKOSGRAPH_RCB_IMPL_HPP

#include "Kokkos_Core.hpp"
#include "KokkosKernels_Utils.hpp"
#include <vector>
#include <algorithm>

namespace KokkosGraph {
namespace Impl {

template <typename perm_view_type>
struct FillOneIncrementFunctor {
  using ordinal_t = typename perm_view_type::value_type;
  perm_view_type A;

  FillOneIncrementFunctor(perm_view_type &A_) : A(A_) {}
  KOKKOS_INLINE_FUNCTION void operator()(ordinal_t i) const { A(i) = i; }
};

template <typename reducer_type, typename view_type, typename ordinal_t>
struct MinMaxReducerFunctor {
  using reducer_value_type = typename reducer_type::value_type;
  view_type A;
  MinMaxReducerFunctor(const view_type &A_) : A(A_) {}
  KOKKOS_INLINE_FUNCTION void operator()(ordinal_t i, reducer_value_type &lminmax) const {
    typename view_type::value_type val = A(i);
    if (val < lminmax.min_val) lminmax.min_val = val;
    if (val > lminmax.max_val) lminmax.max_val = val;
  }
};

template <typename perm_view_type1, typename perm_view_type2, typename coors_view_type>
struct UpdatePermAndMeshFunctor {
  using ordinal_t = typename perm_view_type1::value_type;
  perm_view_type1 reverse_perm_bisect;  // a subview
  perm_view_type1 prev_reverse_perm;    // a subview
  perm_view_type2 perm;                 // a full-length view
  perm_view_type2 reverse_perm;         // a full-length view
  coors_view_type coordinates_orig;     // a subview
  coors_view_type coordinates_new;      // a subview
  ordinal_t p1_size;
  ordinal_t offset;
  ordinal_t N;  // length of reverse_perm_bisect
  ordinal_t ndim;

  UpdatePermAndMeshFunctor(const perm_view_type1 &reverse_perm_bisect_, const perm_view_type1 &prev_reverse_perm_,
                           perm_view_type2 &perm_, perm_view_type2 &reverse_perm_,
                           const coors_view_type &coordinates_orig_, coors_view_type &coordinates_new_,
                           const ordinal_t &p1_size_, const ordinal_t &offset_)
      : reverse_perm_bisect(reverse_perm_bisect_),
        prev_reverse_perm(prev_reverse_perm_),
        perm(perm_),
        reverse_perm(reverse_perm_),
        coordinates_orig(coordinates_orig_),
        coordinates_new(coordinates_new_),
        p1_size(p1_size_),
        offset(offset_) {
    N    = static_cast<ordinal_t>(reverse_perm_bisect.extent(0));
    ndim = static_cast<ordinal_t>(coordinates_orig.extent(1));
  }
  KOKKOS_INLINE_FUNCTION void operator()(ordinal_t i) const {
    // orig_lcl_idx: 0 --> (N-1)
    ordinal_t orig_lcl_idx = reverse_perm_bisect(i);

    // new_lcl_idx: 0 --> (N-1)
    ordinal_t new_lcl_idx;
    if (i < p1_size) {
      new_lcl_idx = i;
    } else {
      new_lcl_idx = (N - 1 - i) + p1_size;
    }
    // Calculate new_gbl_idx by adding an offset
    ordinal_t new_gbl_idx = new_lcl_idx + offset;

    // Retrieve gbl_orig_idx
    ordinal_t gbl_orig_idx = prev_reverse_perm(orig_lcl_idx);

    // Update perm at gbl_orig_idx location
    perm(gbl_orig_idx) = new_gbl_idx;

    // Update reverse_perm at new_gbl_idx location
    reverse_perm(new_gbl_idx) = gbl_orig_idx;

    // Update coordinates with bisecting results
    for (ordinal_t j = 0; j < ndim; j++) {
      coordinates_new(new_lcl_idx, j) = coordinates_orig(orig_lcl_idx, j);
    }
  }
};

template <typename view_type, typename value_type>
void find_min_max(const view_type &A, value_type &min_val, value_type &max_val) {
  using execution_space = typename view_type::device_type::execution_space;
  using reducer_type    = Kokkos::MinMax<value_type>;
  size_t n_elements     = static_cast<size_t>(A.extent(0));
  if (n_elements == 0) {
    min_val = static_cast<value_type>(0);
    max_val = static_cast<value_type>(0);
    return;
  }
  typename reducer_type::value_type result;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<execution_space>(0, n_elements),
                          MinMaxReducerFunctor<reducer_type, view_type, size_t>(A), reducer_type(result));
  min_val = result.min_val;
  max_val = result.max_val;
}

/**
 * @brief Bisect and assign partition indices to coordinate list
 */
template <typename coors_view_type, typename value_type, typename index_view_type, typename ordinal_type>
inline void bisect(const coors_view_type &coors_1d, const value_type &init_min_val, const value_type &init_max_val,
                   index_view_type &reverse_perm_bisect, ordinal_type &p1_size, ordinal_type &p2_size) {
  const ordinal_type N = static_cast<ordinal_type>(coors_1d.extent(0));
  value_type min_val   = init_min_val;
  value_type max_val   = init_max_val;
  value_type p1_weight, p2_weight;
  value_type weight_ratio;
  const int max_bisection_steps = 11;
  // For now, limit the number of times finding mid point to ten to make RCB work with Sierra T/F coordinates
  // TODO: - allow users to pass max_bisection_steps as an input parameter
  //       - or switch to use the median, which is probably a better solution
  for (int bisection_step = 0; bisection_step < max_bisection_steps; ++bisection_step) {
    value_type mid_point = (max_val + min_val) / 2.0;
    p1_size              = 0;
    p2_size              = 0;
    // Use one array "reverse_perm_bisect" to store indices of both partitions
    for (ordinal_type i = 0; i < N; i++) {
      if (coors_1d(i) > mid_point) {  // partition 1: store forward
        reverse_perm_bisect(p1_size) = i;
        p1_size++;
      } else {  // partition 2: store backward
        reverse_perm_bisect(N - 1 - p2_size) = i;
        p2_size++;
      }
    }
    p1_weight = static_cast<value_type>(p1_size);
    p2_weight = static_cast<value_type>(p2_size);

    weight_ratio = std::max(p1_weight, p2_weight) / std::min(p1_weight, p2_weight);

    if (weight_ratio < 1.1)
      break;
    else {
      // Update min_val or max_val to calculate a new mid_point
      // Idea: shift mid_point to the heavier partition
      if (p1_weight > p2_weight)
        min_val = mid_point;
      else
        max_val = mid_point;
    }
  }
}

/**
 * @brief Recursive coordinate bisection on the coordinate list
 */
template <typename coors_view_type, typename perm_view_type>
std::vector<typename perm_view_type::value_type> rcb(coors_view_type &coordinates, perm_view_type &perm,
                                                     perm_view_type &reverse_perm, const int &n_levels) {
  using execution_space = typename coors_view_type::device_type::execution_space;
  using scalar_t        = typename coors_view_type::value_type;
  using ordinal_type    = typename perm_view_type::value_type;

  const ordinal_type N    = static_cast<ordinal_type>(coordinates.extent(0));
  const ordinal_type ndim = static_cast<ordinal_type>(coordinates.extent(1));

  // Allocate coordinates views on device memory
  coors_view_type coordinates_bisect(Kokkos::view_alloc(Kokkos::WithoutInitializing, "coordinates_bisect"), N, ndim);
  perm_view_type reverse_perm_bisect(Kokkos::view_alloc(Kokkos::WithoutInitializing, "reverse_perm_bisect"), N);
  perm_view_type prev_reverse_perm(Kokkos::view_alloc(Kokkos::WithoutInitializing, "prev_reverse_perm"), N);

  // Create host mirrors of device views
  typename coors_view_type::host_mirror_type h_coordinates        = Kokkos::create_mirror_view(coordinates);
  typename perm_view_type::host_mirror_type h_perm                = Kokkos::create_mirror_view(perm);
  typename perm_view_type::host_mirror_type h_reverse_perm        = Kokkos::create_mirror_view(reverse_perm);
  typename perm_view_type::host_mirror_type h_reverse_perm_bisect = Kokkos::create_mirror_view(reverse_perm_bisect);

  // Copy coordinates from device memory to host memory because bisecting is currently executed on host
  Kokkos::deep_copy(h_coordinates, coordinates);

  // Initialize
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space, ordinal_type>(0, N),
                       FillOneIncrementFunctor<perm_view_type>(perm));
  Kokkos::deep_copy(reverse_perm, perm);

  ordinal_type n_partitions =
      1;  // number of partitions at the previous level (initial value is 1, i.e., starting with the entire mesh points)
  const ordinal_type max_n_partitions = static_cast<ordinal_type>(std::pow(2, n_levels - 1));
  std::vector<ordinal_type> partition_sizes(
      max_n_partitions);   // contain the number of basis functions (or elements) per partition in the previous level
  partition_sizes[0] = N;  // starting with the entire mesh points
  std::vector<ordinal_type> partition_sizes_tmp(max_n_partitions);

  // Start RCB
  for (int lvl = 1; lvl < n_levels; lvl++) {  // skip level 0 and start from level 1
    ordinal_type coordinates_offset = 0;      // always start from beginning of the mesh points
    ordinal_type cnt_partitions     = 0;

    // Keep a copy of reverse permutation list
    Kokkos::deep_copy(prev_reverse_perm, reverse_perm);

    for (ordinal_type p = 0; p < n_partitions; p++) {  // go through each partition of previous level and do bisecting
      if (p > 0) {
        // Calculate coordinates offset of the current partition based on the previous partition
        coordinates_offset += partition_sizes[p - 1];
      }

      ordinal_type N0      = partition_sizes[p];  // partition size (or length)
      ordinal_type p1_size = 0;
      ordinal_type p2_size = 0;
      auto sub_coordinates =
          Kokkos::subview(coordinates, Kokkos::make_pair(coordinates_offset, coordinates_offset + N0), Kokkos::ALL());
      auto sub_h_coordinates =
          Kokkos::subview(h_coordinates, Kokkos::make_pair(coordinates_offset, coordinates_offset + N0), Kokkos::ALL());
      auto sub_coordinates_bisect = Kokkos::subview(
          coordinates_bisect, Kokkos::make_pair(coordinates_offset, coordinates_offset + N0), Kokkos::ALL());
      auto sub_reverse_perm_bisect =
          Kokkos::subview(reverse_perm_bisect, Kokkos::make_pair(coordinates_offset, coordinates_offset + N0));
      auto sub_prev_reverse_perm =
          Kokkos::subview(prev_reverse_perm, Kokkos::make_pair(coordinates_offset, coordinates_offset + N0));
      auto sub_h_reverse_perm_bisect =
          Kokkos::subview(h_reverse_perm_bisect, Kokkos::make_pair(coordinates_offset, coordinates_offset + N0));

      // Find min, max, and span of each dimension
      scalar_t x_min, x_max, y_min, y_max, z_min, z_max;
      scalar_t x_span = 0.0;
      scalar_t y_span = 0.0;
      scalar_t z_span = 0.0;

      auto x_coors = Kokkos::subview(sub_coordinates, Kokkos::ALL(), 0);
      find_min_max(x_coors, x_min, x_max);
      x_span = x_max - x_min;

      if (ndim > 1) {
        auto y_coors = Kokkos::subview(sub_coordinates, Kokkos::ALL(), 1);
        find_min_max(y_coors, y_min, y_max);
        y_span = y_max - y_min;
      }

      if (ndim > 2) {
        auto z_coors = Kokkos::subview(sub_coordinates, Kokkos::ALL(), 2);
        find_min_max(z_coors, z_min, z_max);
        z_span = z_max - z_min;
      }

      // Bisect partition on the most elongated dimension (host execution, for now)
      if ((x_span >= y_span) && (x_span >= z_span)) {
        auto h_x_coors = Kokkos::subview(sub_h_coordinates, Kokkos::ALL(), 0);
        bisect(h_x_coors, x_min, x_max, sub_h_reverse_perm_bisect, p1_size, p2_size);
      } else if ((y_span >= x_span) && (y_span >= z_span)) {
        auto h_y_coors = Kokkos::subview(sub_h_coordinates, Kokkos::ALL(), 1);
        bisect(h_y_coors, y_min, y_max, sub_h_reverse_perm_bisect, p1_size, p2_size);
      } else {
        auto h_z_coors = Kokkos::subview(sub_h_coordinates, Kokkos::ALL(), 2);
        bisect(h_z_coors, z_min, z_max, sub_h_reverse_perm_bisect, p1_size, p2_size);
      }

      Kokkos::deep_copy(sub_reverse_perm_bisect, sub_h_reverse_perm_bisect);

      // Update global permutation and reverse permutation lists and shuffle coordinates using bisecting results
      UpdatePermAndMeshFunctor<decltype(sub_reverse_perm_bisect), perm_view_type, decltype(sub_coordinates)> func(
          sub_reverse_perm_bisect, sub_prev_reverse_perm, perm, reverse_perm, sub_coordinates, sub_coordinates_bisect,
          p1_size, coordinates_offset);
      Kokkos::RangePolicy<execution_space, ordinal_type> policy(0, N0);
      Kokkos::parallel_for(policy, func);

      if (p1_size != 0) {
        partition_sizes_tmp[cnt_partitions] = p1_size;
        cnt_partitions++;
      }
      if (p2_size != 0) {
        partition_sizes_tmp[cnt_partitions] = p2_size;
        cnt_partitions++;
      }
    }  // end Partition loop

    // Update coordinates
    Kokkos::deep_copy(coordinates, coordinates_bisect);
    Kokkos::deep_copy(h_coordinates, coordinates_bisect);

    // Update the number of partitions of this level (used for bisections in the next level)
    n_partitions = cnt_partitions;

    // Update partition sizes of this level (used for bisections in the next level)
    std::copy(partition_sizes_tmp.begin(), partition_sizes_tmp.end(), partition_sizes.begin());
  }  // end Level loop

  if (n_partitions < max_n_partitions) partition_sizes.resize(n_partitions);

  return partition_sizes;
}

}  // namespace Impl
}  // namespace KokkosGraph
#endif
