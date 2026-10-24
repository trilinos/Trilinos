// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSGRAPH_RCB_IMPL_HPP
#define KOKKOSGRAPH_RCB_IMPL_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_Sort.hpp"
#include "KokkosKernels_Utils.hpp"
#include <vector>
#include <algorithm>

namespace KokkosGraph {
namespace Impl {

template <typename value_type>
struct MinMaxPair {
  value_type min_val;
  value_type max_val;
};

template <typename perm_view_type>
struct FillOneIncrementFunctor {
  using ordinal_t = typename perm_view_type::value_type;
  perm_view_type A;

  FillOneIncrementFunctor(perm_view_type &A_) : A(A_) {}
  KOKKOS_INLINE_FUNCTION void operator()(ordinal_t i) const { A(i) = i; }
};

template <typename view_type, typename ordinal_t>
struct ColumnMinMaxFunctor {
  using scalar_type = typename view_type::value_type;
  using layout_type = typename view_type::array_layout;
  using device_type = typename view_type::device_type;
  using internal_view_type =
      Kokkos::View<const scalar_type **, layout_type, device_type, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
  using minmax_type = MinMaxPair<scalar_type>[];
  using value_type  = minmax_type;
  internal_view_type A;
  ordinal_t value_count;

  ColumnMinMaxFunctor(const view_type &A_, const ordinal_t &num_cols_) : A(A_), value_count(num_cols_) {}

  KOKKOS_INLINE_FUNCTION void init(minmax_type dst) const {
    for (ordinal_t c = 0; c < value_count; c++) {
      dst[c].min_val = Kokkos::reduction_identity<scalar_type>::min();
      dst[c].max_val = Kokkos::reduction_identity<scalar_type>::max();
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(ordinal_t i, minmax_type dst) const {
    for (ordinal_t c = 0; c < value_count; c++) {
      scalar_type val = A(i, c);
      if (val < dst[c].min_val) dst[c].min_val = val;
      if (val > dst[c].max_val) dst[c].max_val = val;
    }
  }

  KOKKOS_INLINE_FUNCTION void join(minmax_type dst, const minmax_type src) const {
    for (ordinal_t c = 0; c < value_count; c++) {
      if (src[c].min_val < dst[c].min_val) dst[c].min_val = src[c].min_val;
      if (src[c].max_val > dst[c].max_val) dst[c].max_val = src[c].max_val;
    }
  }
};

template <typename perm_view_type>
struct UpdatePermFunctor {
  using ordinal_t = typename perm_view_type::value_type;
  perm_view_type reverse_perm;
  perm_view_type perm;

  UpdatePermFunctor(const perm_view_type &reverse_perm_, perm_view_type &perm_)
      : reverse_perm(reverse_perm_), perm(perm_) {}
  KOKKOS_INLINE_FUNCTION void operator()(ordinal_t i) const {
    ordinal_t orig_idx = reverse_perm(i);
    perm(orig_idx)     = i;
  }
};

template <typename coors_view_type, typename dim_view_type, typename minmax_view_type, typename perm_view_type>
struct ScanLeftRightPartitionsFunctor {
  using value_t   = typename coors_view_type::value_type;
  using ordinal_t = typename perm_view_type::value_type;

  coors_view_type coordinates;
  dim_view_type partitioned_dim;
  minmax_view_type minmax;
  perm_view_type dest_indices;

  ScanLeftRightPartitionsFunctor(const coors_view_type &coordinates_, const dim_view_type &partitioned_dim_,
                                 const minmax_view_type &minmax_, perm_view_type &dest_indices_)
      : coordinates(coordinates_), partitioned_dim(partitioned_dim_), minmax(minmax_), dest_indices(dest_indices_) {}

  // Kokkos scan operator: 'update' keeps track of the exclusive prefix sum,
  // and 'final_pass' is true during the second pass when writing to global memory.
  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_t i, ordinal_t &update, const bool final_pass) const {
    value_t mid_value = (minmax(1) + minmax(0)) / 2.0;
    // Condition: True if the point belongs to the LEFT chunk
    bool is_left = (coordinates(i, partitioned_dim()) > mid_value);
    if (final_pass) {
      if (is_left) {
        // 'update' holds the number of left points seen before index i
        dest_indices(i) = update;
      } else {
        // Total left points = 'update' at the final element + its contribution.
        // However, during final_pass, the total left points is not globally known inside
        // the stream until the end. We handle the Right index layout in the next step.
        dest_indices(i) = i - update;
      }
    }
    // Increment the running prefix sum if this point goes left
    if (is_left) {
      update += 1;
    }
  }
};

template <typename coors_view_type, typename perm_view_type, typename dim_view_type, typename minmax_view_type>
struct ShuffleCoordsAndGenerateReversePermFunctor {
  using value_t   = typename coors_view_type::value_type;
  using ordinal_t = typename perm_view_type::value_type;

  coors_view_type src_coordinates;
  perm_view_type src_reverse_perm;
  perm_view_type dest_indices;
  dim_view_type partitioned_dim;
  minmax_view_type minmax;
  coors_view_type dst_coordinates;
  perm_view_type dst_reverse_perm;
  ordinal_t left_size;
  ordinal_t ndim;

  ShuffleCoordsAndGenerateReversePermFunctor(const coors_view_type &src_coordinates_,
                                             const perm_view_type &src_reverse_perm_,
                                             const perm_view_type &dest_indices_, const dim_view_type &partitioned_dim_,
                                             const minmax_view_type &minmax_, coors_view_type &dst_coordinates_,
                                             perm_view_type &dst_reverse_perm_, const ordinal_t &left_size_)
      : src_coordinates(src_coordinates_),
        src_reverse_perm(src_reverse_perm_),
        dest_indices(dest_indices_),
        partitioned_dim(partitioned_dim_),
        minmax(minmax_),
        dst_coordinates(dst_coordinates_),
        dst_reverse_perm(dst_reverse_perm_),
        left_size(left_size_) {
    ndim = static_cast<ordinal_t>(src_coordinates.extent(1));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_t i) const {
    int dest_idx;
    value_t mid_value = (minmax(1) + minmax(0)) / 2.0;
    if (src_coordinates(i, partitioned_dim()) > mid_value) {
      dest_idx = dest_indices(i);
    } else {
      dest_idx = left_size + dest_indices(i);
    }
    for (ordinal_t j = 0; j < ndim; j++) {
      dst_coordinates(dest_idx, j) = src_coordinates(i, j);
    }
    dst_reverse_perm(dest_idx) = src_reverse_perm(i);
  }
};

template <typename view_type, typename minmax_view_type>
void find_min_max(const view_type &A, minmax_view_type &minmax) {
  using execution_space = typename view_type::device_type::execution_space;
  size_t n_elements     = static_cast<size_t>(A.extent(0));
  size_t n_dims         = static_cast<size_t>(A.extent(1));

  Kokkos::parallel_reduce(Kokkos::RangePolicy<execution_space>(0, n_elements),
                          ColumnMinMaxFunctor<view_type, size_t>(A, n_dims), minmax);
}

/**
 * @brief Bisect and assign partition indices to coordinate list
 */
template <typename coors_view_type, typename dim_view_type, typename minmax_view_type, typename index_view_type,
          typename ordinal_type>
inline void bisect(const coors_view_type &coors, const dim_view_type &partitioned_dim,
                   const minmax_view_type &minmax_bisect, index_view_type &partition_dests, ordinal_type &p1_size,
                   ordinal_type &p2_size, const int &max_bisection_steps) {
  using execution_space = typename coors_view_type::device_type::execution_space;
  using value_type      = typename coors_view_type::value_type;

  const ordinal_type N = static_cast<ordinal_type>(coors.extent(0));

  Kokkos::RangePolicy<execution_space, ordinal_type> policy(0, N);
  ScanLeftRightPartitionsFunctor scan_functor(coors, partitioned_dim, minmax_bisect, partition_dests);

  for (int bisection_step = 0; bisection_step < max_bisection_steps; ++bisection_step) {
    Kokkos::parallel_scan("RCB_PrefixScan", policy, scan_functor, p1_size);

    p2_size = N - p1_size;

    value_type p1_weight = static_cast<value_type>(p1_size);
    value_type p2_weight = static_cast<value_type>(p2_size);

    value_type weight_ratio = std::max(p1_weight, p2_weight) / std::min(p1_weight, p2_weight);

    if ((weight_ratio < 1.1) || (bisection_step == (max_bisection_steps - 1)))
      break;
    else {
      // Update min_val or max_val to calculate a new mid_point
      // Idea: shift mid_point to the heavier partition
      Kokkos::parallel_for(
          Kokkos::RangePolicy<execution_space>(0, 1), KOKKOS_LAMBDA(const ordinal_type /* unused */) {
            if (p1_weight > p2_weight)
              minmax_bisect(0) = (minmax_bisect(1) + minmax_bisect(0)) / 2.0;
            else
              minmax_bisect(1) = (minmax_bisect(1) + minmax_bisect(0)) / 2.0;
          });
    }
  }  // End for
}

/**
 * @brief Recursive coordinate bisection on the coordinate list
 */
template <typename coors_view_type, typename perm_view_type>
std::vector<typename perm_view_type::value_type> rcb(coors_view_type &coordinates, perm_view_type &perm,
                                                     perm_view_type &reverse_perm, const int &n_levels,
                                                     const int &max_bisection_steps = 11) {
  using execution_space = typename coors_view_type::device_type::execution_space;
  using memory_space    = typename coors_view_type::device_type::memory_space;
  using scalar_t        = typename coors_view_type::value_type;
  using ordinal_type    = typename perm_view_type::value_type;

  const ordinal_type N    = static_cast<ordinal_type>(coordinates.extent(0));
  const ordinal_type ndim = static_cast<ordinal_type>(coordinates.extent(1));

  // Allocate working buffers
  std::vector<coors_view_type> coordinates_buf(2);
  std::vector<perm_view_type> reverse_perm_buf(2);
  coordinates_buf[0]  = coordinates;
  coordinates_buf[1]  = coors_view_type(Kokkos::view_alloc(Kokkos::WithoutInitializing, "coordinates_buf"), N, ndim);
  reverse_perm_buf[0] = reverse_perm;
  reverse_perm_buf[1] = perm_view_type(Kokkos::view_alloc(Kokkos::WithoutInitializing, "reverse_perm_buf"), N);
  Kokkos::View<MinMaxPair<scalar_t> *, Kokkos::LayoutLeft, memory_space> minmax(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "minmax"), ndim);
  Kokkos::View<ordinal_type, memory_space> partitioned_dim(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "partitioned_dim"));
  Kokkos::View<scalar_t *, Kokkos::LayoutLeft, memory_space> minmax_bisect(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "minmax_bisect"), 2);

  // Initialize reverse_perm
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space, ordinal_type>(0, N),
                       FillOneIncrementFunctor<perm_view_type>(reverse_perm_buf[0]));

  // Allocating a temporary view for tracking local partitioning offsets
  perm_view_type partition_dests(Kokkos::view_alloc(Kokkos::WithoutInitializing, "partition_dests"), N);

  ordinal_type n_partitions =
      1;  // number of partitions at the previous level (initial value is 1, i.e., starting with the entire mesh points)
  std::vector<std::vector<ordinal_type>> partition_sizes(
      n_levels);                    // contain the numbers of basis functions (or elements) per partition in levels
  partition_sizes[0].push_back(N);  // starting with the entire mesh points

  ordinal_type src = 0;
  ordinal_type dst = 1;
  // Start RCB
  for (int lvl = 1; lvl < n_levels; lvl++) {  // skip level 0 and start from level 1
    ordinal_type coordinates_offset = 0;      // always start from beginning of the mesh points

    for (ordinal_type p = 0; p < n_partitions; p++) {  // go through each partition of previous level and do bisecting
      if (p > 0) {
        // Calculate coordinates offset of the current partition based on the previous partition
        coordinates_offset += partition_sizes[lvl - 1][p - 1];
      }

      ordinal_type N0          = partition_sizes[lvl - 1][p];  // partition size (or length)
      ordinal_type p1_size     = 0;
      ordinal_type p2_size     = 0;
      auto sub_coordinates_src = Kokkos::subview(
          coordinates_buf[src], Kokkos::make_pair(coordinates_offset, coordinates_offset + N0), Kokkos::ALL());
      auto sub_reverse_perm_src =
          Kokkos::subview(reverse_perm_buf[src], Kokkos::make_pair(coordinates_offset, coordinates_offset + N0));
      auto sub_coordinates_dst = Kokkos::subview(
          coordinates_buf[dst], Kokkos::make_pair(coordinates_offset, coordinates_offset + N0), Kokkos::ALL());
      auto sub_reverse_perm_dst =
          Kokkos::subview(reverse_perm_buf[dst], Kokkos::make_pair(coordinates_offset, coordinates_offset + N0));
      auto sub_partition_dests =
          Kokkos::subview(partition_dests, Kokkos::make_pair(coordinates_offset, coordinates_offset + N0));

      // Find min, max, span of each dimension and select the most elongated dimension for partitioning
      find_min_max(sub_coordinates_src, minmax);

      Kokkos::parallel_for(
          Kokkos::RangePolicy<execution_space>(0, 1), KOKKOS_LAMBDA(const ordinal_type /* unused */) {
            scalar_t max_span       = minmax(0).max_val - minmax(0).min_val;
            ordinal_type dim_select = 0;
            for (ordinal_type j = 1; j < ndim; j++) {
              scalar_t span = minmax(j).max_val - minmax(j).min_val;
              if (max_span < span) {
                max_span   = span;
                dim_select = j;
              }
            }
            partitioned_dim() = dim_select;
            minmax_bisect(0)  = minmax(dim_select).min_val;
            minmax_bisect(1)  = minmax(dim_select).max_val;
          });

      // Bisect partition on the most elongated dimension
      bisect(sub_coordinates_src, partitioned_dim, minmax_bisect, sub_partition_dests, p1_size, p2_size,
             max_bisection_steps);

      // Shuffle coordinates using bisecting results and generate reverse_perm
      ShuffleCoordsAndGenerateReversePermFunctor<decltype(sub_coordinates_src), decltype(sub_reverse_perm_src),
                                                 decltype(partitioned_dim), decltype(minmax_bisect)>
          func(sub_coordinates_src, sub_reverse_perm_src, sub_partition_dests, partitioned_dim, minmax_bisect,
               sub_coordinates_dst, sub_reverse_perm_dst, p1_size);
      Kokkos::RangePolicy<execution_space, ordinal_type> policy(0, N0);
      Kokkos::parallel_for(policy, func);

      partition_sizes[lvl].push_back(p1_size);
      partition_sizes[lvl].push_back(p2_size);
    }  // end Partition loop

    // Update the number of partitions of this level (used for bisections in the next level)
    n_partitions = n_partitions * 2;

    // Swap working buffers
    std::swap(src, dst);
  }  // end Level loop

  // Perform deep_copy to obtain partitioned coordinates and reverse_perm (if needed)
  if (n_levels % 2 == 0) {
    Kokkos::deep_copy(coordinates_buf[0], coordinates_buf[1]);
    Kokkos::deep_copy(reverse_perm_buf[0], reverse_perm_buf[1]);
  }

  // Create perm from reverse_perm
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space, ordinal_type>(0, N),
                       UpdatePermFunctor<perm_view_type>(reverse_perm, perm));

  return partition_sizes[n_levels - 1];
}

}  // namespace Impl
}  // namespace KokkosGraph
#endif
