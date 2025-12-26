// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include "KokkosGraph_RCB.hpp"

#include <cmath>
#include <vector>

// Generate 1-D coordinates on a line
template <typename coors_view_t>
coors_view_t generate_coordinates_1d() {
  using scalar_t = typename coors_view_t::value_type;

  scalar_t x_max = 1.1;
  scalar_t x_min = 0;

  int n_pts_per_dim = 121;
  int n_spaces      = n_pts_per_dim - 1;

  int n_pts = n_pts_per_dim;

  scalar_t dx = (x_max - x_min) / n_spaces;

  coors_view_t coordinates(Kokkos::view_alloc(Kokkos::WithoutInitializing, "coordinates"), n_pts, 1);
  auto h_coordinates = Kokkos::create_mirror_view(coordinates);

  int cnt = 0;
  for (int kk = 0; kk < n_pts_per_dim; kk++) {  // x
    scalar_t x            = x_min + kk * dx;
    h_coordinates(cnt, 0) = x;
    cnt++;
  }

  Kokkos::deep_copy(coordinates, h_coordinates);

  return coordinates;
}

// Generate 2-D coordinates in a rectangular
template <typename coors_view_t>
coors_view_t generate_coordinates_2d() {
  using scalar_t = typename coors_view_t::value_type;

  scalar_t x_max = 1.1;
  scalar_t x_min = 0;
  scalar_t y_max = 1.05;
  scalar_t y_min = 0;

  int n_pts_per_dim = 11;
  int n_spaces      = n_pts_per_dim - 1;

  int n_pts = n_pts_per_dim * n_pts_per_dim;

  scalar_t dx = (x_max - x_min) / n_spaces;
  scalar_t dy = (y_max - y_min) / n_spaces;

  coors_view_t coordinates(Kokkos::view_alloc(Kokkos::WithoutInitializing, "coordinates"), n_pts, 2);
  auto h_coordinates = Kokkos::create_mirror_view(coordinates);

  int cnt = 0;
  for (int jj = 0; jj < n_pts_per_dim; jj++) {  // y
    scalar_t y = y_min + jj * dy;
    for (int kk = 0; kk < n_pts_per_dim; kk++) {  // x
      scalar_t x            = x_min + kk * dx;
      h_coordinates(cnt, 0) = x;
      h_coordinates(cnt, 1) = y;
      cnt++;
    }
  }

  Kokkos::deep_copy(coordinates, h_coordinates);

  return coordinates;
}

// Generate 3-D coordinates in a cuboid
template <typename coors_view_t>
coors_view_t generate_coordinates_3d() {
  using scalar_t = typename coors_view_t::value_type;

  scalar_t x_max = 1.1;
  scalar_t x_min = 0;
  scalar_t y_max = 1.05;
  scalar_t y_min = 0;
  scalar_t z_max = 1.0;
  scalar_t z_min = 0;

  int n_pts_per_dim = 5;
  int n_spaces      = n_pts_per_dim - 1;

  int n_pts = n_pts_per_dim * n_pts_per_dim * n_pts_per_dim;

  scalar_t dx = (x_max - x_min) / n_spaces;
  scalar_t dy = (y_max - y_min) / n_spaces;
  scalar_t dz = (z_max - z_min) / n_spaces;

  coors_view_t coordinates(Kokkos::view_alloc(Kokkos::WithoutInitializing, "coordinates"), n_pts, 3);
  auto h_coordinates = Kokkos::create_mirror_view(coordinates);

  int cnt = 0;
  for (int ii = 0; ii < n_pts_per_dim; ii++) {  // z
    scalar_t z = z_min + ii * dz;
    for (int jj = 0; jj < n_pts_per_dim; jj++) {  // y
      scalar_t y = y_min + jj * dy;
      for (int kk = 0; kk < n_pts_per_dim; kk++) {  // x
        scalar_t x            = x_min + kk * dx;
        h_coordinates(cnt, 0) = x;
        h_coordinates(cnt, 1) = y;
        h_coordinates(cnt, 2) = z;
        cnt++;
      }
    }
  }

  Kokkos::deep_copy(coordinates, h_coordinates);

  return coordinates;
}

template <typename scalar_t, typename lno_t, typename device>
void test_rcb(lno_t ndim, lno_t np) {
  using coors_view_t = Kokkos::View<scalar_t**, Kokkos::LayoutLeft>;
  using perm_view_t  = Kokkos::View<lno_t*, Kokkos::LayoutLeft>;

  coors_view_t coordinates;
  if (ndim == 1)
    coordinates = generate_coordinates_1d<coors_view_t>();
  else if (ndim == 2)
    coordinates = generate_coordinates_2d<coors_view_t>();
  else
    coordinates = generate_coordinates_3d<coors_view_t>();

  lno_t n_coordinates = static_cast<lno_t>(coordinates.extent(0));
  perm_view_t perm_rcb("perm_rcb", n_coordinates);
  perm_view_t reverse_perm_rcb("reverse_perm_rcb", n_coordinates);
  auto h_coordinates = Kokkos::create_mirror(coordinates);
  Kokkos::deep_copy(h_coordinates, coordinates);

  lno_t n_levels = static_cast<lno_t>(std::log2(static_cast<double>(np)) + 1);

  // Run RCB
  std::vector<lno_t> partition_sizes =
      KokkosGraph::Experimental::recursive_coordinate_bisection(coordinates, perm_rcb, reverse_perm_rcb, n_levels);

  // Check partition sizes
  lno_t sum_partition_sizes = 0;
  for (lno_t i = 0; i < static_cast<lno_t>(partition_sizes.size()); i++) {
    sum_partition_sizes += partition_sizes[i];
  }
  ASSERT_EQ(sum_partition_sizes, n_coordinates);

  // Check permutations
  auto h_perm_rcb         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), perm_rcb);
  auto h_reverse_perm_rcb = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), reverse_perm_rcb);
  auto h_coordinates_perm = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), coordinates);

  // Output permutation and reverse permutation should be different
  bool dif_flag = false;
  for (lno_t i = 0; i < n_coordinates; i++) {
    if (h_perm_rcb(i) != h_reverse_perm_rcb(i)) dif_flag = true;
  }
  ASSERT_TRUE(dif_flag);

  for (lno_t i = 0; i < n_coordinates; i++) {
    for (lno_t j = 0; j < ndim; j++) {
      ASSERT_EQ(h_coordinates(i, j), h_coordinates_perm(h_perm_rcb(i), j));
    }
  }

  for (lno_t i = 0; i < n_coordinates; i++) {
    for (lno_t j = 0; j < ndim; j++) {
      ASSERT_EQ(h_coordinates_perm(i, j), h_coordinates(h_reverse_perm_rcb(i), j));
    }
  }

  // Run RCB on RCB-reordered coordinates. Output permutation and reverse permutation are supposed to be identical
  KokkosGraph::Experimental::recursive_coordinate_bisection(coordinates, perm_rcb, reverse_perm_rcb, n_levels);

  Kokkos::deep_copy(h_perm_rcb, perm_rcb);
  Kokkos::deep_copy(h_reverse_perm_rcb, reverse_perm_rcb);

  for (lno_t i = 0; i < n_coordinates; i++) {
    ASSERT_EQ(h_perm_rcb(i), h_reverse_perm_rcb(i));
  }
}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                                \
  TEST_F(TestCategory, graph##_##rcb##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_rcb<SCALAR, ORDINAL, DEVICE>(1, 2);                                         \
    test_rcb<SCALAR, ORDINAL, DEVICE>(1, 4);                                         \
    test_rcb<SCALAR, ORDINAL, DEVICE>(1, 8);                                         \
    test_rcb<SCALAR, ORDINAL, DEVICE>(1, 16);                                        \
    test_rcb<SCALAR, ORDINAL, DEVICE>(2, 2);                                         \
    test_rcb<SCALAR, ORDINAL, DEVICE>(2, 4);                                         \
    test_rcb<SCALAR, ORDINAL, DEVICE>(2, 8);                                         \
    test_rcb<SCALAR, ORDINAL, DEVICE>(2, 16);                                        \
    test_rcb<SCALAR, ORDINAL, DEVICE>(3, 2);                                         \
    test_rcb<SCALAR, ORDINAL, DEVICE>(3, 4);                                         \
    test_rcb<SCALAR, ORDINAL, DEVICE>(3, 8);                                         \
    test_rcb<SCALAR, ORDINAL, DEVICE>(3, 16);                                        \
  }

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_OFFSET_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT) && defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int, size_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T) && defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, int64_t, size_t, TestDevice)
#endif

#undef EXECUTE_TEST
