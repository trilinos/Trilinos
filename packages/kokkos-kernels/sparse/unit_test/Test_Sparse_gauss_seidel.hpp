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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosSparse_IOUtils.hpp"
// #include <Kokkos_Sparse_CrsMatrix.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <cstdlib>
#include <iostream>
#include <complex>
#include <map>
#include <random>
#include <vector>
#include "KokkosSparse_gauss_seidel.hpp"
#include "KokkosSparse_partitioning_impl.hpp"
#include "KokkosSparse_sor_sequential_impl.hpp"
#include "KokkosSparse_SortCrs.hpp"
#include "KokkosKernels_TestUtils.hpp"
#include "Test_Sparse_Utils.hpp"

// #ifndef kokkos_complex_double
// #define kokkos_complex_double Kokkos::complex<double>
// #define kokkos_complex_float Kokkos::complex<float>
// #endif

typedef Kokkos::complex<double> kokkos_complex_double;
typedef Kokkos::complex<float> kokkos_complex_float;

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;
using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;

namespace Test {

// Run GS on the given vectors, where the handle is already set up.
template <typename Handle, typename crsMat_t, typename vec_t>
void run_gauss_seidel(Handle &kh, crsMat_t input_mat, vec_t x_vector, vec_t y_vector, bool is_symmetric_graph,
                      typename crsMat_t::value_type omega,
                      int apply_type = 0  // 0 for symmetric, 1 for forward, 2 for backward.
) {
  const size_t num_rows = input_mat.numRows();
  const size_t num_cols = input_mat.numCols();
  const int apply_count = 2;

  gauss_seidel_symbolic(&kh, num_rows, num_cols, input_mat.graph.row_map, input_mat.graph.entries, is_symmetric_graph);
  gauss_seidel_numeric(&kh, num_rows, num_cols, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values,
                       is_symmetric_graph);

  switch (apply_type) {
    case 0:
      symmetric_gauss_seidel_apply(&kh, num_rows, num_cols, input_mat.graph.row_map, input_mat.graph.entries,
                                   input_mat.values, x_vector, y_vector, false, true, omega, apply_count);
      break;
    case 1:
      forward_sweep_gauss_seidel_apply(&kh, num_rows, num_cols, input_mat.graph.row_map, input_mat.graph.entries,
                                       input_mat.values, x_vector, y_vector, false, true, omega, apply_count);
      break;
    case 2:
      backward_sweep_gauss_seidel_apply(&kh, num_rows, num_cols, input_mat.graph.row_map, input_mat.graph.entries,
                                        input_mat.values, x_vector, y_vector, false, true, omega, apply_count);
      break;
    default:
      symmetric_gauss_seidel_apply(&kh, num_rows, num_cols, input_mat.graph.row_map, input_mat.graph.entries,
                                   input_mat.values, x_vector, y_vector, false, true, omega, apply_count);
      break;
  }
}

template <typename crsMat_t, typename vec_t>
void run_gauss_seidel(crsMat_t input_mat, GSAlgorithm gs_algorithm, vec_t x_vector, vec_t y_vector,
                      bool is_symmetric_graph,
                      int apply_type   = 0,  // 0 for symmetric, 1 for forward, 2 for backward.
                      int cluster_size = 1,
                      bool classic     = false,  // only with two-stage, true for sptrsv instead of richardson
                      ClusteringAlgorithm clusterAlgo             = CLUSTER_DEFAULT,
                      KokkosGraph::ColoringAlgorithm coloringAlgo = KokkosGraph::COLORING_DEFAULT) {
  using size_type = typename crsMat_t::size_type;
  using lno_t     = typename crsMat_t::ordinal_type;
  using scalar_t  = typename crsMat_t::value_type;
  using device    = typename crsMat_t::device_type;

  typedef KokkosKernelsHandle<size_type, lno_t, scalar_t, typename device::execution_space,
                              typename device::memory_space, typename device::memory_space>
      KernelHandle;

  scalar_t omega(0.9);

  KernelHandle kh;
  if (gs_algorithm == GS_CLUSTER)
    kh.create_gs_handle(clusterAlgo, cluster_size, coloringAlgo);
  else if (gs_algorithm == GS_TWOSTAGE) {
    // test for two-stage/classical gs
    kh.create_gs_handle(gs_algorithm);
    kh.set_gs_twostage(!classic, input_mat.numRows());
    if (classic) {
      // two-stage with SpTRSV supports only omega = one
      omega = Kokkos::ArithTraits<scalar_t>::one();
    }
  } else {
    kh.create_gs_handle(GS_DEFAULT, coloringAlgo);
  }

  run_gauss_seidel(kh, input_mat, x_vector, y_vector, is_symmetric_graph, omega, apply_type);

  kh.destroy_gs_handle();
}

template <typename ExecSpace, typename Handle, typename crsMat_t, typename vec_t>
void run_gauss_seidel_streams(std::vector<ExecSpace> &instances, std::vector<Handle> &kh,
                              std::vector<crsMat_t> &input_mat, std::vector<vec_t> &x_vector,
                              std::vector<vec_t> &y_vector, bool is_symmetric_graph,
                              typename crsMat_t::value_type omega,
                              int apply_type,  // 0 for symmetric, 1 for forward, 2 for backward.
                              int nstreams = 1) {
  for (int i = 0; i < nstreams; i++) {
    gauss_seidel_symbolic(instances[i], &kh[i], input_mat[i].numRows(), input_mat[i].numCols(),
                          input_mat[i].graph.row_map, input_mat[i].graph.entries, is_symmetric_graph);
    gauss_seidel_numeric(instances[i], &kh[i], input_mat[i].numRows(), input_mat[i].numCols(),
                         input_mat[i].graph.row_map, input_mat[i].graph.entries, input_mat[i].values,
                         is_symmetric_graph);
  }

  const int apply_count = 2;
  for (int i = 0; i < nstreams; i++) {
    switch (apply_type) {
      case 0:
        symmetric_gauss_seidel_apply(instances[i], &kh[i], input_mat[i].numRows(), input_mat[i].numCols(),
                                     input_mat[i].graph.row_map, input_mat[i].graph.entries, input_mat[i].values,
                                     x_vector[i], y_vector[i], false, true, omega, apply_count);
        break;
      case 1:
        forward_sweep_gauss_seidel_apply(instances[i], &kh[i], input_mat[i].numRows(), input_mat[i].numCols(),
                                         input_mat[i].graph.row_map, input_mat[i].graph.entries, input_mat[i].values,
                                         x_vector[i], y_vector[i], false, true, omega, apply_count);
        break;
      case 2:
        backward_sweep_gauss_seidel_apply(instances[i], &kh[i], input_mat[i].numRows(), input_mat[i].numCols(),
                                          input_mat[i].graph.row_map, input_mat[i].graph.entries, input_mat[i].values,
                                          x_vector[i], y_vector[i], false, true, omega, apply_count);
        break;
      default:
        symmetric_gauss_seidel_apply(instances[i], &kh[i], input_mat[i].numRows(), input_mat[i].numCols(),
                                     input_mat[i].graph.row_map, input_mat[i].graph.entries, input_mat[i].values,
                                     x_vector[i], y_vector[i], false, true, omega, apply_count);
        break;
    }
  }
}
}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_gauss_seidel_rank1(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance, bool symmetric) {
  using namespace Test;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename Kokkos::ArithTraits<scalar_t>::mag_type mag_t;
  srand(245);
  lno_t numCols      = numRows;
  crsMat_t input_mat = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<crsMat_t>(
      numRows, numCols, nnz, row_size_variance, bandwidth);
  if (symmetric) {
    // Symmetrize on host, rather than relying on the parallel versions (those
    // can be tested for symmetric=false)
    input_mat = Test::symmetrize<scalar_t, lno_t, size_type, device, crsMat_t>(input_mat);
  }
  lno_t nv = input_mat.numRows();
  scalar_view_t solution_x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "X (correct)"), nv);
  create_random_x_vector(solution_x);
  mag_t initial_norm_res = KokkosBlas::nrm2(solution_x);
  scalar_view_t y_vector = create_random_y_vector(input_mat, solution_x);
  // GS_DEFAULT is GS_TEAM on CUDA and GS_PERMUTED on other spaces, and the
  // behavior of each algorithm _should be_ the same on every execution space,
  // which is why we just test GS_DEFAULT.
  int apply_count = 3;  // test symmetric, forward, backward
  scalar_view_t x_vector(Kokkos::view_alloc(Kokkos::WithoutInitializing, "x vector"), nv);
  const scalar_t one  = Kokkos::ArithTraits<scalar_t>::one();
  const scalar_t zero = Kokkos::ArithTraits<scalar_t>::zero();
  //*** Point-coloring version ****
  for (int apply_type = 0; apply_type < apply_count; ++apply_type) {
    Kokkos::Timer timer1;
    Kokkos::deep_copy(x_vector, zero);
    run_gauss_seidel(input_mat, GS_DEFAULT, x_vector, y_vector, symmetric, apply_type);
    // double gs = timer1.seconds();
    // KokkosKernels::Impl::print_1Dview(x_vector);
    KokkosBlas::axpby(one, solution_x, -one, x_vector);
    mag_t result_norm_res = KokkosBlas::nrm2(x_vector);
    EXPECT_LT(result_norm_res, initial_norm_res);
  }
  //*** Cluster-coloring version ****
  int clusterSizes[3]                              = {2, 5, 34};
  std::vector<ClusteringAlgorithm> clusteringAlgos = {CLUSTER_MIS2, CLUSTER_BALLOON};
  for (int csize = 0; csize < 3; csize++) {
    for (auto clusterAlgo : clusteringAlgos) {
      for (int apply_type = 0; apply_type < apply_count; ++apply_type) {
        Kokkos::Timer timer1;
        // Zero out X before solving
        Kokkos::deep_copy(x_vector, zero);
        run_gauss_seidel(input_mat, GS_CLUSTER, x_vector, y_vector, symmetric, apply_type, clusterSizes[csize], false,
                         clusterAlgo);
        KokkosBlas::axpby(one, solution_x, -one, x_vector);
        mag_t result_norm_res = KokkosBlas::nrm2(x_vector);
        EXPECT_LT(result_norm_res, initial_norm_res);
      }
    }
  }
  //*** Two-stage version ****
  for (int apply_type = 0; apply_type < apply_count; ++apply_type) {
    Kokkos::deep_copy(x_vector, zero);
    run_gauss_seidel(input_mat, GS_TWOSTAGE, x_vector, y_vector, symmetric, apply_type);
    KokkosBlas::axpby(one, solution_x, -one, x_vector);
    mag_t result_norm_res = KokkosBlas::nrm2(x_vector);
    EXPECT_LT(result_norm_res, initial_norm_res);
  }
  //*** Two-stage version (classic) ****
  for (int apply_type = 0; apply_type < apply_count; ++apply_type) {
    Kokkos::deep_copy(x_vector, zero);
    run_gauss_seidel(input_mat, GS_TWOSTAGE, x_vector, y_vector, symmetric, apply_type, 0, true);
    KokkosBlas::axpby(one, solution_x, -one, x_vector);
    mag_t result_norm_res = KokkosBlas::nrm2(x_vector);
    EXPECT_LT(result_norm_res, initial_norm_res);
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_gauss_seidel_rank2(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance, lno_t numVecs,
                             bool symmetric) {
  using namespace Test;
  srand(245);
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef Kokkos::View<scalar_t **, default_layout, device> scalar_view2d_t;
  typedef Kokkos::View<scalar_t **, default_layout, Kokkos::HostSpace> host_scalar_view2d_t;
  typedef typename Kokkos::ArithTraits<scalar_t>::mag_type mag_t;

  lno_t numCols      = numRows;
  crsMat_t input_mat = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<crsMat_t>(
      numRows, numCols, nnz, row_size_variance, bandwidth);
  if (symmetric) {
    // Symmetrize on host, rather than relying on the parallel versions (those
    // can be tested for symmetric=false)
    input_mat = Test::symmetrize<scalar_t, lno_t, size_type, device, crsMat_t>(input_mat);
  }
  lno_t nv = input_mat.numRows();
  host_scalar_view2d_t solution_x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "X (correct)"), nv, numVecs);
  create_random_x_vector(solution_x);
  scalar_view2d_t x_vector(Kokkos::view_alloc(Kokkos::WithoutInitializing, "X"), nv, numVecs);
  Kokkos::deep_copy(x_vector, solution_x);
  scalar_view2d_t y_vector = create_random_y_vector_mv(input_mat, x_vector);
  auto x_host              = Kokkos::create_mirror_view(x_vector);
  std::vector<mag_t> initial_norms(numVecs);
  for (lno_t i = 0; i < numVecs; i++) {
    scalar_t sum = 0;
    for (lno_t j = 0; j < nv; j++) {
      sum += solution_x(j, i) * solution_x(j, i);
    }
    initial_norms[i] = Kokkos::ArithTraits<mag_t>::sqrt(Kokkos::ArithTraits<scalar_t>::abs(sum));
  }
  int apply_count     = 3;  // test symmetric, forward, backward
  const scalar_t zero = Kokkos::ArithTraits<scalar_t>::zero();
  //*** Point-coloring version ****
  for (int apply_type = 0; apply_type < apply_count; ++apply_type) {
    Kokkos::Timer timer1;
    // Zero out X before solving
    Kokkos::deep_copy(x_vector, zero);
    run_gauss_seidel(input_mat, GS_DEFAULT, x_vector, y_vector, symmetric, apply_type);
    Kokkos::deep_copy(x_host, x_vector);
    for (lno_t i = 0; i < numVecs; i++) {
      scalar_t diffDot = 0;
      for (lno_t j = 0; j < numRows; j++) {
        scalar_t diff = x_host(j, i) - solution_x(j, i);
        diffDot += diff * diff;
      }
      mag_t res = Kokkos::ArithTraits<mag_t>::sqrt(Kokkos::ArithTraits<scalar_t>::abs(diffDot));
      EXPECT_LT(res, initial_norms[i]);
    }
  }
  //*** Cluster-coloring version ****
  int clusterSizes[3] = {2, 5, 34};
  for (int csize = 0; csize < 3; csize++) {
    for (int algo = 0; algo < (int)NUM_CLUSTERING_ALGORITHMS; algo++) {
      for (int apply_type = 0; apply_type < apply_count; ++apply_type) {
        Kokkos::Timer timer1;
        // Zero out X before solving
        Kokkos::deep_copy(x_vector, zero);
        run_gauss_seidel(input_mat, GS_CLUSTER, x_vector, y_vector, symmetric, apply_type, clusterSizes[csize], false,
                         (ClusteringAlgorithm)algo);
        Kokkos::deep_copy(x_host, x_vector);
        for (lno_t i = 0; i < numVecs; i++) {
          scalar_t diffDot = 0;
          for (lno_t j = 0; j < numRows; j++) {
            scalar_t diff = x_host(j, i) - solution_x(j, i);
            diffDot += diff * diff;
          }
          mag_t res = Kokkos::ArithTraits<mag_t>::sqrt(Kokkos::ArithTraits<scalar_t>::abs(diffDot));
          EXPECT_LT(res, initial_norms[i]);
        }
      }
    }
  }
  //*** Two-stage version ****
  for (int apply_type = 0; apply_type < apply_count; ++apply_type) {
    // Zero out X before solving
    Kokkos::deep_copy(x_vector, zero);
    run_gauss_seidel(input_mat, GS_TWOSTAGE, x_vector, y_vector, symmetric, apply_type);
    Kokkos::deep_copy(x_host, x_vector);
    for (lno_t i = 0; i < numVecs; i++) {
      scalar_t diffDot = 0;
      for (lno_t j = 0; j < numRows; j++) {
        scalar_t diff = x_host(j, i) - solution_x(j, i);
        diffDot += diff * diff;
      }
      mag_t res = Kokkos::ArithTraits<mag_t>::sqrt(Kokkos::ArithTraits<scalar_t>::abs(diffDot));
      EXPECT_LT(res, initial_norms[i]);
    }
  }
  //*** Two-stage version (classic) ****
  for (int apply_type = 0; apply_type < apply_count; ++apply_type) {
    // Zero out X before solving
    Kokkos::deep_copy(x_vector, zero);
    run_gauss_seidel(input_mat, GS_TWOSTAGE, x_vector, y_vector, symmetric, apply_type, 0, true);
    Kokkos::deep_copy(x_host, x_vector);
    for (lno_t i = 0; i < numVecs; i++) {
      scalar_t diffDot = 0;
      for (lno_t j = 0; j < numRows; j++) {
        scalar_t diff = x_host(j, i) - solution_x(j, i);
        diffDot += diff * diff;
      }
      mag_t res = Kokkos::ArithTraits<mag_t>::sqrt(Kokkos::ArithTraits<scalar_t>::abs(diffDot));
      EXPECT_LT(res, initial_norms[i]);
    }
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_sequential_sor(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance) {
  const scalar_t zero = Kokkos::ArithTraits<scalar_t>::zero();
  const scalar_t one  = Kokkos::ArithTraits<scalar_t>::one();
  srand(245);
  typedef typename device::execution_space exec_space;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  lno_t numCols      = numRows;
  crsMat_t input_mat = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<crsMat_t>(
      numRows, numCols, nnz, row_size_variance, bandwidth);
  auto rowmap  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), input_mat.graph.row_map);
  auto entries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), input_mat.graph.entries);
  auto values  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), input_mat.values);
  // create raw x (unkown), y (rhs) vectors
  using vector_t = typename crsMat_t::values_type::non_const_type;
  // Create random x
  vector_t x("X", numRows);
  auto x_host = Kokkos::create_mirror_view(x);
  for (lno_t i = 0; i < numRows; i++) {
    x_host(i) = one * scalar_t(10.0 * rand() / RAND_MAX);
  }
  Kokkos::deep_copy(x, x_host);
  // record the correct solution, to compare against at the end
  vector_t xgold("X gold", numRows);
  Kokkos::deep_copy(xgold, x);
  vector_t y = Test::create_random_y_vector(input_mat, x);
  exec_space().fence();
  auto y_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y);
  // initial solution is zero
  Kokkos::deep_copy(x_host, zero);
  // get the inverse diagonal (only needed on host)
  Kokkos::View<scalar_t *, Kokkos::HostSpace> invDiag("diag^-1", numRows);
  for (lno_t i = 0; i < numRows; i++) {
    for (size_type j = rowmap(i); j < rowmap(i + 1); j++) {
      if (entries(j) == i) invDiag(i) = one / values(j);
    }
  }
  for (int i = 0; i < 1; i++) {
    KokkosSparse::Impl::Sequential::gaussSeidel<lno_t, size_type, scalar_t, scalar_t, scalar_t>(
        numRows, 1, rowmap.data(), entries.data(), values.data(), y_host.data(), numRows, x_host.data(), numRows,
        invDiag.data(),
        one,  // omega
        "F");
    KokkosSparse::Impl::Sequential::gaussSeidel<lno_t, size_type, scalar_t, scalar_t, scalar_t>(
        numRows, 1, rowmap.data(), entries.data(), values.data(), y_host.data(), numRows, x_host.data(), numRows,
        invDiag.data(),
        one,  // omega
        "B");
  }
  // Copy solution back
  Kokkos::deep_copy(x, x_host);
  // Check against gold solution
  scalar_t xSq             = KokkosBlas::dot(x, x);
  scalar_t solnDot         = KokkosBlas::dot(x, xgold);
  double scaledSolutionDot = Kokkos::ArithTraits<scalar_t>::abs(solnDot / xSq);
  EXPECT_TRUE(0.99 < scaledSolutionDot);
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_balloon_clustering(lno_t numRows, size_type nnzPerRow, lno_t bandwidth) {
  using namespace Test;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type const_lno_row_view_t;
  typedef typename graph_t::entries_type const_lno_nnz_view_t;
  typedef typename graph_t::row_map_type::non_const_type lno_row_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef KokkosKernelsHandle<size_type, lno_t, scalar_t, typename device::execution_space,
                              typename device::memory_space, typename device::memory_space>
      KernelHandle;
  srand(245);
  size_type nnzTotal = nnzPerRow * numRows;
  lno_t nnzVariance  = nnzPerRow / 4;
  crsMat_t A =
      KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows, numRows, nnzTotal, nnzVariance, bandwidth);
  lno_row_view_t symRowmap;
  lno_nnz_view_t symEntries;
  KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap<const_lno_row_view_t, const_lno_nnz_view_t, lno_row_view_t,
                                                         lno_nnz_view_t, typename device::execution_space>(
      numRows, A.graph.row_map, A.graph.entries, symRowmap, symEntries);
  KokkosSparse::Impl::BalloonClustering<KernelHandle, lno_row_view_t, lno_nnz_view_t> balloon(numRows, symRowmap,
                                                                                              symEntries);
  for (int clusterSize = 1; clusterSize <= numRows / 16; clusterSize = std::ceil(clusterSize * 1.3)) {
    auto vertClusters = balloon.run(clusterSize);
    // validate results: make sure cluster labels are in bounds, and that the
    // number of clusters is correct
    auto vertClustersHost = Kokkos::create_mirror_view(vertClusters);
    Kokkos::deep_copy(vertClustersHost, vertClusters);
    lno_t numClusters = (numRows + clusterSize - 1) / clusterSize;
    // check the hard constraints of the clustering
    std::set<lno_t> uniqueClusterIDs;
    for (lno_t i = 0; i < numRows; i++) {
      EXPECT_TRUE(vertClustersHost(i) >= 0);
      EXPECT_TRUE(vertClustersHost(i) < numClusters);
      uniqueClusterIDs.insert(vertClustersHost(i));
    }
    EXPECT_TRUE(uniqueClusterIDs.size() == static_cast<size_t>(numClusters));
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_gauss_seidel_empty() {
  using namespace Test;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type row_map_type;
  typedef typename graph_t::entries_type::non_const_type entries_type;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef KokkosKernelsHandle<size_type, lno_t, scalar_t, typename device::execution_space,
                              typename device::memory_space, typename device::memory_space>
      KernelHandle;
  // The rowmap of a zero-row matrix can be length 0 or 1, so Gauss-Seidel
  // should work with both (the setup and apply are essentially no-ops but they
  // shouldn't crash or throw exceptions) For this test, create size-0 and
  // size-1 rowmaps separately, and make sure each work with both point and
  // cluster. Check also 5x5 matrix with empty rows (0-nnz), which can trigger
  // different bugs.
  for (int doingCluster = 0; doingCluster < 2; doingCluster++) {
    for (const int rowmapLen : {0, 1, 5}) {
      KernelHandle kh;
      if (doingCluster)
        kh.create_gs_handle(CLUSTER_DEFAULT, 10);
      else
        kh.create_gs_handle(GS_DEFAULT);
      const auto nRows = KOKKOSKERNELS_MACRO_MAX(0, rowmapLen - 1);
      // initialized to 0
      row_map_type rowmap("Rowmap", rowmapLen);
      entries_type entries("Entries", 0);
      scalar_view_t values("Values", 0);
      // also, make sure graph symmetrization doesn't crash on zero rows
      gauss_seidel_symbolic(&kh, nRows, nRows, rowmap, entries, false);
      gauss_seidel_numeric(&kh, nRows, nRows, rowmap, entries, values, false);
      scalar_view_t x("X", nRows);
      scalar_view_t y("Y", nRows);
      scalar_t omega(0.9);
      symmetric_gauss_seidel_apply(&kh, nRows, nRows, rowmap, entries, values, x, y, false, true, omega, 3);
      kh.destroy_gs_handle();
    }
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_gauss_seidel_long_rows(lno_t numRows, lno_t numLongRows, lno_t nnzPerShortRow, bool symmetric) {
  using namespace Test;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::index_type::non_const_type entries_view_t;
  typedef typename crsMat_t::row_map_type::non_const_type rowmap_view_t;
  typedef typename Kokkos::ArithTraits<scalar_t>::mag_type mag_t;
  const scalar_t one = Kokkos::ArithTraits<scalar_t>::one();
  srand(245);
  std::vector<size_type> rowmap = {0};
  std::vector<lno_t> entries;
  std::vector<scalar_t> values;
  std::vector<lno_t> rowLengths;
  for (lno_t i = 0; i < numRows; i++) {
    if (i < numLongRows)
      rowLengths.push_back(numRows);
    else
      rowLengths.push_back(nnzPerShortRow);
  }
  std::shuffle(rowLengths.begin(), rowLengths.end(), std::mt19937(std::random_device()()));
  size_type totalEntries = 0;
  int randSteps          = 1000000;
  scalar_t offDiagBase;
  {
    scalar_t unused;
    Test::getRandomBounds(0.6, unused, offDiagBase);
  }
  for (lno_t i = 0; i < numRows; i++) {
    for (lno_t ent = 0; ent < rowLengths[i]; ent++) {
      if (ent == 0) {
        entries.push_back(i);
        values.push_back(2.5 * one);
      } else {
        entries.push_back(rand() % numRows);
        values.push_back((-0.3 + (0.6 * (rand() % randSteps) / randSteps)) * offDiagBase);
      }
    }
    totalEntries += rowLengths[i];
    rowmap.push_back(totalEntries);
  }
  scalar_view_t valuesView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values"), totalEntries);
  entries_view_t entriesView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries"), totalEntries);
  rowmap_view_t rowmapView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Rowmap"), numRows + 1);
  Kokkos::deep_copy(valuesView, Kokkos::View<scalar_t *, Kokkos::HostSpace>(values.data(), totalEntries));
  Kokkos::deep_copy(entriesView, Kokkos::View<lno_t *, Kokkos::HostSpace>(entries.data(), totalEntries));
  Kokkos::deep_copy(rowmapView, Kokkos::View<size_type *, Kokkos::HostSpace>(rowmap.data(), numRows + 1));
  crsMat_t input_mat("A", numRows, numRows, totalEntries, valuesView, rowmapView, entriesView);
  input_mat = KokkosSparse::sort_and_merge_matrix(input_mat);
  if (symmetric) {
    // Symmetrize on host, rather than relying on the parallel versions (those
    // can be tested for symmetric=false)
    input_mat = Test::symmetrize<scalar_t, lno_t, size_type, device, crsMat_t>(input_mat);
  }
  lno_t nv = input_mat.numRows();
  scalar_view_t solution_x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "X (correct)"), nv);
  create_random_x_vector(solution_x);
  mag_t initial_norm_res = KokkosBlas::nrm2(solution_x);
  scalar_view_t y_vector = create_random_y_vector(input_mat, solution_x);
  // GS_DEFAULT is GS_TEAM on CUDA and GS_PERMUTED on other spaces, and the
  // behavior of each algorithm _should be_ the same on every execution space,
  // which is why we just test GS_DEFAULT.
  int apply_count = 1;  // test symmetric, forward, backward
  scalar_view_t x_vector(Kokkos::view_alloc(Kokkos::WithoutInitializing, "x vector"), nv);
  for (int apply_type = 0; apply_type < apply_count; ++apply_type) {
    typedef KokkosKernelsHandle<size_type, lno_t, scalar_t, typename device::execution_space,
                                typename device::memory_space, typename device::memory_space>
        KernelHandle;

    KernelHandle kh;
    kh.create_gs_handle(GS_DEFAULT);
    auto gsHandle = kh.get_point_gs_handle();
    gsHandle->set_long_row_threshold(3 * nnzPerShortRow);
    // Reset x vector to 0
    Kokkos::deep_copy(x_vector, scalar_t());
    run_gauss_seidel(kh, input_mat, x_vector, y_vector, symmetric, 0.9, apply_type);
    KokkosBlas::axpby(one, solution_x, -one, x_vector);
    mag_t result_norm_res = KokkosBlas::nrm2(x_vector);
    EXPECT_LT(result_norm_res, 0.25 * initial_norm_res);
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_gauss_seidel_custom_coloring(lno_t numRows, lno_t nnzPerRow) {
  using namespace Test;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename Kokkos::ArithTraits<scalar_t>::mag_type mag_t;
  const scalar_t one = Kokkos::ArithTraits<scalar_t>::one();
  size_type nnz      = nnzPerRow * numRows;
  crsMat_t input_mat = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<crsMat_t>(
      numRows, numRows, nnz, 0, numRows / 10, 2.0 * one);
  input_mat = Test::symmetrize<scalar_t, lno_t, size_type, device, crsMat_t>(input_mat);
  input_mat = KokkosSparse::sort_and_merge_matrix(input_mat);
  scalar_view_t solution_x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "X (correct)"), numRows);
  create_random_x_vector(solution_x);
  mag_t initial_norm_res = KokkosBlas::nrm2(solution_x);
  scalar_view_t y_vector = create_random_y_vector(input_mat, solution_x);
  scalar_view_t x_vector(Kokkos::view_alloc(Kokkos::WithoutInitializing, "x vector"), numRows);
  typedef KokkosKernelsHandle<size_type, lno_t, scalar_t, typename device::execution_space,
                              typename device::memory_space, typename device::memory_space>
      KernelHandle;

  KernelHandle kh;
  kh.create_gs_handle(GS_DEFAULT, KokkosGraph::COLORING_VBBIT);
  EXPECT_EQ(kh.get_point_gs_handle()->get_coloring_algorithm(), KokkosGraph::COLORING_VBBIT);
  // Reset x vector to 0
  Kokkos::deep_copy(x_vector, scalar_t());
  run_gauss_seidel(kh, input_mat, x_vector, y_vector, true, 0.9, 0);
  KokkosBlas::axpby(one, solution_x, -one, x_vector);
  mag_t result_norm_res = KokkosBlas::nrm2(x_vector);
  EXPECT_LT(result_norm_res, 0.25 * initial_norm_res);
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_gauss_seidel_streams_rank1(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance,
                                     bool symmetric, double omega,
                                     KokkosGraph::ColoringAlgorithm coloringAlgo = KokkosGraph::COLORING_DEFAULT,
                                     int nstreams                                = 1) {
  using namespace Test;
  using crsMat_t        = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type>;
  using scalar_view_t   = typename crsMat_t::values_type::non_const_type;
  using mag_t           = typename Kokkos::ArithTraits<scalar_t>::mag_type;
  using execution_space = typename device::execution_space;

  using const_size_type = const size_type;
  using const_lno_t     = const lno_t;
  using const_scalar_t  = const scalar_t;
  using KernelHandle    = KokkosKernelsHandle<const_size_type, const_lno_t, const_scalar_t, execution_space,
                                           typename device::memory_space, typename device::memory_space>;
  srand(245);
  lno_t numCols                         = numRows;
  typename crsMat_t::value_type m_omega = omega;

#ifdef KOKKOS_ENABLE_OPENMP
  if (std::is_same_v<execution_space, Kokkos::OpenMP>) {
    int exec_concurrency = execution_space().concurrency();
    if (exec_concurrency < nstreams) {
      std::cerr << "TEST SKIPPED: Not enough concurrency to partition "
                   "execution space. exec_concurrency: "
                << exec_concurrency << std::endl;
      return;
    }
  }
#endif  // KOKKOS_ENABLE_OPENMP

  auto instances = Kokkos::Experimental::partition_space(execution_space(), std::vector<int>(nstreams, 1));

  std::vector<KernelHandle> kh_v(nstreams);
  std::vector<crsMat_t> input_mat_v(nstreams);
  std::vector<scalar_view_t> solution_x_v(nstreams);
  std::vector<scalar_view_t> x_vector_v(nstreams);
  std::vector<scalar_view_t> y_vector_v(nstreams);
  std::vector<mag_t> initial_norm_res_v(nstreams);

  const scalar_t one  = Kokkos::ArithTraits<scalar_t>::one();
  const scalar_t zero = Kokkos::ArithTraits<scalar_t>::zero();

  for (int i = 0; i < nstreams; i++) {
    input_mat_v[i] = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<crsMat_t>(
        numRows, numCols, nnz, row_size_variance, bandwidth);

    if (symmetric) {
      // Symmetrize on host, rather than relying on the parallel versions (those
      // can be tested for symmetric=false)
      input_mat_v[i] = Test::symmetrize<scalar_t, lno_t, size_type, device, crsMat_t>(input_mat_v[i]);
    }
    lno_t nv = input_mat_v[i].numRows();
    scalar_view_t solution_x_tmp(Kokkos::view_alloc(Kokkos::WithoutInitializing, "X (correct)"), nv);
    solution_x_v[i] = solution_x_tmp;
    create_random_x_vector(solution_x_v[i]);
    initial_norm_res_v[i] = KokkosBlas::nrm2(solution_x_v[i]);
    y_vector_v[i]         = create_random_y_vector(input_mat_v[i], solution_x_v[i]);
    // GS_DEFAULT is GS_TEAM on CUDA and GS_PERMUTED on other spaces, and the
    // behavior of each algorithm _should be_ the same on every execution space,
    // which is why we just test GS_DEFAULT.

    scalar_view_t x_vector_tmp(Kokkos::view_alloc(Kokkos::WithoutInitializing, "x vector"), nv);
    x_vector_v[i] = x_vector_tmp;

    kh_v[i] = KernelHandle();  // Initialize KokkosKernelsHandle defaults.
    kh_v[i].create_gs_handle(instances[i], nstreams, GS_DEFAULT, coloringAlgo);
  }

  int apply_count = 3;  // test symmetric, forward, backward
  //*** Point-coloring version ****
  for (int apply_type = 0; apply_type < apply_count; ++apply_type) {
    for (int i = 0; i < nstreams; i++) Kokkos::deep_copy(instances[i], x_vector_v[i], zero);

    run_gauss_seidel_streams(instances, kh_v, input_mat_v, x_vector_v, y_vector_v, symmetric, m_omega, apply_type,
                             nstreams);
    for (int i = 0; i < nstreams; i++) {
      KokkosBlas::axpby(instances[i], one, solution_x_v[i], -one, x_vector_v[i]);
      mag_t result_norm_res = KokkosBlas::nrm2(instances[i], x_vector_v[i]);
      EXPECT_LT(result_norm_res, initial_norm_res_v[i]) << "on stream_idx: " << i;
    }
  }

  for (int i = 0; i < nstreams; i++) kh_v[i].destroy_gs_handle();
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                                                    \
  TEST_F(TestCategory, sparse##_##gauss_seidel_asymmetric_rank1##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {        \
    test_gauss_seidel_rank1<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 2000 * 20, 200, 10, false);                         \
  }                                                                                                                    \
  TEST_F(TestCategory,                                                                                                 \
         sparse##_##gauss_seidel_asymmetric_streams_rank1##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {              \
    test_gauss_seidel_streams_rank1<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 2000 * 20, 200, 10, false, 0.9,             \
                                                                     KokkosGraph::COLORING_DEFAULT, 1);                \
    test_gauss_seidel_streams_rank1<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 2000 * 20, 200, 10, false, 0.9,             \
                                                                     KokkosGraph::COLORING_DEFAULT, 2);                \
    test_gauss_seidel_streams_rank1<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 2000 * 20, 200, 10, false, 0.9,             \
                                                                     KokkosGraph::COLORING_DEFAULT, 3);                \
    test_gauss_seidel_streams_rank1<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 2000 * 20, 200, 10, false, 0.9,             \
                                                                     KokkosGraph::COLORING_DEFAULT, 4);                \
  }                                                                                                                    \
  TEST_F(TestCategory, sparse##_##gauss_seidel_asymmetric_rank2##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {        \
    test_gauss_seidel_rank2<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 2000 * 20, 200, 10, 3, false);                      \
  }                                                                                                                    \
  TEST_F(TestCategory, sparse##_##gauss_seidel_symmetric_rank1##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {         \
    test_gauss_seidel_rank1<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 2000 * 20, 200, 10, true);                          \
  }                                                                                                                    \
  TEST_F(TestCategory, sparse##_##gauss_seidel_symmetric_streams_rank1##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_gauss_seidel_streams_rank1<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 2000 * 20, 200, 10, true, 0.9,              \
                                                                     KokkosGraph::COLORING_DEFAULT, 1);                \
    test_gauss_seidel_streams_rank1<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 2000 * 20, 200, 10, true, 0.9,              \
                                                                     KokkosGraph::COLORING_DEFAULT, 2);                \
    test_gauss_seidel_streams_rank1<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 2000 * 20, 200, 10, true, 0.9,              \
                                                                     KokkosGraph::COLORING_DEFAULT, 3);                \
    test_gauss_seidel_streams_rank1<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 2000 * 20, 200, 10, true, 0.9,              \
                                                                     KokkosGraph::COLORING_DEFAULT, 4);                \
  }                                                                                                                    \
  TEST_F(TestCategory, sparse##_##gauss_seidel_symmetric_rank2##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {         \
    test_gauss_seidel_rank2<SCALAR, ORDINAL, OFFSET, DEVICE>(2000, 2000 * 20, 200, 10, 3, true);                       \
  }                                                                                                                    \
  TEST_F(TestCategory, sparse##_##gauss_seidel_empty##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {                   \
    test_gauss_seidel_empty<SCALAR, ORDINAL, OFFSET, DEVICE>();                                                        \
  }                                                                                                                    \
  TEST_F(TestCategory, sparse##_##balloon_clustering##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {                   \
    test_balloon_clustering<SCALAR, ORDINAL, OFFSET, DEVICE>(5000, 100, 2000);                                         \
  }                                                                                                                    \
  TEST_F(TestCategory, sparse##_##sequential_sor##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {                       \
    test_sequential_sor<SCALAR, ORDINAL, OFFSET, DEVICE>(1000, 1000 * 15, 50, 10);                                     \
  }                                                                                                                    \
  TEST_F(TestCategory, sparse##_##gauss_seidel_long_rows##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {               \
    test_gauss_seidel_long_rows<SCALAR, ORDINAL, OFFSET, DEVICE>(500, 10, 20, true);                                   \
  }                                                                                                                    \
  TEST_F(TestCategory, sparse##_##gauss_seidel_custom_coloring##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {         \
    test_gauss_seidel_custom_coloring<SCALAR, ORDINAL, OFFSET, DEVICE>(500, 10);                                       \
  }

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
