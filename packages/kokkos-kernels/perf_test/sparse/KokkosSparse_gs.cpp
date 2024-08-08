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

#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>
#include <KokkosKernels_Handle.hpp>
#include <KokkosKernels_TestUtils.hpp>
#include <KokkosSparse_gauss_seidel.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <KokkosKernels_config.h>
#include "KokkosKernels_default_types.hpp"
#include "KokkosSparse_IOUtils.hpp"
#include <iostream>
#include <random>
#include <vector>
#include <string>
#include <unordered_set>

using std::cout;
using std::string;
using namespace KokkosSparse;

static char* getNextArg(int& i, int argc, char** argv) {
  i++;
  if (i >= argc) {
    std::cerr << "Error: expected additional command-line argument!\n";
    exit(1);
  }
  return argv[i];
}

struct GS_Parameters {
  const char* matrix_path = nullptr;
  int n                   = 10000;
  int nnzPerRow           = 27;
  int numLongRows         = 0;
  int minNnzPerLongRow    = 1000;
  int maxNnzPerLongRow    = 2000;
  bool graph_symmetric    = false;
  int sweeps              = 1;
  int nstreams            = 1;
  GSAlgorithm algo        = GS_DEFAULT;
  GSDirection direction   = GS_FORWARD;
  // Point:
  int longRowThreshold = 0;
  // Cluster:
  ClusteringAlgorithm coarse_algo = CLUSTER_DEFAULT;
  int cluster_size                = 10;
  // Two stage:
  bool classic = false;
};

template <typename crsMat_t>
crsMat_t generateLongRowMatrix(const GS_Parameters& params) {
  typedef typename crsMat_t::value_type scalar_t;
  typedef typename crsMat_t::ordinal_type lno_t;
  typedef typename crsMat_t::size_type size_type;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::index_type::non_const_type entries_view_t;
  typedef typename crsMat_t::row_map_type::non_const_type rowmap_view_t;
  typedef typename crsMat_t::device_type device;
  // Generate random diag. dominant matrix
  srand(245);
  std::vector<size_type> rowmap = {0};
  std::vector<lno_t> entries;
  std::vector<scalar_t> values;
  std::vector<lno_t> rowLengths;
  lno_t numRows = params.n;
  for (lno_t i = 0; i < numRows; i++) {
    if (i < params.numLongRows) {
      lno_t interval = params.maxNnzPerLongRow - params.minNnzPerLongRow;
      lno_t rowLen;
      if (interval == 0)
        rowLen = params.maxNnzPerLongRow;
      else
        rowLen = params.minNnzPerLongRow + rand() % interval;
      if (rowLen > numRows) rowLen = numRows;
      rowLengths.push_back(rowLen);
    } else
      rowLengths.push_back(params.nnzPerRow);
  }
  std::shuffle(rowLengths.begin(), rowLengths.end(), std::mt19937(std::random_device()()));
  size_type totalEntries = 0;
  int randSteps          = 1000000;
  // Set of columns inserted so far into current short row
  std::unordered_set<lno_t> shortRowEntries;
  // Set of all possible rows, randomly permuted (select long row entries from
  // head)
  std::vector<lno_t> longRowEntries(numRows);
  for (lno_t i = 0; i < numRows; i++) longRowEntries[i] = i;
  const scalar_t one = Kokkos::reduction_identity<scalar_t>::prod();
  for (lno_t i = 0; i < numRows; i++) {
    shortRowEntries.clear();
    bool rowIsLong = rowLengths[i] > params.nnzPerRow;
    if (rowIsLong) std::shuffle(longRowEntries.begin(), longRowEntries.end(), std::mt19937(std::random_device()()));
    for (lno_t ent = 0; ent < rowLengths[i]; ent++) {
      if (ent == 0) {
        entries.push_back(i);
        values.push_back(5.0 + 3.0 * (rand() % randSteps) / randSteps * one);
      } else {
        if (rowIsLong)
          entries.push_back(longRowEntries[ent]);
        else {
          // re-roll random column until one is found that isn't already in row
          lno_t col;
          while (true) {
            col = rand() % numRows;
            if (shortRowEntries.find(col) == shortRowEntries.end()) {
              shortRowEntries.insert(col);
              break;
            }
          }
          entries.push_back(col);
        }
        values.push_back((-0.1 + (0.2 * (rand() % randSteps) / randSteps)) * one);
      }
    }
    totalEntries += rowLengths[i];
    rowmap.push_back(totalEntries);
  }
  scalar_view_t valuesView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values"), totalEntries);
  entries_view_t entriesView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries"), totalEntries);
  rowmap_view_t rowmapView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Rowmap"), numRows + 1);
  Kokkos::deep_copy(valuesView, Kokkos::View<scalar_t*, Kokkos::HostSpace>(values.data(), totalEntries));
  Kokkos::deep_copy(entriesView, Kokkos::View<lno_t*, Kokkos::HostSpace>(entries.data(), totalEntries));
  Kokkos::deep_copy(rowmapView, Kokkos::View<size_type*, Kokkos::HostSpace>(rowmap.data(), numRows + 1));
  crsMat_t A("A", numRows, numRows, totalEntries, valuesView, rowmapView, entriesView);
  A = KokkosSparse::sort_and_merge_matrix(A);
  if (params.graph_symmetric) {
    // Symmetrize on host, rather than relying on the parallel versions (those
    // can be tested for symmetric=false)
    A = Test::symmetrize<scalar_t, lno_t, size_type, device, crsMat_t>(A);
  }
  return A;
}

template <typename device_t>
void runGS(const GS_Parameters& params) {
  typedef default_scalar scalar_t;
  typedef default_lno_t lno_t;
  typedef default_size_type size_type;
  typedef typename device_t::execution_space exec_space;
  typedef typename device_t::memory_space mem_space;
  typedef KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, exec_space, mem_space, mem_space>
      KernelHandle;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device_t, void, size_type> crsMat_t;
  // typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  crsMat_t A;
  if (params.matrix_path)
    A = KokkosSparse::Impl::read_kokkos_crst_matrix<crsMat_t>(params.matrix_path);
  else
    A = generateLongRowMatrix<crsMat_t>(params);
  lno_t nrows = A.numRows();
  lno_t ncols = A.numCols();
  if (nrows != ncols) {
    cout << "ERROR: Gauss-Seidel only works for square matrices\n";
    Kokkos::finalize();
    exit(1);
  }
  std::vector<exec_space> instances;
  // size_type nnz = A.nnz();
  std::vector<KernelHandle> kh(params.nstreams);
  // use a random RHS - uniformly distributed over (-5, 5)
  std::vector<scalar_view_t> b(params.nstreams);
  // initial LHS is 0
  std::vector<scalar_view_t> x(params.nstreams);
  // Extract diagonal blocks of CRS matrix
  std::vector<crsMat_t> DiagBlks(params.nstreams);
  // how long symbolic/numeric phases take (the graph reuse case isn't that
  // interesting since numeric doesn't do much)
  Kokkos::Timer timer;

  {
    namespace KE = Kokkos::Experimental;
    auto ns      = params.nstreams;
    auto es      = exec_space();
    std::vector<int> weights(ns);
    std::fill(weights.begin(), weights.end(), 1);
    instances = KE::partition_space(es, weights);
  }

  double blockExtractionTime = 0, symbolicLaunchTimeTotal = 0, symbolicComputeTimeTotal = 0, numericLaunchTimeTotal = 0,
         numericComputeTimeTotal = 0, applyLaunchTimeTotal = 0, applyComputeTimeTotal = 0;

  timer.reset();
  KokkosSparse::Impl::kk_extract_diagonal_blocks_crsmatrix_sequential(A, DiagBlks);
  Kokkos::fence();
  blockExtractionTime = timer.seconds();

  /////////////////// Handle creation ///////////////////
  for (int i = 0; i < params.nstreams; i++) {
    auto blk_A     = DiagBlks[i];
    auto blk_nrows = blk_A.numRows();
    auto blk_ncols = blk_A.numCols();
    if (blk_nrows != blk_ncols) {
      cout << "ERROR: Gauss-Seidel only works for square matrices\n";
      Kokkos::finalize();
      exit(1);
    }
    b[i] = scalar_view_t("b[" + std::to_string(i) + "]", blk_nrows);
    x[i] = scalar_view_t("x[" + std::to_string(i) + "]", blk_nrows);
    {
      srand(54321 + i);
      auto bhost = Kokkos::create_mirror_view(b[i]);
      for (lno_t row_id = 0; row_id < blk_nrows; row_id++) {
        bhost(row_id) = 10.0 * rand() / RAND_MAX - 5.0;
      }
      Kokkos::deep_copy(instances[i], b[i], bhost);
    }
    // cluster size of 1 is standard multicolor GS
    if (params.algo == GS_DEFAULT) {
      kh[i].create_gs_handle(instances[i], params.nstreams);
      kh[i].get_point_gs_handle()->set_long_row_threshold(params.longRowThreshold);
    } else if (params.algo == GS_CLUSTER) {
      kh[i].create_gs_handle(params.coarse_algo, params.cluster_size);
    } else {
      kh[i].create_gs_handle(params.algo);
      if (params.algo == GS_TWOSTAGE) kh[i].set_gs_twostage(!params.classic, blk_nrows);
    }
  }

  /////////////////// Symbolic /////////////////
  timer.reset();
  for (int i = 0; i < params.nstreams; i++) {
    auto blk_A     = DiagBlks[i];
    auto blk_nrows = blk_A.numRows();
    KokkosSparse::Experimental::gauss_seidel_symbolic(instances[i], &kh[i], blk_nrows, blk_nrows, blk_A.graph.row_map,
                                                      blk_A.graph.entries, params.graph_symmetric);
  }
  symbolicLaunchTimeTotal = timer.seconds();
  timer.reset();
  Kokkos::fence();
  symbolicComputeTimeTotal = timer.seconds();

  /////////////////// Numeric /////////////////
  timer.reset();
  for (int i = 0; i < params.nstreams; i++) {
    auto blk_A     = DiagBlks[i];
    auto blk_nrows = blk_A.numRows();
    KokkosSparse::Experimental::gauss_seidel_numeric(instances[i], &kh[i], blk_nrows, blk_nrows, blk_A.graph.row_map,
                                                     blk_A.graph.entries, blk_A.values, params.graph_symmetric);
  }
  numericLaunchTimeTotal = timer.seconds();
  timer.reset();
  Kokkos::fence();
  numericComputeTimeTotal = timer.seconds();

  /////////////////// Apply /////////////////
  timer.reset();
  for (int i = 0; i < params.nstreams; i++) {
    auto blk_A     = DiagBlks[i];
    auto blk_nrows = blk_A.numRows();
    // Last two parameters are damping factor (should be 1) and sweeps
    switch (params.direction) {
      case GS_SYMMETRIC:
        KokkosSparse::Experimental::symmetric_gauss_seidel_apply(instances[i], &kh[i], blk_nrows, blk_nrows,
                                                                 blk_A.graph.row_map, blk_A.graph.entries, blk_A.values,
                                                                 x[i], b[i], true, true, 1.0, params.sweeps);
        break;
      case GS_FORWARD:
        KokkosSparse::Experimental::forward_sweep_gauss_seidel_apply(
            instances[i], &kh[i], blk_nrows, blk_nrows, blk_A.graph.row_map, blk_A.graph.entries, blk_A.values, x[i],
            b[i], true, true, 1.0, params.sweeps);
        break;
      case GS_BACKWARD:
        KokkosSparse::Experimental::backward_sweep_gauss_seidel_apply(
            instances[i], &kh[i], blk_nrows, blk_nrows, blk_A.graph.row_map, blk_A.graph.entries, blk_A.values, x[i],
            b[i], true, true, 1.0, params.sweeps);
        break;
    }
  }
  applyLaunchTimeTotal = timer.seconds();
  timer.reset();
  Kokkos::fence();
  applyComputeTimeTotal = timer.seconds();
  timer.reset();

  for (int i = 0; i < params.nstreams; i++) {
    auto blk_A     = DiagBlks[i];
    auto blk_nrows = blk_A.numRows();
    kh[i].destroy_gs_handle();
    // Now, compute the 2-norm of residual
    scalar_view_t res("Ax-b", blk_nrows);
    Kokkos::deep_copy(instances[i], res, b[i]);
    double bnorm   = KokkosBlas::nrm2(instances[i], b[i]);
    scalar_t alpha = Kokkos::reduction_identity<scalar_t>::prod();
    scalar_t beta  = -alpha;
    KokkosSparse::spmv<exec_space, scalar_t, crsMat_t, scalar_view_t, scalar_t, scalar_view_t>(instances[i], "N", alpha,
                                                                                               blk_A, x[i], beta, res);
    double resnorm = KokkosBlas::nrm2(instances[i], res);
    // note: this still works if the solution diverges
    std::cout << "StreamID(" << i << "): Relative res norm: " << resnorm / bnorm << '\n';
  }
  std::cout << "\n*** Total block extraction time: " << blockExtractionTime << '\n';
  std::cout << "\n*** Total Symbolic launch time: " << symbolicLaunchTimeTotal << '\n';
  std::cout << "*** Total Symbolic compute time: " << symbolicComputeTimeTotal << '\n';
  std::cout << "\n*** Total Numeric launch time: " << numericLaunchTimeTotal << '\n';
  std::cout << "*** Total Numeric compute time: " << numericComputeTimeTotal << '\n';
  std::cout << "\n*** Total Apply launch time: " << applyLaunchTimeTotal << '\n';
  std::cout << "*** Total Apply compute time: " << applyComputeTimeTotal << '\n';
  double launchTimeTotal = symbolicLaunchTimeTotal + numericLaunchTimeTotal + applyLaunchTimeTotal;
  std::cout << "\n*** Total launch time: " << launchTimeTotal << '\n';
  double computeTimeTotal = symbolicComputeTimeTotal + numericComputeTimeTotal + applyComputeTimeTotal;
  std::cout << "*** Total compute time: " << computeTimeTotal << '\n';
  std::cout << "\n*** Total compute and launch time: " << launchTimeTotal + computeTimeTotal << '\n';
}

int main(int argc, char** argv) {
  // Expect two args: matrix name and device flag.
  if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
    cout << "Usage: ./sparse_gs [--device] --amtx matrix.mtx [other args]\n\n";
    cout << "\"--device-type\" flag can be \"--serial\", \"--openmp\", "
            "\"--cuda\" or \"--threads\".\n";
    cout << "If device is not given, the default device for this build is "
            "used.\n";
    cout << "\nOther flags:\n";
    cout << "--sym-graph : pass if matrix is known to be structurally "
            "symmetric.\n";
    cout << "            : if generating matrix randomly, it is symmetrized\n";
    cout << "--sweeps S: run S times (default 1)\n";
    cout << "--streams N: partition matrix and run across N streams (default "
            "1)\n";
    cout << "Randomized matrix settings, if not reading from file:\n";
    cout << "  --n <N> : number of rows/columns\n";
    cout << "  --nnz <N> : number of nonzeros in each regular row\n";
    cout << "  --long-rows <N> : number of long rows\n";
    cout << "  --min-long-row-nnz <N> : min number of nonzeros in each long "
            "row (default 1000)\n";
    cout << "  --max-long-row-nnz <N> : max number of nonzeros in each long "
            "row (default 2000)\n";
    cout << "Randomized matrix settings, if not reading from file:\n";
    cout << "Randomized matrix settings, if not reading from file:\n";
    cout << "4 main algorithms (required, choose one):\n";
    cout << "  --point\n";
    cout << "  --cluster\n";
    cout << "  --twostage\n";
    cout << "  --classic\n\n";
    cout << "Apply direction (default is forward)\n";
    cout << "  --forward\n";
    cout << "  --backward\n";
    cout << "  --symmetric\n";
    cout << "Options for point:\n";
    cout << "  --long-row-threshold <N> : rows with at least this many entries "
            "are processed separately.\n";
    cout << "Options for cluster:\n";
    cout << "  --cluster-size N (default: 10)\n";
    cout << "  --coarse-algo ALGO\n";
    cout << "     ALGO may be: \"balloon\" or \"mis2\"\n";
    cout << "     Default is chosen by the library. If using mis2, "
            "--cluster-size option has no effect.\n";
    return 0;
  }
  Kokkos::initialize(argc, argv);
  // device is just the name of the execution space, lowercase
  string deviceName;
  GS_Parameters params;
  int i = 1;
  for (; i < argc; i++) {
    if (!strcmp(argv[i], "--amtx"))
      params.matrix_path = getNextArg(i, argc, argv);
    else if (!strcmp(argv[i], "--n"))
      params.n = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--nnz"))
      params.nnzPerRow = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--long-rows"))
      params.numLongRows = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--min-long-row-nnz"))
      params.minNnzPerLongRow = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--max-long-row-nnz"))
      params.maxNnzPerLongRow = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--serial"))
      deviceName = "serial";
    else if (!strcmp(argv[i], "--openmp"))
      deviceName = "openmp";
    else if (!strcmp(argv[i], "--threads"))
      deviceName = "threads";
    else if (!strcmp(argv[i], "--cuda"))
      deviceName = "cuda";
    else if (!strcmp(argv[i], "--sym-graph"))
      params.graph_symmetric = true;
    else if (!strcmp(argv[i], "--symmetric"))
      params.direction = GS_SYMMETRIC;
    else if (!strcmp(argv[i], "--forward"))
      params.direction = GS_FORWARD;
    else if (!strcmp(argv[i], "--backward"))
      params.direction = GS_BACKWARD;
    else if (!strcmp(argv[i], "--sweeps"))
      params.sweeps = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--streams"))
      params.nstreams = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--point"))
      params.algo = GS_DEFAULT;
    else if (!strcmp(argv[i], "--cluster"))
      params.algo = GS_CLUSTER;
    else if (!strcmp(argv[i], "--twostage"))
      params.algo = GS_TWOSTAGE;
    else if (!strcmp(argv[i], "--classic")) {
      params.algo    = GS_TWOSTAGE;
      params.classic = true;
    } else if (!strcmp(argv[i], "--long-row-threshold"))
      params.longRowThreshold = atoi(getNextArg(i, argc, argv));
    else if (!strcmp(argv[i], "--coarse-algo")) {
      const char* algo = getNextArg(i, argc, argv);
      if (!strcmp(algo, "balloon"))
        params.coarse_algo = CLUSTER_BALLOON;
      else if (!strcmp(algo, "mis2"))
        params.coarse_algo = CLUSTER_MIS2;
      else {
        std::cout << "Error: invalid coarsening algorithm. Options are balloon "
                     "and mis2.\n";
        Kokkos::finalize();
        exit(1);
      }
    } else if (!strcmp(argv[i], "--cluster-size"))
      params.cluster_size = atoi(getNextArg(i, argc, argv));
    else {
      cout << "Error: unknown argument " << argv[i] << '\n';
      Kokkos::finalize();
      exit(1);
    }
  }
  bool run = false;
  if (!deviceName.length()) {
    runGS<Kokkos::DefaultExecutionSpace>(params);
    run = true;
  }
#ifdef KOKKOS_ENABLE_SERIAL
  if (deviceName == "serial") {
    runGS<Kokkos::Serial>(params);
    run = true;
  }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
  if (deviceName == "openmp") {
    runGS<Kokkos::OpenMP>(params);
    run = true;
  }
#endif
#ifdef KOKKOS_ENABLE_THREADS
  if (deviceName == "threads") {
    runGS<Kokkos::Threads>(params);
    run = true;
  }
#endif
#ifdef KOKKOS_ENABLE_CUDA
  if (deviceName == "cuda") {
    runGS<Kokkos::Cuda>(params);
    run = true;
  }
#endif
  if (!run) {
    std::cerr << "Error: device " << deviceName << " was requested but it's not enabled in this build.\n";
    return 1;
  }
  Kokkos::finalize();
  return 0;
}
