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

#include <cstdio>

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <limits.h>
#include <cmath>
#include <unordered_map>

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_IOUtils.hpp>
#include <KokkosSparse_Utils.hpp>  // for graph_max_degree
#include <KokkosSparse_spmv.hpp>
#include "KokkosKernels_default_types.hpp"

using Scalar  = default_scalar;
using Ordinal = default_lno_t;
using Offset  = default_size_type;
using KAT     = Kokkos::ArithTraits<Scalar>;

struct SPMVBenchmarking {
  // note: CLI currently only allows square matrices to be randomly generated
  // and nz/row is fixed at 10
  Ordinal num_rows     = 110503;
  Ordinal num_cols     = 110503;
  char mode            = 'N';
  int loop             = 100;
  int num_vecs         = 1;
  Scalar beta          = KAT::zero();
  std::string filename = "";
  bool flush_cache     = false;
  bool non_reuse       = false;

  // Using the parameters above, run and time spmv where x and y use the given
  // memory layout.
  template <typename Layout>
  void run() {
    using matrix_type = KokkosSparse::CrsMatrix<Scalar, Ordinal, Kokkos::DefaultExecutionSpace, void, Offset>;
    using mv_type     = Kokkos::View<Scalar**, Layout>;
    using h_mv_type   = typename mv_type::HostMirror;

    srand(17312837);
    matrix_type A;
    if (filename != "") {
      std::cout << "Reading A from file \"" << filename << "\"...\n";
      A        = KokkosSparse::Impl::read_kokkos_crst_matrix<matrix_type>(filename.c_str());
      num_rows = A.numRows();
      num_cols = A.numCols();
    } else {
      std::cout << "Randomly generating A...\n";
      Offset nnz = 10 * num_rows;
      // note: the help text says the bandwidth is fixed at 0.01 * numRows
      A = KokkosSparse::Impl::kk_generate_sparse_matrix<matrix_type>(num_rows, num_cols, nnz, 0, 0.01 * num_rows);
    }

    std::cout << "A is " << A.numRows() << "x" << A.numCols() << ", with " << A.nnz() << " nonzeros\n";
    std::cout << "Mean nnz/row: " << (double)A.nnz() / A.numRows() << '\n';
    std::cout << "Max nnz/row: "
              << KokkosSparse::Impl::graph_max_degree<Kokkos::DefaultExecutionSpace, Ordinal>(A.graph.row_map) << '\n';
    std::cout << "SpMV mode " << mode << ", " << num_vecs << " vectors, beta = " << beta << ", multivectors are ";
    std::cout << (std::is_same_v<Layout, Kokkos::LayoutLeft> ? "LayoutLeft" : "LayoutRight");
    std::cout << '\n';

    bool transpose_like = (mode == 'T') || (mode == 'H');

    Ordinal xlen = transpose_like ? A.numRows() : A.numCols();
    Ordinal ylen = transpose_like ? A.numCols() : A.numRows();

    mv_type x("X", xlen, num_vecs);
    mv_type y("Y", ylen, num_vecs);

    h_mv_type h_x         = Kokkos::create_mirror_view(x);
    h_mv_type h_y         = Kokkos::create_mirror_view(y);
    h_mv_type h_y_compare = Kokkos::create_mirror(y);

    for (int v = 0; v < num_vecs; v++) {
      for (Ordinal i = 0; i < xlen; i++) {
        h_x(i, v) = (Scalar)(1.0 * (rand() % 40) - 20.);
      }
    }

    Kokkos::deep_copy(x, h_x);

    // Benchmark
    auto x0 = Kokkos::subview(x, Kokkos::ALL(), 0);
    auto y0 = Kokkos::subview(y, Kokkos::ALL(), 0);

    // Create handles for both rank-1 and rank-2 cases,
    // even though only 1 will get used below (depending on num_vecs)

    KokkosSparse::SPMVHandle<Kokkos::DefaultExecutionSpace, matrix_type, decltype(x0), decltype(y0)> handle_rank1;
    KokkosSparse::SPMVHandle<Kokkos::DefaultExecutionSpace, matrix_type, mv_type, mv_type> handle_rank2;
    // Assuming that 1GB is enough to fully clear the L3 cache of a CPU, or the
    // L2 of a GPU. (Some AMD EPYC chips have 768 MB L3)
    Kokkos::View<char*, Kokkos::DefaultExecutionSpace> cacheFlushData;
    if (flush_cache) {
      Kokkos::resize(cacheFlushData, 1024 * 1024 * 1024);
    }

    Kokkos::DefaultExecutionSpace space;

    // Do 5 warm up calls (not timed). This will also initialize the handle.
    for (int i = 0; i < 5; i++) {
      if (num_vecs == 1) {
        // run the rank-1 version
        if (non_reuse)
          KokkosSparse::spmv(space, &mode, 1.0, A, x0, beta, y0);
        else
          KokkosSparse::spmv(space, &handle_rank1, &mode, 1.0, A, x0, beta, y0);
      } else {
        // rank-2
        if (non_reuse)
          KokkosSparse::spmv(space, &mode, 1.0, A, x, beta, y);
        else
          KokkosSparse::spmv(space, &handle_rank2, &mode, 1.0, A, x, beta, y);
      }
      space.fence();
    }

    double totalTime = 0;
    Kokkos::Timer timer;
    for (int i = 0; i < loop; i++) {
      if (flush_cache) {
        // Copy some non-zero data to the view multiple times to flush the
        // cache. (nonzero in case the system has an optimized path for zero
        // pages)
        for (int rep = 0; rep < 4; rep++) Kokkos::deep_copy(space, cacheFlushData, char(rep + 1));
      }
      space.fence();
      timer.reset();
      if (num_vecs == 1) {
        // run the rank-1 version
        if (non_reuse)
          KokkosSparse::spmv(space, &mode, 1.0, A, x0, beta, y0);
        else
          KokkosSparse::spmv(space, &handle_rank1, &mode, 1.0, A, x0, beta, y0);
      } else {
        // rank-2
        if (non_reuse)
          KokkosSparse::spmv(space, &mode, 1.0, A, x, beta, y);
        else
          KokkosSparse::spmv(space, &handle_rank2, &mode, 1.0, A, x, beta, y);
      }
      space.fence();
      totalTime += timer.seconds();
    }
    double avg_time = totalTime / loop;
    std::cout << avg_time << " s\n";
  }
};

void print_help() {
  printf("  -s [nrows]            : matrix dimension (square)\n");
  printf(
      "  --nv n                : number of columns in x/y multivector (default "
      "1).\n");
  printf(
      "  --layout left|right   : memory layout of x/y. Default depends on "
      "build's default execution space\n");
  printf(
      "  -m N|T|H|C            : matrix apply mode:\n"
      "                          N - normal, default\n"
      "                          T - transpose\n"
      "                          H - conjugate transpose\n"
      "                          C - conjugate\n");
  printf(
      "  -f [file],-fb [file]  : Read in Matrix Market (.mtx), or binary "
      "(.bin) matrix file.\n");
  printf(
      "  -l [LOOP]             : How many spmv to run to aggregate average "
      "time. \n");
  printf("  -b beta               : beta, as in y := Ax + (beta)y\n");
  printf(
      "  --flush               : Flush the cache between each spmv call "
      "(slow!)\n");
  printf(
      "  --non-reuse           : Use non-reuse interface (without "
      "SPMVHandle)\n");
}

int main(int argc, char** argv) {
  SPMVBenchmarking sb;
  char layout;
  if (std::is_same<default_layout, Kokkos::LayoutLeft>::value)
    layout = 'L';
  else
    layout = 'R';

  if (argc == 1) {
    print_help();
    return 0;
  }

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "-s") == 0)) {
      // only square matrices supported now
      sb.num_rows = atoi(argv[++i]);
      sb.num_cols = sb.num_rows;
      continue;
    }
    if ((strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "-fb") == 0)) {
      sb.filename = argv[++i];
      continue;
    }
    if ((strcmp(argv[i], "-l") == 0)) {
      sb.loop = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-m") == 0)) {
      sb.mode = toupper(argv[++i][0]);
      if (sb.mode != 'N' && sb.mode != 'T' && sb.mode != 'C' && sb.mode != 'H')
        throw std::invalid_argument("Mode must be one of N, T, C or H.");
      continue;
    }
    if ((strcmp(argv[i], "--nv") == 0)) {
      sb.num_vecs = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-b") == 0)) {
      sb.beta = atof(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "--layout") == 0)) {
      i++;
      if (toupper(argv[i][0]) == 'L')
        layout = 'L';
      else if (toupper(argv[i][0]) == 'R')
        layout = 'R';
      else
        throw std::runtime_error("Invalid layout");
      continue;
    }
    if ((strcmp(argv[i], "--flush") == 0)) {
      sb.flush_cache = true;
      continue;
    }
    if ((strcmp(argv[i], "--non-reuse") == 0)) {
      sb.non_reuse = true;
      continue;
    }
    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      print_help();
      return 0;
    }
  }

  Kokkos::initialize(argc, argv);

  if (layout == 'L')
    sb.template run<Kokkos::LayoutLeft>();
  else
    sb.template run<Kokkos::LayoutRight>();

  Kokkos::finalize();
}
