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

#include <iostream>
#include <algorithm>
#include "KokkosKernels_config.h"
#include "KokkosSparse_IOUtils.hpp"
#include "KokkosKernels_perf_test_utilities.hpp"

#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_SortCrs.hpp"

using perf_test::CommonInputParams;

struct LocalParams {
  std::string mtxFile;
};

void print_options() {
  std::cerr << "Options\n" << std::endl;

  std::cerr << perf_test::list_common_options();

  std::cerr << "\t[Required] --mtx <path> :: matrix to sort\n";
  std::cerr << "\t[Optional] --repeat      :: how many times to repeat sorting\n";
}

int parse_inputs(LocalParams& params, int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    if (perf_test::check_arg_str(i, argc, argv, "--mtx", params.mtxFile)) {
      ++i;
    } else {
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
      print_options();
      return 1;
    }
  }
  return 0;
}

template <typename exec_space>
void run_experiment(int argc, char** argv, const CommonInputParams& common_params) {
  using namespace KokkosSparse;

  using mem_space = typename exec_space::memory_space;
  using device_t  = typename Kokkos::Device<exec_space, mem_space>;
  using size_type = default_size_type;
  using lno_t     = default_lno_t;
  using scalar_t  = default_scalar;
  using crsMat_t  = KokkosSparse::CrsMatrix<scalar_t, lno_t, device_t, void, size_type>;

  using graph_t = typename crsMat_t::StaticCrsGraphType;

  LocalParams params;
  if (parse_inputs(params, argc, argv)) return;

  crsMat_t A = KokkosSparse::Impl::read_kokkos_crst_matrix<crsMat_t>(params.mtxFile.c_str());
  std::cout << "Loaded matrix: " << A.numRows() << "x" << A.numCols() << " with " << A.nnz() << " entries.\n";
  // This first sort call serves as a warm-up
  KokkosSparse::sort_crs_matrix(A);
  lno_t m          = A.numRows();
  lno_t n          = A.numCols();
  auto rowmapHost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
  auto entriesHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);
  typename crsMat_t::index_type shuffledEntries("shuffled entries", A.nnz());
  // Randomly shuffle the entries within each row, so that the rows aren't
  // already sorted. Leave the values alone; this changes the matrix numerically
  // but this doesn't affect sorting.
  for (lno_t i = 0; i < m; i++) {
    std::random_shuffle(entriesHost.data() + i, entriesHost.data() + i + 1);
  }
  Kokkos::deep_copy(shuffledEntries, entriesHost);
  exec_space exec;
  Kokkos::Timer timer;
  double totalTime = 0;
  for (int rep = 0; rep < common_params.repeat; rep++) {
    Kokkos::deep_copy(exec, A.graph.entries, shuffledEntries);
    exec.fence();
    timer.reset();
    KokkosSparse::sort_crs_matrix(exec, A);
    exec.fence();
    totalTime += timer.seconds();
  }
  std::cout << "Mean sort_crs_matrix time over " << common_params.repeat << " trials: ";
  std::cout << totalTime / common_params.repeat << "\n";
}

#define KOKKOSKERNELS_PERF_TEST_NAME run_experiment
#include "KokkosKernels_perf_test_instantiation.hpp"
int main(int argc, char** argv) { return main_instantiation(argc, argv); }  // main
