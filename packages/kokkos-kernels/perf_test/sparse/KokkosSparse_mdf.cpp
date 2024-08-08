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
#include "KokkosKernels_config.h"
#include "KokkosKernels_Handle.hpp"
#include "KokkosSparse_IOUtils.hpp"
#include "KokkosSparse_Utils_cusparse.hpp"
#include "KokkosKernels_TestUtils.hpp"
#include "KokkosKernels_perf_test_utilities.hpp"

#include "KokkosSparse_mdf.hpp"

using perf_test::CommonInputParams;

struct LocalParams {
  std::string amtx;
  int m         = 10000;
  int n         = 10000;
  int nnzPerRow = 30;
  bool diag     = false;  // Whether B should be diagonal only (requires A square)
  bool verbose  = false;
  int repeat    = 1;
};

template <class row_map_t, class entries_t>
struct diag_generator_functor {
  using size_type = typename row_map_t::non_const_value_type;

  row_map_t row_map;
  entries_t entries;

  diag_generator_functor(row_map_t row_map_, entries_t entries_) : row_map(row_map_), entries(entries_){};

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type rowIdx) const {
    row_map(rowIdx + 1) = rowIdx + 1;
    entries(rowIdx)     = rowIdx;
  }
};

void print_options() {
  std::cerr << "Options\n" << std::endl;

  std::cerr << perf_test::list_common_options();

  std::cerr << "\t[Optional] --amtx <path> :: input matrix" << std::endl;
  std::cerr << "\t[Optional] --repeat      :: how many times to repeat overall "
               "MDF"
            << std::endl;
  std::cerr << "\t[Optional] --verbose     :: enable verbose output" << std::endl;
  std::cerr << "\nSettings for randomly generated A matrix" << std::endl;
  std::cerr << "\t[Optional] --m           :: number of rows to generate" << std::endl;
  std::cerr << "\t[Optional] --n           :: number of cols to generate" << std::endl;
  std::cerr << "\t[Optional] --nnz         :: number of entries per row to generate" << std::endl;
  std::cerr << "\t[Optional] --diag        :: generate a diagonal matrix" << std::endl;
}  // print_options

int parse_inputs(LocalParams& params, int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    if (perf_test::check_arg_str(i, argc, argv, "--amtx", params.amtx)) {
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--m", params.m)) {
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--n", params.n)) {
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--nnz", params.nnzPerRow)) {
      ++i;
    } else if (perf_test::check_arg_bool(i, argc, argv, "--diag", params.diag)) {
    } else if (perf_test::check_arg_int(i, argc, argv, "--repeat", params.repeat)) {
      ++i;
    } else if (perf_test::check_arg_bool(i, argc, argv, "--verbose", params.verbose)) {
    } else {
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
      print_options();
      return 1;
    }
  }
  return 0;
}  // parse_inputs

template <typename execution_space>
void run_experiment(int argc, char** argv, CommonInputParams /*params*/) {
  using crsMat_t   = KokkosSparse::CrsMatrix<double, int, execution_space, void, int>;
  using size_type  = typename crsMat_t::size_type;
  using lno_t      = typename crsMat_t::ordinal_type;
  using scalar_t   = typename crsMat_t::value_type;
  using device_t   = typename crsMat_t::device_type;
  using exec_space = typename device_t::execution_space;

  using graph_t   = typename crsMat_t::StaticCrsGraphType;
  using rowmap_t  = typename graph_t::row_map_type::non_const_type;
  using entries_t = typename graph_t::entries_type::non_const_type;
  using values_t  = typename crsMat_t::values_type::non_const_type;

  LocalParams localParams;
  parse_inputs(localParams, argc, argv);

  std::cout << "************************************* \n";
  std::cout << "************************************* \n";
  crsMat_t A;
  lno_t m = localParams.m;
  lno_t n = localParams.n;
  if (localParams.amtx.length()) {
    std::cout << "Loading A from " << localParams.amtx << '\n';
    A = KokkosSparse::Impl::read_kokkos_crst_matrix<crsMat_t>(localParams.amtx.c_str());
    m = A.numRows();
    n = A.numCols();
  } else {
    if (localParams.diag) {
      std::cout << "Randomly generating diag matrix\n";
      rowmap_t rowmapA("A row map", m + 1);
      entries_t entriesA("A entries", m);
      values_t valuesA("A values", m);

      // Generate the graph of A
      diag_generator_functor diag_generator(rowmapA, entriesA);
      Kokkos::parallel_for(Kokkos::RangePolicy<size_type, exec_space>(0, m), diag_generator);

      // Generate the values of A
      Kokkos::Random_XorShift64_Pool<exec_space> rand_pool(13718);
      Kokkos::fill_random(valuesA, rand_pool, 10 * Kokkos::ArithTraits<scalar_t>::one());

      // Actually put A together
      graph_t graph(entriesA, rowmapA);
      A = crsMat_t("A matrix", m, valuesA, graph);
    } else {
      std::cout << "Randomly generating matrix\n";
      size_type nnzUnused = m * localParams.nnzPerRow;
      A                   = KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(m, n, nnzUnused, 0, (n + 3) / 3);
    }
  }

  if (localParams.verbose) {
    std::cout << "Matrix A" << std::endl;
    std::cout << "  row_map A:" << std::endl;
    KokkosKernels::Impl::print_1Dview(A.graph.row_map);
    std::cout << "  entries A:" << std::endl;
    KokkosKernels::Impl::print_1Dview(A.graph.entries);
    std::cout << "  values A:" << std::endl;
    KokkosKernels::Impl::print_1Dview(A.values);
    std::cout << std::endl;
  }

  Kokkos::Timer timer;
  double handleTime   = 0;
  double symbolicTime = 0;
  double numericTime  = 0;

  timer.reset();
  KokkosSparse::Experimental::MDF_handle<crsMat_t> handle(A);
  handle.set_verbosity(0);
  if (localParams.verbose) {
    handle.set_verbosity(1);
  }
  handleTime += timer.seconds();

  for (int sumRep = 0; sumRep < localParams.repeat; sumRep++) {
    timer.reset();
    KokkosSparse::Experimental::mdf_symbolic(A, handle);
    Kokkos::fence();
    symbolicTime += timer.seconds();

    timer.reset();
    KokkosSparse::Experimental::mdf_numeric(A, handle);
    Kokkos::fence();
    numericTime += timer.seconds();
  }

  std::cout << "Mean total time:    "
            << handleTime + (symbolicTime / localParams.repeat) + (numericTime / localParams.repeat) << std::endl
            << "Handle time: " << handleTime << std::endl
            << "Mean symbolic time: " << (symbolicTime / localParams.repeat) << std::endl
            << "Mean numeric time:  " << (numericTime / localParams.repeat) << std::endl;

  if (localParams.verbose) {
    entries_t permutation = handle.get_permutation();

    std::cout << "MDF permutation:" << std::endl;
    KokkosKernels::Impl::print_1Dview(permutation);
  }
}  // run_experiment

#define KOKKOSKERNELS_PERF_TEST_NAME run_experiment
#include "KokkosKernels_perf_test_instantiation.hpp"
int main(int argc, char** argv) { return main_instantiation(argc, argv); }  // main
