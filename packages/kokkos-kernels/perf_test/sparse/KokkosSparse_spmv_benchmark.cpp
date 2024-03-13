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

// Headers needed to create initial data
// and to check results at the end
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_IOUtils.hpp>
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_TestUtils.hpp"
#include "KokkosKernels_perf_test_utilities.hpp"

// Headers for benchmark library
#include <benchmark/benchmark.h>
#include "Benchmark_Context.hpp"

// Headers for spmv
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosSparse_spmv.hpp>

namespace {

struct spmv_parameters {
  int N, offset;
  std::string filename;
  std::string alg;
  std::string tpl;

  spmv_parameters(const int N_)
      : N(N_), offset(0), filename(""), alg(""), tpl("") {}
};

void print_options() {
  std::cerr << "Options\n" << std::endl;

  std::cerr << perf_test::list_common_options();

  std::cerr
      << "\t[Optional] --repeat      :: how many times to repeat overall test"
      << std::endl;
  std::cerr << "  -n [N]          :: generate a semi-random banded (band size "
               "0.01xN)\n"
               "NxN matrix with average of 10 entries per row."
            << std::endl;
  std::cerr << "\t[Optional] --alg           :: the algorithm to run (default, "
               "native, merge)"
            << std::endl;
  std::cerr
      << "\t[Optional] --alg           :: the algorithm to run (classic, merge)"
      << std::endl;
  std::cerr << "\t[Optional] --TPL       :: when available and compatible with "
               "alg, a TPL can be used (cusparse, rocsparse, MKL)"
            << std::endl;
  std::cerr
      << "  -f [file]       : Read in Matrix Market formatted text file 'file'."
      << std::endl;
  std::cerr << "  --offset [O]    : Subtract O from every index.\n"
            << "                    Useful in case the matrix market file is "
               "not 0 based."
            << std::endl;
}  // print_options

void parse_inputs(int argc, char** argv, spmv_parameters& params) {
  for (int i = 1; i < argc; ++i) {
    if (perf_test::check_arg_int(i, argc, argv, "-n", params.N)) {
      ++i;
    } else if (perf_test::check_arg_str(i, argc, argv, "--alg", params.alg)) {
      if ((params.alg != "") && (params.alg != "default") &&
          (params.alg != "native") && (params.alg != "merge")) {
        throw std::runtime_error(
            "--alg can only be an empty string, `default`, `native` or "
            "`merge`!");
      }
      ++i;
    } else if (perf_test::check_arg_str(i, argc, argv, "--TPL", params.tpl)) {
      ++i;
    } else if (perf_test::check_arg_str(i, argc, argv, "-f", params.filename)) {
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--offset",
                                        params.offset)) {
      ++i;
    } else {
      print_options();
      KK_USER_REQUIRE_MSG(false, "Unrecognized command line argument #"
                                     << i << ": " << argv[i]);
    }
  }
}  // parse_inputs

template <class execution_space>
void run_spmv(benchmark::State& state, const spmv_parameters& inputs) {
  using matrix_type =
      KokkosSparse::CrsMatrix<double, int, execution_space, void, int>;
  using mv_type = Kokkos::View<double*, execution_space>;

  KokkosKernels::Experimental::Controls controls;
  if ((inputs.alg == "default") || (inputs.alg == "native") ||
      (inputs.alg == "merge")) {
    controls.setParameter("algorithm", inputs.alg);
  }

  // Create test matrix
  srand(17312837);
  matrix_type A;
  if (inputs.filename == "") {
    int nnz = 10 * inputs.N;
    A       = KokkosSparse::Impl::kk_generate_sparse_matrix<matrix_type>(
        inputs.N, inputs.N, nnz, 0, 0.01 * inputs.N);
  } else {
    A = KokkosSparse::Impl::read_kokkos_crst_matrix<matrix_type>(
        inputs.filename.c_str());
  }

  // Create input vectors
  mv_type x("X", A.numRows());
  mv_type y("Y", A.numCols());

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  Kokkos::fill_random(x, rand_pool, 10);
  Kokkos::fill_random(y, rand_pool, 10);

  // Run the actual experiments
  for (auto _ : state) {
    KokkosSparse::spmv(controls, KokkosSparse::NoTranspose, 1.0, A, x, 0.0, y);
    Kokkos::fence();
  }
}

}  // namespace

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);

  benchmark::Initialize(&argc, argv);
  benchmark::SetDefaultTimeUnit(benchmark::kMillisecond);
  KokkosKernelsBenchmark::add_benchmark_context(true);

  perf_test::CommonInputParams common_params;
  perf_test::parse_common_options(argc, argv, common_params);

  std::string bench_name = "KokkosSparse_spmv";

  // Set input parameters, default to random 100000x100000
  spmv_parameters inputs(100000);
  parse_inputs(argc, argv, inputs);

  // Google benchmark will report the wrong n if an input file matrix is used.
  KokkosKernelsBenchmark::register_benchmark_real_time(
      bench_name.c_str(), run_spmv<Kokkos::DefaultExecutionSpace>, {"n"},
      {inputs.N}, common_params.repeat, inputs);
  benchmark::RunSpecifiedBenchmarks();

  benchmark::Shutdown();
  Kokkos::finalize();

  return 0;
}
