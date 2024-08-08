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
#include <cmath>
#include <unordered_map>
#include <iomanip>  // std::setprecision

#include <Kokkos_Core.hpp>

#include "KokkosSparse_Utils.hpp"
#include "KokkosSparse_spiluk.hpp"
#include "KokkosSparse_par_ilut.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_default_types.hpp"
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_IOUtils.hpp>
#include "KokkosKernels_perf_test_utilities.hpp"

#include "Benchmark_Context.hpp"
#include <benchmark/benchmark.h>

#ifdef USE_GINKGO
#include <ginkgo/ginkgo.hpp>
#endif

namespace {

using KokkosSparse::Experimental::par_ilut_numeric;
using KokkosSparse::Experimental::par_ilut_symbolic;

using KokkosSparse::Experimental::spiluk_numeric;
using KokkosSparse::Experimental::spiluk_symbolic;
using KokkosSparse::Experimental::SPILUKAlgorithm;

// Build up useful types
using scalar_t  = default_scalar;
using lno_t     = default_lno_t;
using size_type = default_size_type;
using exe_space = Kokkos::DefaultExecutionSpace;
using mem_space = typename exe_space::memory_space;
using device    = Kokkos::Device<exe_space, mem_space>;

using RowMapType  = Kokkos::View<size_type*, device>;
using EntriesType = Kokkos::View<lno_t*, device>;
using ValuesType  = Kokkos::View<scalar_t*, device>;

using sp_matrix_type = KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type>;
using KernelHandle =
    KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, exe_space, mem_space, mem_space>;
using float_t = typename Kokkos::ArithTraits<scalar_t>::mag_type;

///////////////////////////////////////////////////////////////////////////////
void run_par_ilut_test(benchmark::State& state, KernelHandle& kh, const sp_matrix_type& A, int& num_iters)
///////////////////////////////////////////////////////////////////////////////
{
  const int rows = state.range(0);

  auto par_ilut_handle = kh.get_par_ilut_handle();

  // Pull out views from CRS
  auto A_row_map = A.graph.row_map;
  auto A_entries = A.graph.entries;
  auto A_values  = A.values;

  // Allocate L and U CRS views as outputs
  RowMapType L_row_map("L_row_map", rows + 1);
  RowMapType U_row_map("U_row_map", rows + 1);

  // Initial L/U approximations for A
  EntriesType L_entries("L_entries", 0);
  ValuesType L_values("L_values", 0);
  EntriesType U_entries("U_entries", 0);
  ValuesType U_values("U_values", 0);

  for (auto _ : state) {
    state.ResumeTiming();
    par_ilut_symbolic(&kh, A_row_map, A_entries, L_row_map, U_row_map);

    size_type nnzL = par_ilut_handle->get_nnzL();
    size_type nnzU = par_ilut_handle->get_nnzU();

    Kokkos::resize(L_entries, nnzL);
    Kokkos::resize(U_entries, nnzU);
    Kokkos::resize(L_values, nnzL);
    Kokkos::resize(U_values, nnzU);
    Kokkos::deep_copy(L_entries, 0);
    Kokkos::deep_copy(U_entries, 0);
    Kokkos::deep_copy(L_values, 0);
    Kokkos::deep_copy(U_values, 0);

    par_ilut_numeric(&kh, A_row_map, A_entries, A_values, L_row_map, L_entries, L_values, U_row_map, U_entries,
                     U_values);
    Kokkos::fence();
    state.PauseTiming();

    // Check worked
    num_iters = par_ilut_handle->get_num_iters();
    KK_REQUIRE_MSG(num_iters < par_ilut_handle->get_max_iter(), "par_ilut hit max iters");

    // Reset inputs
    Kokkos::deep_copy(L_row_map, 0);
    Kokkos::deep_copy(U_row_map, 0);
  }
}

#ifdef USE_GINKGO
///////////////////////////////////////////////////////////////////////////////
static constexpr bool IS_GPU = KokkosKernels::Impl::kk_is_gpu_exec_space<exe_space>();

using ginkgo_exec = std::conditional_t<IS_GPU, gko::CudaExecutor, gko::OmpExecutor>;

template <typename GinkgoT>
std::shared_ptr<GinkgoT> get_ginkgo_exec() {
  return GinkgoT::create();
}

#ifdef KOKKOS_ENABLE_CUDA
template <>
std::shared_ptr<gko::CudaExecutor> get_ginkgo_exec<gko::CudaExecutor>() {
  auto ref_exec = gko::ReferenceExecutor::create();
  return gko::CudaExecutor::create(0 /*device id*/, ref_exec);
}
#endif

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
void run_par_ilut_test_ginkgo(benchmark::State& state, KernelHandle& kh, const sp_matrix_type& A, const int& num_iters)
///////////////////////////////////////////////////////////////////////////////
{
  const int rows = state.range(0);

  auto par_ilut_handle = kh.get_par_ilut_handle();

  // Pull out views from CRS
  auto A_row_map = A.graph.row_map;
  auto A_entries = A.graph.entries;
  auto A_values  = A.values;

  using mtx = gko::matrix::Csr<scalar_t, lno_t>;

  auto exec = get_ginkgo_exec<ginkgo_exec>();

  // ginkgo does not differentiate between index type and size type. We need
  // to convert A_row_map to lno_t.
  EntriesType A_row_map_cp("A_row_map_cp", rows + 1);
  Kokkos::deep_copy(A_row_map_cp, A_row_map);

  // Populate mtx
  auto a_mtx_uniq = mtx::create_const(exec, gko::dim<2>(rows, rows),
                                      gko::array<scalar_t>::const_view(exec, A_values.extent(0), A_values.data()),
                                      gko::array<lno_t>::const_view(exec, A_entries.extent(0), A_entries.data()),
                                      gko::array<lno_t>::const_view(exec, A_row_map_cp.extent(0), A_row_map_cp.data()));

  std::shared_ptr<const mtx> a_mtx = std::move(a_mtx_uniq);

  for (auto _ : state) {
    auto fact = gko::factorization::ParIlut<scalar_t, lno_t>::build()
                    .with_fill_in_limit(par_ilut_handle->get_fill_in_limit())
                    .with_approximate_select(false)
                    .with_iterations(num_iters)
                    .on(exec)
                    ->generate(a_mtx);
  }
}
#endif

///////////////////////////////////////////////////////////////////////////////
void run_spiluk_test(benchmark::State& state, KernelHandle& kh, const sp_matrix_type& A, const int& team_size,
                     const bool measure_symbolic)
///////////////////////////////////////////////////////////////////////////////
{
  const int rows = state.range(0);

  constexpr int EXPAND_FACT  = 10;
  const lno_t fill_lev       = 2;
  const size_type handle_nnz = EXPAND_FACT * A.nnz() * (fill_lev + 1);
  kh.create_spiluk_handle(SPILUKAlgorithm::SEQLVLSCHD_TP1, rows, handle_nnz, handle_nnz);
  auto spiluk_handle = kh.get_spiluk_handle();
  spiluk_handle->set_team_size(team_size);

  // Pull out views from CRS
  auto A_row_map = A.graph.row_map;
  auto A_entries = A.graph.entries;
  auto A_values  = A.values;

  // Allocate L and U CRS views as outputs
  RowMapType L_row_map("L_row_map", rows + 1);
  RowMapType U_row_map("U_row_map", rows + 1);

  // Initial L/U approximations for A
  EntriesType L_entries("L_entries", handle_nnz);
  ValuesType L_values("L_values", handle_nnz);
  EntriesType U_entries("U_entries", handle_nnz);
  ValuesType U_values("U_values", handle_nnz);

  for (auto _ : state) {
    if (measure_symbolic) {
      state.ResumeTiming();
    }
    spiluk_symbolic(&kh, fill_lev, A_row_map, A_entries, L_row_map, L_entries, U_row_map, U_entries);
    Kokkos::fence();
    state.PauseTiming();

    const size_type nnzL = spiluk_handle->get_nnzL();
    const size_type nnzU = spiluk_handle->get_nnzU();

    Kokkos::resize(L_entries, nnzL);
    Kokkos::resize(U_entries, nnzU);
    Kokkos::resize(L_values, nnzL);
    Kokkos::resize(U_values, nnzU);

    if (!measure_symbolic) {
      state.ResumeTiming();
      spiluk_numeric(&kh, fill_lev, A_row_map, A_entries, A_values, L_row_map, L_entries, L_values, U_row_map,
                     U_entries, U_values);
      Kokkos::fence();
      state.PauseTiming();
    }

    // Reset inputs
    Kokkos::deep_copy(L_row_map, 0);
    Kokkos::deep_copy(U_row_map, 0);
    Kokkos::deep_copy(L_entries, 0);
    Kokkos::deep_copy(U_entries, 0);
    Kokkos::deep_copy(L_values, 0);
    Kokkos::deep_copy(U_values, 0);
    Kokkos::resize(L_entries, handle_nnz);
    Kokkos::resize(U_entries, handle_nnz);

    spiluk_handle->reset_handle(rows, handle_nnz, handle_nnz);
  }
}

///////////////////////////////////////////////////////////////////////////////
int test_par_ilut_perf(const std::string& matrix_file, int rows, int nnz_per_row, const int bandwidth, int team_size,
                       const int loop, const int test)
///////////////////////////////////////////////////////////////////////////////
{
  KernelHandle kh;
  kh.create_par_ilut_handle();

  // Generate or read A
  sp_matrix_type A;
  if (matrix_file == "") {
    size_type nnz                 = rows * nnz_per_row;
    const lno_t row_size_variance = 0;
    const scalar_t diag_dominance = 1;
    A                             = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<sp_matrix_type>(
        rows, rows, nnz, row_size_variance, bandwidth, diag_dominance);
  } else {
    A           = KokkosSparse::Impl::read_kokkos_crst_matrix<sp_matrix_type>(matrix_file.c_str());
    rows        = A.numRows();
    nnz_per_row = A.nnz() / rows;
  }

  // Now that we have A, we can set team_size
  if (team_size == -1) {
    team_size = KokkosKernels::Impl::kk_is_gpu_exec_space<exe_space>() ? nnz_per_row : 1;
  }

  KokkosSparse::sort_crs_matrix(A);

  // Make handles
  auto par_ilut_handle = kh.get_par_ilut_handle();
  par_ilut_handle->set_team_size(team_size);
  par_ilut_handle->set_nrows(rows);

  const auto default_policy = par_ilut_handle->get_default_team_policy();

  // Report test config to user
  if (matrix_file == "") {
    std::cout << "Testing par_ilut with rows=" << rows << "\n  nnz_per_row=" << nnz_per_row
              << "\n  bandwidth=" << bandwidth;
  } else {
    std::cout << "Testing par_ilut with input matrix=" << matrix_file;
  }
  std::cout << "\n  total nnz=" << A.nnz() << "\n  league_size=" << default_policy.league_size()
            << "\n  team_size=" << default_policy.team_size()
            << "\n  concurrent teams=" << exe_space().concurrency() / default_policy.team_size() << "\n  loop=" << loop
            << std::endl;

  std::string name     = "KokkosSparse_par_ilut";
  int num_iters        = 6;
  const auto arg_names = std::vector<std::string>{"rows"};
  const auto args      = std::vector<int64_t>{rows};

  if (test & 1) {
    auto plambda = [&](benchmark::State& state) { run_par_ilut_test(state, kh, A, num_iters); };
    KokkosKernelsBenchmark::register_benchmark_real_time((name + "_par_ilut").c_str(), plambda, arg_names, args, loop);
  }

#ifdef USE_GINKGO
  if (test & 2) {
    auto glambda = [&](benchmark::State& state) { run_par_ilut_test_ginkgo(state, kh, A, num_iters); };
    KokkosKernelsBenchmark::register_benchmark_real_time((name + "_gingko").c_str(), glambda, arg_names, args, loop);
  }
#endif

  if (test & 4) {
    auto s1lambda = [&](benchmark::State& state) { run_spiluk_test(state, kh, A, team_size, true); };
    auto s2lambda = [&](benchmark::State& state) { run_spiluk_test(state, kh, A, team_size, false); };
    KokkosKernelsBenchmark::register_benchmark_real_time((name + "_spiluk_symbolic").c_str(), s1lambda, arg_names, args,
                                                         loop);

    KokkosKernelsBenchmark::register_benchmark_real_time((name + "_spiluk_numeric").c_str(), s2lambda, arg_names, args,
                                                         loop);
  }

  // Need to run before vars used by lambdas go out of scope
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
void print_help_par_ilut()
///////////////////////////////////////////////////////////////////////////////
{
  printf("Options:\n");
  printf("  -f [F]  : Read in Matrix Market formatted text file.\n");
  printf("  -n [N]  : generate a semi-random banded NxN matrix.\n");
  printf("  -z [Z]  : number nnz per row. Default is min(1%% of N, 50).\n");
  printf("  -b [B]  : bandwidth per row. Default is max(2 * n^(1/2), nnz).\n");
  printf(
      "  -ts [T] : Number of threads per team. Default is 1 on OpenMP, "
      "nnz_per_row on CUDA\n");
  // printf("  -vl [V] : Vector-length (i.e. how many Cuda threads are a Kokkos
  // 'thread').\n");
  printf("  -l [L]  : How many runs to aggregate average time. Default is 4\n\n");
  printf(
      "  -t [T]  : Which tests to run. Bitwise. e.g. 7 => run all, 1 => "
      "par_ilut, 2 => ginkgo, 4 => spiluk,. Default is 7\n\n");
}

///////////////////////////////////////////////////////////////////////////////
void handle_int_arg(int argc, char** argv, int& i, std::map<std::string, int*> option_map)
///////////////////////////////////////////////////////////////////////////////
{
  std::string arg = argv[i];
  auto it         = option_map.find(arg);
  KK_USER_REQUIRE_MSG(it != option_map.end(), "Unknown option: " << arg);
  KK_USER_REQUIRE_MSG(i + 1 < argc, "Missing option value for option: " << arg);
  *(it->second) = atoi(argv[++i]);
}

}  // namespace

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
///////////////////////////////////////////////////////////////////////////////
{
  std::string mfile = "";
  int rows          = -1;
  int nnz_per_row   = -1;  // depends on other options, so don't set to default yet
  int bandwidth     = -1;
  int team_size     = -1;
  int test          = 7;

  std::map<std::string, int*> option_map = {
      {"-n", &rows}, {"-z", &nnz_per_row}, {"-b", &bandwidth}, {"-ts", &team_size}, {"-t", &test}};

  if (argc == 1) {
    print_help_par_ilut();
    return 0;
  }

  // Handle common params
  perf_test::CommonInputParams common_params;
  perf_test::parse_common_options(argc, argv, common_params);

  // Handle user options
  for (int i = 1; i < argc; i++) {
    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      print_help_par_ilut();
      return 0;
    } else if ((strcmp(argv[i], "-f") == 0)) {
      mfile = argv[++i];
    } else {
      handle_int_arg(argc, argv, i, option_map);
    }
  }

  // Determine where A is coming from
  if (rows != -1) {
    // We are randomly generating the input A
    KK_USER_REQUIRE_MSG(rows >= 100, "Need to have at least 100 rows");

    KK_USER_REQUIRE_MSG(mfile == "", "Need provide either -n or -f argument to this program, not both");
  } else {
    // We are reading A from a file
    KK_USER_REQUIRE_MSG(mfile != "", "Need provide either -n or -f argument to this program, not both");
  }

  // Set dependent defaults. Default team_size cannot be set
  // until we know more about A
  if (nnz_per_row == -1) {
    nnz_per_row = std::min(rows / 100, 50);
  }
  if (bandwidth == -1) {
    bandwidth = std::max(2 * (int)std::sqrt(rows), 2 * nnz_per_row);
  }

  Kokkos::initialize(argc, argv);
  {
    benchmark::Initialize(&argc, argv);
    benchmark::SetDefaultTimeUnit(benchmark::kSecond);
    KokkosKernelsBenchmark::add_benchmark_context(true);

    test_par_ilut_perf(mfile, rows, nnz_per_row, bandwidth, team_size, common_params.repeat, test);

    benchmark::Shutdown();
  }
  Kokkos::finalize();
  return 0;
}
