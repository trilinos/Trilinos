// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <cstdio>
#include <chrono>
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
#include "KokkosSparse_gmres.hpp"
#include "KokkosSparse_LUPrec.hpp"
#include "KokkosSparse_par_ilut.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_default_types.hpp"
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_IOUtils.hpp>
#include "KokkosKernels_benchmark_utilities.hpp"

#include "Benchmark_Context.hpp"
#include <benchmark/benchmark.h>

#ifdef USE_GINKGO
#include <ginkgo/ginkgo.hpp>
#endif

namespace {

using KokkosSparse::Experimental::gmres;
using KokkosSparse::Experimental::par_ilut_numeric;
using KokkosSparse::Experimental::par_ilut_symbolic;

using KokkosSparse::spiluk_numeric;
using KokkosSparse::spiluk_symbolic;
using KokkosSparse::Experimental::SPILUKAlgorithm;

// Build up useful types
using scalar_t    = KokkosKernels::default_scalar;
using lno_t       = KokkosKernels::default_lno_t;
using size_type_t = lno_t;  // ginkgo does not have the concept of size_type separate from ordinal type
using exe_space_t = Kokkos::DefaultExecutionSpace;
using rpolicy_t   = Kokkos::RangePolicy<exe_space_t>;
using mem_space_t = typename exe_space_t::memory_space;
using device_t    = Kokkos::Device<exe_space_t, mem_space_t>;
using rows_t      = Kokkos::View<size_type_t*, device_t>;
using entries_t   = Kokkos::View<lno_t*, device_t>;
using values_t    = Kokkos::View<scalar_t*, device_t>;
using result_t    = std::tuple<rows_t, entries_t, values_t, rows_t, entries_t, values_t>;  // full result
using lu_result_t = std::tuple<rows_t, entries_t, values_t>;                               // result for one of L or U
using sp_matrix_t = KokkosSparse::CrsMatrix<scalar_t, lno_t, device_t, void, size_type_t>;
using float_t     = typename KokkosKernels::ArithTraits<scalar_t>::mag_type;
using kkhandle_t  = KokkosKernels::Experimental::KokkosKernelsHandle<size_type_t, lno_t, scalar_t, exe_space_t,
                                                                    mem_space_t, mem_space_t>;

///////////////////////////////////////////////////////////////////////////////
auto try_lu_prec(kkhandle_t& kh, const sp_matrix_t& A, const result_t& results, const int max_subspace,
                 const bool do_baseline = false)
///////////////////////////////////////////////////////////////////////////////
{
  // Unpack
  auto [l_row_map, l_entries, l_values, u_row_map, u_entries, u_values] = results;

  const int rows     = l_row_map.size() - 1;
  const auto m       = max_subspace;
  constexpr auto tol = 1e-8;

  // Create gmres handle
  kh.create_gmres_handle(m, tol);
  auto gmres_handle = kh.get_gmres_handle();
  gmres_handle->set_verbose(false);
  using GMRESHandle    = typename std::remove_reference<decltype(*gmres_handle)>::type;
  using ViewVectorType = typename GMRESHandle::nnz_value_view_t;

  // Create CRSs
  sp_matrix_t L("L", rows, rows, l_values.extent(0), l_values, l_row_map, l_entries),
      U("U", rows, rows, u_values.extent(0), u_values, u_row_map, u_entries);

  // Set initial vectors:
  ViewVectorType X("X", rows);    // Solution and initial guess
  ViewVectorType Wj("Wj", rows);  // For checking residuals at end.
  ViewVectorType B(Kokkos::view_alloc(Kokkos::WithoutInitializing, "B"),
                   rows);  // right-hand side vec
  // Make rhs ones so that results are repeatable:
  Kokkos::deep_copy(B, 1.0);

  int num_iters_plain(9999), num_iters_precond(0);

  // Solve Ax = b
  if (do_baseline) {
    auto start = std::chrono::steady_clock::now();
    gmres(&kh, A, B, X);

    // Double check residuals at end of solve:
    float_t nrmB = KokkosBlas::nrm2(B);
    KokkosSparse::spmv("N", 1.0, A, X, 0.0, Wj);  // wj = Ax
    KokkosBlas::axpy(-1.0, Wj, B);                // b = b-Ax.
    float_t endRes = KokkosBlas::nrm2(B) / nrmB;

    const auto conv_flag = gmres_handle->get_conv_flag_val();
    num_iters_plain      = gmres_handle->get_num_iters();

    KK_REQUIRE(num_iters_plain > 0);
    if (conv_flag == GMRESHandle::Flag::Conv || endRes < gmres_handle->get_tol()) {
      std::cout << "WARNING: baseline did not converge" << std::endl;
    }
    auto end                                       = std::chrono::steady_clock::now();
    std::chrono::duration<double> seconds_duration = end - start;
    std::cout << "Baseline gmres took " << seconds_duration.count() << " seconds" << std::endl;
  }

  // Solve Ax = b with LU preconditioner.
  {
    auto start = std::chrono::steady_clock::now();
    gmres_handle->reset_handle(m, tol);
    gmres_handle->set_verbose(false);

    // Make precond
    KokkosSparse::Experimental::LUPrec<sp_matrix_t, kkhandle_t> myPrec(L, U);

    // reset X for next gmres call
    Kokkos::deep_copy(X, 0.0);

    gmres(&kh, A, B, X, &myPrec);

    // Double check residuals at end of solve:
    float_t nrmB = KokkosBlas::nrm2(B);
    KokkosSparse::spmv("N", 1.0, A, X, 0.0, Wj);  // wj = Ax
    KokkosBlas::axpy(-1.0, Wj, B);                // b = b-Ax.
    float_t endRes = KokkosBlas::nrm2(B) / nrmB;

    const auto conv_flag = gmres_handle->get_conv_flag_val();
    num_iters_precond    = gmres_handle->get_num_iters();

    KK_REQUIRE(endRes < gmres_handle->get_tol());
    KK_REQUIRE(conv_flag == GMRESHandle::Flag::Conv);
    KK_REQUIRE(num_iters_precond < num_iters_plain);

    auto end                                       = std::chrono::steady_clock::now();
    std::chrono::duration<double> seconds_duration = end - start;
    std::cout << "LUPrec gmres took " << seconds_duration.count() << " seconds with max subspace " << max_subspace
              << std::endl;
  }

  return std::tuple<int, int>{num_iters_plain, num_iters_precond};
}

///////////////////////////////////////////////////////////////////////////////
template <typename GkoMtxT>
lu_result_t copy_to_kokkos_and_dealloc(GkoMtxT& factor, const int nrows)
///////////////////////////////////////////////////////////////////////////////
{
  lno_t* rows      = const_cast<lno_t*>(factor->get_const_row_ptrs());
  lno_t* entries   = const_cast<lno_t*>(factor->get_const_col_idxs());
  scalar_t* values = const_cast<scalar_t*>(factor->get_const_values());

  const auto nnz = factor->get_num_stored_elements();

  // Allocate kokkos views
  rows_t vrow_map("row_map", nrows + 1);
  entries_t ventries("entries", nnz);
  values_t vvalues("values", nnz);

  Kokkos::parallel_for(
      "cp rows", rpolicy_t(0, nrows + 1), KOKKOS_LAMBDA(const int row) { vrow_map(row) = rows[row]; });

  Kokkos::parallel_for(
      "cp entries", rpolicy_t(0, nnz), KOKKOS_LAMBDA(const int i) {
        ventries(i) = entries[i];
        vvalues(i)  = values[i];
      });

  // Deallocate factor to save memory
  factor = GkoMtxT();

  return lu_result_t{vrow_map, ventries, vvalues};
}

///////////////////////////////////////////////////////////////////////////////
void run_par_ilut_test(benchmark::State& state, kkhandle_t& kh, const sp_matrix_t& A, int& num_iters, result_t& result,
                       const bool validate, const int gmres_max_subspace, bool compare = false)
///////////////////////////////////////////////////////////////////////////////
{
  const int rows = state.range(0);

  auto par_ilut_handle = kh.get_par_ilut_handle();

  // Pull out views from CRS
  auto A_row_map = A.graph.row_map;
  auto A_entries = A.graph.entries;
  auto A_values  = A.values;

  // Allocate L and U CRS views as outputs
  rows_t L_row_map("L_row_map", rows + 1);
  rows_t U_row_map("U_row_map", rows + 1);

  // Initial L/U approximations for A
  entries_t L_entries("L_entries", 0);
  values_t L_values("L_values", 0);
  entries_t U_entries("U_entries", 0);
  values_t U_values("U_values", 0);

  for (auto _ : state) {
    // Reset inputs
    Kokkos::deep_copy(L_row_map, 0);
    Kokkos::deep_copy(U_row_map, 0);

    state.ResumeTiming();
    par_ilut_symbolic(&kh, A_row_map, A_entries, L_row_map, U_row_map);

    size_type_t nnzL = par_ilut_handle->get_nnzL();
    size_type_t nnzU = par_ilut_handle->get_nnzU();

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
    std::cout << "PAR_ILUT finished in " << num_iters << " iters" << std::endl;
    KK_REQUIRE_MSG(num_iters <= par_ilut_handle->get_max_iter(), "par_ilut hit max iters");
  }

  result = {L_row_map, L_entries, L_values, U_row_map, U_entries, U_values};

  if (validate) {
    auto [plain, pr] = try_lu_prec(kh, A, result, gmres_max_subspace, true /*do a baseline (no prec) run*/);
    std::cout << "LUPrec results: " << std::endl;
    std::cout << "  num iters no prec: " << plain << std::endl;
    std::cout << "  num iters par_ilut: " << pr << std::endl;
  }

  // If we want to compare results, we store the results. This can dramatically increase the memory
  // requirement.
  if (!compare) {
    result = result_t();
  }
}

#ifdef USE_GINKGO
///////////////////////////////////////////////////////////////////////////////
static constexpr bool IS_GPU = KokkosKernels::Impl::is_gpu_exec_space_v<exe_space_t>;

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
void run_par_ilut_test_ginkgo(benchmark::State& state, kkhandle_t& kh, sp_matrix_t& A, const int& num_iters,
                              result_t& result, const bool validate, const int gmres_max_subspace, bool compare = false)
///////////////////////////////////////////////////////////////////////////////
{
  using gko_par_ilut_t   = typename gko::factorization::ParIlut<scalar_t, lno_t>;
  using gko_matrix_t     = typename gko_par_ilut_t::matrix_type;
  using gko_matrix_ptr_t = typename std::shared_ptr<const gko_matrix_t>;

  const int rows = state.range(0);

  auto par_ilut_handle = kh.get_par_ilut_handle();

  std::cout << "Running ginkgo par_ilut with " << num_iters << " iters" << std::endl;

  // Pull out views from CRS
  auto A_row_map = A.graph.row_map;
  auto A_entries = A.graph.entries;
  auto A_values  = A.values;

  using mtx = gko::matrix::Csr<scalar_t, lno_t>;

  auto exec = get_ginkgo_exec<ginkgo_exec>();

  gko_matrix_ptr_t l_factor, u_factor;

  {
    // Populate mtx
    auto a_mtx_uniq = gko_matrix_t::create_const(
        exec, gko::dim<2>(rows, rows), gko::array<scalar_t>::const_view(exec, A_values.extent(0), A_values.data()),
        gko::array<lno_t>::const_view(exec, A_entries.extent(0), A_entries.data()),
        gko::array<lno_t>::const_view(exec, A_row_map.extent(0), A_row_map.data()));

    gko_matrix_ptr_t a_mtx = std::move(a_mtx_uniq);

    for (auto _ : state) {
      state.ResumeTiming();
      auto fact = gko::factorization::ParIlut<scalar_t, lno_t>::build()
                      .with_fill_in_limit(par_ilut_handle->get_fill_in_limit())
                      .with_approximate_select(false)
                      .with_iterations(num_iters)
                      .on(exec)
                      ->generate(a_mtx);
      state.PauseTiming();

      // Store for later use
      l_factor = fact->get_l_factor();
      u_factor = fact->get_u_factor();
    }
  }  // This should deallocate all the gko stuff except for l and u factors

  // As we copy data over to kokkos, we need to be very careful not to use too much memory
  auto [L_row_map, L_entries, L_values] = copy_to_kokkos_and_dealloc(l_factor, rows);
  auto [U_row_map, U_entries, U_values] = copy_to_kokkos_and_dealloc(u_factor, rows);

  result = {L_row_map, L_entries, L_values, U_row_map, U_entries, U_values};

  if (validate) {
    auto [plain, pr] = try_lu_prec(kh, A, result, gmres_max_subspace);
    std::cout << "LUPrec results: " << std::endl;
    std::cout << "  num iters ginkgo: " << pr << std::endl;
  }

  // If we want to compare results, we store the results. This can dramatically increase the memory
  // requirement.
  if (!compare) {
    result = result_t();
  }
}
#endif

///////////////////////////////////////////////////////////////////////////////
void run_spiluk_test(benchmark::State& state, kkhandle_t& kh, const sp_matrix_t& A, const int& team_size,
                     const bool measure_symbolic, result_t& result, const bool validate, const int gmres_max_subspace,
                     bool compare = false)
///////////////////////////////////////////////////////////////////////////////
{
  const int rows = state.range(0);

  constexpr int EXPAND_FACT    = 4;  // Be careful with this. Too high of a value can make you run out of mem
  const lno_t fill_lev         = 2;
  const size_type_t handle_nnz = EXPAND_FACT * A.nnz() * (fill_lev + 1);
  kh.create_spiluk_handle(SPILUKAlgorithm::SEQLVLSCHD_TP1, rows, handle_nnz, handle_nnz);
  auto spiluk_handle = kh.get_spiluk_handle();
  spiluk_handle->set_team_size(team_size);

  // Pull out views from CRS
  auto A_row_map = A.graph.row_map;
  auto A_entries = A.graph.entries;
  auto A_values  = A.values;

  // Allocate L and U CRS views as outputs
  rows_t L_row_map("L_row_map", rows + 1);
  rows_t U_row_map("U_row_map", rows + 1);

  // Initial L/U approximations for A
  entries_t L_entries("L_entries", handle_nnz);
  values_t L_values("L_values", handle_nnz);
  entries_t U_entries("U_entries", handle_nnz);
  values_t U_values("U_values", handle_nnz);

  for (auto _ : state) {
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

    if (measure_symbolic) {
      state.ResumeTiming();
    }
    KokkosSparse::spiluk_symbolic(&kh, fill_lev, A_row_map, A_entries, L_row_map, L_entries, U_row_map, U_entries);
    Kokkos::fence();
    state.PauseTiming();

    const size_type_t nnzL = spiluk_handle->get_nnzL();
    const size_type_t nnzU = spiluk_handle->get_nnzU();

    Kokkos::resize(L_entries, nnzL);
    Kokkos::resize(U_entries, nnzU);
    Kokkos::resize(L_values, nnzL);
    Kokkos::resize(U_values, nnzU);

    if (!measure_symbolic) {
      state.ResumeTiming();
      KokkosSparse::spiluk_numeric(&kh, fill_lev, A_row_map, A_entries, A_values, L_row_map, L_entries, L_values,
                                   U_row_map, U_entries, U_values);
      Kokkos::fence();
      state.PauseTiming();
    }
  }
  result = {L_row_map, L_entries, L_values, U_row_map, U_entries, U_values};

  if (validate && !measure_symbolic) {
    auto [plain, pr] = try_lu_prec(kh, A, result, gmres_max_subspace);
    std::cout << "LUPrec results: " << std::endl;
    std::cout << "  num iters spiluk: " << pr << std::endl;
  }

  // If we want to compare results, we store the results. This can dramatically increase the memory
  // requirement.
  if (!compare) {
    result = result_t();
  }
}

///////////////////////////////////////////////////////////////////////////////
int run_ilu_perf_tests(const std::string& matrix_file, int rows, int nnz_per_row, const int bandwidth, int team_size,
                       const int loop, const int test, const bool validate, const int gmres_max_subspace,
                       const int max_iter, const double fill_in_limit, const double residual_norm_delta_stop)
///////////////////////////////////////////////////////////////////////////////
{
  kkhandle_t kh;
  kh.create_par_ilut_handle(max_iter, residual_norm_delta_stop, fill_in_limit);

  std::cout << "Running par_ilut with max_iter=" << max_iter
            << ", residual_norm_delta_stop=" << residual_norm_delta_stop << ", fill_in_limit=" << fill_in_limit
            << std::endl;

  // Generate or read A
  auto start = std::chrono::steady_clock::now();
  sp_matrix_t A;
  if (matrix_file == "") {
    size_type_t nnz               = rows * nnz_per_row;
    const lno_t row_size_variance = 0;
    const scalar_t diag_dominance = 1;
    A                             = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<sp_matrix_t>(
        rows, rows, nnz, row_size_variance, bandwidth, diag_dominance);
  } else {
    A           = KokkosSparse::Impl::read_kokkos_crst_matrix<sp_matrix_t>(matrix_file.c_str());
    rows        = A.numRows();
    nnz_per_row = A.nnz() / rows;
  }
  KokkosSparse::sort_crs_matrix(A);

  auto end                                       = std::chrono::steady_clock::now();
  std::chrono::duration<double> seconds_duration = end - start;

  // Now that we have A, we can set team_size
  if (team_size == -1) {
    team_size = KokkosKernels::Impl::is_gpu_exec_space_v<exe_space_t> ? nnz_per_row : 1;
  }

  // Make handles
  auto par_ilut_handle = kh.get_par_ilut_handle();
  par_ilut_handle->set_team_size(team_size);
  par_ilut_handle->set_nrows(rows);

  const auto default_policy = par_ilut_handle->get_default_team_policy();

  // Report test config to user
  std::cout << "A generated/read in " << seconds_duration.count() << " seconds" << std::endl;
  if (matrix_file == "") {
    std::cout << "Testing ILU with rows=" << rows << "\n  nnz_per_row=" << nnz_per_row << "\n  bandwidth=" << bandwidth;
  } else {
    std::cout << "Testing ILU with input matrix=" << matrix_file;
  }
  std::cout << "\n  total nnz=" << A.nnz() << "\n  league_size=" << default_policy.league_size()
            << "\n  team_size=" << default_policy.team_size()
            << "\n  concurrent teams=" << exe_space_t().concurrency() / default_policy.team_size()
            << "\n  loop=" << loop << std::endl;

  std::string name     = "KokkosSparse_par_ilut";
  int num_iters        = max_iter;
  const auto arg_names = std::vector<std::string>{"rows"};
  const auto args      = std::vector<int64_t>{rows};

  result_t parilut_results, ginkgo_results, spiluk_results;

  if (test & 1) {
    auto plambda = [&](benchmark::State& state) {
      run_par_ilut_test(state, kh, A, num_iters, parilut_results, validate, gmres_max_subspace);
    };
    KokkosKernelsBenchmark::register_benchmark_real_time((name + "_par_ilut").c_str(), plambda, arg_names, args, loop);
  }

  if (test & 2) {
#ifdef USE_GINKGO
    auto glambda = [&](benchmark::State& state) {
      run_par_ilut_test_ginkgo(state, kh, A, num_iters, ginkgo_results, validate, gmres_max_subspace);
    };
    KokkosKernelsBenchmark::register_benchmark_real_time((name + "_gingko").c_str(), glambda, arg_names, args, loop);
#else
    KK_REQUIRE_MSG(false, "Requested ginkgo perf test but it is not enabled");
#endif
  }

  if (test & 4) {
    auto s1lambda = [&](benchmark::State& state) {
      run_spiluk_test(state, kh, A, team_size, true, spiluk_results, false, gmres_max_subspace);
    };
    auto s2lambda = [&](benchmark::State& state) {
      run_spiluk_test(state, kh, A, team_size, false, spiluk_results, validate, gmres_max_subspace);
    };
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
  printf("  -v      : Validate ILU results by doing luprec+gmres.\n");
  printf("  -m [M]  : For GMRES, set the max subspace size.\n");
  printf("  -i [I]  : For par_ilut, set max iters.\n");
  printf("  -l [L]  : For par_ilut, set fill in limit.\n");
  printf("  -r [R]  : For par_ilut, set residual norm delta stop. \n");
  printf("  -C      : Do not print context.\n");
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
void handle_arg(int argc, char** argv, int& i, std::map<std::string, int*> option_map,
                std::map<std::string, double*> foption_map)
///////////////////////////////////////////////////////////////////////////////
{
  std::string arg = argv[i];
  auto it         = option_map.find(arg);
  auto fit        = foption_map.find(arg);
  KK_USER_REQUIRE_MSG(it != option_map.end() || fit != foption_map.end(), "Unknown option: " << arg);
  KK_USER_REQUIRE_MSG(i + 1 < argc, "Missing option value for option: " << arg);
  if (it != option_map.end()) {
    *(it->second) = atoi(argv[++i]);
  } else {
    *(fit->second) = atof(argv[++i]);
  }
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
#ifdef USE_GINKGO
  int test = 7;  // Default to all if ginkgo is enabled
#else
  int test = 5;  // Default to par_ilut and spiluk
#endif
  bool validate                   = false;
  bool no_context                 = false;
  int gmres_max_subspace          = 50;
  int max_iter                    = 20;
  double fill_in_limit            = 0.75;
  double residual_norm_delta_stop = 1e-2;

  std::map<std::string, int*> option_map = {{"-n", &rows},       {"-z", &nnz_per_row}, {"-b", &bandwidth},
                                            {"-ts", &team_size}, {"-t", &test},        {"-m", &gmres_max_subspace},
                                            {"-i", &max_iter}};

  std::map<std::string, double*> foption_map = {{"-l", &fill_in_limit}, {"-r", &residual_norm_delta_stop}};

  if (argc == 1) {
    print_help_par_ilut();
    return 0;
  }

  // Handle common params
  benchmark::CommonInputParams common_params;
  benchmark::parse_common_options(argc, argv, common_params);

  // Handle user options
  for (int i = 1; i < argc; i++) {
    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      print_help_par_ilut();
      return 0;
    } else if ((strcmp(argv[i], "-f") == 0)) {
      mfile = argv[++i];
    } else if ((strcmp(argv[i], "-v") == 0)) {
      validate = true;
    } else if ((strcmp(argv[i], "-C") == 0)) {
      no_context = true;
    } else {
      handle_arg(argc, argv, i, option_map, foption_map);
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
    // bandwidth = std::max(2 * (int)std::sqrt(rows), 2 * nnz_per_row) / 2;
    bandwidth = nnz_per_row * 8;
  }

  Kokkos::initialize(argc, argv);
  {
    benchmark::Initialize(&argc, argv);
    benchmark::SetDefaultTimeUnit(benchmark::kSecond);
    if (!no_context) {
      KokkosKernelsBenchmark::add_benchmark_context(true);
    }

    run_ilu_perf_tests(mfile, rows, nnz_per_row, bandwidth, team_size, common_params.repeat, test, validate,
                       gmres_max_subspace, max_iter, fill_in_limit, residual_norm_delta_stop);

    benchmark::Shutdown();
  }
  Kokkos::finalize();
  return 0;
}
