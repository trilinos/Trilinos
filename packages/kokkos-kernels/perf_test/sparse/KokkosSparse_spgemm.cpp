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
#include "KokkosKernels_default_types.hpp"
#include "KokkosSparse_IOUtils.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spgemm.hpp"
#include "KokkosSparse_SortCrs.hpp"
#include "KokkosBlas1_nrminf.hpp"
#include "KokkosBlas1_axpby.hpp"
#include "KokkosKernels_TestParameters.hpp"
#include "KokkosKernels_TestUtils.hpp"
#include "KokkosKernels_perf_test_utilities.hpp"

#define TRANSPOSEFIRST false
#define TRANSPOSESECOND false

template <typename crsMat_t, typename device>
bool is_same_matrix(crsMat_t output_mat_actual, crsMat_t output_mat_reference) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  size_t nrows1    = output_mat_actual.graph.row_map.extent(0);
  size_t ncols1    = output_mat_actual.graph.row_map.extent(0);
  size_t nentries1 = output_mat_actual.graph.entries.extent(0);
  size_t nvals1    = output_mat_actual.values.extent(0);

  size_t nrows2    = output_mat_reference.graph.row_map.extent(0);
  size_t ncols2    = output_mat_reference.graph.row_map.extent(0);
  size_t nentries2 = output_mat_reference.graph.entries.extent(0);
  size_t nvals2    = output_mat_reference.values.extent(0);

  if (nrows1 != nrows2 || ncols1 != ncols2) {
    std::cerr << "Wrong dimensions: is " << nrows1 << 'x' << ncols1 << " but should be " << nrows2 << 'x' << ncols2
              << '\n';
    return false;
  }
  if (nentries1 != nentries2) {
    std::cerr << "Wrong number of entries: " << nentries1 << ", but should have " << nentries2 << '\n';
    return false;
  }
  if (nvals1 != nvals2) {
    std::cerr << "Wrong number of values: " << nvals1 << ", but should have " << nvals2 << '\n';
    return false;
  }

  bool is_identical = true;
  is_identical =
      KokkosKernels::Impl::kk_is_identical_view<typename graph_t::row_map_type, typename graph_t::row_map_type,
                                                typename lno_view_t::value_type, typename device::execution_space>(
          output_mat_actual.graph.row_map, output_mat_reference.graph.row_map, 0);
  if (!is_identical) {
    std::cerr << "Wrong rowmap:\n";
    KokkosKernels::Impl::print_1Dview(std::cerr, output_mat_actual.graph.row_map);
    std::cerr << "but should be:\n";
    KokkosKernels::Impl::print_1Dview(std::cerr, output_mat_reference.graph.row_map);
    return false;
  }

  is_identical =
      KokkosKernels::Impl::kk_is_identical_view<lno_nnz_view_t, lno_nnz_view_t, typename lno_nnz_view_t::value_type,
                                                typename device::execution_space>(
          output_mat_actual.graph.entries, output_mat_reference.graph.entries, 0);
  if (!is_identical) {
    for (size_t i = 0; i < nrows1; ++i) {
      size_t rb      = output_mat_actual.graph.row_map(i);
      size_t re      = output_mat_actual.graph.row_map(i + 1);
      bool incorrect = false;
      for (size_t j = rb; j < re; ++j) {
        if (output_mat_actual.graph.entries(j) != output_mat_reference.graph.entries(j)) {
          incorrect = true;
          break;
        }
      }
      if (incorrect) {
        for (size_t j = rb; j < re; ++j) {
          std::cerr << "row:" << i << " j:" << j << " h_ent1(j):" << output_mat_actual.graph.entries(j)
                    << " h_ent2(j):" << output_mat_reference.graph.entries(j) << " rb:" << rb << " re:" << re
                    << std::endl;
        }
      }
    }
    std::cerr << "Wrong entries, see above." << std::endl;
    return false;
  }

  scalar_view_t valueDiff(Kokkos::view_alloc(Kokkos::WithoutInitializing, "spgemm values diff"),
                          output_mat_actual.values.extent(0));
  Kokkos::deep_copy(valueDiff, output_mat_actual.values);
  KokkosBlas::axpy(-1.0, output_mat_reference.values, valueDiff);
  auto maxDiff = KokkosBlas::nrminf(valueDiff);

  std::cout << "Absolute maximum difference between actual and reference C values: " << maxDiff << '\n';

  return true;
}

void print_options() {
  std::cerr << "Options\n" << std::endl;

  std::cerr << perf_test::list_common_options();

  std::cerr << "\t[Required] INPUT MATRIX: '--amtx [left_hand_side.mtx]' -- for C=AxA" << std::endl;

  std::cerr << "\t[Optional] '--algorithm "
               "[DEFAULT=KKDEFAULT=KKSPGEMM|KKMEM|KKDENSE]' --> to choose algorithm. "
               "KKMEM is outdated, use KKSPGEMM instead."
            << std::endl;
  std::cerr << "\t[Optional] --bmtx [righ_hand_side.mtx]' for C = AxB" << std::endl;
  std::cerr << "\t[Optional] OUTPUT MATRICES: '--cmtx [output_matrix.mtx]' --> "
               "to write output C=AxB"
            << std::endl;
  std::cerr << "\t[Optional] --DENSEACCMAX: on CPUs default algorithm may "
               "choose to use dense accumulators. This parameter defaults to "
               "250k, which is max k value to choose dense accumulators. This "
               "can be increased with more memory bandwidth."
            << std::endl;
  std::cerr << "\t[Optional] '--dynamic': Use this for dynamic "
               "loop scheduling. (Better performance most of the time)"
            << std::endl;
  std::cerr << "\t[Optional] '--verbose': detailed output about SpGEMM and the "
               "output matrix"
            << std::endl;
  std::cerr << "\t[Optional] '--checkoutput': verify result against serial "
               "reference implementation"
            << std::endl;
}

int parse_inputs(KokkosKernels::Experiment::Parameters& params, int argc, char** argv) {
  std::string algoStr;
  for (int i = 1; i < argc; ++i) {
    if (perf_test::check_arg_int(i, argc, argv, "--repeat", params.repeat)) {
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--hashscale", params.minhashscale)) {
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--chunksize", params.chunk_size)) {
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--teamsize", params.team_size)) {
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--vectorsize", params.vector_size)) {
      ++i;
    } else if (perf_test::check_arg_bool(i, argc, argv, "--compression2step", params.compression2step)) {
    } else if (perf_test::check_arg_int(i, argc, argv, "--shmem", params.shmemsize)) {
      ++i;
    } else if (perf_test::check_arg_bool(i, argc, argv, "--CRWC", params.calculate_read_write_cost)) {
    } else if (perf_test::check_arg_str(i, argc, argv, "--CIF", params.coloring_input_file)) {
      ++i;
    } else if (perf_test::check_arg_str(i, argc, argv, "--COF", params.coloring_output_file)) {
      ++i;
    } else if (perf_test::check_arg_double(i, argc, argv, "--CCO", params.compression_cut_off)) {
      // if 0.85 set, if compression does not reduce flops by at least 15%
      // symbolic will run on original matrix. otherwise, it will compress the
      // graph and run symbolic on compressed one.
      ++i;
    } else if (perf_test::check_arg_double(i, argc, argv, "--FLHCO", params.first_level_hash_cut_off)) {
      // if linear probing is used as hash, what is the max occupancy percantage
      // we allow in the hash.
      ++i;
    } else if (perf_test::check_arg_bool(i, argc, argv, "--flop", params.calculate_read_write_cost)) {
      // print flop statistics. only for the first repeat.
      // note: if either --CRWC or --flop is passed, this parameter is set to
      // true
    } else if (perf_test::check_arg_int(i, argc, argv, "--mklsort", params.mkl_sort_option)) {
      // when mkl2 is run, the sort option to use.
      // 7:not to sort the output
      // 8:to sort the output
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--mklkeepout", params.mkl_keep_output)) {
      // mkl output is not kept.
      ++i;
    } else if (perf_test::check_arg_bool(i, argc, argv, "--checkoutput", params.check_output)) {
      // check correctness
    } else if (perf_test::check_arg_str(i, argc, argv, "--amtx", params.a_mtx_bin_file)) {
      // A at C=AxB
      ++i;
    } else if (perf_test::check_arg_str(i, argc, argv, "--bmtx", params.b_mtx_bin_file)) {
      // B at C=AxB.
      // if not provided, C = AxA will be performed.
      ++i;
    } else if (perf_test::check_arg_str(i, argc, argv, "--cmtx", params.c_mtx_bin_file)) {
      // if provided, C will be written to given file.
      // has to have ".bin", or ".crs" extension.
      ++i;
    } else if (perf_test::check_arg_bool(i, argc, argv, "--dynamic", params.use_dynamic_scheduling)) {
      // dynamic scheduling will be used for loops.
      // currently it is default already.
      // so has to use the dynamic schedulin.
    } else if (perf_test::check_arg_int(i, argc, argv, "--DENSEACCMAX", params.MaxColDenseAcc)) {
      // on CPUs and KNLs if DEFAULT algorithm or KKSPGEMM is chosen,
      // it uses dense accumulators for smaller matrices based on the size of
      // column (k) in B. Max column size is 250,000 for k to use dense
      // accumulators. this parameter overwrites this. with cache mode, or CPUs
      // with smaller thread count, where memory bandwidth is not an issue, this
      // cut-off can be increased to be more than 250,000
      ++i;
    } else if (perf_test::check_arg_bool(i, argc, argv, "--verbose", params.verbose)) {
      // print the timing and information about the inner steps.
      // if you are timing TPL libraries, for correct timing use verbose option,
      // because there are pre- post processing in these TPL kernel wraps.
    } else if (perf_test::check_arg_str(i, argc, argv, "--algorithm", algoStr)) {
      if (0 == Test::string_compare_no_case(algoStr, "DEFAULT")) {
        params.algorithm = KokkosSparse::SPGEMM_KK;
      } else if (0 == Test::string_compare_no_case(algoStr, "KKDEFAULT")) {
        params.algorithm = KokkosSparse::SPGEMM_KK;
      } else if (0 == Test::string_compare_no_case(algoStr, "KKSPGEMM")) {
        params.algorithm = KokkosSparse::SPGEMM_KK;
      }

      else if (0 == Test::string_compare_no_case(algoStr, "KKMEM")) {
        params.algorithm = KokkosSparse::SPGEMM_KK_MEMORY;
      } else if (0 == Test::string_compare_no_case(algoStr, "KKDENSE")) {
        params.algorithm = KokkosSparse::SPGEMM_KK_DENSE;
      } else if (0 == Test::string_compare_no_case(algoStr, "KKLP")) {
        params.algorithm = KokkosSparse::SPGEMM_KK_LP;
      } else if (0 == Test::string_compare_no_case(algoStr, "KKDEBUG")) {
        params.algorithm = KokkosSparse::SPGEMM_DEBUG;
      }

      else {
        std::cerr << "Unrecognized value for --algorithm (argument #" << i << "): " << argv[i] << std::endl;
        print_options();
        return 1;
      }
      ++i;
    } else {
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
      print_options();
      return 1;
    }
  }
  return 0;
}

template <typename ExecSpace>
void run_spgemm(int argc, char** argv, perf_test::CommonInputParams) {
  using namespace KokkosSparse;
  using namespace KokkosSparse::Experimental;

  using MemSpace  = typename ExecSpace::memory_space;
  using size_type = default_size_type;
  using lno_t     = default_lno_t;
  using scalar_t  = default_scalar;
  using device_t  = Kokkos::Device<ExecSpace, MemSpace>;
  using crsMat_t  = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device_t, void, size_type>;
  using KernelHandle =
      KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, ExecSpace, MemSpace, MemSpace>;

  KokkosKernels::Experiment::Parameters params;

  if (parse_inputs(params, argc, argv)) {
    return;
  }
  if (params.a_mtx_bin_file == "") {
    std::cerr << "Provide a and b matrix files" << std::endl;
    print_options();
    return;
  }

  crsMat_t A, B, C;

  // read a and b matrices

  A = KokkosSparse::Impl::read_kokkos_crst_matrix<crsMat_t>(params.a_mtx_bin_file.c_str());

  if ((params.b_mtx_bin_file == "" || params.a_mtx_bin_file == params.b_mtx_bin_file)) {
    std::cout << "B is not provided or is the same as A. Multiplying AxA." << std::endl;
    B = A;
  } else {
    B = KokkosSparse::Impl::read_kokkos_crst_matrix<crsMat_t>(params.b_mtx_bin_file.c_str());
  }

  int algorithm  = params.algorithm;
  int repeat     = params.repeat;
  int chunk_size = params.chunk_size;

  int shmemsize                 = params.shmemsize;
  int team_size                 = params.team_size;
  int use_dynamic_scheduling    = params.use_dynamic_scheduling;
  int verbose                   = params.verbose;
  int calculate_read_write_cost = params.calculate_read_write_cost;
  // char spgemm_step = params.spgemm_step;
  int vector_size     = params.vector_size;
  int check_output    = params.check_output;
  int mkl_keep_output = params.mkl_keep_output;
  // spgemm_step++;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::row_map_type::non_const_type lno_view_t;
  typedef typename crsMat_t::index_type::non_const_type lno_nnz_view_t;

  lno_view_t row_mapC;
  lno_nnz_view_t entriesC;
  scalar_view_t valuesC;

  KernelHandle kh;
  kh.set_team_work_size(chunk_size);
  kh.set_shmem_size(shmemsize);
  kh.set_suggested_team_size(team_size);
  kh.set_suggested_vector_size(vector_size);

  if (use_dynamic_scheduling) {
    kh.set_dynamic_scheduling(true);
  }
  if (verbose) {
    kh.set_verbose(true);
  }

  const lno_t m = A.numRows();
  const lno_t n = B.numRows();
  const lno_t k = B.numCols();

  if (verbose) std::cout << "m:" << m << " n:" << n << " k:" << k << std::endl;
  if (n < A.numCols()) {
    std::cerr << "left.numCols():" << A.numCols() << " right.numRows():" << B.numRows() << std::endl;
    exit(1);
  }

  // The reference product (for verifying correctness)
  // Don't allocate them if they won't be used, but they must be declared here.
  lno_view_t row_mapC_ref;
  lno_nnz_view_t entriesC_ref;
  scalar_view_t valuesC_ref;
  // Reference output has same type as actual output
  crsMat_t C_ref;

  if (check_output) {
    if (verbose) std::cout << "Running a reference algorithm" << std::endl;
    row_mapC_ref = lno_view_t("non_const_lnow_row", m + 1);
    KernelHandle sequential_kh;
    sequential_kh.set_team_work_size(chunk_size);
    sequential_kh.set_shmem_size(shmemsize);
    sequential_kh.set_suggested_team_size(team_size);
    sequential_kh.create_spgemm_handle(KokkosSparse::SPGEMM_SERIAL);

    if (use_dynamic_scheduling) {
      sequential_kh.set_dynamic_scheduling(true);
    }

    spgemm_symbolic(&sequential_kh, m, n, k, A.graph.row_map, A.graph.entries, TRANSPOSEFIRST, B.graph.row_map,
                    B.graph.entries, TRANSPOSESECOND, row_mapC_ref);

    ExecSpace().fence();

    size_type c_nnz_size = sequential_kh.get_spgemm_handle()->get_c_nnz();
    entriesC_ref         = lno_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC"), c_nnz_size);
    valuesC_ref          = scalar_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC"), c_nnz_size);

    spgemm_numeric(&sequential_kh, m, n, k, A.graph.row_map, A.graph.entries, A.values, TRANSPOSEFIRST,

                   B.graph.row_map, B.graph.entries, B.values, TRANSPOSESECOND, row_mapC_ref, entriesC_ref,
                   valuesC_ref);
    ExecSpace().fence();

    C_ref = crsMat_t("CorrectC", m, k, valuesC_ref.extent(0), valuesC_ref, row_mapC_ref, entriesC_ref);
  }

  for (int i = 0; i < repeat; ++i) {
    kh.create_spgemm_handle(KokkosSparse::SPGEMMAlgorithm(algorithm));

    kh.get_spgemm_handle()->mkl_keep_output = mkl_keep_output;
    kh.get_spgemm_handle()->set_mkl_sort_option(params.mkl_sort_option);

    // if mkl2 input needs to be converted to 1base.
    kh.get_spgemm_handle()->mkl_convert_to_1base = true;

    // 250000 default. if cache-mode is used on KNL can increase to 1M.
    kh.get_spgemm_handle()->MaxColDenseAcc = params.MaxColDenseAcc;

    if (i == 0) {
      kh.get_spgemm_handle()->set_read_write_cost_calc(calculate_read_write_cost);
    }
    // do the compression whether in 2 step, or 1 step.
    kh.get_spgemm_handle()->set_compression_steps(!params.compression2step);
    // whether to scale the hash more. default is 1, so no scale.
    kh.get_spgemm_handle()->set_min_hash_size_scale(params.minhashscale);
    // max occupancy in 1-level LP hashes. LL hashes can be 100%
    kh.get_spgemm_handle()->set_first_level_hash_cut_off(params.first_level_hash_cut_off);
    // min reduction on FLOPs to run compression
    kh.get_spgemm_handle()->set_compression_cut_off(params.compression_cut_off);

    row_mapC = lno_view_t("non_const_lnow_row", m + 1);
    entriesC = lno_nnz_view_t("entriesC (empty)", 0);
    valuesC  = scalar_view_t("valuesC (empty)", 0);

    Kokkos::Timer timer1;
    spgemm_symbolic(&kh, m, n, k, A.graph.row_map, A.graph.entries, TRANSPOSEFIRST, B.graph.row_map, B.graph.entries,
                    TRANSPOSESECOND, row_mapC);

    ExecSpace().fence();
    double symbolic_time = timer1.seconds();

    Kokkos::Timer timer3;
    size_type c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
    if (verbose) std::cout << "C SIZE:" << c_nnz_size << std::endl;
    if (c_nnz_size) {
      entriesC = lno_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC"), c_nnz_size);
      valuesC  = scalar_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC"), c_nnz_size);
    }
    spgemm_numeric(&kh, m, n, k, A.graph.row_map, A.graph.entries, A.values, TRANSPOSEFIRST, B.graph.row_map,
                   B.graph.entries, B.values, TRANSPOSESECOND, row_mapC, entriesC, valuesC);

    ExecSpace().fence();
    double numeric_time = timer3.seconds();

    std::cout << "mm_time:" << symbolic_time + numeric_time << " symbolic_time:" << symbolic_time
              << " numeric_time:" << numeric_time << std::endl;
  }
  if (verbose) {
    std::cout << "row_mapC:" << row_mapC.extent(0) << std::endl;
    std::cout << "entriesC:" << entriesC.extent(0) << std::endl;
    std::cout << "valuesC:" << valuesC.extent(0) << std::endl;
    KokkosKernels::Impl::print_1Dview(valuesC);
    KokkosKernels::Impl::print_1Dview(entriesC);
    KokkosKernels::Impl::print_1Dview(row_mapC);
  }
  crsMat_t C_result("CrsMatrixC", m, k, valuesC.extent(0), valuesC, row_mapC, entriesC);
  if (check_output) {
    bool is_identical = is_same_matrix<crsMat_t, device_t>(C_result, C_ref);
    if (!is_identical) {
      std::cerr << "SpGEMM result differs with reference implementation.\n";
      exit(1);
    } else {
      std::cerr << "SpGEMM result matches reference implementation.\n";
    }
  }

  if (params.c_mtx_bin_file != "") {
    KokkosSparse::sort_crs_matrix(C_result);

    KokkosSparse::Impl::write_graph_bin((lno_t)(C_result.numRows()), (size_type)(C_result.nnz()),
                                        C_result.graph.row_map.data(), C_result.graph.entries.data(),
                                        C_result.values.data(), params.c_mtx_bin_file.c_str());
  }
}

#define KOKKOSKERNELS_PERF_TEST_NAME run_spgemm
#include "KokkosKernels_perf_test_instantiation.hpp"
int main(int argc, char** argv) { return main_instantiation(argc, argv); }  // main
