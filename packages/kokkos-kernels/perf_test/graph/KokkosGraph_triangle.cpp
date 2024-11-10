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
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosGraph_Triangle.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_IOUtils.hpp"  //for read_kokkos_crst_graph
#include "KokkosKernels_TestUtils.hpp"
#include "KokkosKernels_TestParameters.hpp"
#include "KokkosKernels_perf_test_utilities.hpp"

template <size_t BufSize, typename SpaceType = Kokkos::DefaultExecutionSpace>
struct Flush {
  typedef double value_type;

  // flush a large host buffer
  Kokkos::View<value_type *, SpaceType> _buf;
  Flush(int flush_option) : _buf("Flush::buf", BufSize) {
    Kokkos::deep_copy(_buf, 1);
    Kokkos::fence();
    if (flush_option == 2) {
      for (size_t i = 0; i < BufSize; ++i) {
        _buf(i) = rand();
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type &update) { update = 0; }

  KOKKOS_INLINE_FUNCTION
  void join(value_type &update, const value_type &input) { update += input; }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type &update) const { update += _buf[i]; }

  void run() {
    double sum = 0;
    Kokkos::parallel_reduce("KokkosGraph::PerfTest::Flush", Kokkos::RangePolicy<SpaceType>(0, BufSize / sizeof(double)),
                            *this, sum);
    SpaceType().fence();
    std::cout << "Flush sum:" << sum << std::endl;
    FILE *fp = fopen("/dev/null", "w");
    fprintf(fp, "%f\n", sum);
    fclose(fp);
  }
};

void print_options() {
  std::cerr << "Options\n" << std::endl;
  std::cerr << perf_test::list_common_options();
  std::cerr << "Input Matrix                       : --amtx [path_to_input_matrix]" << std::endl;
  std::cerr << "\tInput Matrix format can be multiple formats. If it ends with:" << std::endl;
  std::cerr << "\t\t.mtx: it will read matrix market format." << std::endl;
  std::cerr << "\t\t.bin: it will read binary crs matrix format." << std::endl;
  std::cerr << "\t\t.crs: it will read text crs matrix format." << std::endl;
  std::cerr << "--algorithm                          :" << std::endl;
  // BMK 3-28-23: these algorithms do not give correct triangle counts
  // std::cerr << "\tTRIANGLEAI: for Adj x Incidence" << std::endl;
  // std::cerr << "\tTRIANGLEIA: for Incidence x Adj -- implementing set "
  //             "intersection (2D) -- 3rd fastest"
  //          << std::endl;
  // std::cerr
  //    << "\tTRIANGLEIAUNION: for Incidence x Adj -- implementing set union "
  //    << std::endl;
  std::cerr << "\tTRIANGLELL: Lower x Lower -- usually fastest " << std::endl;
  std::cerr << "\tTRIANGLELU: Lower x Upper -- usually 2nd fastest " << std::endl;
  std::cerr << "--FLOP                               : Calculate and print the "
               "number of operations. This will be calculated on the first run."
            << std::endl;
  std::cerr << "--COMPRESSION [0|1]                   : Enable disable "
               "compression. Default:1."
            << std::endl;
  std::cerr << "--RS [0|1|2]                         : Whether to sort lower "
               "triangular matrix. 0 - no sort, 1 - sort, 2 - algorithm "
               "decides based on max row size (default)"
            << std::endl;
  std::cerr << "--accumulator [default|dense|sparse] : what type of "
               "accumulator to use."
            << std::endl;
  std::cerr << "--RLT                                : If given, lower triangle will "
               "be used for AdjxIncidence or Incidence x Adj algorithms."
            << std::endl;
  std::cerr << "--dynamic                            : If set, dynamic schedule will "
               "be used. Currently default is dynamic scheduling as well."
            << std::endl;
  std::cerr << "--verbose                            : If set, the inner timer "
               "stats will be printed."
            << std::endl;
  std::cerr << "--repeat [repeatnum]                 : how many repeats will be run." << std::endl;
  std::cerr << "--chunksize [chunksize]              : how many vertices are "
               "executed with in a loop index. Default is 16."
            << std::endl;
  std::cerr << "--sort_option [0|1|2]                : How lower triangle will "
               "be sorted. 0: for largest to bottom, 1 for largest to top, 2 "
               "for interleaved."
            << std::endl;
  std::cerr << "--cache_flush [0|1|2]                : Flush between repetitions. 0 "
               "- no flush, 1 - soft flush, 2 - hard flush with random numbers."
            << std::endl;

  std::cerr << "\nSuggested use of LL: executable --amtx path_to_file.bin "
               "--algorithm TRIANGLELL --repeat 6 --verbose --chunksize [4|16]"
            << std::endl;
  std::cerr << "Suggested use of LU: executable --amtx path_to_file.bin "
               "--algorithm TRIANGLELU --repeat 6 --verbose --chunksize [4|16]"
            << std::endl;
  // std::cerr
  //    << "Suggested use of AI: executable --amtx path_to_file.bin --algorithm
  //    "
  //       "TRIANGLEIA --repeat 6 --verbose --chunksize [4|16] rlt"
  //    << std::endl;
}

int parse_inputs(KokkosKernels::Experiment::Parameters &params, int argc, char **argv) {
  for (int i = 1; i < argc; ++i) {
    if (0 == Test::string_compare_no_case(argv[i], "--repeat")) {
      params.repeat = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--triangle_operation")) {
      params.triangle_options = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--chunksize")) {
      params.chunk_size = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--teamsize")) {
      params.team_size = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--vectorsize")) {
      params.vector_size = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--compression")) {
      params.apply_compression = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--flop")) {
      params.calculate_read_write_cost = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--CIF")) {
      params.coloring_input_file = argv[++i];
    } else if (0 == Test::string_compare_no_case(argv[i], "--COF")) {
      params.coloring_output_file = argv[++i];
    } else if (0 == Test::string_compare_no_case(argv[i], "--mhscale")) {
      params.minhashscale = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--mcscale")) {
      params.multi_color_scale = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--shmem")) {
      params.shmemsize = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--compression2step")) {
      params.compression2step = true;
    } else if (0 == Test::string_compare_no_case(argv[i], "--mklsort")) {
      params.mkl_sort_option = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--mklkeepout")) {
      params.mkl_keep_output = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--checkoutput")) {
      params.check_output = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--amtx")) {
      params.a_mtx_bin_file = argv[++i];
    } else if (0 == Test::string_compare_no_case(argv[i], "--dynamic")) {
      params.use_dynamic_scheduling = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--cache_flush")) {
      params.cache_flush = atoi(argv[++i]);
    }

    else if (0 == Test::string_compare_no_case(argv[i], "--RLT")) {
      params.right_lower_triangle = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--RS")) {
      params.right_sort = atoi(argv[++i]);
    }

    else if (0 == Test::string_compare_no_case(argv[i], "--verbose")) {
      params.verbose = 1;
    }

    else if (0 == Test::string_compare_no_case(argv[i], "--accumulator")) {
      ++i;
      if (0 == Test::string_compare_no_case(argv[i], "default")) {
        params.accumulator = 0;
      } else if (0 == Test::string_compare_no_case(argv[i], "dense")) {
        params.accumulator = 1;
      } else if (0 == Test::string_compare_no_case(argv[i], "sparse")) {
        params.accumulator = 2;
      } else {
        std::cerr << "1-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
        print_options();
        return 1;
      }
    } else if (0 == Test::string_compare_no_case(argv[i], "--algorithm")) {
      ++i;
      if (0 == Test::string_compare_no_case(argv[i], "TRIANGLEAI")) {
        params.algorithm = 16;
        std::cerr << "\nAlgorithm TRIANGLEAI is disabled (produces incorrect "
                     "triangle count)\n";
        return 1;
      } else if (0 == Test::string_compare_no_case(argv[i], "TRIANGLEIA")) {
        params.algorithm = 17;
        std::cerr << "\nAlgorithm TRIANGLEIA is disabled (produces incorrect "
                     "triangle count)\n";
        return 1;
      } else if (0 == Test::string_compare_no_case(argv[i], "TRIANGLEIAUNION")) {
        params.algorithm = 18;
        std::cerr << "\nAlgorithm TRIANGLEIAUNION is disabled (produces "
                     "incorrect triangle count)\n";
        return 1;
      } else if (0 == Test::string_compare_no_case(argv[i], "TRIANGLELL")) {
        params.algorithm = 19;
      } else if (0 == Test::string_compare_no_case(argv[i], "TRIANGLELU")) {
        params.algorithm = 20;
      } else {
        std::cerr << "2-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
        print_options();
        return 1;
      }
    } else {
      std::cerr << "3-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
      print_options();
      return 1;
    }
  }
  return 0;
}

template <typename exec_space>
void run_experiment(int argc, char **argv, perf_test::CommonInputParams) {
  using namespace KokkosSparse;
  using mem_space = typename exec_space::memory_space;
  using device_t  = Kokkos::Device<exec_space, mem_space>;
  using lno_t     = default_lno_t;
  using size_type = default_size_type;
  using graph_t   = Kokkos::StaticCrsGraph<lno_t, default_layout, device_t, void, size_type>;
  using KernelHandle =
      KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, lno_t, exec_space, mem_space, mem_space>;

  if (KokkosKernels::Impl::kk_is_gpu_exec_space<exec_space>()) {
    std::cerr << "** Triangle counting is currently not supported on GPU backends.\n";
    return;
  }

  KokkosKernels::Experiment::Parameters params;

  if (parse_inputs(params, argc, argv)) {
    return;
  }
  if (params.a_mtx_bin_file == "") {
    std::cerr << "Provide a graph file" << std::endl;
    print_options();
    return;
  }

  std::cout << "Sizeof(idx):" << sizeof(lno_t) << " sizeof(size_type):" << sizeof(size_type) << std::endl;

  // read graph
  graph_t crsGraph = KokkosSparse::Impl::read_kokkos_crst_graph<graph_t>(params.a_mtx_bin_file.c_str());

  int algorithm  = params.algorithm;
  int repeat     = params.repeat;
  int chunk_size = params.chunk_size;

  int shmemsize              = params.shmemsize;
  int team_size              = params.team_size;
  int use_dynamic_scheduling = params.use_dynamic_scheduling;
  int verbose                = params.verbose;

  int accumulator = params.accumulator;
  int vector_size = params.vector_size;

  Kokkos::View<size_t *, exec_space> row_mapC;

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
  const lno_t m = crsGraph.numRows();

  for (int i = 0; i < repeat; ++i) {
    size_type rowmap_size = crsGraph.entries.extent(0);
    switch (algorithm) {
      case 16:
        kh.create_spgemm_handle(SPGEMM_KK_TRIANGLE_AI);
        rowmap_size = m;
        break;
      case 17:
        kh.create_spgemm_handle(SPGEMM_KK_TRIANGLE_IA);
        std::cout << "IA" << std::endl;
        break;
      case 18: kh.create_spgemm_handle(SPGEMM_KK_TRIANGLE_IA_UNION); break;
      case 19:
        kh.create_spgemm_handle(SPGEMM_KK_TRIANGLE_LL);
        rowmap_size = m;
        break;
      case 20:
        kh.create_spgemm_handle(SPGEMM_KK_TRIANGLE_LU);
        rowmap_size = m;
        break;
      default: kh.create_spgemm_handle(SPGEMM_KK_TRIANGLE_IA); break;
    }

    kh.get_spgemm_handle()->set_compression_steps(!params.compression2step);

    kh.get_spgemm_handle()->set_sort_lower_triangular(params.right_sort);
    kh.get_spgemm_handle()->set_create_lower_triangular(params.right_lower_triangle);
    kh.get_spgemm_handle()->set_compression(params.apply_compression);
    kh.get_spgemm_handle()->set_min_hash_size_scale(params.minhashscale);

    switch (accumulator) {
      case 0:
      default: kh.get_spgemm_handle()->set_accumulator_type(SPGEMM_ACC_DEFAULT); break;
      case 1: kh.get_spgemm_handle()->set_accumulator_type(SPGEMM_ACC_DENSE); break;
      case 2: kh.get_spgemm_handle()->set_accumulator_type(SPGEMM_ACC_SPARSE); break;
    }

    constexpr size_t LLC_CAPACITY = 128 * 1024 * 1024;
    if (params.cache_flush) {
      std::cout << "Flushing cache with option:" << params.cache_flush << std::endl;
      Flush<LLC_CAPACITY, exec_space> flush(params.cache_flush);
      flush.run();
    }
    if (i == 0) {
      kh.get_spgemm_handle()->set_read_write_cost_calc(params.calculate_read_write_cost);
    }

    Kokkos::Timer timer1;

    row_mapC = Kokkos::View<size_t *, exec_space>("non_const_lnow_row", rowmap_size);

    double symbolic_time = 0;
    if (params.triangle_options == 0) {
      if (params.apply_compression) {
        KokkosGraph::Experimental::triangle_generic(
            &kh, m, crsGraph.row_map, crsGraph.entries,
            KOKKOS_LAMBDA(const lno_t &row, const lno_t & /* col_set_index */, const lno_t &col_set,
                          const lno_t & /* thread_id */) { row_mapC(row) += KokkosKernels::Impl::pop_count(col_set); });
      } else {
        KokkosGraph::Experimental::triangle_generic(
            &kh, m, crsGraph.row_map, crsGraph.entries,
            KOKKOS_LAMBDA(const lno_t &row, const lno_t & /*col_set_index*/, const lno_t & /*col_set*/,
                          const lno_t & /*thread_id*/) { row_mapC(row)++; });
      }

      size_t num_triangles = 0;
      KokkosKernels::Impl::kk_reduce_view<Kokkos::View<size_t *, exec_space>, exec_space>(rowmap_size, row_mapC,
                                                                                          num_triangles);
      symbolic_time = timer1.seconds();
      std::cout << "num_triangles:" << num_triangles << std::endl;
    }
    kh.destroy_spgemm_handle();
    std::cout << "mm_time:" << symbolic_time << std::endl;
  }
}

#define KOKKOSKERNELS_PERF_TEST_NAME run_experiment
#include "KokkosKernels_perf_test_instantiation.hpp"
int main(int argc, char **argv) { return main_instantiation(argc, argv); }  // main
