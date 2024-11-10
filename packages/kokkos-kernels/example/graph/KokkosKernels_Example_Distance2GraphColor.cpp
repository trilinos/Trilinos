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

// Standard Library Headers
#include <iomanip>
#include <iostream>
#include <string>
#include <unistd.h>

// Kokkos Headers
#include <Kokkos_Core.hpp>
#include <Kokkos_StaticCrsGraph.hpp>

// Kokkos-Kernels Headers
#include <KokkosGraph_Distance2Color.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

using namespace KokkosGraph;

#ifdef KOKKOSKERNELS_INST_DOUBLE
using kk_scalar_type = double;
#else
#ifdef KOKKOSKERNELS_INST_FLOAT
using kk_scalar_type = float;
#endif
#endif

#ifdef KOKKOSKERNELS_INST_OFFSET_INT
using kk_size_type = int;
#else
#ifdef KOKKOSKERNELS_INST_OFFSET_SIZE_T
using kk_size_type   = size_t;
#endif
#endif

#ifdef KOKKOSKERNELS_INST_ORDINAL_INT
using kk_lno_type = int;
#else
#ifdef KOKKOSKERNELS_INST_ORDINAL_INT64_T
using kk_lno_type    = int64_t;
#endif
#endif

using namespace KokkosGraph;

namespace KokkosKernels {
namespace Example {

struct Parameters {
  int algorithm;
  int repeat;
  int chunk_size;
  int output_graphviz_vert_max;
  int output_graphviz;
  int shmemsize;
  int verbose_level;
  int check_output;
  char* coloring_input_file;
  char* coloring_output_file;
  int output_histogram;
  int use_threads;
  int use_openmp;
  int use_cuda;
  int use_serial;
  int validate;
  char* mtx_bin_file;

  Parameters() {
    algorithm                = 0;
    repeat                   = 6;
    chunk_size               = -1;
    shmemsize                = 16128;
    verbose_level            = 0;
    check_output             = 0;
    coloring_input_file      = NULL;
    coloring_output_file     = NULL;
    output_histogram         = 0;
    output_graphviz          = 0;
    output_graphviz_vert_max = 1500;
    use_threads              = 0;
    use_openmp               = 0;
    use_cuda                 = 0;
    use_serial               = 0;
    validate                 = 0;
    mtx_bin_file             = NULL;
  }
};

void print_options(std::ostream& os, const char* app_name, unsigned int indent = 0) {
  std::string spaces(indent, ' ');
  os << "Usage:" << std::endl
     << spaces << "  " << app_name << " [parameters]" << std::endl
     << std::endl
     << spaces << "Parameters:" << std::endl
     << spaces << "  Parallelism (select one of the following):" << std::endl
     << spaces << "      --serial <N>        Execute serially." << std::endl
     << spaces << "      --threads <N>       Use N posix threads." << std::endl
     << spaces << "      --openmp <N>        Use OpenMP with N threads." << std::endl
     << spaces << "      --cuda              Use CUDA" << std::endl
     << std::endl
     << spaces << "  Required Parameters:" << std::endl
     << spaces << "      --amtx <filename>   Input file in Matrix Market format (.mtx)." << std::endl
     << std::endl
     << spaces
     << "      --algorithm <algorithm_name>   Set the algorithm to use.  "
        "Allowable values are:"
     << std::endl
     << spaces
     << "                 COLORING_D2_MATRIX_SQUARED  - Matrix-squared + "
        "Distance-1 method."
     << std::endl
     << spaces
     << "                 COLORING_D2_SERIAL          - Serial algorithm (must "
        "use with 'serial' mode)"
     << std::endl
     << spaces
     << "                 COLORING_D2_VB              - Vertex Based method "
        "using boolean forbidden array (Default)."
     << std::endl
     << spaces
     << "                 COLORING_D2_VB_BIT          - VB with Bitvector "
        "Forbidden Array"
     << std::endl
     << spaces
     << "                 COLORING_D2_VB_BIT_EF       - VB_BIT with Edge "
        "Filtering"
     << std::endl
     << std::endl
     << spaces << "  Optional Parameters:" << std::endl
     << spaces
     << "      --output-histogram              Print out a histogram of the "
        "colors."
     << std::endl
     << spaces
     << "      --output-graphviz               Write the output to a graphviz "
        "file (G.dot)."
     << std::endl
     << spaces
     << "                                      Note: Vertices with color 0 "
        "will be filled in and colored"
     << std::endl
     << spaces
     << "      --output-graphviz-vert-max <N>  Upper limit of vertices in G to "
        "allow graphviz output. Default=1500."
     << std::endl
     << spaces
     << "                                      Requires --output-graphviz to "
        "also be enabled."
     << std::endl
     << spaces
     << "      --validate                      Check that the coloring is a "
        "valid distance-2 graph coloring"
     << std::endl
     << spaces
     << "      --verbose-level <N>             Set verbosity level [0..5] "
        "where N > 0 means print verbose messags."
     << std::endl
     << spaces << "                                      Default: 0" << std::endl
     << spaces << "      --help                          Print out command line help." << std::endl
     << spaces << " " << std::endl;
}

int parse_inputs(KokkosKernels::Example::Parameters& params, int argc, char** argv) {
  bool got_required_param_amtx      = false;
  bool got_required_param_algorithm = false;

  for (int i = 1; i < argc; ++i) {
    if (0 == Test::string_compare_no_case(argv[i], "--threads")) {
      params.use_threads = atoi(argv[++i]);
      // std::cout << "use_threads = " << params.use_threads << std::endl;
    } else if (0 == Test::string_compare_no_case(argv[i], "--serial")) {
      params.use_serial = atoi(argv[++i]);
      // std::cout << "use_serial = " << params.use_serial << std::endl;
    } else if (0 == Test::string_compare_no_case(argv[i], "--openmp")) {
      params.use_openmp = atoi(argv[++i]);
      // std::cout << "use_openmp = " << params.use_openmp << std::endl;
    } else if (0 == Test::string_compare_no_case(argv[i], "--cuda")) {
      params.use_cuda = 1;
      // std::cout << "use_cuda = " << params.use_cuda << std::endl;
    } else if (0 == Test::string_compare_no_case(argv[i], "--amtx")) {
      got_required_param_amtx = true;
      params.mtx_bin_file     = argv[++i];
    } else if (0 == Test::string_compare_no_case(argv[i], "--validate")) {
      params.validate = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--verbose-level")) {
      params.verbose_level = atoi(argv[++i]);
      params.verbose_level = std::min(5, params.verbose_level);
      params.verbose_level = std::max(0, params.verbose_level);
    } else if (0 == Test::string_compare_no_case(argv[i], "--output-histogram")) {
      params.output_histogram = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--output-graphviz")) {
      params.output_graphviz = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--output-graphviz-vert-max")) {
      params.output_graphviz_vert_max = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--algorithm")) {
      ++i;
      if (0 == Test::string_compare_no_case(argv[i], "COLORING_D2_MATRIX_SQUARED")) {
        params.algorithm             = 1;
        got_required_param_algorithm = true;
      } else if (0 == Test::string_compare_no_case(argv[i], "COLORING_D2_SERIAL")) {
        params.algorithm             = 2;
        got_required_param_algorithm = true;
      } else if (0 == Test::string_compare_no_case(argv[i], "COLORING_D2_VB") ||
                 0 == Test::string_compare_no_case(argv[i], "COLORING_D2")) {
        params.algorithm             = 3;
        got_required_param_algorithm = true;
      } else if (0 == Test::string_compare_no_case(argv[i], "COLORING_D2_VB_BIT")) {
        params.algorithm             = 4;
        got_required_param_algorithm = true;
      } else if (0 == Test::string_compare_no_case(argv[i], "COLORING_D2_VB_BIT_EF")) {
        params.algorithm             = 5;
        got_required_param_algorithm = true;
      } else {
        std::cerr << "2-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
        print_options(std::cout, argv[0]);
        return 1;
      }
    } else if (0 == Test::string_compare_no_case(argv[i], "--help") ||
               0 == Test::string_compare_no_case(argv[i], "-h")) {
      print_options(std::cout, argv[0]);
      return 1;
    } else {
      std::cerr << "3-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
      print_options(std::cout, argv[0]);
      return 1;
    }
  }

  if (!got_required_param_amtx) {
    std::cout << "Missing required parameter amtx" << std::endl << std::endl;
    print_options(std::cout, argv[0]);
    return 1;
  }
  if (!got_required_param_algorithm) {
    std::cout << "Missing required parameter algorithm" << std::endl << std::endl;
    print_options(std::cout, argv[0]);
    return 1;
  }
  if (!params.use_serial && !params.use_threads && !params.use_openmp && !params.use_cuda) {
    print_options(std::cout, argv[0]);
    return 1;
  }
  return 0;
}

template <typename ExecSpace, typename DataType, typename CrsGraph_type, typename TempMemSpace,
          typename PersistentMemSpace>
void run_example(CrsGraph_type crsGraph, DataType num_cols, Parameters params) {
  using namespace KokkosGraph;
  using namespace KokkosGraph::Experimental;

  int algorithm = params.algorithm;
  int shmemsize = params.shmemsize;

  using lno_view_type     = typename CrsGraph_type::row_map_type::non_const_type;
  using lno_nnz_view_type = typename CrsGraph_type::entries_type::non_const_type;
  using size_type         = typename lno_view_type::non_const_value_type;
  using lno_type          = typename lno_nnz_view_type::non_const_value_type;
  using KernelHandle_type =
      KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_type, kk_scalar_type, ExecSpace, TempMemSpace,
                                                       PersistentMemSpace>;

  // Create a kernel handle
  KernelHandle_type kh;
  kh.set_shmem_size(shmemsize);

  if (params.verbose_level > 0) {
    kh.set_verbose(true);
  }

  // ------------------------------------------
  // Set up the D2 coloring kernel handle
  // ------------------------------------------
  std::string label_algorithm;
  switch (algorithm) {
    case 1:
      kh.create_distance2_graph_coloring_handle(COLORING_D2_MATRIX_SQUARED);
      label_algorithm = "COLORING_D2_MATRIX_SQUARED";
      break;
    case 2:
      kh.create_distance2_graph_coloring_handle(COLORING_D2_SERIAL);
      label_algorithm = "COLORING_D2_SERIAL";
      break;
    case 3:
      kh.create_distance2_graph_coloring_handle(COLORING_D2_VB);
      label_algorithm = "COLORING_D2_VB";
      break;
    case 4:
      kh.create_distance2_graph_coloring_handle(COLORING_D2_VB_BIT);
      label_algorithm = "COLORING_D2_VB_BIT";
      break;
    case 5:
      kh.create_distance2_graph_coloring_handle(COLORING_D2_VB_BIT_EF);
      label_algorithm = "COLORING_D2_VB_BIT_EF";
      break;
    default:
      kh.create_distance2_graph_coloring_handle(COLORING_D2_VB);
      label_algorithm = "COLORING_D2_VB";
      break;
  }

  std::cout << std::endl << "Run Graph Color D2 (" << label_algorithm << ")" << std::endl;

  // ------------------------------------------
  // Call the distance-2 graph coloring routine
  // ------------------------------------------
  graph_compute_distance2_color(&kh, crsGraph.numRows(), num_cols, crsGraph.row_map, crsGraph.entries, crsGraph.row_map,
                                crsGraph.entries);

  // ------------------------------------------
  // Get the results
  // ------------------------------------------
  size_t num_colors = kh.get_distance2_graph_coloring_handle()->get_num_colors();
  size_t num_phases = kh.get_distance2_graph_coloring_handle()->get_num_phases();

  if (params.verbose_level > 0) {
    std::cout << "Total Time: " << kh.get_distance2_graph_coloring_handle()->get_overall_coloring_time() << std::endl
              << "Num colors: " << kh.get_distance2_graph_coloring_handle()->get_num_colors() << std::endl
              << "Num Phases: " << kh.get_distance2_graph_coloring_handle()->get_num_phases() << std::endl
              << "Colors:\n\t";
    KokkosKernels::Impl::print_1Dview(kh.get_distance2_graph_coloring_handle()->get_vertex_colors());
    std::cout << std::endl;
  }

  // ------------------------------------------
  // Save coloring to a GraphViz file
  // ------------------------------------------
  if (params.output_graphviz && crsGraph.numRows() <= params.output_graphviz_vert_max) {
    auto colors = kh.get_distance2_graph_coloring_handle()->get_vertex_colors();

    std::ofstream os("G.dot", std::ofstream::out);

    kh.get_distance2_graph_coloring_handle()->dump_graphviz(os, crsGraph.numRows(), crsGraph.row_map, crsGraph.entries,
                                                            colors);
  }

  // ------------------------------------------
  // Verify correctness
  // ------------------------------------------
  std::string str_color_is_valid = "UNKNOWN";
  if (0 != params.validate) {
    str_color_is_valid = "VALID";

    bool d2_coloring_is_valid            = false;
    bool d2_coloring_validation_flags[4] = {false};

    d2_coloring_is_valid = KokkosGraph::Impl::graph_verify_distance2_color(
        &kh, crsGraph.numRows(),
        // crsGraph.numCols(),
        num_cols, crsGraph.row_map, crsGraph.entries, crsGraph.row_map, crsGraph.entries, d2_coloring_validation_flags);

    // Print out messages based on coloring validation check.
    if (d2_coloring_is_valid) {
      std::cout << std::endl << "Distance-2 Graph Coloring is VALID" << std::endl << std::endl;
    } else {
      str_color_is_valid = "INVALID";
      std::cout << std::endl
                << "Distance-2 Graph Coloring is NOT VALID" << std::endl
                << "  - Vert(s) left uncolored : " << d2_coloring_validation_flags[1] << std::endl
                << "  - Invalid D2 Coloring    : " << d2_coloring_validation_flags[2] << std::endl
                << std::endl;
    }
    if (d2_coloring_validation_flags[3]) {
      std::cout << "Distance-2 Graph Coloring may have poor quality." << std::endl
                << "  - Vert(s) have high color value : " << d2_coloring_validation_flags[3] << std::endl
                << std::endl;
    }
  }

  // ------------------------------------------
  // Print out a histogram of the colors
  // ------------------------------------------
  if (0 != params.output_histogram) {
    KokkosGraph::Impl::graph_print_distance2_color_histogram(&kh, crsGraph.numRows(), num_cols, crsGraph.row_map,
                                                             crsGraph.entries, crsGraph.row_map, crsGraph.entries,
                                                             false);
  }

  // ------------------------------------------
  // Print out a summary
  // ------------------------------------------
  std::string mtx_bin_file = params.mtx_bin_file;
  mtx_bin_file             = mtx_bin_file.substr(mtx_bin_file.find_last_of("/\\") + 1);

  std::cout << "Summary" << std::endl
            << "-------" << std::endl
            << "    KExecSName     : " << Kokkos::DefaultExecutionSpace::name() << std::endl
            << "    Filename       : " << mtx_bin_file << std::endl
            << "    Num Verts      : " << crsGraph.numRows() << std::endl
            << "    Num Edges      : " << crsGraph.entries.extent(0) << std::endl
            << "    Concurrency    : " << Kokkos::DefaultExecutionSpace().concurrency() << std::endl
            << "    Algorithm      : " << label_algorithm << std::endl
            << "Coloring Stats" << std::endl
            << "    Num colors     : " << num_colors << std::endl
            << "    Num Phases     : " << num_phases << std::endl
            << "    Validation     : " << str_color_is_valid << std::endl
            << std::endl;

}  // run_example()

template <typename size_type, typename lno_type, typename exec_space, typename hbm_mem_space>
void driver(Parameters params) {
  using myExecSpace  = exec_space;
  using myFastDevice = Kokkos::Device<exec_space, hbm_mem_space>;
  using crstmat_type = typename KokkosSparse::CrsMatrix<double, lno_type, myFastDevice, void, size_type>;
  using graph_type   = typename crstmat_type::StaticCrsGraphType;
  using data_type    = typename graph_type::data_type;

  char* mat_file = params.mtx_bin_file;

  crstmat_type crsmat = KokkosKernels::Impl::read_kokkos_crst_matrix<crstmat_type>(mat_file);
  graph_type crsgraph = crsmat.graph;
  data_type num_cols  = crsmat.numCols();

  KokkosKernels::Example::run_example<myExecSpace, data_type, graph_type, hbm_mem_space, hbm_mem_space>(
      crsgraph, num_cols, params);

}  // driver()

}  // namespace Example
}  // namespace KokkosKernels

int main(int argc, char* argv[]) {
  KokkosKernels::Example::Parameters params;

  if (parse_inputs(params, argc, argv)) {
    return 1;
  }

  if (params.mtx_bin_file == NULL) {
    std::cerr << "Provide a matrix file" << std::endl;
    return 0;
  }

  const int num_threads = params.use_openmp;  // Assumption is that use_openmp variable is provided
                                              // as number of threads
  const int device_id = 0;
  Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(num_threads).set_device_id(device_id));

  // Print out information about the configuration of the run if verbose_level
  // >= 5
  if (params.verbose_level >= 5) {
    Kokkos::print_configuration(std::cout);
  }

#if defined(KOKKOS_ENABLE_OPENMP)
  if (params.use_openmp) {
    KokkosKernels::Example::driver<kk_size_type, kk_lno_type, Kokkos::OpenMP, Kokkos::OpenMP::memory_space>(params);
  }
#endif

#if defined(KOKKOS_ENABLE_CUDA)
  if (params.use_cuda) {
    KokkosKernels::Example::driver<kk_size_type, kk_lno_type, Kokkos::Cuda, Kokkos::Cuda::memory_space>(params);
  }
#endif

#if defined(KOKKOS_ENABLE_SERIAL)
  if (params.use_serial) {
    KokkosKernels::Example::driver<kk_size_type, kk_lno_type, Kokkos::Serial, Kokkos::Serial::memory_space>(params);
  }
#endif

  Kokkos::finalize();

  return 0;
}
