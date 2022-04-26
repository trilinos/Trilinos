/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <stdlib.h>
#include <string>
#include <unistd.h>

#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <limits>
#include <string>
#include <sys/time.h>

#include <Kokkos_Core.hpp>

#include <KokkosKernels_IOUtils.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include <KokkosKernels_TestParameters.hpp>
#include <KokkosGraph_Distance2Color.hpp>
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosGraph;

enum ColoringMode {
  MODE_D2_SYMMETRIC,
  MODE_BIPARTITE_ROWS,
  MODE_BIPARTITE_COLS
};

struct D2Parameters {
  GraphColoringAlgorithmDistance2 algorithm;
  int repeat;
  int verbose;
  int use_threads;
  int use_openmp;
  int use_cuda;
  int use_hip;
  int use_serial;
  const char* mtx_file;
  ColoringMode d2_color_type;

  D2Parameters() {
    algorithm     = COLORING_D2_DEFAULT;
    repeat        = 6;
    verbose       = 0;
    use_threads   = 0;
    use_openmp    = 0;
    use_cuda      = 0;
    use_hip       = 0;
    use_serial    = 0;
    mtx_file      = NULL;
    d2_color_type = MODE_D2_SYMMETRIC;
  }
};

typedef default_scalar kk_scalar_t;
typedef default_size_type kk_size_type;
typedef default_lno_t kk_lno_t;

using namespace KokkosGraph;

void print_options(std::ostream& os, const char* app_name,
                   unsigned int indent = 0) {
  std::string spaces(indent, ' ');
  os << "Usage:" << std::endl
     << spaces << "  " << app_name << " [parameters]" << std::endl
     << std::endl
     << spaces << "Parameters:" << std::endl
     << spaces << "  Required Parameters:" << std::endl
     << spaces
     << "      --amtx <filename>   Input file in Matrix Market format (.mtx)."
     << std::endl
     << std::endl
     << spaces << "      Device type (the following are enabled in this build):"
     << std::endl
#ifdef KOKKOS_ENABLE_SERIAL
     << spaces << "          --serial            Execute serially." << std::endl
#endif
#ifdef KOKKOS_ENABLE_THREADS
     << spaces << "          --threads <N>       Use N posix threads."
     << std::endl
#endif
#ifdef KOKKOS_ENABLE_OPENMP
     << spaces << "          --openmp <N>        Use OpenMP with N threads."
     << std::endl
#endif
#ifdef KOKKOS_ENABLE_CUDA
     << spaces << "          --cuda <device id>  Use given CUDA device"
     << std::endl
#endif
#ifdef KOKKOS_ENABLE_HIP
     << spaces << "          --hip <device id>  Use given HIP device"
     << std::endl
#endif
     << std::endl
     << spaces << "  Coloring modes:" << std::endl
     << spaces
     << "      --symmetric_d2  (default): distance-2 on undirected/symmmetric "
        "graph"
     << std::endl
     << spaces
     << "      --bipartite_rows: color rows (left side of bipartite graph)"
     << std::endl
     << spaces
     << "      --bipartite_cols: color columns (right side of bipartite graph)"
     << std::endl
     << std::endl
     << spaces << "  Optional Parameters:" << std::endl
     << spaces
     << "      --algorithm <algorithm_name>   Set the algorithm to use.  "
        "Allowable values are:"
     << std::endl
     << spaces
     << "          COLORING_D2_SERIAL          - Serial net-based algorithm"
     << std::endl
     << spaces
     << "          COLORING_D2_VB              - Vertex Based method using "
        "boolean forbidden array (Default)."
     << std::endl
     << spaces
     << "          COLORING_D2_VB_BIT          - VB with Bitvector Forbidden "
        "Array"
     << std::endl
     << spaces
     << "          COLORING_D2_VB_BIT_EF       - VB_BIT with Edge Filtering"
     << std::endl
     << spaces
     << "          COLORING_D2_NB_BIT          - Net-based (fastest parallel "
        "algorithm)"
     << std::endl
     << spaces
     << "      --repeat <N>        Set number of test repetitions (Default: 1) "
     << std::endl
     << spaces
     << "      --verbose           Enable verbose mode (record and print "
        "timing + extra information)"
     << std::endl
     << spaces << "      --help              Print out command line help."
     << std::endl
     << spaces << " " << std::endl;
}

static char* getNextArg(int& i, int argc, char** argv) {
  i++;
  if (i >= argc) {
    std::cerr << "Error: expected additional command-line argument!\n";
    exit(1);
  }
  return argv[i];
}

int parse_inputs(D2Parameters& params, int argc, char** argv) {
  bool got_required_param_amtx = false;
  for (int i = 1; i < argc; ++i) {
    if (0 == Test::string_compare_no_case(argv[i], "--threads")) {
      params.use_threads = atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--serial")) {
      params.use_serial = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--openmp")) {
      params.use_openmp = atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--cuda")) {
      params.use_cuda = 1 + atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--hip")) {
      params.use_hip = 1 + atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--repeat")) {
      params.repeat = atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--amtx")) {
      got_required_param_amtx = true;
      params.mtx_file         = getNextArg(i, argc, argv);
    } else if (0 == Test::string_compare_no_case(argv[i], "--verbose")) {
      params.verbose = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--algorithm")) {
      ++i;
      if (0 == Test::string_compare_no_case(argv[i], "COLORING_D2_SERIAL")) {
        params.algorithm = COLORING_D2_SERIAL;
      } else if (0 == Test::string_compare_no_case(argv[i], "COLORING_D2_VB")) {
        params.algorithm = COLORING_D2_VB;
      } else if (0 ==
                 Test::string_compare_no_case(argv[i], "COLORING_D2_VB_BIT")) {
        params.algorithm = COLORING_D2_VB_BIT;
      } else if (0 == Test::string_compare_no_case(argv[i],
                                                   "COLORING_D2_VB_BIT_EF")) {
        params.algorithm = COLORING_D2_VB_BIT_EF;
      } else if (0 ==
                 Test::string_compare_no_case(argv[i], "COLORING_D2_NB_BIT")) {
        params.algorithm = COLORING_D2_NB_BIT;
      } else {
        std::cerr << "2-Unrecognized command line argument #" << i << ": "
                  << argv[i] << std::endl;
        print_options(std::cout, argv[0]);
        return 1;
      }
    } else if (0 == Test::string_compare_no_case(argv[i], "--symmetric_d2")) {
      params.d2_color_type = MODE_D2_SYMMETRIC;
    } else if (0 == Test::string_compare_no_case(argv[i], "--bipartite_rows")) {
      params.d2_color_type = MODE_BIPARTITE_ROWS;
    } else if (0 == Test::string_compare_no_case(argv[i], "--bipartite_cols")) {
      params.d2_color_type = MODE_BIPARTITE_COLS;
    } else if (0 == Test::string_compare_no_case(argv[i], "--help") ||
               0 == Test::string_compare_no_case(argv[i], "-h")) {
      print_options(std::cout, argv[0]);
      return 1;
    } else {
      std::cerr << "3-Unrecognized command line argument #" << i << ": "
                << argv[i] << std::endl;
      print_options(std::cout, argv[0]);
      return 1;
    }
  }

  if (!got_required_param_amtx) {
    std::cout << "Missing required parameter amtx" << std::endl << std::endl;
    print_options(std::cout, argv[0]);
    return 1;
  }
  if (!params.use_serial && !params.use_threads && !params.use_openmp &&
      !params.use_cuda && !params.use_hip) {
    print_options(std::cout, argv[0]);
    return 1;
  }
  return 0;
}

namespace KokkosKernels {
namespace Experiment {

std::string getCurrentDateTimeStr() {
  // Note: This could be replaced with `std::put_time(&tm, "%FT%T%z")` but
  // std::put_time isn't
  //       supported on the intel C++ compilers as of v. 17.0.x
  time_t now = time(0);
  char output[100];
  std::strftime(output, sizeof(output), "%FT%T%Z", std::localtime(&now));
  return output;
}

template <typename crsGraph_t>
void run_experiment(crsGraph_t crsGraph, int num_cols,
                    const D2Parameters& params) {
  using namespace KokkosGraph;
  using namespace KokkosGraph::Experimental;

  using device_t       = typename crsGraph_t::device_type;
  using exec_space     = typename device_t::execution_space;
  using mem_space      = typename device_t::memory_space;
  using lno_view_t     = typename crsGraph_t::row_map_type::non_const_type;
  using lno_nnz_view_t = typename crsGraph_t::entries_type::non_const_type;
  using size_type      = typename lno_view_t::non_const_value_type;
  using lno_t          = typename lno_nnz_view_t::non_const_value_type;

  int repeat = params.repeat;

  int verbose = params.verbose;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle<
      size_type, lno_t, kk_scalar_t, exec_space, mem_space, mem_space>
      KernelHandle;

  std::cout << "Num verts: " << crsGraph.numRows() << std::endl
            << "Num edges: " << crsGraph.entries.extent(0) << std::endl;

  KernelHandle kh;

  if (verbose) {
    kh.set_verbose(true);
  }

  // accumulators for average stats
  size_t total_colors = 0;
  size_t total_phases = 0;

  kh.create_distance2_graph_coloring_handle(params.algorithm);
  std::string label_algorithm =
      kh.get_distance2_graph_coloring_handle()->getD2AlgorithmName();
  std::cout << std::endl
            << "Run Graph Color D2 (" << label_algorithm << ")" << std::endl;

  // If any of the runs have an invalid result, this will be set to false.
  bool all_results_valid = true;

  // Loop over # of experiments to run
  for (int i = 0; i < repeat; ++i) {
    switch (params.d2_color_type) {
      case MODE_D2_SYMMETRIC:
        graph_color_distance2(&kh, crsGraph.numRows(), crsGraph.row_map,
                              crsGraph.entries);
        break;
      case MODE_BIPARTITE_ROWS:
        bipartite_color_rows(&kh, crsGraph.numRows(), num_cols,
                             crsGraph.row_map, crsGraph.entries);
        break;
      case MODE_BIPARTITE_COLS:
        bipartite_color_columns(&kh, crsGraph.numRows(), num_cols,
                                crsGraph.row_map, crsGraph.entries);
        break;
    }
    total_colors += kh.get_distance2_graph_coloring_handle()->get_num_colors();
    total_phases += kh.get_distance2_graph_coloring_handle()->get_num_phases();

    std::cout
        << "Total Time: "
        << kh.get_distance2_graph_coloring_handle()->get_overall_coloring_time()
        << std::endl
        << "Num colors: "
        << kh.get_distance2_graph_coloring_handle()->get_num_colors()
        << std::endl
        << "Num Phases: "
        << kh.get_distance2_graph_coloring_handle()->get_num_phases()
        << std::endl;

    std::cout << "\t";
    KokkosKernels::Impl::print_1Dview(
        kh.get_distance2_graph_coloring_handle()->get_vertex_colors());
    std::cout << std::endl;

    // If verbose mode is on and there the graph has fewer than 1500 verts, dump
    // a GraphVIZ DOT file.
    if (verbose && repeat == i + 1 && crsGraph.numRows() < 1500) {
      auto colors =
          kh.get_distance2_graph_coloring_handle()->get_vertex_colors();
      std::ofstream os("G.dot", std::ofstream::out);
      kh.get_distance2_graph_coloring_handle()->dump_graphviz(
          os, crsGraph.numRows(), crsGraph.row_map, crsGraph.entries, colors);
    }

    // ------------------------------------------
    // Verify correctness
    // ------------------------------------------
    // TODO bmk: write a faster color verification
    /*
    bool d2_coloring_is_valid            = false;
    bool d2_coloring_validation_flags[4] = { false };

    d2_coloring_is_valid = KokkosGraph::Impl::graph_verify_distance2_color(&kh,
    crsGraph.numRows(), num_cols, crsGraph.row_map, crsGraph.entries,
    crsGraph.row_map, crsGraph.entries, d2_coloring_validation_flags);

    // Print out messages based on coloring validation check.
    if(d2_coloring_is_valid)
    {
        std::cout << std::endl << "Distance-2 Graph Coloring is VALID" <<
    std::endl << std::endl;
    }
    else
    {
        all_results_valid = false;
        std::cout << std::endl
                  << "Distance-2 Graph Coloring is NOT VALID" << std::endl
                  << "  - Vert(s) left uncolored : " <<
    d2_coloring_validation_flags[1] << std::endl
                  << "  - Invalid D2 Coloring    : " <<
    d2_coloring_validation_flags[2] << std::endl
                  << std::endl;
    }
    if(d2_coloring_validation_flags[3])
    {
        std::cout << "Distance-2 Graph Coloring may have poor quality." <<
    std::endl
                  << "  - Vert(s) have high color value : " <<
    d2_coloring_validation_flags[3] << std::endl
                  << std::endl;
    }
    */

    // ------------------------------------------
    // Print out the colors histogram
    // ------------------------------------------
    KokkosGraph::Impl::graph_print_distance2_color_histogram(&kh, false);

  }  // for i...

  // ------------------------------------------
  // Compute Distance 2 Degree Stats
  // ------------------------------------------
  std::cout << "Compute Distance-2 Degree " << std::endl;

  Kokkos::Timer timer;

  double total_time =
      kh.get_distance2_graph_coloring_handle()->get_overall_coloring_time();
  double total_time_color_greedy = kh.get_distance2_graph_coloring_handle()
                                       ->get_overall_coloring_time_phase1();
  double total_time_find_conflicts = kh.get_distance2_graph_coloring_handle()
                                         ->get_overall_coloring_time_phase2();
  double total_time_resolve_conflicts =
      kh.get_distance2_graph_coloring_handle()
          ->get_overall_coloring_time_phase3();
  double total_time_matrix_squared = kh.get_distance2_graph_coloring_handle()
                                         ->get_overall_coloring_time_phase4();
  double total_time_matrix_squared_d1 =
      kh.get_distance2_graph_coloring_handle()
          ->get_overall_coloring_time_phase5();

  double avg_time                = total_time / (double)repeat;
  double avg_time_color_greedy   = total_time_color_greedy / (double)repeat;
  double avg_time_find_conflicts = total_time_find_conflicts / (double)repeat;
  double avg_time_resolve_conflicts =
      total_time_resolve_conflicts / (double)repeat;
  double avg_colors              = total_colors / (double)repeat;
  double avg_phases              = total_phases / (double)repeat;
  double avg_time_matrix_squared = total_time_matrix_squared / (double)repeat;
  double avg_time_matrix_squared_d1 =
      total_time_matrix_squared_d1 / (double)repeat;

  std::string short_mtx_file(params.mtx_file);
  short_mtx_file =
      short_mtx_file.substr(short_mtx_file.find_last_of("/\\") + 1);

  int result;
  char hostname[100];
  char username[100];

  result = gethostname(hostname, 100);
  if (result) {
    perror("gethostname");
  }

  result = getlogin_r(username, 100);
  if (result) {
    perror("getlogin_r");
  }

  std::string all_results_valid_str = "PASSED";
  if (!all_results_valid) all_results_valid_str = "FAILED";

  std::string currentDateTimeStr = getCurrentDateTimeStr();

  std::cout << "Summary" << std::endl
            << "-------" << std::endl
            << "    Date/Time      : " << currentDateTimeStr << std::endl
            << "    KExecSName     : " << Kokkos::DefaultExecutionSpace::name()
            << std::endl
            << "    Filename       : " << short_mtx_file << std::endl
            << "    Num Verts      : " << crsGraph.numRows() << std::endl
            << "    Num Edges      : " << crsGraph.entries.extent(0)
            << std::endl
            << "    Concurrency    : "
            << Kokkos::DefaultExecutionSpace::concurrency() << std::endl
            << "    Algorithm      : " << label_algorithm << std::endl
            << "Overall Time/Stats" << std::endl
            << "    Total Time     : " << total_time << std::endl
            << "    Avg Time       : " << avg_time << std::endl
            << "VB Distance[1|2] Stats " << std::endl
            << "    Avg Time CG    : " << avg_time_color_greedy << std::endl
            << "    Avg Time FC    : " << avg_time_find_conflicts << std::endl
            << "    Avg Time RC    : " << avg_time_resolve_conflicts
            << std::endl
            << "Matrix-Squared + D1 Stats" << std::endl
            << "    Avg Time to M^2: " << avg_time_matrix_squared << std::endl
            << "    Avg Time to D1 : " << avg_time_matrix_squared_d1
            << std::endl
            << "Coloring Stats" << std::endl
            << "    Avg colors     : " << avg_colors << std::endl
            << "    Avg Phases     : " << avg_phases << std::endl
            << "    Validation     : " << all_results_valid_str << std::endl
            << std::endl;

  std::cout << "CSVTIMEHDR"
            << ","
            << "Filename"
            << ","
            << "Host"
            << ","
            << "DateTime"
            << ","
            << "Num Rows"
            << ","
            << "Num Edges"
            << ","
            << "Execution Space"
            << ","
            << "Algorithm"
            << ","
            << "Concurrency"
            << ","
            << "Repetitions"
            << ","
            << "Total Time"
            << ","
            << "Total Time to M^2"
            << ","
            << "Total Time D1(M^2)"
            << ","
            << "Total Time CG"
            << ","
            << "Total Time FC"
            << ","
            << "Total Time RC"
            << ","
            << "Avg Colors"
            << ","
            << "Avg Num Phases"
            << ","
            << "Validation" << std::endl;

  std::cout << "CSVTIMEDATA"
            << "," << short_mtx_file << "," << hostname << ","
            << currentDateTimeStr << "," << crsGraph.numRows() << ","
            << crsGraph.entries.extent(0) << ","
            << Kokkos::DefaultExecutionSpace::name() << "," << label_algorithm
            << "," << Kokkos::DefaultExecutionSpace::concurrency() << ","
            << repeat << "," << total_time << "," << total_time_matrix_squared
            << "," << total_time_matrix_squared_d1 << ","
            << total_time_color_greedy << "," << total_time_find_conflicts
            << "," << total_time_resolve_conflicts

            << "," << avg_colors << "," << avg_phases << ","
            << all_results_valid_str << std::endl;

  std::cout << "CSVHISTHDR"
            << ","
            << "Filename"
            << ","
            << "Host"
            << ","
            << "DateTime"
            << ","
            << "Num Rows"
            << ","
            << "Num Edges"
            << ","
            << "Execution Space"
            << ","
            << "Algorithm"
            << ","
            << "Concurrency"
            << ","
            << "Histogram: 1 .. N" << std::endl;

  std::cout << "CSVHISTDATA"
            << "," << short_mtx_file << "," << hostname << ","
            << currentDateTimeStr << "," << crsGraph.numRows() << ","
            << crsGraph.entries.extent(0) << ","
            << Kokkos::DefaultExecutionSpace::name() << "," << label_algorithm
            << "," << Kokkos::DefaultExecutionSpace::concurrency() << ",";
  KokkosGraph::Impl::graph_print_distance2_color_histogram(&kh, true);
  std::cout << std::endl;

  // Kokkos::print_configuration(std::cout);
}

template <typename size_type, typename lno_t, typename exec_space,
          typename mem_space>
void experiment_driver(const D2Parameters& params) {
  using device_t = Kokkos::Device<exec_space, mem_space>;
  using crsMat_t = typename KokkosSparse::CrsMatrix<double, lno_t, device_t,
                                                    void, size_type>;
  using graph_t  = typename crsMat_t::StaticCrsGraphType;

  crsMat_t A =
      KokkosKernels::Impl::read_kokkos_crst_matrix<crsMat_t>(params.mtx_file);
  graph_t Agraph = A.graph;
  int num_cols   = A.numCols();

  KokkosKernels::Experiment::run_experiment<graph_t>(Agraph, num_cols, params);
}

}  // namespace Experiment
}  // namespace KokkosKernels

int main(int argc, char* argv[]) {
  D2Parameters params;

  // Override default repeats (default is 6)
  params.repeat = 1;

  if (parse_inputs(params, argc, argv)) {
    return 1;
  }

  if (params.mtx_file == NULL) {
    std::cerr << "Provide a matrix file" << std::endl;
    return 0;
  }

  std::cout << "Sizeof(kk_lno_t) : " << sizeof(kk_lno_t) << std::endl
            << "Sizeof(size_type): " << sizeof(kk_size_type) << std::endl;

  const int num_threads =
      params.use_openmp;  // Assumption is that use_openmp variable is provided
                          // as number of threads
  int device_id = 0;
  if (params.use_cuda)
    device_id = params.use_cuda - 1;
  else if (params.use_hip)
    device_id = params.use_hip - 1;
  Kokkos::initialize(Kokkos::InitArguments(num_threads, -1, device_id));

  // Print out verbose information about the configuration of the run.
  // Kokkos::print_configuration(std::cout);

#if defined(KOKKOS_MULTI_MEM)
  const bool use_multi_mem = true;
// todo: Throw an error or print a message if KOKKOS_MULTI_MEM is enabled for
// this test?  (WCMCLEN--SCAFFOLDING)
#else
  const bool use_multi_mem = false;
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
  if (params.use_openmp) {
    if (!use_multi_mem) {
      KokkosKernels::Experiment::experiment_driver<
          kk_size_type, kk_lno_t, Kokkos::OpenMP, Kokkos::OpenMP::memory_space>(
          params);
    }
  }
#endif

#if defined(KOKKOS_ENABLE_THREADS)
  if (params.use_threads) {
    if (!use_multi_mem) {
      KokkosKernels::Experiment::experiment_driver<
          kk_size_type, kk_lno_t, Kokkos::Threads,
          Kokkos::Threads::memory_space>(params);
    }
  }
#endif

#if defined(KOKKOS_ENABLE_CUDA)
  if (params.use_cuda) {
    if (!use_multi_mem) {
      KokkosKernels::Experiment::experiment_driver<
          kk_size_type, kk_lno_t, Kokkos::Cuda, Kokkos::Cuda::memory_space>(
          params);
    }
  }
#endif

#if defined(KOKKOS_ENABLE_HIP)
  if (params.use_hip) {
    if (!use_multi_mem) {
      KokkosKernels::Experiment::experiment_driver<
          kk_size_type, kk_lno_t, Kokkos::Experimental::HIP,
          Kokkos::Experimental::HIPSpace>(params);
    }
  }
#endif

#if defined(KOKKOS_ENABLE_SERIAL)
  if (params.use_serial) {
    if (!use_multi_mem) {
      KokkosKernels::Experiment::experiment_driver<
          kk_size_type, kk_lno_t, Kokkos::Serial, Kokkos::Serial::memory_space>(
          params);
    }
  }
#endif

  Kokkos::finalize();

  return 0;
}
