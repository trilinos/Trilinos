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
#include <KokkosKernels_Handle.hpp>

#include <cstdlib>
#include <iostream>

#include <random>     // std::default_random_engine
#include <algorithm>  // std::shuffle
#include <vector>

#include "KokkosKernels_IOUtils.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_TestParameters.hpp"
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosKernels_TestUtils.hpp"
#include "KokkosSparse_IOUtils.hpp"

void print_options(std::ostream &os, const char *app_name, unsigned int indent = 0) {
  std::string spaces(indent, ' ');
  os << "Usage:" << std::endl
     << spaces << "  " << app_name << " [parameters]" << std::endl
     << std::endl
     << spaces << "Parameters:" << std::endl
     << spaces << "  Parallelism (select one of the following):" << std::endl
#if defined(KOKKOS_ENABLE_SERIAL)
     << spaces << "      --serial            Execute serially." << std::endl
#endif
#if defined(KOKKOS_ENABLE_THREADS)
     << spaces << "      --threads <N>       Use N posix threads." << std::endl
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
     << spaces << "      --openmp <N>        Use OpenMP with N threads." << std::endl
#endif
#if defined(KOKKOS_ENABLE_CUDA)
     << spaces << "      --cuda <id>         Use CUDA (device $id)" << std::endl
#endif
#if defined(KOKKOS_ENABLE_HIP)
     << spaces << "      --hip <id>          Use HIP (device $id)" << std::endl
#endif
     << std::endl
     << spaces << "  Required Parameters:" << std::endl
     << spaces << "      --amtx <filename>   Input file in Matrix Market format (.mtx)." << std::endl
     << std::endl
     << spaces
     << "      --algorithm <algorithm_name>   Set the algorithm to use.  "
        "Allowable values are:"
     << std::endl
     << spaces
     << "                 COLORING_DEFAULT  - Use the default coloring method, "
        "architecture dependent."
     << std::endl
     << spaces << "                 COLORING_SERIAL   - Use the serial algorithm." << std::endl
     << spaces
     << "                 COLORING_VB       - Use the parallel vertex-based "
        "method."
     << std::endl
     << spaces
     << "                 COLORING_VBBIT    - Use the parallel vertex-based "
        "with bit vectors method."
     << std::endl
     << spaces << "                 COLORING_EB       - Use edge based method." << std::endl
     << spaces
     << "                 COLORING_VBD      - Use the vertex-based "
        "deterministic method."
     << std::endl
     << spaces
     << "                 COLORING_VBDBIT   - Use the vertex-based "
        "deterministic with bit vectors method."
     << std::endl
     << std::endl
     << spaces << "  Optional Parameters:" << std::endl
     << spaces << "      --chunksize <N>     Set the chunk size." << std::endl
     << spaces << "      --dynamic           Use dynamic scheduling." << std::endl
     << spaces << "      --outputfile <FILE> Output the colors of the nodes to the file." << std::endl
     << spaces << "      --repeat <N>        Set number of test repetitions (Default: 1) " << std::endl
     << spaces << "      --teamsize  <N>     Set the team size." << std::endl
     << spaces << "      --vectorsize <N>    Set the vector size." << std::endl
     << spaces
     << "      --verbose           Enable verbose mode (record and print "
        "timing + extra information)"
     << std::endl
     << spaces << "      --help              Print out command line help." << std::endl
     << spaces << " " << std::endl;
}

static char *getNextArg(int &i, int argc, char **argv) {
  i++;
  if (i >= argc) {
    std::cerr << "Error: expected additional command-line argument!\n";
    exit(1);
  }
  return argv[i];
}

int parse_inputs(KokkosKernels::Experiment::Parameters &params, int argc, char **argv) {
  bool got_required_param_amtx      = false;
  bool got_required_param_algorithm = false;

  for (int i = 1; i < argc; ++i) {
    if (0 == Test::string_compare_no_case(argv[i], "--threads")) {
      params.use_threads = atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--serial")) {
      params.use_serial = atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--openmp")) {
      params.use_openmp = atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--cuda")) {
      params.use_cuda = 1 + atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--hip")) {
      params.use_hip = 1 + atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--repeat")) {
      params.repeat = atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--chunksize")) {
      params.chunk_size = atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--teamsize")) {
      params.team_size = atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--vectorsize")) {
      params.vector_size = atoi(getNextArg(i, argc, argv));
    } else if (0 == Test::string_compare_no_case(argv[i], "--amtx")) {
      got_required_param_amtx = true;
      params.a_mtx_bin_file   = getNextArg(i, argc, argv);
    } else if (0 == Test::string_compare_no_case(argv[i], "--dynamic")) {
      params.use_dynamic_scheduling = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--verbose")) {
      params.verbose = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--outputfile") ||
               0 == Test::string_compare_no_case(argv[i], "-o")) {
      params.coloring_output_file = getNextArg(i, argc, argv);
    } else if (0 == Test::string_compare_no_case(argv[i], "--algorithm")) {
      got_required_param_algorithm = true;
      ++i;
      if (0 == Test::string_compare_no_case(argv[i], "COLORING_DEFAULT")) {
        params.algorithm = 1;
      } else if (0 == Test::string_compare_no_case(argv[i], "COLORING_SERIAL")) {
        params.algorithm = 2;
      } else if (0 == Test::string_compare_no_case(argv[i], "COLORING_VB")) {
        params.algorithm = 3;
      } else if (0 == Test::string_compare_no_case(argv[i], "COLORING_VBBIT")) {
        params.algorithm = 4;
      } else if (0 == Test::string_compare_no_case(argv[i], "COLORING_VBCS")) {
        params.algorithm = 5;
      } else if (0 == Test::string_compare_no_case(argv[i], "COLORING_EB")) {
        params.algorithm = 6;
      } else if (0 == Test::string_compare_no_case(argv[i], "COLORING_VBD")) {
        params.algorithm = 7;
      } else if (0 == Test::string_compare_no_case(argv[i], "COLORING_VBDBIT")) {
        params.algorithm = 8;
      } else if (0 == Test::string_compare_no_case(argv[i], "--help") ||
                 0 == Test::string_compare_no_case(argv[i], "-h")) {
        print_options(std::cout, argv[0]);
        return 1;
      } else {
        std::cerr << "2-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
        print_options(std::cout, argv[0]);
        return 1;
      }
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
  if (!params.use_serial && !params.use_threads && !params.use_openmp && !params.use_cuda && !params.use_hip) {
    print_options(std::cout, argv[0]);
    return 1;
  }

  return 0;
}

using KokkosKernels::Impl::xorshiftHash;

template <typename lno_t, typename size_type, typename rowmap_t, typename entries_t>
bool verifySymmetric(lno_t numVerts, const rowmap_t &d_rowmap, const entries_t &d_entries) {
  auto rowmap  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d_rowmap);
  auto entries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d_entries);
  size_t hash  = 0;
  for (lno_t v = 0; v < numVerts; v++) {
    size_type rowBegin = rowmap(v);
    size_type rowEnd   = rowmap(v + 1);
    for (size_type i = rowBegin; i < rowEnd; i++) {
      lno_t nei = entries(i);
      if (nei < numVerts && nei != v) {
        hash ^= xorshiftHash<size_t>(xorshiftHash<size_t>(v) ^ xorshiftHash<size_t>(nei));
      }
    }
  }
  return hash == 0U;
}

template <typename lno_t, typename size_type, typename rowmap_t, typename entries_t, typename colors_t>
bool verifyColoring(lno_t numVerts, const rowmap_t &d_rowmap, const entries_t &d_entries, const colors_t &d_colors) {
  auto rowmap  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d_rowmap);
  auto entries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d_entries);
  auto colors  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d_colors);
  // Just do the simplest possible neighbors-of-neighbors loop to find conflicts
  for (lno_t v = 0; v < numVerts; v++) {
    if (colors(v) == 0) {
      std::cout << "Vertex " << v << " is uncolored.\n";
      return false;
    }
    size_type rowBegin = rowmap(v);
    size_type rowEnd   = rowmap(v + 1);
    for (size_type i = rowBegin; i < rowEnd; i++) {
      lno_t nei = entries(i);
      if (nei < numVerts && nei != v) {
        // check for dist-1 conflict
        if (colors(v) == colors(nei)) {
          std::cout << "Dist-1 conflict between " << v << " and " << nei << '\n';
          return false;
        }
      }
    }
  }
  return true;
}

namespace KokkosKernels {

namespace Experiment {

template <typename ExecSpace, typename crsGraph_t, typename crsGraph_t2, typename crsGraph_t3, typename TempMemSpace,
          typename PersistentMemSpace>
void run_experiment(crsGraph_t crsGraph, int num_cols, Parameters params) {
  // using namespace KokkosSparse;
  using namespace KokkosGraph;
  using namespace KokkosGraph::Experimental;
  // using namespace KokkosSparse::Experimental;

  int algorithm  = params.algorithm;
  int repeat     = params.repeat;
  int chunk_size = params.chunk_size;

  int shmemsize              = params.shmemsize;
  int team_size              = params.team_size;
  int use_dynamic_scheduling = params.use_dynamic_scheduling;
  int verbose                = params.verbose;

  // char spgemm_step = params.spgemm_step;
  int vector_size = params.vector_size;

  typedef typename crsGraph_t3::row_map_type::non_const_type lno_view_t;
  typedef typename crsGraph_t3::entries_type::non_const_type lno_nnz_view_t;

  typedef typename lno_view_t::non_const_value_type size_type;
  typedef typename lno_nnz_view_t::non_const_value_type lno_t;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, lno_t, ExecSpace, TempMemSpace,
                                                           PersistentMemSpace>
      KernelHandle;

  if (verbose) {
    if (verifySymmetric<lno_t, size_type, decltype(crsGraph.row_map), decltype(crsGraph.entries)>(
            crsGraph.numRows(), crsGraph.row_map, crsGraph.entries)) {
      std::cout << std::endl << "Graph is symmetric (valid input)" << std::endl;
    } else {
      std::cout << std::endl << "Graph is nonsymmetric (INVALID INPUT)" << std::endl;
      // Don't attempt coloring when input is invalid
      return;
    }
  }

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

  std::cout << "algorithm: " << algorithm << std::endl;

  double totalTime = 0.0;
  for (int i = 0; i < repeat; ++i) {
    switch (algorithm) {
      case 1: kh.create_graph_coloring_handle(COLORING_DEFAULT); break;
      case 2: kh.create_graph_coloring_handle(COLORING_SERIAL); break;
      case 3: kh.create_graph_coloring_handle(COLORING_VB); break;
      case 4: kh.create_graph_coloring_handle(COLORING_VBBIT); break;
      case 5: kh.create_graph_coloring_handle(COLORING_VBCS); break;
      case 6: kh.create_graph_coloring_handle(COLORING_EB); break;

      case 7: kh.create_graph_coloring_handle(COLORING_VBD); break;

      case 8: kh.create_graph_coloring_handle(COLORING_VBDBIT); break;

      default: kh.create_graph_coloring_handle(COLORING_DEFAULT);
    }

    graph_color_symbolic(&kh, crsGraph.numRows(), num_cols, crsGraph.row_map, crsGraph.entries);

    std::cout << std::endl
              << "Time:" << kh.get_graph_coloring_handle()->get_overall_coloring_time()
              << " sec. "
                 "Num colors:"
              << kh.get_graph_coloring_handle()->get_num_colors()
              << " "
                 "Num Phases:"
              << kh.get_graph_coloring_handle()->get_num_phases() << std::endl;
    std::cout << "\t";

    auto colors = kh.get_graph_coloring_handle()->get_vertex_colors();
    KokkosKernels::Impl::print_1Dview(colors);

    if (verbose) {
      if (verifyColoring<lno_t, size_type, decltype(crsGraph.row_map), decltype(crsGraph.entries), decltype(colors)>(
              crsGraph.numRows(), crsGraph.row_map, crsGraph.entries, colors)) {
        std::cout << std::endl << "Graph Coloring is VALID" << std::endl << std::endl;
      } else {
        std::cout << std::endl << "Graph Coloring is NOT VALID" << std::endl;
        break;
      }
    }

    if (params.coloring_output_file != "") {
      std::ofstream os(params.coloring_output_file, std::ofstream::out);
      KokkosKernels::Impl::print_1Dview(os, colors, true, "\n");
    }
    totalTime += kh.get_graph_coloring_handle()->get_overall_coloring_time();
  }
  std::cout << "Average time over " << repeat << " trials: " << totalTime / repeat << " sec.\n";
}

template <typename size_type, typename lno_t, typename exec_space, typename hbm_mem_space, typename sbm_mem_space>
void run_multi_mem_experiment(Parameters params) {
  typedef exec_space myExecSpace;
  typedef Kokkos::Device<exec_space, hbm_mem_space> myFastDevice;
  typedef Kokkos::Device<exec_space, sbm_mem_space> mySlowExecSpace;

  typedef typename KokkosSparse::CrsMatrix<double, lno_t, myFastDevice, void, size_type> fast_crstmat_t;
  typedef typename fast_crstmat_t::StaticCrsGraphType fast_graph_t;
  // typedef typename fast_graph_t::row_map_type::non_const_type
  // fast_row_map_view_t; typedef typename
  // fast_graph_t::entries_type::non_const_type   fast_cols_view_t;

  // typedef typename fast_graph_t::row_map_type::const_type
  // const_fast_row_map_view_t; typedef typename
  // fast_graph_t::entries_type::const_type   const_fast_cols_view_t;

  typedef typename KokkosSparse::CrsMatrix<double, lno_t, mySlowExecSpace, void, size_type> slow_crstmat_t;
  typedef typename slow_crstmat_t::StaticCrsGraphType slow_graph_t;

  // typedef typename slow_graph_t::row_map_type::non_const_type
  // slow_row_map_view_t; typedef typename
  // slow_graph_t::entries_type::non_const_type   slow_cols_view_t; typedef
  // typename slow_graph_t::row_map_type::const_type const_slow_row_map_view_t;
  // typedef typename slow_graph_t::entries_type::const_type
  // const_slow_cols_view_t;

  const char *a_mat_file = params.a_mtx_bin_file.c_str();
  // char *b_mat_file = params.b_mtx_bin_file;
  // char *c_mat_file = params.c_mtx_bin_file;

  slow_graph_t a_slow_crsgraph, /*b_slow_crsgraph,*/ c_slow_crsgraph;
  fast_graph_t a_fast_crsgraph, /*b_fast_crsgraph,*/ c_fast_crsgraph;

  int num_cols = 0;

  // read a and b matrices and store them on slow or fast memory.
  if (params.a_mem_space == 1) {
    fast_crstmat_t a_fast_crsmat;
    a_fast_crsmat   = KokkosSparse::Impl::read_kokkos_crst_matrix<fast_crstmat_t>(a_mat_file);
    a_fast_crsgraph = a_fast_crsmat.graph;
    num_cols        = a_fast_crsmat.numCols();

  } else {
    slow_crstmat_t a_slow_crsmat;
    a_slow_crsmat   = KokkosSparse::Impl::read_kokkos_crst_matrix<slow_crstmat_t>(a_mat_file);
    a_slow_crsgraph = a_slow_crsmat.graph;
    num_cols        = a_slow_crsmat.numCols();
  }

  if (params.a_mem_space == 1) {
    if (params.b_mem_space == 1) {
      if (params.c_mem_space == 1) {
        if (params.work_mem_space == 1) {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment<myExecSpace, fast_graph_t, fast_graph_t, fast_graph_t,
                                                    hbm_mem_space, hbm_mem_space>(a_fast_crsgraph, num_cols, params);
        } else {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment<myExecSpace, fast_graph_t, fast_graph_t, fast_graph_t,
                                                    sbm_mem_space, sbm_mem_space>(a_fast_crsgraph, num_cols, params);
        }

      } else {
        // C is in slow memory.
        if (params.work_mem_space == 1) {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment<myExecSpace, fast_graph_t, fast_graph_t, slow_graph_t,
                                                    hbm_mem_space, hbm_mem_space>(a_fast_crsgraph, num_cols, params);
        } else {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment<myExecSpace, fast_graph_t, fast_graph_t, slow_graph_t,
                                                    sbm_mem_space, sbm_mem_space>(a_fast_crsgraph, num_cols, params);
        }
      }
    } else {
      // B is in slow memory
      if (params.c_mem_space == 1) {
        if (params.work_mem_space == 1) {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment<myExecSpace, fast_graph_t, slow_graph_t, fast_graph_t,
                                                    hbm_mem_space, hbm_mem_space>(a_fast_crsgraph, num_cols, params);
        } else {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment<myExecSpace, fast_graph_t, slow_graph_t, fast_graph_t,
                                                    sbm_mem_space, sbm_mem_space>(a_fast_crsgraph, num_cols, params);
        }

      } else {
        // C is in slow memory.
        if (params.work_mem_space == 1) {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment<myExecSpace, fast_graph_t, slow_graph_t, slow_graph_t,
                                                    hbm_mem_space, hbm_mem_space>(a_fast_crsgraph, num_cols, params);
        } else {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment<myExecSpace, fast_graph_t, slow_graph_t, slow_graph_t,
                                                    sbm_mem_space, sbm_mem_space>(a_fast_crsgraph, num_cols, params);
        }
      }
    }
  } else {
    // A is in slow memory
    if (params.b_mem_space == 1) {
      if (params.c_mem_space == 1) {
        if (params.work_mem_space == 1) {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment<myExecSpace, slow_graph_t, fast_graph_t, fast_graph_t,
                                                    hbm_mem_space, hbm_mem_space>(a_slow_crsgraph, num_cols, params);
        } else {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment<myExecSpace, slow_graph_t, fast_graph_t, fast_graph_t,
                                                    sbm_mem_space, sbm_mem_space>(a_slow_crsgraph, num_cols, params);
        }

      } else {
        // C is in slow memory.
        if (params.work_mem_space == 1) {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment<myExecSpace, slow_graph_t, fast_graph_t, slow_graph_t,
                                                    hbm_mem_space, hbm_mem_space>(a_slow_crsgraph, num_cols, params);
        } else {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment<myExecSpace, slow_graph_t, fast_graph_t, slow_graph_t,
                                                    sbm_mem_space, sbm_mem_space>(a_slow_crsgraph, num_cols, params);
        }
      }
    } else {
      // B is in slow memory
      if (params.c_mem_space == 1) {
        if (params.work_mem_space == 1) {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment<myExecSpace, slow_graph_t, slow_graph_t, fast_graph_t,
                                                    hbm_mem_space, hbm_mem_space>(a_slow_crsgraph, num_cols, params);
        } else {
          /* c_fast_crsgraph = */
          KokkosKernels::Experiment::run_experiment<myExecSpace, slow_graph_t, slow_graph_t, fast_graph_t,
                                                    sbm_mem_space, sbm_mem_space>(a_slow_crsgraph, num_cols, params);
        }

      } else {
        // C is in slow memory.
        if (params.work_mem_space == 1) {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment<myExecSpace, slow_graph_t, slow_graph_t, slow_graph_t,
                                                    hbm_mem_space, hbm_mem_space>(a_slow_crsgraph, num_cols, params);
        } else {
          /*c_slow_crsgraph =*/
          KokkosKernels::Experiment::run_experiment<myExecSpace, slow_graph_t, slow_graph_t, slow_graph_t,
                                                    sbm_mem_space, sbm_mem_space>(a_slow_crsgraph, num_cols, params);
        }
      }
    }
  }
}

}  // namespace Experiment
}  // namespace KokkosKernels

int main(int argc, char **argv) {
  typedef unsigned size_type;
  typedef int idx;
  // typedef int size_type;
  // typedef int idx;

  KokkosKernels::Experiment::Parameters params;

  if (parse_inputs(params, argc, argv)) {
    return 1;
  }
  if (params.a_mtx_bin_file == "") {
    std::cerr << "Provide a matrix file" << std::endl;
    return 0;
  }
  std::cout << "Sizeof(idx):" << sizeof(idx) << " sizeof(size_type):" << sizeof(size_type) << std::endl;

  const int num_threads = params.use_openmp;  // Assumption is that use_openmp variable is provided
                                              // as number of threads
  const int device_id = std::max(params.use_cuda, params.use_hip) - 1;
  Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(num_threads).set_device_id(device_id));
  Kokkos::print_configuration(std::cout);

#if defined(KOKKOS_ENABLE_OPENMP)

  if (params.use_openmp) {
#ifdef KOKKOSKERNELS_MULTI_MEM
    KokkosKernels::Experiment::run_multi_mem_experiment<size_type, idx, Kokkos::OpenMP, Kokkos::OpenMP::memory_space,
                                                        Kokkos::HostSpace>(params);
#else

    KokkosKernels::Experiment::run_multi_mem_experiment<size_type, idx, Kokkos::OpenMP, Kokkos::OpenMP::memory_space,
                                                        Kokkos::OpenMP::memory_space>(params);
#endif
  }
#endif

#if defined(KOKKOS_ENABLE_CUDA)
  if (params.use_cuda) {
#ifdef KOKKOSKERNELS_MULTI_MEM
    KokkosKernels::Experiment::run_multi_mem_experiment<size_type, idx, Kokkos::Cuda, Kokkos::Cuda::memory_space,
                                                        Kokkos::CudaHostPinnedSpace>(params);
#else
    KokkosKernels::Experiment::run_multi_mem_experiment<size_type, idx, Kokkos::Cuda, Kokkos::Cuda::memory_space,
                                                        Kokkos::Cuda::memory_space>(params);

#endif
  }

#endif

#if defined(KOKKOS_ENABLE_HIP)
  if (params.use_hip) {
    KokkosKernels::Experiment::run_multi_mem_experiment<size_type, idx, Kokkos::HIP, Kokkos::HIPSpace,
                                                        Kokkos::HIPSpace>(params);
  }
#endif

#if defined(KOKKOS_ENABLE_SERIAL)
  if (params.use_serial) {
#ifdef KOKKOSKERNELS_MULTI_MEM
    KokkosKernels::Experiment::run_multi_mem_experiment<size_type, idx, Kokkos::Serial, Kokkos::Serial::memory_space,
                                                        Kokkos::HostSpace>(params);
#else

    KokkosKernels::Experiment::run_multi_mem_experiment<size_type, idx, Kokkos::Serial, Kokkos::Serial::memory_space,
                                                        Kokkos::Serial::memory_space>(params);
#endif
  }
#endif

  Kokkos::finalize();

  return 0;
}
