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

#include <stdlib.h>
#include <string>
#include <set>
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

#include "KokkosKernels_Utils.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spadd.hpp"
#include "KokkosGraph_MIS2.hpp"
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_TestUtils.hpp"
#include "KokkosSparse_IOUtils.hpp"

using namespace KokkosGraph;

struct MIS2Parameters {
  int repeat           = 1;
  bool verbose         = false;
  int use_threads      = 0;
  int use_openmp       = 0;
  int use_cuda         = 0;
  int use_hip          = 0;
  int use_serial       = 0;
  const char* mtx_file = NULL;
  MIS2_Algorithm algo  = MIS2_FAST;
};

template <typename lno_t, typename size_type, typename rowmap_t, typename entries_t, typename mis_t>
bool verifyD2MIS(lno_t numVerts, const rowmap_t& rowmap, const entries_t& entries, const mis_t& misArray) {
  // set a std::set of the mis, for fast membership test
  std::set<lno_t> mis;
  for (size_t i = 0; i < misArray.extent(0); i++) mis.insert(misArray(i));
  for (lno_t i = 0; i < numVerts; i++) {
    // determine whether another vertex in the set is
    // within 2 hops of i.
    bool misIn2Hops = false;
    for (size_type j = rowmap(i); j < rowmap(i + 1); j++) {
      lno_t nei1 = entries(j);
      if (nei1 == i || nei1 >= numVerts) continue;
      if (mis.find(nei1) != mis.end()) {
        misIn2Hops = true;
        break;
      }
      for (size_type k = rowmap(nei1); k < rowmap(nei1 + 1); k++) {
        lno_t nei2 = entries(k);
        if (nei2 == i || nei2 >= numVerts) continue;
        if (mis.find(nei2) != mis.end()) {
          misIn2Hops = true;
          break;
        }
      }
    }
    if (mis.find(i) == mis.end()) {
      // i is not in the set
      if (!misIn2Hops) {
        std::cout << "INVALID D2 MIS: vertex " << i << " is not in the set,\n";
        std::cout << "but there are no vertices in the set within 2 hops.\n";
        return false;
      }
    } else {
      // i is in the set
      if (misIn2Hops) {
        std::cout << "INVALID D2 MIS: vertex " << i << " is in the set,\n";
        std::cout << "but there is another vertex within 2 hops which is also "
                     "in the set.\n";
        return false;
      }
    }
  }
  return true;
}

void print_options(std::ostream& os, const char* app_name, unsigned int indent = 0) {
  std::string spaces(indent, ' ');
  os << "Usage:" << std::endl
     << spaces << "  " << app_name << " [parameters]" << std::endl
     << std::endl
     << spaces << "Parameters:" << std::endl
     << spaces << "  Required Parameters:" << std::endl
     << spaces << "      --amtx <filename>   Input file in Matrix Market format (.mtx)." << std::endl
     << std::endl
     << spaces << "      Device type (the following are enabled in this build):" << std::endl
#ifdef KOKKOS_ENABLE_SERIAL
     << spaces << "          --serial            Execute serially." << std::endl
#endif
#ifdef KOKKOS_ENABLE_THREADS
     << spaces << "          --threads           Use posix threads.\n"
#endif
#ifdef KOKKOS_ENABLE_OPENMP
     << spaces << "          --openmp            Use OpenMP.\n"
#endif
#ifdef KOKKOS_ENABLE_CUDA
     << spaces << "          --cuda              Use CUDA.\n"
#endif
#ifdef KOKKOS_ENABLE_HIP
     << spaces << "          --hip               Use HIP.\n"
#endif
     << std::endl
     << spaces << "  Optional Parameters:" << std::endl
     << spaces << "      --algo alg          alg: fast, quality" << std::endl
     << spaces << "      --repeat <N>        Set number of test repetitions (Default: 1) " << std::endl
     << spaces
     << "      --verbose           Enable verbose mode (record and print "
        "timing + extra information)"
     << std::endl
     << spaces << "      --help              Print out command line help." << std::endl
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

int parse_inputs(MIS2Parameters& params, int argc, char** argv) {
  bool got_required_param_amtx = false;
  for (int i = 1; i < argc; ++i) {
    if (0 == Test::string_compare_no_case(argv[i], "--threads")) {
      params.use_threads = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--serial")) {
      params.use_serial = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--openmp")) {
      params.use_openmp = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--cuda")) {
      params.use_cuda = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--hip")) {
      params.use_hip = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--repeat")) {
      params.repeat = atoi(getNextArg(i, argc, argv));
      if (params.repeat <= 0) {
        std::cout << "*** Repeat count must be positive, defaulting to 1.\n";
        params.repeat = 1;
      }
    } else if (0 == Test::string_compare_no_case(argv[i], "--amtx")) {
      got_required_param_amtx = true;
      params.mtx_file         = getNextArg(i, argc, argv);
    } else if (0 == Test::string_compare_no_case(argv[i], "--algo")) {
      const char* algName = getNextArg(i, argc, argv);
      if (!Test::string_compare_no_case(algName, "fast"))
        params.algo = MIS2_FAST;
      else if (!Test::string_compare_no_case(algName, "quality"))
        params.algo = MIS2_QUALITY;
      else
        throw std::invalid_argument("Algorithm not valid: must be 'fast' or 'quality'");
    } else if (0 == Test::string_compare_no_case(argv[i], "--verbose")) {
      params.verbose = true;
    } else if (0 == Test::string_compare_no_case(argv[i], "--help") ||
               0 == Test::string_compare_no_case(argv[i], "-h")) {
      print_options(std::cout, argv[0]);
      return 1;
    } else {
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
      print_options(std::cout, argv[0]);
      return 1;
    }
  }

  if (!got_required_param_amtx) {
    std::cout << "Missing required parameter amtx" << std::endl << std::endl;
    print_options(std::cout, argv[0]);
    return 1;
  }
  if (!params.use_serial && !params.use_threads && !params.use_openmp && !params.use_cuda && !params.use_hip) {
    print_options(std::cout, argv[0]);
    return 1;
  }
  return 0;
}

template <typename device_t>
void run_mis2(const MIS2Parameters& params) {
  using size_type  = default_size_type;
  using lno_t      = default_lno_t;
  using exec_space = typename device_t::execution_space;
  using mem_space  = typename device_t::memory_space;
  using crsMat_t   = typename KokkosSparse::CrsMatrix<default_scalar, lno_t, device_t, void, size_type>;
  using lno_view_t = typename crsMat_t::index_type::non_const_type;
  using KKH = KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, default_scalar, exec_space, mem_space,
                                                               mem_space>;

  Kokkos::Timer t;
  crsMat_t A_in = KokkosSparse::Impl::read_kokkos_crst_matrix<crsMat_t>(params.mtx_file);
  std::cout << "I/O time: " << t.seconds() << " s\n";
  t.reset();
  // Symmetrize the matrix just in case
  crsMat_t At_in = KokkosSparse::Impl::transpose_matrix(A_in);
  crsMat_t A;
  KKH kkh;
  const default_scalar one = Kokkos::ArithTraits<default_scalar>::one();
  kkh.create_spadd_handle(false);
  KokkosSparse::spadd_symbolic(&kkh, A_in, At_in, A);
  KokkosSparse::spadd_numeric(&kkh, one, A_in, one, At_in, A);
  kkh.destroy_spadd_handle();
  std::cout << "Time to symmetrize: " << t.seconds() << " s\n";
  auto rowmap    = A.graph.row_map;
  auto entries   = A.graph.entries;
  lno_t numVerts = A.numRows();

  std::cout << "Num verts: " << numVerts << '\n' << "Num edges: " << A.nnz() << '\n';

  lno_view_t mis;

  t.reset();
  for (int rep = 0; rep < params.repeat; rep++) {
    mis = KokkosGraph::graph_d2_mis<device_t, decltype(rowmap), decltype(entries)>(rowmap, entries, params.algo);
    exec_space().fence();
  }
  double totalTime = t.seconds();
  std::cout << "MIS-2 average time: " << totalTime / params.repeat << '\n';
  std::cout << "MIS size: " << mis.extent(0) << '\n';

  if (params.verbose) {
    std::cout << "Vertices in independent set:\n";
    KokkosKernels::Impl::print_1Dview(mis);
    auto rowmapHost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowmap);
    auto entriesHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), entries);
    auto misHost     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), mis);
    if (verifyD2MIS<lno_t, size_type, decltype(rowmapHost), decltype(entriesHost), decltype(misHost)>(
            numVerts, rowmapHost, entriesHost, misHost))
      std::cout << "MIS-2 is correct.\n";
    else
      std::cout << "*** MIS-2 not correct! ***\n";
  }
}

int main(int argc, char* argv[]) {
  MIS2Parameters params;

  if (parse_inputs(params, argc, argv)) {
    return 1;
  }

  if (params.mtx_file == NULL) {
    std::cerr << "Provide a matrix file" << std::endl;
    return 0;
  }

  Kokkos::initialize();

  bool run = false;

#if defined(KOKKOS_ENABLE_OPENMP)
  if (params.use_openmp) {
    run_mis2<Kokkos::OpenMP>(params);
    run = true;
  }
#endif

#if defined(KOKKOS_ENABLE_THREADS)
  if (params.use_threads) {
    run_mis2<Kokkos::Threads>(params);
    run = true;
  }
#endif

#if defined(KOKKOS_ENABLE_CUDA)
  if (params.use_cuda) {
    run_mis2<Kokkos::Cuda>(params);
    run = true;
  }
#endif

#if defined(KOKKOS_ENABLE_HIP)
  if (params.use_hip) {
    run_mis2<Kokkos::HIP>(params);
    run = true;
  }
#endif

#if defined(KOKKOS_ENABLE_SERIAL)
  if (params.use_serial) {
    run_mis2<Kokkos::Serial>(params);
    run = true;
  }
#endif

  if (!run) {
    std::cerr << "*** ERROR: did not run, none of the supported device types "
                 "were selected.\n";
  }

  Kokkos::finalize();

  return 0;
}
