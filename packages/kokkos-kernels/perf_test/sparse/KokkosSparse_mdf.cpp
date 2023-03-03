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
#include "KokkosSparse_mdf.hpp"
#include "KokkosKernels_TestUtils.hpp"

struct Params {
  int use_cuda    = 0;
  int use_hip     = 0;
  int use_sycl    = 0;
  int use_openmp  = 0;
  int use_threads = 0;
  std::string amtx;
  int m         = 10000;
  int n         = 10000;
  int nnzPerRow = 30;
  bool diag = false;  // Whether B should be diagonal only (requires A square)
  bool verbose = false;
  int repeat   = 1;
};

template <class row_map_t, class entries_t>
struct diag_generator_functor {
  using size_type = typename row_map_t::non_const_value_type;

  row_map_t row_map;
  entries_t entries;

  diag_generator_functor(row_map_t row_map_, entries_t entries_)
      : row_map(row_map_), entries(entries_){};

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type rowIdx) const {
    row_map(rowIdx + 1) = rowIdx + 1;
    entries(rowIdx)     = rowIdx;
  }
};

template <typename crsMat_t>
void run_experiment(const Params& params) {
  using size_type  = typename crsMat_t::size_type;
  using lno_t      = typename crsMat_t::ordinal_type;
  using scalar_t   = typename crsMat_t::value_type;
  using device_t   = typename crsMat_t::device_type;
  using exec_space = typename device_t::execution_space;

  using graph_t   = typename crsMat_t::StaticCrsGraphType;
  using rowmap_t  = typename graph_t::row_map_type::non_const_type;
  using entries_t = typename graph_t::entries_type::non_const_type;
  using values_t  = typename crsMat_t::values_type::non_const_type;

  std::cout << "************************************* \n";
  std::cout << "************************************* \n";
  crsMat_t A;
  lno_t m = params.m;
  lno_t n = params.n;
  if (params.amtx.length()) {
    std::cout << "Loading A from " << params.amtx << '\n';
    A = KokkosSparse::Impl::read_kokkos_crst_matrix<crsMat_t>(
        params.amtx.c_str());
    m = A.numRows();
    n = A.numCols();
  } else {
    if (params.diag) {
      std::cout << "Randomly generating diag matrix\n";
      rowmap_t rowmapA("A row map", m + 1);
      entries_t entriesA("A entries", m);
      values_t valuesA("A values", m);

      // Generate the graph of A
      diag_generator_functor diag_generator(rowmapA, entriesA);
      Kokkos::parallel_for(Kokkos::RangePolicy<size_type, exec_space>(0, m),
                           diag_generator);

      // Generate the values of A
      Kokkos::Random_XorShift64_Pool<exec_space> rand_pool(13718);
      Kokkos::fill_random(valuesA, rand_pool,
                          10 * Kokkos::ArithTraits<scalar_t>::one());

      // Actually put A together
      graph_t graph(entriesA, rowmapA);
      A = crsMat_t("A matrix", m, valuesA, graph);
    } else {
      std::cout << "Randomly generating matrix\n";
      size_type nnzUnused = m * params.nnzPerRow;
      A = KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(
          m, n, nnzUnused, 0, (n + 3) / 3);
    }
  }

  if (params.verbose) {
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
  handleTime += timer.seconds();

  for (int sumRep = 0; sumRep < params.repeat; sumRep++) {
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
            << handleTime + (symbolicTime / params.repeat) +
                   (numericTime / params.repeat)
            << std::endl
            << "Handle time: " << handleTime << std::endl
            << "Mean symbolic time: " << (symbolicTime / params.repeat)
            << std::endl
            << "Mean numeric time:  " << (numericTime / params.repeat)
            << std::endl;

  if (params.verbose) {
    entries_t permutation = handle.get_permutation();

    std::cout << "MDF permutation:" << std::endl;
    KokkosKernels::Impl::print_1Dview(permutation);
  }
}  // run_experiment

void print_options() {
  std::cerr << "Options\n" << std::endl;

  std::cerr
      << "\t[Required] BACKEND: '--threads[numThreads]' | '--openmp "
         "[numThreads]' | '--cuda [cudaDeviceIndex]' | '--hip [hipDeviceIndex]'"
         " | '--sycl [syclDeviceIndex]'"
      << std::endl;

  std::cerr << "\t[Optional] --amtx <path> :: input matrix" << std::endl;
  std::cerr << "\t[Optional] --repeat      :: how many times to repeat overall "
               "MDF"
            << std::endl;
  std::cerr << "\t[Optional] --verbose     :: enable verbose output"
            << std::endl;
  std::cerr << "\nSettings for randomly generated A matrix" << std::endl;
  std::cerr << "\t[Optional] --m           :: number of rows to generate"
            << std::endl;
  std::cerr << "\t[Optional] --n           :: number of cols to generate"
            << std::endl;
  std::cerr
      << "\t[Optional] --nnz         :: number of entries per row to generate"
      << std::endl;
  std::cerr << "\t[Optional] --diag        :: generate a diagonal matrix"
            << std::endl;
}  // print_options

int parse_inputs(Params& params, int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    if (0 == Test::string_compare_no_case(argv[i], "--threads")) {
      params.use_threads = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--openmp")) {
      params.use_openmp = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--cuda")) {
      params.use_cuda = atoi(argv[++i]) + 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--hip")) {
      params.use_hip = atoi(argv[++i]) + 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--sycl")) {
      params.use_sycl = atoi(argv[++i]) + 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--amtx")) {
      params.amtx = argv[++i];
    } else if (0 == Test::string_compare_no_case(argv[i], "--m")) {
      params.m = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--n")) {
      params.n = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--nnz")) {
      params.nnzPerRow = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--diag")) {
      params.diag = true;
    } else if (0 == Test::string_compare_no_case(argv[i], "--repeat")) {
      params.repeat = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--verbose")) {
      params.verbose = true;
    } else {
      std::cerr << "Unrecognized command line argument #" << i << ": "
                << argv[i] << std::endl;
      print_options();
      return 1;
    }
  }
  return 0;
}  // parse_inputs

int main(int argc, char** argv) {
  Params params;

  if (parse_inputs(params, argc, argv)) {
    return 1;
  }
  const int num_threads =
      std::max(params.use_openmp,
               params.use_threads);  // Assumption is that use_openmp variable
                                     // is provided as number of threads

  // If cuda, hip or sycl is used, set device_id
  int device_id = 0;
  if (params.use_cuda > 0) {
    device_id = params.use_cuda - 1;
  }
  if (params.use_hip > 0) {
    device_id = params.use_hip - 1;
  }
  if (params.use_sycl > 0) {
    device_id = params.use_sycl - 1;
  }

  Kokkos::initialize(Kokkos::InitializationSettings()
                         .set_num_threads(num_threads)
                         .set_device_id(device_id));

  bool useOMP     = params.use_openmp != 0;
  bool useThreads = params.use_threads != 0;
  bool useCUDA    = params.use_cuda != 0;
  bool useHIP     = params.use_hip != 0;
  bool useSYCL    = params.use_sycl != 0;
  bool useSerial  = !useOMP && !useCUDA && !useHIP && !useSYCL;

  if (useOMP) {
#if defined(KOKKOS_ENABLE_OPENMP)
    using crsMat_t =
        KokkosSparse::CrsMatrix<double, int, Kokkos::OpenMP, void, int>;
    run_experiment<crsMat_t>(params);
#else
    std::cout << "ERROR: OpenMP requested, but not available.\n";
    return 1;
#endif
  }
  if (useThreads) {
#if defined(KOKKOS_ENABLE_THREADS)
    using crsMat_t =
        KokkosSparse::CrsMatrix<double, int, Kokkos::Threads, void, int>;
    run_experiment<crsMat_t>(params);
#else
    std::cout << "ERROR: OpenMP requested, but not available.\n";
    return 1;
#endif
  }
  if (useCUDA) {
#if defined(KOKKOS_ENABLE_CUDA)
    using crsMat_t =
        KokkosSparse::CrsMatrix<double, int, Kokkos::Cuda, void, int>;
    run_experiment<crsMat_t>(params);
#else
    std::cout << "ERROR: CUDA requested, but not available.\n";
    return 1;
#endif
  }
  if (useHIP) {
#if defined(KOKKOS_ENABLE_HIP)
    using crsMat_t =
        KokkosSparse::CrsMatrix<double, int, Kokkos::HIP, void, int>;
    run_experiment<crsMat_t>(params);
#else
    std::cout << "ERROR: HIP requested, but not available.\n";
    return 1;
#endif
  }
  if (useSYCL) {
#if defined(KOKKOS_ENABLE_SYCL)
    using crsMat_t =
        KokkosSparse::CrsMatrix<double, int, Kokkos::Experimental::SYCL, void,
                                int>;
    run_experiment<crsMat_t>(params);
#else
    std::cout << "ERROR: SYCL requested, but not available.\n";
    return 1;
#endif
  }
  if (useSerial) {
#if defined(KOKKOS_ENABLE_SERIAL)
    using crsMat_t =
        KokkosSparse::CrsMatrix<double, int, Kokkos::Serial, void, int>;
    run_experiment<crsMat_t>(params);
#else
    std::cout << "ERROR: Serial device requested, but not available.\n";
    return 1;
#endif
  }
  Kokkos::finalize();
  return 0;
}  // main
