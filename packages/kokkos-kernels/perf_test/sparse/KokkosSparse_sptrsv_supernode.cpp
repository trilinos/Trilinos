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

#include "Kokkos_Random.hpp"
#include "KokkosSparse_Utils.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_IOUtils.hpp"

#include "KokkosSparse_sptrsv.hpp"
#include "KokkosSparse_sptrsv_supernode.hpp"

#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA) && (!defined(KOKKOS_ENABLE_CUDA) || (8000 <= CUDA_VERSION)) && \
    defined(KOKKOSKERNELS_INST_DOUBLE)

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)

#include "KokkosSparse_sptrsv_aux.hpp"

namespace KSExp = KokkosSparse::Experimental;
namespace KSPTE = KokkosSparse::PerfTest::Experimental;

enum { CUSPARSE, SUPERNODAL_NAIVE, SUPERNODAL_ETREE, SUPERNODAL_DAG, SUPERNODAL_SPMV, SUPERNODAL_SPMV_DAG };

/* =========================================================================================
 */
template <typename scalar_type>
int test_sptrsv_perf(std::vector<int> tests, bool verbose, std::string& lower_filename, std::string& upper_filename,
                     std::string& supernode_filename, bool merge, bool invert_offdiag, bool u_in_csr, int loop) {
  using ordinal_type = int;
  using size_type    = int;
  using STS          = Kokkos::ArithTraits<scalar_type>;
  using mag_type     = typename STS::mag_type;

  // Default spaces
  // using execution_space = Kokkos::OpenMP;
  using execution_space = Kokkos::DefaultExecutionSpace;
  using memory_space    = typename execution_space::memory_space;

  // Host spaces
  using host_execution_space = Kokkos::DefaultHostExecutionSpace;
  using host_memory_space    = typename host_execution_space::memory_space;

  //
  using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<size_type, ordinal_type, scalar_type,
                                                                        execution_space, memory_space, memory_space>;

  //
  using host_crsmat_t = typename KernelHandle::SPTRSVHandleType::host_crsmat_t;

  //
  using host_scalar_view_t = Kokkos::View<scalar_type*, host_memory_space>;
  using scalar_view_t      = Kokkos::View<scalar_type*, memory_space>;

  const scalar_type ZERO(0.0);
  const scalar_type ONE(1.0);

  // tolerance
  mag_type tol = STS::epsilon();

  int num_failed = 0;
  std::cout << std::endl;
  std::cout << "Execution space: " << execution_space::name() << std::endl;
  std::cout << "Memory space   : " << memory_space::name() << std::endl;
  std::cout << std::endl;
  if ((!lower_filename.empty() || !upper_filename.empty()) && !supernode_filename.empty()) {
    // ==============================================
    // read the CRS matrix ** on host **
    // it stores the supernodal triangular matrix, stored by blocks with
    // explicit zeros
    std::cout << " Supernode Tester Begin:" << std::endl;
    std::string matrix_filename = (lower_filename.empty() ? upper_filename : lower_filename);
    std::cout << " > Read a triangular-matrix filename " << matrix_filename << std::endl;
    host_crsmat_t M       = KokkosSparse::Impl::read_kokkos_crst_matrix<host_crsmat_t>(matrix_filename.c_str());
    const size_type nrows = M.graph.numRows();
    // transpose the matrix to be stored in CCS
    //{
    auto graphM   = M.graph;  // in_graph
    auto row_mapM = graphM.row_map;
    auto entriesM = graphM.entries;
    auto valuesM  = M.values;

    using host_graph_t      = typename host_crsmat_t::StaticCrsGraphType;
    using row_map_view_t    = typename host_graph_t::row_map_type::non_const_type;
    using cols_view_t       = typename host_graph_t::entries_type::non_const_type;
    using values_view_t     = typename host_crsmat_t::values_type::non_const_type;
    using in_row_map_view_t = typename host_graph_t::row_map_type;
    using in_cols_view_t    = typename host_graph_t::entries_type;
    using in_values_view_t  = typename host_crsmat_t::values_type;

    size_type nnzL = row_mapM(nrows);
    row_map_view_t row_map("rowmap_view", nrows + 1);
    cols_view_t entries("colmap_view", nnzL);
    values_view_t values("values_view", nnzL);
    // transpose L
    KokkosSparse::Impl::transpose_matrix<in_row_map_view_t, in_cols_view_t, in_values_view_t, row_map_view_t,
                                         cols_view_t, values_view_t, row_map_view_t, host_execution_space>(
        nrows, nrows, row_mapM, entriesM, valuesM, row_map, entries, values);

    // store L in CSC
    host_graph_t static_graph(entries, row_map);
    host_crsmat_t T("CrsMatrix", nrows, values, static_graph);
    //}
    host_crsmat_t A, L;
    if (!lower_filename.empty()) {
      A = M;  // original L in CSR
      L = T;  // transposed L in CSC
    } else {
      A = T;  // original L was in CSC, so transposed L in CSR
      L = M;  // original L in CSC
    }
    // read supernode sizes from a file
    // first entry in the file specifies the number of supernode
    // the rest of the entries specifies the column offsets to the beginning of
    // supernodes
    std::cout << " > Read a supernode filename " << supernode_filename << std::endl;
    std::ifstream fp(supernode_filename.c_str(), std::ios::in);
    if (!fp.is_open()) {
      std::cout << std::endl << " failed to open " << supernode_filename << std::endl << std::endl;
      return 0;
    }
    int nsuper;
    fp >> nsuper;
    Kokkos::View<int*, Kokkos::HostSpace> supercols("supercols", 1 + nsuper);
    int* etree = NULL;
    for (int i = 0; i <= nsuper; i++) {
      fp >> supercols(i);
    }
    fp.close();

    Kokkos::Timer timer;
    // ==============================================
    // Run all requested algorithms
    for (auto test : tests) {
      std::cout << "\ntest = " << test << std::endl;

      KernelHandle khL;
      KernelHandle khU;  // TODO: can I just make a copy later (khU = khL)?
      switch (test) {
        case SUPERNODAL_NAIVE:
        case SUPERNODAL_ETREE:
        case SUPERNODAL_DAG:
        case SUPERNODAL_SPMV:
        case SUPERNODAL_SPMV_DAG: {
          // ==============================================
          // create an handle
          if (test == SUPERNODAL_NAIVE) {
            std::cout << " > create handle for SUPERNODAL_NAIVE" << std::endl << std::endl;
            khL.create_sptrsv_handle(KSExp::SPTRSVAlgorithm::SUPERNODAL_NAIVE, nrows, true);
            khU.create_sptrsv_handle(KSExp::SPTRSVAlgorithm::SUPERNODAL_NAIVE, nrows, true);
          } else if (test == SUPERNODAL_DAG) {
            std::cout << " > create handle for SUPERNODAL_DAG" << std::endl << std::endl;
            khL.create_sptrsv_handle(KSExp::SPTRSVAlgorithm::SUPERNODAL_DAG, nrows, true);
            khU.create_sptrsv_handle(KSExp::SPTRSVAlgorithm::SUPERNODAL_DAG, nrows, true);
          } else if (test == SUPERNODAL_SPMV_DAG) {
            std::cout << " > create handle for SUPERNODAL_SPMV_DAG" << std::endl << std::endl;
            khL.create_sptrsv_handle(KSExp::SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG, nrows, true);
            khU.create_sptrsv_handle(KSExp::SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG, nrows, true);
          }
          // verbose (optional, default is false)
          khL.set_sptrsv_verbose(verbose);
          khU.set_sptrsv_verbose(verbose);

          // specify if U is stored in CSR or CSC
          khU.set_sptrsv_column_major(!u_in_csr);

          // specify wheather to merge supernodes (optional, default merge is
          // false)
          khL.set_sptrsv_merge_supernodes(merge);
          khU.set_sptrsv_merge_supernodes(merge);

          // specify wheather to apply diagonal-inversion to off-diagonal blocks
          // (optional, default is false)
          khL.set_sptrsv_invert_offdiagonal(invert_offdiag);
          khU.set_sptrsv_invert_offdiagonal(invert_offdiag);

          // ==============================================
          // do symbolic analysis (preprocssing, e.g., merging supernodes,
          // inverting diagonal/offdiagonal blocks, and scheduling based on
          // graph/dag)
          khU.get_sptrsv_handle()->set_column_major(!khL.get_sptrsv_handle()->is_column_major());
          KSExp::sptrsv_supernodal_symbolic(nsuper, supercols.data(), etree, L.graph, &khL, L.graph, &khU);

          // ==============================================
          // do numeric compute (copy numerical values from SuperLU data
          // structure to our sptrsv data structure)
          KSExp::sptrsv_compute(&khL, L);

          // ==============================================
          // Preaparing for the first solve
          //> create the known solution and set to all 1's ** on host **
          host_scalar_view_t sol_host("sol_host", nrows);
          // Kokkos::deep_copy (sol_host, ONE);
          Kokkos::Random_XorShift64_Pool<host_execution_space> random(13718);
          Kokkos::fill_random(sol_host, random, scalar_type(1));

          // > create the rhs ** on host **
          // A*sol generates rhs: rhs is dense, use spmv
          host_scalar_view_t rhs_host("rhs_host", nrows);
          KokkosSparse::spmv("N", ONE, A, sol_host, ZERO, rhs_host);

          // ==============================================
          // copy rhs to the default host/device
          scalar_view_t rhs("rhs", nrows);
          scalar_view_t sol("sol", nrows);
          Kokkos::deep_copy(rhs, rhs_host);

          // ==============================================
          // do L solve
          timer.reset();
          KSExp::sptrsv_solve(&khL, sol, rhs);
          Kokkos::fence();
          std::cout << " > Lower-TRI: " << std::endl;
          std::cout << "   Solve Time   : " << timer.seconds() << std::endl;

          // copy solution to host
          Kokkos::deep_copy(sol_host, sol);

          // ==============================================
          // Error Check ** on host **
          Kokkos::fence();
          std::cout << std::endl;
          if (!KSPTE::check_errors(tol, A, rhs_host, sol_host)) {
            num_failed++;
          }

          // Benchmark
          // L-solve
          double min_time = 1.0e32;
          double max_time = 0.0;
          double ave_time = 0.0;
          Kokkos::fence();
          for (int i = 0; i < loop; i++) {
            timer.reset();
            KSExp::sptrsv_solve(&khL, sol, rhs);
            Kokkos::fence();
            double time = timer.seconds();
            ave_time += time;
            if (time > max_time) max_time = time;
            if (time < min_time) min_time = time;
            // std::cout << time << std::endl;
          }
          std::cout << " L-solve: loop = " << loop << std::endl;
          std::cout << "  LOOP_AVG_TIME:  " << ave_time / loop << std::endl;
          std::cout << "  LOOP_MAX_TIME:  " << max_time << std::endl;
          std::cout << "  LOOP_MIN_TIME:  " << min_time << std::endl << std::endl;
        } break;

        default: std::cout << " > Invalid test ID < " << std::endl; exit(0);
      }
    }
  }
  std::cout << std::endl << std::endl;

  return num_failed;
}

void print_help_sptrsv() {
  printf("Options:\n");
  printf("  --test [OPTION] : Use different kernel implementations\n");
  printf("                    Options:\n");
  printf("                    superlu-naive, superlu-dag, superlu-spmv-dag\n\n");
  printf(
      "  -lf [file]      : Read in Lower-triangular matrix in Matrix Market "
      "formatted text file 'file'.\n");
  printf(
      "  -uf [file]      : Read in Upper-triangular matrix in Matrix Market "
      "formatted text file 'file'.\n");
  printf("  -sf [file]      : Read in Supernode sizes from 'file'.\n");
  printf("  --loop [LOOP]   : How many spmv to run to aggregate average time. \n");
}

int main(int argc, char** argv) {
  std::vector<int> tests;
  std::string lower_filename;
  std::string upper_filename;
  std::string supernode_filename;

  int loop = 1;
  // merge supernodes
  bool merge = false;
  // invert off-diagonal of L-factor
  bool invert_offdiag = false;
  // store U in CSR, or CSC
  bool u_in_csr = true;
  // verbose
  bool verbose = true;

  if (argc == 1) {
    print_help_sptrsv();
    return 0;
  }

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "--test") == 0)) {
      i++;
      if ((strcmp(argv[i], "superlu-naive") == 0)) {
        tests.push_back(SUPERNODAL_NAIVE);
      }
      if ((strcmp(argv[i], "superlu-dag") == 0)) {
        tests.push_back(SUPERNODAL_DAG);
      }
      if ((strcmp(argv[i], "superlu-spmv-dag") == 0)) {
        tests.push_back(SUPERNODAL_SPMV_DAG);
      }
      if ((strcmp(argv[i], "cusparse") == 0)) {
        tests.push_back(CUSPARSE);
      }
      continue;
    }
    if ((strcmp(argv[i], "-lf") == 0)) {
      lower_filename = argv[++i];
      continue;
    }
    if ((strcmp(argv[i], "-uf") == 0)) {
      upper_filename = argv[++i];
      continue;
    }
    if ((strcmp(argv[i], "-sf") == 0)) {
      supernode_filename = argv[++i];
      continue;
    }
    if ((strcmp(argv[i], "--quiet") == 0)) {
      verbose = false;
      continue;
    }
    if ((strcmp(argv[i], "--loop") == 0)) {
      loop = atoi(argv[++i]);
      continue;
    }
    /* not supported through this interface, yet
     *  if((strcmp(argv[i],"--merge")==0)) {
      merge = true;
      continue;
    }*/
    if ((strcmp(argv[i], "--invert-offdiag") == 0)) {
      invert_offdiag = true;
      continue;
    }
    if ((strcmp(argv[i], "--u-in-csc") == 0)) {
      u_in_csr = false;
      continue;
    }
    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      print_help_sptrsv();
      return 0;
    }
  }

  std::cout << std::endl;
  for (size_t i = 0; i < tests.size(); ++i) {
    std::cout << "tests[" << i << "] = " << tests[i] << std::endl;
  }

  {
    using scalar_t = double;
    // using scalar_t = Kokkos::complex<double>;
    Kokkos::ScopeGuard kokkosScope(argc, argv);
    int total_errors = test_sptrsv_perf<scalar_t>(tests, verbose, lower_filename, upper_filename, supernode_filename,
                                                  merge, invert_offdiag, u_in_csr, loop);
    if (total_errors == 0)
      std::cout << "Kokkos::SPTRSV Test: Passed" << std::endl << std::endl;
    else
      std::cout << "Kokkos::SPTRSV Test: Failed (" << total_errors << " / " << 2 * tests.size() << " failed)"
                << std::endl
                << std::endl;
  }
  return 0;
}
#else  // defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
int main() {
  std::cout << std::endl << " ** SUPERNODAL NOT ENABLED **" << std::endl << std::endl;
  exit(0);
  return 1;
}
#endif

#else  // defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA ) &&
       // (!defined(KOKKOS_ENABLE_CUDA) || ( 8000 <= CUDA_VERSION ))

int main() {
#if !defined(KOKKOSKERNELS_INST_DOUBLE)
  std::cout << " Only supported with double precision" << std::endl;
#endif
#if !defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
  std::cout << " KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA **not** defined" << std::endl;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
  std::cout << " KOKKOS_ENABLE_CUDA defined" << std::endl;
#if !defined(KOKKOS_ENABLE_CUDA_LAMBDA)
  std::cout << " KOKKOS_ENABLE_CUDA_LAMBDA **not** defined" << std::endl;
#endif
  std::cout << " CUDA_VERSION = " << CUDA_VERSION << std::endl;
#endif
  return 1;
}
#endif
