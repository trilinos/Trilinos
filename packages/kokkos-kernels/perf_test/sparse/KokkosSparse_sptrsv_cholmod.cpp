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
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_IOUtils.hpp"

#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_sptrsv.hpp"
#include "KokkosSparse_sptrsv_cholmod.hpp"

#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA) && (!defined(KOKKOS_ENABLE_CUDA) || (8000 <= CUDA_VERSION)) && \
    defined(KOKKOSKERNELS_INST_DOUBLE)

#if defined(KOKKOSKERNELS_ENABLE_TPL_CHOLMOD) && defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)

#include "cholmod.h"
// auxiliary functions in perf_test (e.g., pivoting, printing)
#include "KokkosSparse_sptrsv_aux.hpp"

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;
using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosSparse::PerfTest::Experimental;

enum { CUSPARSE, SUPERNODAL_NAIVE, SUPERNODAL_ETREE, SUPERNODAL_SPMV };

/* =========================================================================================
 */
template <typename cholmod_int_type, typename scalar_type, typename size_type, typename ordinal_type>
cholmod_factor *factor_cholmod(const size_type nrow, const size_type nnz, scalar_type *nzvals, size_type *rowptr,
                               ordinal_type *colind, cholmod_common *Comm, int **etree) {
  // Start Cholmod
  cholmod_common *cm = Comm;
  if (std::is_same<cholmod_int_type, long>::value == true) {
    std::cout << " > calling with long " << std::endl;
    cholmod_l_start(cm);
  } else if (std::is_same<cholmod_int_type, int>::value == true) {
    std::cout << " > calling with int " << std::endl;
    cholmod_start(cm);
  }
  cm->supernodal = CHOLMOD_SUPERNODAL;
  cm->print      = 5;

  // Manually, initialize the matrix
  cholmod_sparse A;
  A.stype  = 1;  // symmetric
  A.sorted = 0;
  A.packed = 1;
  if (std::is_same<cholmod_int_type, long>::value == true) {
    A.itype = CHOLMOD_LONG;
  } else if (std::is_same<cholmod_int_type, int>::value == true) {
    A.itype = CHOLMOD_INT;
  }
  if (std::is_same<scalar_type, double>::value == true) {
    A.xtype = CHOLMOD_REAL;
    A.dtype = CHOLMOD_DOUBLE;
  } else if (std::is_same<scalar_type, Kokkos::complex<double>>::value == true) {
    A.xtype = CHOLMOD_COMPLEX;
    A.dtype = CHOLMOD_DOUBLE;
  }
  A.nrow  = nrow;
  A.ncol  = nrow;
  A.nzmax = nnz;
  // Make a copy of Crs's integer pointers with Cholmod int type
  A.p = new cholmod_int_type[nrow + 1];
  A.i = new cholmod_int_type[nnz];
  for (size_type i = 0; i <= nrow; i++) {
    ((cholmod_int_type *)A.p)[i] = rowptr[i];
  }
  for (size_type i = 0; i < nnz; i++) {
    ((cholmod_int_type *)A.i)[i] = colind[i];
  }
  //
  A.x = nzvals;

  // Symbolic factorization
  cholmod_factor *L;
  printf(" ** calling cholmod_analyze **\n");
  if (std::is_same<cholmod_int_type, long>::value == true) {
    L = cholmod_l_analyze(&A, cm);
  } else if (std::is_same<cholmod_int_type, int>::value == true) {
    L = cholmod_analyze(&A, cm);
  }
  if (cm->status != CHOLMOD_OK) {
    printf(" ** cholmod_analyze returned with status = %d **", cm->status);
  }

  // Numerical factorization
  int cholmod_stat = 0;
  printf(" ** calling cholmod_factorize **\n");
  if (std::is_same<cholmod_int_type, long>::value == true) {
    cholmod_stat = cholmod_l_factorize(&A, L, cm);
  } else if (std::is_same<cholmod_int_type, int>::value == true) {
    cholmod_stat = cholmod_factorize(&A, L, cm);
  }
  if (!cholmod_stat) {
    printf(" ** cholmod_factorize returned FALSE **\n");
  }
  if (cm->status != CHOLMOD_OK) {
    printf(" ** cholmod_factorize returned with status = %d, minor = %ld **", cm->status, L->minor);
  }
  switch (cm->selected) {
    case CHOLMOD_NATURAL: printf("  > NATURAL ordering (%d)\n", CHOLMOD_NATURAL); break;
    case CHOLMOD_AMD: printf("  > AMD ordering (%d)\n", CHOLMOD_AMD); break;
    case CHOLMOD_METIS: printf("  > METIS ordering (%d)\n", CHOLMOD_METIS); break;
    case CHOLMOD_NESDIS: printf("  > NESDIS ordering (%d)\n", CHOLMOD_NESDIS); break;
  }
  // print_factor_cholmod<scalar_type>(L, cm);
  compute_etree_cholmod<cholmod_int_type>(&A, cm, etree);

  return L;
}

/* =========================================================================================
 */
void free_cholmod(cholmod_factor *L, cholmod_common *cm) {
  /* free matrices */
  cholmod_free_factor(&L, cm);
  cholmod_finish(cm);
}

/* =========================================================================================
 */
template <typename scalar_type>
int test_sptrsv_perf(std::vector<int> tests, std::string &filename, bool u_in_csr, bool invert_diag,
                     bool invert_offdiag, int block_size, int loop) {
  using STS      = Kokkos::ArithTraits<scalar_type>;
  using mag_type = typename STS::mag_type;

  // using cholmod_int_type = long;
  using cholmod_int_type = int;

  // using int (for CuSparse on GPU)
  using ordinal_type = int;
  using size_type    = int;
  // using ordinal_type = long;
  // using size_type    = long;

  // Default spaces
  using execution_space = Kokkos::DefaultExecutionSpace;
  using memory_space    = typename execution_space::memory_space;

  // Host spaces
  using host_execution_space = Kokkos::DefaultHostExecutionSpace;
  using host_memory_space    = typename host_execution_space::memory_space;

  //
  using host_crsmat_t = KokkosSparse::CrsMatrix<scalar_type, ordinal_type, host_execution_space, void, size_type>;
  using crsmat_t      = KokkosSparse::CrsMatrix<scalar_type, ordinal_type, execution_space, void, size_type>;

  //
  using graph_t = typename crsmat_t::StaticCrsGraphType;

  //
  using host_scalar_view_t = Kokkos::View<scalar_type *, host_memory_space>;
  using scalar_view_t      = Kokkos::View<scalar_type *, memory_space>;

  //
  using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<size_type, ordinal_type, scalar_type,
                                                                        execution_space, memory_space, memory_space>;

  const scalar_type ZERO(0.0);
  const scalar_type ONE(1.0);

  // tolerance
  mag_type tol = STS::epsilon();

  int num_failed = 0;
  std::cout << std::endl;
  std::cout << "Execution space: " << execution_space::name() << std::endl;
  std::cout << "Memory space   : " << memory_space::name() << std::endl;
  std::cout << std::endl;
  if (!filename.empty()) {
    // ==============================================
    // read the matrix ** on host **
    std::cout << " CHOLMOD Tester Begin: Read matrix filename " << filename << std::endl;
    host_crsmat_t Mtx     = KokkosSparse::Impl::read_kokkos_crst_matrix<host_crsmat_t>(filename.c_str());  // in_matrix
    auto graph_host       = Mtx.graph;                                                                     // in_graph
    const size_type nrows = graph_host.numRows();
    auto row_map_host     = graph_host.row_map;
    auto entries_host     = graph_host.entries;
    auto values_host      = Mtx.values;
    // print_factor_cholmod(&Mtx);

    Kokkos::Timer timer;
    // ==============================================
    // call CHOLMOD on the host
    cholmod_common cm;
    cholmod_factor *L = NULL;
    int *etree;
    timer.reset();
    std::cout << " > call CHOLMOD for factorization" << std::endl;
    L                  = factor_cholmod<cholmod_int_type, scalar_type>(nrows, Mtx.nnz(), values_host.data(),
                                                      const_cast<size_type *>(row_map_host.data()), entries_host.data(),
                                                      &cm, &etree);
    double factor_time = timer.seconds();
    std::cout << "   Factorization Time: " << factor_time << std::endl << std::endl;
    using integer_view_host_t = typename KernelHandle::SPTRSVHandleType::integer_view_host_t;
    integer_view_host_t iperm_view("perm view", nrows);
    integer_view_host_t perm_view("iperm view", nrows);

    int *iperm = iperm_view.data();
    int *perm  = perm_view.data();
    for (int i = 0; i < nrows; i++) {
      iperm[i]       = ((cholmod_int_type *)(L->Perm))[i];
      perm[iperm[i]] = i;
    }

    // ==============================================
    // Run all requested algorithms
    for (auto test : tests) {
      std::cout << "\ntest = " << test << std::endl;

      KernelHandle khL, khU;
      switch (test) {
        case SUPERNODAL_NAIVE:
        case SUPERNODAL_ETREE:
        case SUPERNODAL_SPMV: {
          // ==============================================
          // Create handles for U and U^T solves
          if (test == SUPERNODAL_ETREE) {
            std::cout << " > create handle for SUPERNODAL_ETREE" << std::endl << std::endl;
            khL.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_ETREE, nrows, true);
            khU.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_ETREE, nrows, false);
          } else if (test == SUPERNODAL_SPMV) {
            std::cout << " > create handle for SUPERNODAL_SPMV" << std::endl << std::endl;
            khL.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_SPMV, nrows, true);
            khU.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_SPMV, nrows, false);
          } else {
            std::cout << " > create handle for SUPERNODAL_NAIVE" << std::endl << std::endl;
            khL.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_NAIVE, nrows, true);
            khU.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_NAIVE, nrows, false);
          }

          // ==============================================
          // set etree (required)
          khL.set_sptrsv_etree(etree);
          khU.set_sptrsv_etree(etree);

          // ==============================================
          // set permutation
          khL.set_sptrsv_perm(perm);
          khU.set_sptrsv_perm(perm);

          // ==============================================
          // specify if U is stored in CSR or CSC
          std::cout << "=============================== " << std::endl;
          std::cout << " U in CSR           : " << u_in_csr << std::endl;
          khU.set_sptrsv_column_major(!u_in_csr);

          // ==============================================
          // specify wheather to invert diagonal blocks
          khL.set_sptrsv_invert_diagonal(invert_diag);
          khU.set_sptrsv_invert_diagonal(invert_diag);
          invert_diag = khU.get_sptrsv_invert_diagonal();
          std::cout << " Invert diagonal    : " << invert_diag << std::endl;

          // ==============================================
          // specify wheather to apply diagonal-inversion to off-diagonal blocks
          // (optional, default is false)
          khU.set_sptrsv_invert_offdiagonal(invert_offdiag);
          // > make sure if the flag is set before setting for L
          invert_offdiag = khU.get_sptrsv_invert_offdiagonal();
          khL.set_sptrsv_invert_offdiagonal(invert_offdiag);
          std::cout << " Invert Off-diagonal: " << invert_offdiag << std::endl;

          // ==============================================
          // block size to switch to device call
          if (block_size >= 0) {
            std::cout << " Block Size         : " << block_size << std::endl;
            khL.set_sptrsv_diag_supernode_sizes(block_size, block_size);
            khU.set_sptrsv_diag_supernode_sizes(block_size, block_size);
          }
          std::cout << std::endl;

          // ==============================================
          // Do symbolic analysis
          timer.reset();
          sptrsv_symbolic<cholmod_int_type>(&khL, &khU, L, &cm);
          double symbolic_time = timer.seconds();
          std::cout << "   Symbolic Time   : " << symbolic_time << std::endl << std::endl;

          // ==============================================
          // Do numerical compute
          timer.reset();
          sptrsv_compute<cholmod_int_type>(&khL, &khU, L, &cm);
          double compute_time = timer.seconds();
          std::cout << "   Numeric Time   : " << compute_time << std::endl << std::endl;

          // ==============================================
          // Create the known solution and set to all 1's ** on host **
          host_scalar_view_t sol_host("sol_host", nrows);
          // Kokkos::deep_copy(sol_host, ONE);
          Kokkos::Random_XorShift64_Pool<host_execution_space> random(13718);
          Kokkos::fill_random(sol_host, random, scalar_type(1));

          // ==============================================
          // Create the rhs ** on host **
          // A*sol
          host_scalar_view_t rhs_host("rhs_host", nrows);
          KokkosSparse::spmv("N", ONE, Mtx, sol_host, ZERO, rhs_host);

          // ==============================================
          // apply forward-pivot on the host
          host_scalar_view_t tmp_host("temp", nrows);
          forwardP_supernode<scalar_type>(nrows, perm, 1, rhs_host.data(), nrows, tmp_host.data(), nrows);

          // ==============================================
          // copy rhs to the default host/device
          scalar_view_t rhs("rhs", nrows);
          scalar_view_t sol("sol", nrows);
          Kokkos::deep_copy(rhs, tmp_host);

          // ==============================================
          // do L solve
          // numeric (only rhs is modified) on the default device/host space
          timer.reset();
          sptrsv_solve(&khL, sol, rhs);
          Kokkos::fence();
          std::cout << "   Solve Time   : " << timer.seconds() << std::endl;

          // ==============================================
          // do L^T solve
          // numeric (only rhs is modified) on the default device/host space
          timer.reset();
          sptrsv_solve(&khU, rhs, sol);
          Kokkos::fence();
          std::cout << "   Solve Time   : " << timer.seconds() << std::endl;

          // ==============================================
          // apply backward-pivot
          // > copy solution to host
          Kokkos::deep_copy(tmp_host, rhs);
          backwardP_supernode<scalar_type>(nrows, perm, 1, tmp_host.data(), nrows, sol_host.data(), nrows);

          // ==============================================
          // Error Check ** on host **
          Kokkos::fence();
          std::cout << std::endl;
          if (!check_errors(tol, Mtx, rhs_host, sol_host)) {
            num_failed++;
          }

          // try again?
          {
            Kokkos::deep_copy(sol_host, ONE);
            KokkosSparse::spmv("N", ONE, Mtx, sol_host, ZERO, rhs_host);
            forwardP_supernode<scalar_type>(nrows, perm, 1, rhs_host.data(), nrows, tmp_host.data(), nrows);
            Kokkos::deep_copy(rhs, tmp_host);
            sptrsv_solve(&khL, sol, rhs);
            sptrsv_solve(&khU, rhs, sol);
            Kokkos::fence();
            Kokkos::deep_copy(tmp_host, rhs);
            backwardP_supernode<scalar_type>(nrows, perm, 1, tmp_host.data(), nrows, sol_host.data(), nrows);

            if (!check_errors(tol, Mtx, rhs_host, sol_host)) {
              num_failed++;
            }
          }
          std::cout << std::endl;

          // Benchmark
          // L-solve
          double min_time = 0.0;
          double max_time = 0.0;
          double ave_time = 0.0;
          Kokkos::fence();
          for (int i = 0; i < loop; i++) {
            timer.reset();
            sptrsv_solve(&khL, sol, rhs);
            Kokkos::fence();
            double time = timer.seconds();
            ave_time += time;
            if (time > max_time || i == 0) max_time = time;
            if (time < min_time || i == 0) min_time = time;
          }
          std::cout << " L-solve: loop = " << loop << std::endl;
          std::cout << "  LOOP_AVG_TIME:  " << ave_time / loop << std::endl;
          std::cout << "  LOOP_MAX_TIME:  " << max_time << std::endl;
          std::cout << "  LOOP_MIN_TIME:  " << min_time << std::endl << std::endl;

          // U-solve
          min_time = 1.0e32;
          max_time = 0.0;
          ave_time = 0.0;
          Kokkos::fence();
          for (int i = 0; i < loop; i++) {
            timer.reset();
            sptrsv_solve(&khU, rhs, sol);
            Kokkos::fence();
            double time = timer.seconds();
            ave_time += time;
            if (time > max_time) max_time = time;
            if (time < min_time) min_time = time;
            // std::cout << time << std::endl;
          }
          std::cout << " U-solve: loop = " << loop << std::endl;
          std::cout << "  LOOP_AVG_TIME:  " << ave_time / loop << std::endl;
          std::cout << "  LOOP_MAX_TIME:  " << max_time << std::endl;
          std::cout << "  LOOP_MIN_TIME:  " << min_time << std::endl << std::endl;
        } break;

        case CUSPARSE: {
          std::cout << " > create handle for CuSparse (SUPERNODAL_NAIVE)" << std::endl << std::endl;
          khL.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_NAIVE, nrows, true);

          // ==============================================
          // read CHOLMOD factor int crsMatrix on the host (cholmodMat_host) and
          // copy to default host/device (cholmodMtx)
          timer.reset();
          std::cout << " > Read Cholmod factor into KokkosSparse::CrsMatrix "
                       "(invert diagonabl, and copy to device) "
                    << std::endl;
          khL.set_sptrsv_invert_diagonal(false);
          auto graph      = read_cholmod_graphL<cholmod_int_type, graph_t>(&khL, L, &cm);
          auto cholmodMtx = read_cholmod_factor<cholmod_int_type, crsmat_t, graph_t>(&khL, L, &cm, graph);
          std::cout << "   Conversion Time: " << timer.seconds() << std::endl << std::endl;

          bool col_majorL = true;
          bool col_majorU = false;
          if (!check_cusparse(Mtx, col_majorL, cholmodMtx, col_majorU, cholmodMtx, perm, perm, tol, loop)) {
            num_failed++;
          }
        } break;

        default: std::cout << " > Testing only Cholmod < " << std::endl; exit(0);
      }
    }

    // free cholmod data structures
    free_cholmod(L, &cm);
  }
  std::cout << std::endl << std::endl;

  return num_failed;
}

void print_help_sptrsv() {
  printf("Options:\n");
  printf("  --test [OPTION]        : Use different kernel implementations\n");
  printf("                           Options:\n");
  printf("                           cholmod_naive, cholmod_etree\n\n");
  printf(
      "  -f [file]              : Read in Matrix Market formatted text file "
      "'file'.\n");
  printf("  --loop [LOOP]          : How many time to run the test.\n");
  printf(
      "  --u-in-csc             : To store U-factor in CSC, needed for "
      "invert.\n");
  printf("  --invert-diag          : To invert diagonal blocks.\n");
  printf("  --invert-offdiag       : To apply inverse to off-diagonal blocks.\n");
  printf(
      "  --block-size [SIZE]    : To specify the threshold to switch device "
      "and bached kernel.\n");
  printf("  --scalar-type [d or z] :\n");
}

int main(int argc, char **argv) {
  std::vector<int> tests;
  std::string filename;

  int loop = 1;
  // store U in CSR, or CSC
  bool u_in_csr = true;
  // invert diagonal of L-factor
  bool invert_diag = false;
  // apply invert of diagonal to offdiagonal
  bool invert_offdiag = false;
  // block size to switch to device call (default is 100)
  int block_size = -1;
  // scalar type
  std::string char_scalar = "d";

  if (argc == 1) {
    print_help_sptrsv();
    return 0;
  }

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "--test") == 0)) {
      i++;
      if ((strcmp(argv[i], "cholmod-naive") == 0)) {
        tests.push_back(SUPERNODAL_NAIVE);
      } else if ((strcmp(argv[i], "cholmod-etree") == 0)) {
        tests.push_back(SUPERNODAL_ETREE);
      } else if ((strcmp(argv[i], "cholmod-spmv") == 0)) {
        tests.push_back(SUPERNODAL_SPMV);
      } else if ((strcmp(argv[i], "cusparse") == 0)) {
        tests.push_back(CUSPARSE);
      }
      continue;
    }
    if ((strcmp(argv[i], "-f") == 0)) {
      filename = argv[++i];
      continue;
    }
    if ((strcmp(argv[i], "--loop") == 0)) {
      loop = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "--u-in-csc") == 0)) {
      u_in_csr = false;
      continue;
    }
    if ((strcmp(argv[i], "--invert-diag") == 0)) {
      invert_diag = true;
      continue;
    }
    if ((strcmp(argv[i], "--invert-offdiag") == 0)) {
      invert_offdiag = true;
      continue;
    }
    if ((strcmp(argv[i], "--block-size") == 0)) {
      block_size = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "--scalar-type") == 0)) {
      char_scalar = argv[++i];
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
    // Cholmod may not support single, yet
    // int total_errors = test_sptrsv_perf<float>(tests, filename, loop);
    int total_errors = 0;
    Kokkos::ScopeGuard kokkosScope(argc, argv);
    if (char_scalar == "z") {
#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
      total_errors = test_sptrsv_perf<Kokkos::complex<double>>(tests, filename, u_in_csr, invert_diag, invert_offdiag,
                                                               block_size, loop);
#else
      std::cout << std::endl << " KOKKOSKERNELS_INST_COMPLEX_DOUBLE  is not enabled ** " << std::endl << std::endl;
#endif
    } else if (char_scalar == "d") {
#if defined(KOKKOSKERNELS_INST_DOUBLE)
      total_errors = test_sptrsv_perf<double>(tests, filename, u_in_csr, invert_diag, invert_offdiag, block_size, loop);
#else
      std::cout << std::endl << " KOKKOSKERNELS_INST_DOUBLE  is not enabled ** " << std::endl << std::endl;
#endif
    }
    if (total_errors == 0)
      std::cout << "Kokkos::SPTRSV Test: Passed" << std::endl << std::endl;
    else
      std::cout << "Kokkos::SPTRSV Test: Failed" << std::endl << std::endl;
  }

  return 0;
}
#else  // defined(KOKKOSKERNELS_ENABLE_TPL_CHOLMOD)
int main() {
  std::cout << std::endl << "** CHOLMOD NOT ENABLED **" << std::endl << std::endl;
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
  std::cout << " KOKKOS_ENABLE_CUDA_LAMBDA not defined\n" << std::endl;
#endif
  std::cout << " CUDA_VERSION = " << CUDA_VERSION << std::endl;
#endif
  return 1;
}
#endif
