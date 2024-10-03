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
#include "KokkosSparse_sptrsv_superlu.hpp"

#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA) && (!defined(KOKKOS_ENABLE_CUDA) || (8000 <= CUDA_VERSION))

#if defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU) && defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)

namespace SLU {
namespace D {
#include "slu_ddefs.h"
}
namespace S {
#include "slu_sdefs.h"
}
namespace C {
#include "slu_cdefs.h"
}
namespace Z {
#include "slu_zdefs.h"
}
}  // namespace SLU
// auxiliary functions from perf_test (e.g., pivoting, printing)
#include "KokkosSparse_sptrsv_aux.hpp"

#ifdef KOKKOSKERNELS_ENABLE_TPL_METIS
// optionally, for matrix ordering before SuperLU
#include "metis.h"
#endif

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;
using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosSparse::PerfTest::Experimental;
using namespace SLU::D;
using namespace SLU::S;
using namespace SLU::C;
using namespace SLU::Z;

enum { CUSPARSE, SUPERNODAL_NAIVE, SUPERNODAL_ETREE, SUPERNODAL_DAG, SUPERNODAL_SPMV, SUPERNODAL_SPMV_DAG };

/* =========================================================================================
 */
template <typename scalar_type>
void factor_superlu(bool symm_mode, bool metis, const int nrow, scalar_type *nzvals, int *rowptr, int *colind,
                    int panel_size, int relax_size, SuperMatrix &L, SuperMatrix &U, int **perm_r, int **perm_c,
                    int **parents) {
  // allocate permutation vectors for SuperLU
  *perm_c = new int[nrow];
  *perm_r = new int[nrow];

  if (metis) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_METIS
    if (symm_mode) {
      idx_t n   = nrow;
      idx_t nnz = rowptr[n];

      // remove diagonal elements (and casting to METIS idx_t)
      idx_t *metis_rowptr = new idx_t[n + 1];
      idx_t *metis_colind = new idx_t[nnz];

      nnz             = 0;
      metis_rowptr[0] = 0;
      for (int i = 0; i < n; i++) {
        for (int k = rowptr[i]; k < rowptr[i + 1]; k++) {
          if (colind[k] != i) {
            metis_colind[nnz] = colind[k];
            nnz++;
          }
        }
        metis_rowptr[i + 1] = nnz;
      }

      // call METIS
      idx_t *metis_perm  = new idx_t[n];
      idx_t *metis_iperm = new idx_t[n];
      std::cout << "  + calling METIS_NodeND: (n=" << n << ", nnz=" << nnz << ") " << std::endl;
      if (METIS_OK != METIS_NodeND(&n, metis_rowptr, metis_colind, NULL, NULL, metis_perm, metis_iperm)) {
        std::cout << std::endl << "METIS_NodeND failed" << std::endl << std::endl;
      }

      // copy permutation to SuperLU
      for (idx_t i = 0; i < n; i++) {
        (*perm_r)[i] = metis_iperm[i];
        (*perm_c)[i] = metis_iperm[i];
      }

      delete[] metis_perm;
      delete[] metis_iperm;
      delete[] metis_rowptr;
      delete[] metis_colind;
    } else {
      std::cout << "   + METIS enabled only for symmetric mode" << std::endl << std::endl;
      metis = false;
    }
#else
    std::cout << std::endl << " ** METIS not ENABLED **" << std::endl << std::endl;
#endif
  }

  SuperMatrix A;
  NCformat *Astore;
  int info;
  superlu_options_t options;
  SuperLUStat_t stat;

  set_default_options(&options);
  if (symm_mode) {
    options.SymmetricMode = YES;
  }
#ifdef KOKKOSKERNELS_ENABLE_TPL_METIS
  if (metis) {
    options.ColPerm = MY_PERMC;
    options.RowPerm = MY_PERMR;
  }
#endif

  int nnz = rowptr[nrow];
  // casting to call d or z version
  scalar_type *nzvals_tran = nzvals;
  int *rowptr_tran         = rowptr;
  int *colind_tran         = colind;
  if (std::is_same<scalar_type, double>::value == true) {
    if (!symm_mode) {
      dCompRow_to_CompCol(nrow, nrow, nnz, reinterpret_cast<double *>(nzvals), colind, rowptr,
                          reinterpret_cast<double **>(&nzvals_tran), &colind_tran, &rowptr_tran);
    }
    dCreate_CompCol_Matrix(&A, nrow, nrow, nnz, reinterpret_cast<double *>(nzvals_tran), colind_tran, rowptr_tran,
                           SLU_NC, SLU_D, SLU_GE);
  } else if (std::is_same<scalar_type, float>::value == true) {
    if (!symm_mode) {
      sCompRow_to_CompCol(nrow, nrow, nnz, reinterpret_cast<float *>(nzvals), colind, rowptr,
                          reinterpret_cast<float **>(&nzvals_tran), &colind_tran, &rowptr_tran);
    }
    sCreate_CompCol_Matrix(&A, nrow, nrow, nnz, reinterpret_cast<float *>(nzvals_tran), colind_tran, rowptr_tran,
                           SLU_NC, SLU_D, SLU_GE);
  } else if (std::is_same<scalar_type, std::complex<float>>::value == true ||
             std::is_same<scalar_type, Kokkos::complex<float>>::value == true) {
    if (!symm_mode) {
      cCompRow_to_CompCol(nrow, nrow, nnz, reinterpret_cast<complex *>(nzvals), colind, rowptr,
                          reinterpret_cast<complex **>(&nzvals_tran), &colind_tran, &rowptr_tran);
    }
    cCreate_CompCol_Matrix(&A, nrow, nrow, nnz, reinterpret_cast<complex *>(nzvals_tran), colind_tran, rowptr_tran,
                           SLU_NC, SLU_Z, SLU_GE);
  } else if (std::is_same<scalar_type, std::complex<double>>::value == true ||
             std::is_same<scalar_type, Kokkos::complex<double>>::value == true) {
    if (!symm_mode) {
      zCompRow_to_CompCol(nrow, nrow, nnz, reinterpret_cast<doublecomplex *>(nzvals), colind, rowptr,
                          reinterpret_cast<doublecomplex **>(&nzvals_tran), &colind_tran, &rowptr_tran);
    }
    zCreate_CompCol_Matrix(&A, nrow, nrow, nnz, reinterpret_cast<doublecomplex *>(nzvals_tran), colind_tran,
                           rowptr_tran, SLU_NC, SLU_Z, SLU_GE);
  }

  /* Initialize the statistics variables. */
  StatInit(&stat);
  int w1 = (sp_ienv(1) > sp_ienv(2) ? sp_ienv(1) : sp_ienv(2));
  int w2 = (panel_size > relax_size ? panel_size : relax_size);
  if (w2 > w1) {
    SUPERLU_FREE(stat.panel_histo);
    stat.panel_histo = intCalloc(w2 + 1);
  }

  /* Call SuperLU to solve the problem. */
  int *etree = new int[A.ncol];
  if (options.ColPerm != MY_PERMC) {
    get_perm_c(options.ColPerm, &A, *perm_c);
  }
  SuperMatrix AC;
  sp_preorder(&options, &A, *perm_c, etree, &AC);

  Astore = (NCformat *)(A.Store);
  printf("  + calling SuperLU dgstrf with panel_size=%d, relax_size=%d..\n", panel_size, relax_size);
  printf("   * Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);

#ifdef HAVE_KOKKOSKERNELS_SUPERLU5_API
  GlobalLU_t Glu;
#endif
  int lwork = 0;
  if (std::is_same<scalar_type, double>::value == true) {
    dgstrf(&options, &AC, relax_size, panel_size, etree, NULL, lwork, *perm_c, *perm_r, &L, &U,
#ifdef HAVE_KOKKOSKERNELS_SUPERLU5_API
           &Glu,
#endif
           &stat, &info);
  } else if (std::is_same<scalar_type, float>::value == true) {
    sgstrf(&options, &AC, relax_size, panel_size, etree, NULL, lwork, *perm_c, *perm_r, &L, &U,
#ifdef HAVE_KOKKOSKERNELS_SUPERLU5_API
           &Glu,
#endif
           &stat, &info);
  } else if (std::is_same<scalar_type, std::complex<float>>::value == true ||
             std::is_same<scalar_type, Kokkos::complex<float>>::value == true) {
    cgstrf(&options, &AC, relax_size, panel_size, etree, NULL, lwork, *perm_c, *perm_r, &L, &U,
#ifdef HAVE_KOKKOSKERNELS_SUPERLU5_API
           &Glu,
#endif
           &stat, &info);
  } else {
    zgstrf(&options, &AC, relax_size, panel_size, etree, NULL, lwork, *perm_c, *perm_r, &L, &U,
#ifdef HAVE_KOKKOSKERNELS_SUPERLU5_API
           &Glu,
#endif
           &stat, &info);
  }
  if (info != 0) printf(" SuperLU failed with info=%d\n", info);
  StatFree(&stat);
  Destroy_SuperMatrix_Store(&A);
  Destroy_CompCol_Permuted(&AC);
  if (!symm_mode) {
    SUPERLU_FREE(nzvals_tran);
    SUPERLU_FREE(rowptr_tran);
    SUPERLU_FREE(colind_tran);
  }

  /* convert etree to parents */
  SCformat *Lstore = (SCformat *)(L.Store);
  int nsuper       = 1 + Lstore->nsuper;  // # of supernodal columns
  *parents         = new int[nsuper];
  for (int s = 0; s < nsuper; s++) {
    int j = Lstore->sup_to_col[s + 1] - 1;  // the last column index of this supernode
    if (etree[j] == nrow) {
      (*parents)[s] = -1;
    } else {
      (*parents)[s] = Lstore->col_to_sup[etree[j]];
    }
  }
  delete[] etree;

  return;
}

/* =========================================================================================
 */
void free_superlu(SuperMatrix &L, SuperMatrix &U, int *perm_r, int *perm_c, int *parents) {
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);

  delete[] perm_r;
  delete[] perm_c;
  delete[] parents;
}

/* =========================================================================================
 */
template <typename scalar_type>
int test_sptrsv_perf(std::vector<int> tests, bool verbose, std::string &filename, bool symm_mode, bool metis,
                     bool merge, bool invert_diag, bool invert_offdiag, bool u_in_csr, bool trmm_on_device,
                     int panel_size, int relax_size, int block_size, int loop) {
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
  using host_crsmat_t = KokkosSparse::CrsMatrix<scalar_type, ordinal_type, host_execution_space, void, size_type>;
  using crsmat_t      = KokkosSparse::CrsMatrix<scalar_type, ordinal_type, execution_space, void, size_type>;

  //
  using graph_t = typename crsmat_t::StaticCrsGraphType;

  //
  using host_scalar_view_t = Kokkos::View<scalar_type *, host_memory_space>;
  using scalar_view_t      = Kokkos::View<scalar_type *, memory_space>;

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
    std::cout << " SuperLU Tester Begin: Read matrix filename " << filename << std::endl;
    host_crsmat_t Mtx = KokkosSparse::Impl::read_kokkos_crst_matrix<host_crsmat_t>(filename.c_str());

    const size_type nrows = Mtx.graph.numRows();

    auto graph_host   = Mtx.graph;  // in_graph
    auto row_map_host = graph_host.row_map;
    auto entries_host = graph_host.entries;
    auto values_host  = Mtx.values;

    // ==============================================
    // call SuperLU on the host
    // > data for SuperLU
    int *etree;
    int *perm_r, *perm_c;
    SuperMatrix L;
    SuperMatrix U;
    // > call SuperLU
    Kokkos::Timer timer;
    std::cout << " > call SuperLU for factorization" << std::endl;
    factor_superlu<scalar_type>(symm_mode, metis, nrows, values_host.data(), const_cast<int *>(row_map_host.data()),
                                entries_host.data(), panel_size, relax_size, L, U, &perm_r, &perm_c, &etree);
    double factor_time = timer.seconds();
    std::cout << "   Factorization Time: " << factor_time << std::endl << std::endl;

    // ==============================================
    // Run all requested algorithms
    for (auto test : tests) {
      std::cout << "\ntest = " << test << std::endl;

      KernelHandle khL, khU;
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
            khL.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_NAIVE, nrows, true);
            khU.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_NAIVE, nrows, false);
          } else if (test == SUPERNODAL_ETREE) {
            std::cout << " > create handle for SUPERNODAL_ETREE" << std::endl << std::endl;
            khL.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_ETREE, nrows, true);
            khU.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_ETREE, nrows, false);
          } else if (test == SUPERNODAL_DAG) {
            std::cout << " > create handle for SUPERNODAL_DAG" << std::endl << std::endl;
            khL.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_DAG, nrows, true);
            khU.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_DAG, nrows, false);
          }
          if (test == SUPERNODAL_SPMV) {
            std::cout << " > create handle for SUPERNODAL_SPMV" << std::endl << std::endl;
            khL.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_SPMV, nrows, true);
            khU.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_SPMV, nrows, false);
          } else if (test == SUPERNODAL_SPMV_DAG) {
            std::cout << " > create handle for SUPERNODAL_SPMV_DAG" << std::endl << std::endl;
            khL.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG, nrows, true);
            khU.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG, nrows, false);
          }
          // verbose (optional, default is false)
          khU.set_sptrsv_verbose(verbose);
          khL.set_sptrsv_verbose(verbose);

          // specify if U is stored in CSR or CSC
          std::cout << "=============================== " << std::endl;
          std::cout << " U in CSR           : " << u_in_csr << std::endl;
          khU.set_sptrsv_column_major(!u_in_csr);

          // specify wheather to merge supernodes (optional, default merge is
          // false)
          std::cout << " Merge Supernode    : " << merge << std::endl;
          khL.set_sptrsv_merge_supernodes(merge);
          khU.set_sptrsv_merge_supernodes(merge);

          // specify wheather to invert diagonal blocks
          std::cout << " Invert diagonal    : " << invert_diag << std::endl;
          khL.set_sptrsv_invert_diagonal(invert_diag);
          khU.set_sptrsv_invert_diagonal(invert_diag);

          // specify wheather to apply diagonal-inversion to off-diagonal blocks
          // (optional, default is false)
          std::cout << " Invert Off-diagonal: " << invert_offdiag << std::endl;
          khL.set_sptrsv_invert_offdiagonal(invert_offdiag);
          khU.set_sptrsv_invert_offdiagonal(invert_offdiag);

          // set etree (required if SUPERNODAL_ETREE)
          if (test == SUPERNODAL_ETREE || test == SUPERNODAL_SPMV) {
            khL.set_sptrsv_etree(etree);
            khU.set_sptrsv_etree(etree);
          }

          // set permutation
          khL.set_sptrsv_perm(perm_r);
          khU.set_sptrsv_perm(perm_c);

          // specify whether to run trmm on device
          std::cout << " TRMM on device     : " << trmm_on_device << std::endl;
          khL.set_sptrsv_trmm_on_device(trmm_on_device);
          khU.set_sptrsv_trmm_on_device(trmm_on_device);

          // block size to switch to device call
          if (block_size >= 0) {
            std::cout << " Block Size         : " << block_size << std::endl;
            khL.set_sptrsv_diag_supernode_sizes(block_size, block_size);
            khU.set_sptrsv_diag_supernode_sizes(block_size, block_size);
          }
          std::cout << std::endl;

          // ==============================================
          // do symbolic analysis (preprocssing, e.g., merging supernodes,
          // inverting diagonal/offdiagonal blocks, and scheduling based on
          // graph/dag)
          timer.reset();
          sptrsv_symbolic(&khL, &khU, L, U);
          double symbolic_time = timer.seconds();
          std::cout << "   Symbolic Time   : " << symbolic_time << std::endl << std::endl;

          // ==============================================
          // do numeric compute (copy numerical values from SuperLU data
          // structure to our sptrsv data structure)
          timer.reset();
          sptrsv_compute(&khL, &khU, L, U);
          double compute_time = timer.seconds();
          std::cout << "   Numeric Time   : " << compute_time << std::endl << std::endl;

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
          KokkosSparse::spmv("N", ONE, Mtx, sol_host, ZERO, rhs_host);

          // ==============================================
          // apply forward-pivot to rhs on the host
          host_scalar_view_t tmp_host("temp", nrows);
          forwardP_supernode<scalar_type>(nrows, perm_r, 1, rhs_host.data(), nrows, tmp_host.data(), nrows);

          // copy rhs to the default host/device
          scalar_view_t rhs("rhs", nrows);
          scalar_view_t sol("sol", nrows);
          Kokkos::deep_copy(rhs, tmp_host);

          // ==============================================
          // do L solve
          timer.reset();
          sptrsv_solve(&khL, sol, rhs);
          Kokkos::fence();
          double solveL_time = timer.seconds();
          std::cout << " > Lower-TRI: " << std::endl;
          std::cout << "   Solve Time   : " << solveL_time << std::endl;

          // ==============================================
          // do U solve
          timer.reset();
          sptrsv_solve(&khU, rhs, sol);
          Kokkos::fence();
          double solveU_time = timer.seconds();
          std::cout << " > Upper-TRI: " << std::endl;
          std::cout << "   Solve Time   : " << solveU_time << std::endl;

          // copy solution to host
          Kokkos::deep_copy(tmp_host, rhs);
          // apply backward-pivot
          backwardP_supernode<scalar_type>(nrows, perm_c, 1, tmp_host.data(), nrows, sol_host.data(), nrows);

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
            forwardP_supernode<scalar_type>(nrows, perm_r, 1, rhs_host.data(), nrows, tmp_host.data(), nrows);
            Kokkos::deep_copy(rhs, tmp_host);

#if 1
            sptrsv_solve(&khL, &khU, sol, rhs);
#else
            sptrsv_solve(&khL, sol, rhs);
            sptrsv_solve(&khU, rhs, sol);
#endif

            Kokkos::fence();
            Kokkos::deep_copy(tmp_host, rhs);
            backwardP_supernode<scalar_type>(nrows, perm_c, 1, tmp_host.data(), nrows, sol_host.data(), nrows);

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
          min_time = 0.0;
          max_time = 0.0;
          ave_time = 0.0;
          Kokkos::fence();
          for (int i = 0; i < loop; i++) {
            timer.reset();
            sptrsv_solve(&khU, rhs, sol);
            Kokkos::fence();
            double time = timer.seconds();
            ave_time += time;
            if (time > max_time || i == 0) max_time = time;
            if (time < min_time || i == 0) min_time = time;
          }
          std::cout << " U-solve: loop = " << loop << std::endl;
          std::cout << "  LOOP_AVG_TIME:  " << ave_time / loop << std::endl;
          std::cout << "  LOOP_MAX_TIME:  " << max_time << std::endl;
          std::cout << "  LOOP_MIN_TIME:  " << min_time << std::endl << std::endl;
        } break;

        case CUSPARSE: {
          std::cout << " > create handle for CuSparse (SUPERNODAL_NAIVE)" << std::endl << std::endl;
          khL.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_NAIVE, nrows, true);
          khU.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_NAIVE, nrows, false);

          khU.set_sptrsv_column_major(!u_in_csr);

          // ==============================================
          // read SuperLU factor on the host (and copy to default host/device)
          timer.reset();
          graph_t graphL;
          crsmat_t superluL;
          khL.set_sptrsv_invert_diagonal(false);
          graphL      = read_superlu_graphL<graph_t>(&khL, &L);
          double time = timer.seconds();
          timer.reset();
          superluL          = read_superlu_valuesL<crsmat_t>(&khL, &L, graphL);
          double readL_time = timer.seconds();
          std::cout << "   Conversion Time for L: " << time << " + " << readL_time << std::endl;

          timer.reset();
          graph_t graphU;
          crsmat_t superluU;
          khU.set_sptrsv_invert_diagonal(false);
          graphU = read_superlu_graphU<graph_t>(&khU, &L, &U);
          time   = timer.seconds();
          timer.reset();
          superluU          = read_superlu_valuesU<crsmat_t, graph_t>(&khU, &L, &U, graphU);
          double readU_time = timer.seconds();
          std::cout << "   Conversion Time for U: " << time << " + " << readU_time << std::endl;

          // remove zeros in L/U
          timer.reset();
          std::cout << "   Compress L-factor: " << std::endl;
          superluL = remove_zeros_crsmat(superluL);
          std::cout << "   Compress U-factor: " << std::endl;
          superluU             = remove_zeros_crsmat(superluU);
          double compress_time = timer.seconds();
          std::cout << "   Compression Time: " << compress_time << std::endl << std::endl;

          bool col_majorL = true;
          bool col_majorU = !u_in_csr;
          if (!check_cusparse(Mtx, col_majorL, superluL, col_majorU, superluU, perm_r, perm_c, tol, loop)) {
            num_failed++;
          }
        } break;

        default: std::cout << " > Invalid test ID < " << std::endl; exit(0);
      }
    }
    // free SuperLU data structures
    free_superlu(L, U, perm_r, perm_c, etree);
  }
  std::cout << std::endl << std::endl;

  return num_failed;
}

void print_help_sptrsv() {
  printf("Options:\n");
  printf("  --test [OPTION]         : Use different kernel implementations\n");
  printf("                            Options:\n");
  printf(
      "                            superlu-naive, superlu-etree, "
      "superlu-dag\n\n");
  printf(
      "  -f [file]               : Read in Matrix Market formatted text file "
      "'file'.\n");
  printf(
      "  --loop [LOOP]           : How many spmv to run to aggregate average "
      "time. \n");
  printf(
      "  --u-in-csc              : To store U-factor in CSC, needed for "
      "invert.\n");
  printf("  --invert-diag           : To invert diagonal blocks.\n");
  printf("  --invert-offdiag        : To apply inverse to off-diagonal blocks.\n");
  printf(
      "  --block-size [SIZE]     : To specify the threshold to switch device "
      "and bached kernel.\n");
  printf("  --scalar-type [s,d,c,z] :\n");
}

int main(int argc, char **argv) {
  std::vector<int> tests;
  std::string filename;

  int loop = 1;
  // use symmetric mode for SuperLU
  bool symm_mode = false;
  // use METIS before SuperLU
  bool metis = false;
  // merge supernodes
  bool merge = false;
  // invert diagonal of L-factor
  bool invert_diag = true;
  // apply invert of diagonal to offdiagonal
  bool invert_offdiag = false;
  // store U in CSR, or CSC
  bool u_in_csr = true;
  // specify whether to run KokkosKernels::trmm on device
  bool trmm_on_device = false;
  // block size to switch to device call (default is 100)
  int block_size = -1;
  // parameters for SuperLU (only affects factorization)
  int panel_size = sp_ienv(1);
  int relax_size = sp_ienv(2);
  // verbose
  bool verbose = true;
  // scalar type
  std::string char_scalar      = "d";
  std::string scalarTypeString = "(scalar_t = double)";

  if (argc == 1) {
    print_help_sptrsv();
    return 0;
  }

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "--test") == 0)) {
      i++;
      if ((strcmp(argv[i], "superlu-naive") == 0)) {
        tests.push_back(SUPERNODAL_NAIVE);
      } else if ((strcmp(argv[i], "superlu-etree") == 0)) {
        tests.push_back(SUPERNODAL_ETREE);
      } else if ((strcmp(argv[i], "superlu-dag") == 0)) {
        tests.push_back(SUPERNODAL_DAG);
      } else if ((strcmp(argv[i], "superlu-spmv") == 0)) {
        tests.push_back(SUPERNODAL_SPMV);
      } else if ((strcmp(argv[i], "superlu-spmv-dag") == 0)) {
        tests.push_back(SUPERNODAL_SPMV_DAG);
      } else if ((strcmp(argv[i], "cusparse") == 0)) {
        tests.push_back(CUSPARSE);
      } else {
        std::cerr << "Invalid --tests option: \"" << argv[i] << "\"" << std::endl;
        return -EINVAL;
      }
      continue;
    }
    if ((strcmp(argv[i], "-f") == 0)) {
      filename = argv[++i];
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
    if ((strcmp(argv[i], "--symm") == 0)) {
      symm_mode = true;
      continue;
    }
    if ((strcmp(argv[i], "--metis") == 0)) {
      metis = true;
      continue;
    }
    if ((strcmp(argv[i], "--merge") == 0)) {
      merge = true;
      continue;
    }
    if ((strcmp(argv[i], "--no-invert-diag") == 0)) {
      invert_diag = false;
      continue;
    }
    if ((strcmp(argv[i], "--invert-offdiag") == 0)) {
      invert_offdiag = true;
      continue;
    }
    if ((strcmp(argv[i], "--u-in-csc") == 0)) {
      u_in_csr = false;
      continue;
    }
    if ((strcmp(argv[i], "--trmm-on-device") == 0)) {
      trmm_on_device = true;
      continue;
    }
    if ((strcmp(argv[i], "--panel-size") == 0)) {
      panel_size = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "--relax-size") == 0)) {
      relax_size = atoi(argv[++i]);
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

  Kokkos::ScopeGuard kokkosScope(argc, argv);

  int total_errors = 0;
  if (char_scalar == "z") {
#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
    scalarTypeString = "(scalar_t = Kokkos::complex<double>)";
    total_errors     = test_sptrsv_perf<Kokkos::complex<double>>(tests, verbose, filename, symm_mode, metis, merge,
                                                             invert_diag, invert_offdiag, u_in_csr, trmm_on_device,
                                                             panel_size, relax_size, block_size, loop);
#else
    std::cout << std::endl << " KOKKOSKERNELS_INST_COMPLEX_DOUBLE  is not enabled ** " << std::endl << std::endl;
#endif
  } else if (char_scalar == "c") {
#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
    scalarTypeString = "(scalar_t = Kokkos::complex<float>)";
    total_errors     = test_sptrsv_perf<Kokkos::complex<float>>(tests, verbose, filename, symm_mode, metis, merge,
                                                            invert_diag, invert_offdiag, u_in_csr, trmm_on_device,
                                                            panel_size, relax_size, block_size, loop);
#else
    std::cout << std::endl << " KOKKOSKERNELS_INST_COMPLEX_FLOAT  is not enabled ** " << std::endl << std::endl;
#endif
  } else if (char_scalar == "d") {
#if defined(KOKKOSKERNELS_INST_DOUBLE)
    scalarTypeString = "(scalar_t = double)";
    total_errors =
        test_sptrsv_perf<double>(tests, verbose, filename, symm_mode, metis, merge, invert_diag, invert_offdiag,
                                 u_in_csr, trmm_on_device, panel_size, relax_size, block_size, loop);
#else
    std::cout << std::endl << " KOKKOSKERNELS_INST_DOUBLE  is not enabled ** " << std::endl << std::endl;
#endif
  } else if (char_scalar == "f") {
#if defined(KOKKOSKERNELS_INST_FLOAT)
    scalarTypeString = "(scalar_t = float)";
    total_errors =
        test_sptrsv_perf<float>(tests, verbose, filename, symm_mode, metis, merge, invert_diag, invert_offdiag,
                                u_in_csr, trmm_on_device, panel_size, relax_size, block_size, loop);
#else
    std::cout << std::endl << " KOKKOSKERNELS_INST_FLOAT  is not enabled ** " << std::endl << std::endl;
#endif
  }
  if (total_errors == 0)
    std::cout << "Kokkos::SPTRSV Test: Passed " << scalarTypeString << std::endl << std::endl;
  else
    std::cout << "Kokkos::SPTRSV Test: Failed (" << total_errors << " / " << 2 * tests.size() << " failed) "
              << scalarTypeString << std::endl
              << std::endl;

  return 0;
}
#else  // defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
int main() {
  std::cout << std::endl;
#if !defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
  std::cout << " ** SUPERLU NOT ENABLED **" << std::endl;
#endif
#if !defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
  std::cout << " ** SUPERNODAL SPTRSV NOT ENABLED **" << std::endl;
#endif
  std::cout << std::endl;

  return 1;
}
#endif

#else  // defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA ) &&
       // (!defined(KOKKOS_ENABLE_CUDA) || ( 8000 <= CUDA_VERSION ))

int main() {
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
