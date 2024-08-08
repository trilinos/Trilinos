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

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include <cusparse_v2.h>
#endif

#include "KokkosBlas1_nrm2.hpp"
#include "KokkosSparse_spmv.hpp"

#ifndef KOKKOSSPARSE_SPTRSV_AUX
#define KOKKOSSPARSE_SPTRSV_AUX

namespace KokkosSparse {
namespace PerfTest {
namespace Experimental {

/* =========================================================================================
 */
template <typename scalar_t>
void forwardP_supernode(int n, int *perm_r, int nrhs, scalar_t *B, int ldb, scalar_t *X, int ldx) {
  /* Permute right hand sides to form Pr*B */
  for (int j = 0; j < nrhs; j++) {
    scalar_t *rhs_work = &B[j * ldb];
    scalar_t *sol_work = &X[j * ldx];
    for (int k = 0; k < n; k++) {
      sol_work[perm_r[k]] = rhs_work[k];
    }
  }
}

template <typename scalar_t>
void backwardP_supernode(int n, int *perm_c, int nrhs, scalar_t *B, int ldb, scalar_t *X, int ldx) {
  /* Compute the final solution X := Pc*X. */
  for (int j = 0; j < nrhs; j++) {
    scalar_t *rhs_work = &B[j * ldb];
    scalar_t *sol_work = &X[j * ldx];
    for (int k = 0; k < n; k++) {
      sol_work[k] = rhs_work[perm_c[k]];
    }
  }
}

/* =========================================================================================
 */
template <typename mag_t, typename crsmat_t, typename scalar_view_t>
bool check_errors(mag_t tol, crsmat_t &Mtx, scalar_view_t rhs, scalar_view_t sol) {
  using graph_t        = typename crsmat_t::StaticCrsGraphType;
  using entries_view_t = typename graph_t::entries_type::non_const_type;
  using lno_t          = typename entries_view_t::non_const_value_type;
  using values_view_t  = typename crsmat_t::values_type::non_const_type;
  using scalar_t       = typename values_view_t::value_type;
  using STS            = Kokkos::ArithTraits<scalar_t>;

  using execution_space = typename scalar_view_t::execution_space;

  const mag_t ZERO(0.0);
  const scalar_t ONE(1.0);

  // normA
  mag_t normA = ZERO;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<execution_space>(0, Mtx.nnz()),
      KOKKOS_LAMBDA(const lno_t i, mag_t &tsum) { tsum += STS::abs(Mtx.values(i)) * STS::abs(Mtx.values(i)); }, normA);
  normA = sqrt(normA);

  // normB
  mag_t normB = KokkosBlas::nrm2(rhs);

  // normX
  mag_t normX = KokkosBlas::nrm2(sol);

  // normR = ||B - AX||
  KokkosSparse::spmv("N", -ONE, Mtx, sol, ONE, rhs);
  mag_t normR = KokkosBlas::nrm2(rhs);

  std::cout << " > check : ||B - AX||/(||B|| + ||A||*||X||) = " << normR << "/(" << normB << " + " << normA << " * "
            << normX << ") = " << normR / (normB + normA * normX) << std::endl;

  const int nrows = Mtx.graph.numRows();
  return (normR / (mag_t(nrows) * (normB + normA * normX)) <= tol);
}

/* =========================================================================================
 */
template <typename crsmat_t>
void print_crsmat(crsmat_t &A) {
  auto graph   = A.graph;  // in_graph
  auto row_map = graph.row_map;
  auto entries = graph.entries;
  auto values  = A.values;

  int n = graph.numRows();
  std::cout << "[";
  for (int i = 0; i < n; i++) {
    for (int k = row_map[i]; k < row_map[i + 1]; k++) {
      std::cout << i << " " << entries[k] << " " << values[k] << " " << k << std::endl;
    }
  }
  std::cout << "];" << std::endl;
}

template <typename graph_t>
void print_graph(graph_t &graph) {
  auto row_map = graph.row_map;
  auto entries = graph.entries;

  int n = graph.numRows();
  std::cout << "[";
  for (int i = 0; i < n; i++) {
    for (int k = row_map[i]; k < row_map[i + 1]; k++) {
      std::cout << i << " " << entries[k] << " " << std::endl;
    }
  }
  std::cout << "];" << std::endl;
}

/* =========================================================================================
 */
template <typename crsmat_t>
crsmat_t remove_zeros_crsmat(crsmat_t &A) {
  using graph_t        = typename crsmat_t::StaticCrsGraphType;
  using row_map_view_t = typename graph_t::row_map_type::non_const_type;
  using cols_view_t    = typename graph_t::entries_type::non_const_type;
  using values_view_t  = typename crsmat_t::values_type::non_const_type;

  using row_map_view_host_t = typename row_map_view_t::HostMirror;
  using cols_view_host_t    = typename cols_view_t::HostMirror;
  using values_view_host_t  = typename values_view_t::HostMirror;
  using scalar_t            = typename values_view_t::value_type;
  using size_type           = typename crsmat_t::size_type;

  using range_type = Kokkos::pair<int, int>;

  const scalar_t zero(0.0);

  auto graph   = A.graph;  // in_graph
  int n        = graph.numRows();
  auto row_map = graph.row_map;
  auto entries = graph.entries;
  auto values  = A.values;

  row_map_view_host_t hr = Kokkos::create_mirror_view(row_map);
  cols_view_host_t hc    = Kokkos::create_mirror_view(entries);
  values_view_host_t hv  = Kokkos::create_mirror_view(values);

  Kokkos::deep_copy(hr, row_map);
  Kokkos::deep_copy(hc, entries);
  Kokkos::deep_copy(hv, values);

  // compress
  size_type nnzA0 = hr(n);
  size_type nnzA  = 0;
  for (int i = 0; i < n; i++) {
    size_type nnz = nnzA;
    for (int k = hr(i); k < hr(i + 1); k++) {
      if (hv(k) != zero) {
        hv(nnzA) = hv(k);
        hc(nnzA) = hc(k);
        nnzA++;
      }
    }
    hr(i) = nnz;
  }
  hr(n) = nnzA;
  std::cout << "   > compressed from " << nnzA0 << " to " << nnzA << " nnzs" << std::endl;

  // allocate & create
  row_map_view_t new_row_map("rowmap_view", n + 1);
  cols_view_t new_entries("colmap_view", nnzA);
  values_view_t new_values("values_view", nnzA);

  Kokkos::deep_copy(new_row_map, hr);
  Kokkos::deep_copy(new_entries, subview(hc, range_type(0, nnzA)));
  Kokkos::deep_copy(new_values, subview(hv, range_type(0, nnzA)));

  graph_t static_graph(new_entries, new_row_map);
  crsmat_t crsmat("CrsMatrix", n, new_values, static_graph);
  return crsmat;
}

/* =========================================================================================
 */
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
std::string getCuSparseErrorString(cusparseStatus_t status) {
#if 0
  return cusparseGetErrorString (status);
#else
  switch (status) {
    case CUSPARSE_STATUS_SUCCESS: return "CUSPARSE_STATUS_SUCCESS";
    case CUSPARSE_STATUS_NOT_INITIALIZED: return "CUSPARSE_STATUS_NOT_INITIALIZED";
    case CUSPARSE_STATUS_ALLOC_FAILED: return "CUSPARSE_STATUS_ALLOC_FAILED";
    case CUSPARSE_STATUS_INVALID_VALUE: return "CUSPARSE_STATUS_INVALID_VALUE";
    case CUSPARSE_STATUS_ARCH_MISMATCH: return "CUSPARSE_STATUS_ARCH_MISMATCH";
    case CUSPARSE_STATUS_EXECUTION_FAILED: return "USPARSE_STATUS_EXECUTION_FAILED";
    case CUSPARSE_STATUS_INTERNAL_ERROR: return "CUSPARSE_STATUS_INTERNAL_ERROR";
    default: return "un-handled error code";
  }
#endif
}
#endif

/* =========================================================================================
 */
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
#if CUSPARSE_VERSION >= 12500
template <typename crsmat_t, typename host_crsmat_t>
bool check_cusparse(host_crsmat_t &, bool, crsmat_t &, bool, crsmat_t &, int *, int *, double, int) {
  // TODO: call KokkosSparse::sptrsv (if hardcoded problem settings below are
  // compatible), or add wrappers for modern interface (cusparseSpSV*)
  throw std::logic_error("Legacy cuSPARSE csrsv interface not available.");
  return false;
}

#else

template <typename crsmat_t, typename host_crsmat_t>
bool check_cusparse(host_crsmat_t &Mtx, bool col_majorL, crsmat_t &L, bool col_majorU, crsmat_t &U, int *perm_r,
                    int *perm_c, double tol, int loop) {
  using values_view_t = typename crsmat_t::values_type::non_const_type;
  using scalar_t      = typename values_view_t::value_type;
  using size_type     = typename crsmat_t::size_type;

  using host_values_view_t = typename host_crsmat_t::values_type::non_const_type;

  using execution_space = typename values_view_t::execution_space;
  using memory_space    = typename execution_space::memory_space;

  using host_execution_space = typename host_values_view_t::execution_space;
  using host_memory_space    = typename host_execution_space::memory_space;

  using host_scalar_view_t = Kokkos::View<scalar_t *, host_memory_space>;
  using scalar_view_t      = Kokkos::View<scalar_t *, memory_space>;

  const scalar_t ZERO(0.0);
  const scalar_t ONE(1.0);

  Kokkos::Timer timer;
  const int nrows = Mtx.graph.numRows();

  // ==============================================
  // > create a handle
  cusparseStatus_t status;
  cusparseHandle_t handle = 0;
  status                  = cusparseCreate(&handle);
  if (CUSPARSE_STATUS_SUCCESS != status) {
    std::cout << " ** cusparseCreate failed with " << getCuSparseErrorString(status) << " ** " << std::endl;
  }
  cusparseSetPointerMode(handle,
                         CUSPARSE_POINTER_MODE_HOST);  // scalars are passed by reference on host

  // > create a empty info structure for L-solve (e.g., analysis results)
  csrsv2Info_t infoL = 0;
  status             = cusparseCreateCsrsv2Info(&infoL);
  if (CUSPARSE_STATUS_SUCCESS != status) {
    std::cout << " ** cusparseCreateCsrsv2Info failed with " << getCuSparseErrorString(status) << " ** " << std::endl;
  }

  // ==============================================
  // Preparing for L-solve
  // step 1: create a descriptor
  size_type nnzL = L.nnz();
  auto graphL    = L.graph;  // in_graph
  auto row_mapL  = graphL.row_map;
  auto entriesL  = graphL.entries;
  auto valuesL   = L.values;

  // NOTE: it is stored in CSC = UPPER + TRANSPOSE
  cusparseMatDescr_t descrL = 0;
  status                    = cusparseCreateMatDescr(&descrL);
  if (CUSPARSE_STATUS_SUCCESS != status) {
    std::cout << " ** cusparseCreateMatDescr failed with " << getCuSparseErrorString(status) << " ** " << std::endl;
  }
  cusparseSetMatIndexBase(descrL, CUSPARSE_INDEX_BASE_ZERO);
  cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_UPPER);
  cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_NON_UNIT);
  cusparseSetMatType(descrL, CUSPARSE_MATRIX_TYPE_GENERAL);

  // ==============================================
  // step 2: query how much memory used in csrsv2, and allocate the buffer
  // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
  int pBufferSize;
  void *pBufferL             = 0;
  cusparseOperation_t transL = (col_majorL ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
  if (std::is_same<scalar_t, double>::value) {
    cusparseDcsrsv2_bufferSize(handle, transL, nrows, nnzL, descrL, reinterpret_cast<double *>(valuesL.data()),
                               row_mapL.data(), entriesL.data(), infoL, &pBufferSize);
  } else {
    cusparseZcsrsv2_bufferSize(handle, transL, nrows, nnzL, descrL, reinterpret_cast<cuDoubleComplex *>(valuesL.data()),
                               row_mapL.data(), entriesL.data(), infoL, &pBufferSize);
  }
  cudaMalloc((void **)&pBufferL, pBufferSize);

  // ==============================================
  // step 3: analysis
  std::cout << "  Lower-Triangular" << std::endl;
  timer.reset();
  const cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
  if (std::is_same<scalar_t, double>::value) {
    status = cusparseDcsrsv2_analysis(handle, transL, nrows, nnzL, descrL, reinterpret_cast<double *>(valuesL.data()),
                                      row_mapL.data(), entriesL.data(), infoL, policy, pBufferL);
  } else {
    status = cusparseZcsrsv2_analysis(handle, transL, nrows, nnzL, descrL,
                                      reinterpret_cast<cuDoubleComplex *>(valuesL.data()), row_mapL.data(),
                                      entriesL.data(), infoL, policy, pBufferL);
  }
  double time_symbolic = timer.seconds();
  std::cout << "  Cusparse Symbolic Time: " << time_symbolic << std::endl;
  if (CUSPARSE_STATUS_SUCCESS != status) {
    std::cout << " ** cusparseZcsrsv2_analysis failed with " << getCuSparseErrorString(status) << " ** " << std::endl;
  }
  // L has unit diagonal, so no structural zero is reported.

  int structural_zero;
  status = cusparseXcsrsv2_zeroPivot(handle, infoL, &structural_zero);
  if (CUSPARSE_STATUS_ZERO_PIVOT == status) {
    printf("L(%d,%d) is missing\n", structural_zero, structural_zero);
  }

  // ==============================================
  // Preaparing for the first solve
  //> create the known solution and set to all 1's ** on host **
  host_scalar_view_t sol_host("sol_host", nrows);
  Kokkos::deep_copy(sol_host, ONE);

  // > create the rhs ** on host **
  // A*sol generates rhs: rhs is dense, use spmv
  host_scalar_view_t rhs_host("rhs_host", nrows);
  KokkosSparse::spmv("N", ONE, Mtx, sol_host, ZERO, rhs_host);

  // ==============================================
  // step 1: apply forward-pivot to rhs on the host
  host_scalar_view_t tmp_host("temp", nrows);
  forwardP_supernode<scalar_t>(nrows, perm_r, 1, rhs_host.data(), nrows, tmp_host.data(), nrows);

  // copy rhs to the default host/device
  scalar_view_t rhs("rhs", nrows);
  scalar_view_t sol("sol", nrows);
  Kokkos::deep_copy(rhs, tmp_host);

  // ==============================================
  // step 2: solve L*y = x
  Kokkos::fence();
  timer.reset();
  if (std::is_same<scalar_t, double>::value) {
    const double alpha = 1.0;
    status =
        cusparseDcsrsv2_solve(handle, transL, nrows, nnzL, &alpha, descrL, reinterpret_cast<double *>(valuesL.data()),
                              row_mapL.data(), entriesL.data(), infoL, reinterpret_cast<double *>(rhs.data()),
                              reinterpret_cast<double *>(sol.data()), policy, pBufferL);
  } else {
    const cuDoubleComplex alpha = make_cuDoubleComplex(1.0, 0.0);
    status                      = cusparseZcsrsv2_solve(handle, transL, nrows, nnzL, &alpha, descrL,
                                                        reinterpret_cast<cuDoubleComplex *>(valuesL.data()), row_mapL.data(),
                                                        entriesL.data(), infoL, reinterpret_cast<cuDoubleComplex *>(rhs.data()),
                                                        reinterpret_cast<cuDoubleComplex *>(sol.data()), policy, pBufferL);
  }
  Kokkos::fence();
  double time_solve = timer.seconds();
  std::cout << "  Cusparse Solve Time   : " << time_solve << std::endl;
  if (CUSPARSE_STATUS_SUCCESS != status) {
    std::cout << " ** cusparseZcsrsv2_solve failed with " << getCuSparseErrorString(status) << " ** " << std::endl;
  }
  // L has unit diagonal, so no numerical zero is reported.
  int numerical_zero;
  status = cusparseXcsrsv2_zeroPivot(handle, infoL, &numerical_zero);
  if (CUSPARSE_STATUS_ZERO_PIVOT == status) {
    printf("L(%d,%d) is zero\n", numerical_zero, numerical_zero);
  }

  // ==============================================
  // Preparing for U-solve
  size_type nnzU = U.nnz();
  auto graphU    = U.graph;  // in_graph
  auto row_mapU  = graphU.row_map;
  auto entriesU  = graphU.entries;
  auto valuesU   = U.values;

  // > create a empty info structure for U-solve (e.g., analysis results)
  csrsv2Info_t infoU = 0;
  status             = cusparseCreateCsrsv2Info(&infoU);
  if (CUSPARSE_STATUS_SUCCESS != status) {
    std::cout << " ** cusparseCreateCsrsv2Info failed with " << getCuSparseErrorString(status) << " ** " << std::endl;
  }

  // ==============================================
  // step 1: create a descriptor
  cusparseMatDescr_t descrU = 0;
  status                    = cusparseCreateMatDescr(&descrU);
  if (CUSPARSE_STATUS_SUCCESS != status) {
    std::cout << " ** cusparseCreateMatDescr create status error name " << getCuSparseErrorString(status) << " ** "
              << std::endl;
  }
  // NOTE: if CSR, UPPER+NO-TRANSPOSE, else LOWER+Trans
  cusparseSetMatIndexBase(descrU, CUSPARSE_INDEX_BASE_ZERO);
  if (col_majorU) {
    cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_LOWER);
  } else {
    cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
  }
  cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
  cusparseSetMatType(descrU, CUSPARSE_MATRIX_TYPE_GENERAL);

  // ==============================================
  // step 2: query how much memory used in csrsv2, and allocate the buffer
  // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
  void *pBufferU             = 0;
  cusparseOperation_t transU = (col_majorU ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE);
  if (std::is_same<scalar_t, double>::value) {
    cusparseDcsrsv2_bufferSize(handle, transU, nrows, nnzU, descrU, reinterpret_cast<double *>(valuesU.data()),
                               row_mapU.data(), entriesU.data(), infoU, &pBufferSize);
  } else {
    cusparseZcsrsv2_bufferSize(handle, transU, nrows, nnzU, descrU, reinterpret_cast<cuDoubleComplex *>(valuesU.data()),
                               row_mapU.data(), entriesU.data(), infoU, &pBufferSize);
  }
  cudaMalloc((void **)&pBufferU, pBufferSize);

  // ==============================================
  // step 3: analysis
  std::cout << std::endl << "  Upper-Triangular" << std::endl;
  timer.reset();
  if (std::is_same<scalar_t, double>::value) {
    status = cusparseDcsrsv2_analysis(handle, transU, nrows, nnzU, descrU, reinterpret_cast<double *>(valuesU.data()),
                                      row_mapU.data(), entriesU.data(), infoU, policy, pBufferU);
  } else {
    status = cusparseZcsrsv2_analysis(handle, transU, nrows, nnzU, descrU,
                                      reinterpret_cast<cuDoubleComplex *>(valuesU.data()), row_mapU.data(),
                                      entriesU.data(), infoU, policy, pBufferU);
  }
  time_symbolic = timer.seconds();
  std::cout << "  Cusparse Symbolic Time: " << time_symbolic << std::endl;
  if (CUSPARSE_STATUS_SUCCESS != status) {
    std::cout << " ** cusparseDcsrsv2_analysis failed with " << getCuSparseErrorString(status) << " ** " << std::endl;
  }
  status = cusparseXcsrsv2_zeroPivot(handle, infoU, &structural_zero);
  if (CUSPARSE_STATUS_ZERO_PIVOT == status) {
    printf("U(%d,%d) is missing\n", structural_zero, structural_zero);
  }

  // ==============================================
  // step 1: solve U*y = x
  timer.reset();
  if (std::is_same<scalar_t, double>::value) {
    const double alpha = 1.0;
    status =
        cusparseDcsrsv2_solve(handle, transU, nrows, nnzU, &alpha, descrU, reinterpret_cast<double *>(valuesU.data()),
                              row_mapU.data(), entriesU.data(), infoU, reinterpret_cast<double *>(sol.data()),
                              reinterpret_cast<double *>(rhs.data()), policy, pBufferU);
  } else {
    const cuDoubleComplex alpha = make_cuDoubleComplex(1.0, 0.0);
    status                      = cusparseZcsrsv2_solve(handle, transU, nrows, nnzU, &alpha, descrU,
                                                        reinterpret_cast<cuDoubleComplex *>(valuesU.data()), row_mapU.data(),
                                                        entriesU.data(), infoU, reinterpret_cast<cuDoubleComplex *>(sol.data()),
                                                        reinterpret_cast<cuDoubleComplex *>(rhs.data()), policy, pBufferU);
  }
  Kokkos::fence();
  time_solve = timer.seconds();
  std::cout << "  Cusparse Solve Time   : " << time_solve << std::endl;
  if (CUSPARSE_STATUS_SUCCESS != status) {
    std::cout << " ** usparseDcsrsv2_solve failed with " << getCuSparseErrorString(status) << " ** " << std::endl;
  }
  status = cusparseXcsrsv2_zeroPivot(handle, infoU, &numerical_zero);
  if (CUSPARSE_STATUS_ZERO_PIVOT == status) {
    printf("U(%d,%d) is zero\n", numerical_zero, numerical_zero);
  }

  // ==============================================
  // copy solution to host
  Kokkos::deep_copy(tmp_host, rhs);
  // apply backward-pivot
  backwardP_supernode<scalar_t>(nrows, perm_c, 1, tmp_host.data(), nrows, sol_host.data(), nrows);

  // ==============================================
  // Error Check ** on host **
  Kokkos::fence();
  std::cout << std::endl;
  bool success = check_errors(tol, Mtx, rhs_host, sol_host);

  // Try again?
  if (success) {
    // reinitialize rhs
    Kokkos::deep_copy(sol_host, ONE);
    KokkosSparse::spmv("N", ONE, Mtx, sol_host, ZERO, rhs_host);

    // forward pivot
    forwardP_supernode<scalar_t>(nrows, perm_r, 1, rhs_host.data(), nrows, tmp_host.data(), nrows);

    // copy & solve & copy back
    Kokkos::deep_copy(rhs, tmp_host);
    if (std::is_same<scalar_t, double>::value) {
      const double alpha = 1.0;
      cusparseDcsrsv2_solve(handle, transL, nrows, nnzL, &alpha, descrL, reinterpret_cast<double *>(valuesL.data()),
                            row_mapL.data(), entriesL.data(), infoL, reinterpret_cast<double *>(rhs.data()),
                            reinterpret_cast<double *>(sol.data()), policy, pBufferL);
      cusparseDcsrsv2_solve(handle, transU, nrows, nnzU, &alpha, descrU, reinterpret_cast<double *>(valuesU.data()),
                            row_mapU.data(), entriesU.data(), infoU, reinterpret_cast<double *>(sol.data()),
                            reinterpret_cast<double *>(rhs.data()), policy, pBufferU);
    } else {
      const cuDoubleComplex alpha = make_cuDoubleComplex(1.0, 0.0);
      cusparseZcsrsv2_solve(handle, transL, nrows, nnzL, &alpha, descrL,
                            reinterpret_cast<cuDoubleComplex *>(valuesL.data()), row_mapL.data(), entriesL.data(),
                            infoL, reinterpret_cast<cuDoubleComplex *>(rhs.data()),
                            reinterpret_cast<cuDoubleComplex *>(sol.data()), policy, pBufferL);
      cusparseZcsrsv2_solve(handle, transU, nrows, nnzU, &alpha, descrU,
                            reinterpret_cast<cuDoubleComplex *>(valuesU.data()), row_mapU.data(), entriesU.data(),
                            infoU, reinterpret_cast<cuDoubleComplex *>(sol.data()),
                            reinterpret_cast<cuDoubleComplex *>(rhs.data()), policy, pBufferU);
    }
    Kokkos::deep_copy(tmp_host, rhs);

    // backward pivot and check
    backwardP_supernode<scalar_t>(nrows, perm_c, 1, tmp_host.data(), nrows, sol_host.data(), nrows);
    success = check_errors(tol, Mtx, rhs_host, sol_host);
  }
  std::cout << std::endl;

  if (success) {
    // ==============================================
    // Benchmark
    // L-solve
    double min_time = 0.0;
    double max_time = 0.0;
    double ave_time = 0.0;
    Kokkos::fence();
    for (int i = 0; i < loop; i++) {
      double time;
      if (std::is_same<scalar_t, double>::value) {
        const double alpha = 1.0;
        timer.reset();
        cusparseDcsrsv2_solve(handle, transL, nrows, nnzL, &alpha, descrL, reinterpret_cast<double *>(valuesL.data()),
                              row_mapL.data(), entriesL.data(), infoL, reinterpret_cast<double *>(rhs.data()),
                              reinterpret_cast<double *>(sol.data()), policy, pBufferL);
        Kokkos::fence();
        time = timer.seconds();
      } else {
        const cuDoubleComplex alpha = make_cuDoubleComplex(1.0, 0.0);
        timer.reset();
        cusparseZcsrsv2_solve(handle, transL, nrows, nnzL, &alpha, descrL,
                              reinterpret_cast<cuDoubleComplex *>(valuesL.data()), row_mapL.data(), entriesL.data(),
                              infoL, reinterpret_cast<cuDoubleComplex *>(rhs.data()),
                              reinterpret_cast<cuDoubleComplex *>(sol.data()), policy, pBufferL);
        Kokkos::fence();
        time = timer.seconds();
      }
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
      double time;
      if (std::is_same<scalar_t, double>::value) {
        double alpha = 1.0;
        timer.reset();
        cusparseDcsrsv2_solve(handle, transU, nrows, nnzU, &alpha, descrU, reinterpret_cast<double *>(valuesU.data()),
                              row_mapU.data(), entriesU.data(), infoU, reinterpret_cast<double *>(sol.data()),
                              reinterpret_cast<double *>(rhs.data()), policy, pBufferU);
        Kokkos::fence();
        time = timer.seconds();
      } else {
        const cuDoubleComplex alpha = make_cuDoubleComplex(1.0, 0.0);
        timer.reset();
        cusparseZcsrsv2_solve(handle, transU, nrows, nnzU, &alpha, descrU,
                              reinterpret_cast<cuDoubleComplex *>(valuesU.data()), row_mapU.data(), entriesU.data(),
                              infoU, reinterpret_cast<cuDoubleComplex *>(sol.data()),
                              reinterpret_cast<cuDoubleComplex *>(rhs.data()), policy, pBufferU);
        Kokkos::fence();
        time = timer.seconds();
      }
      ave_time += time;
      if (time > max_time || i == 0) max_time = time;
      if (time < min_time || i == 0) min_time = time;
    }
    std::cout << " U-solve: loop = " << loop << std::endl;
    std::cout << "  LOOP_AVG_TIME:  " << ave_time / loop << std::endl;
    std::cout << "  LOOP_MAX_TIME:  " << max_time << std::endl;
    std::cout << "  LOOP_MIN_TIME:  " << min_time << std::endl << std::endl;
  }
  return success;
}
#endif

#else
template <typename crsmat_t, typename host_crsmat_t>
bool check_cusparse(host_crsmat_t & /*Mtx*/, bool /*col_majorL*/, crsmat_t & /*L*/, bool /*col_majorU*/,
                    crsmat_t & /*U*/, int * /*perm_r*/, int * /*perm_c*/, double /*tol*/, int /*loop*/) {
  printf(" KOKKOSKERNELS_ENABLE_TPL_CUSPARSE **not** enabled\n");
  return false;
}
#endif

}  // namespace Experimental
}  // namespace PerfTest
}  // namespace KokkosSparse
#endif
