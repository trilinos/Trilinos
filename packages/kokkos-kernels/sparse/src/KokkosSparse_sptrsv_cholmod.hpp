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

/// \file KokkosSparse_sptrsv.hpp
/// \brief Parallel sparse triangular solve
///
/// This file provides KokkosSparse::sptrsv.  This function performs a
/// local (no MPI) sparse triangular solve on matrices stored in
/// compressed row sparse ("Crs") format.

#ifndef KOKKOSSPARSE_SPTRSV_CHOLMOD_HPP_
#define KOKKOSSPARSE_SPTRSV_CHOLMOD_HPP_

#if defined(KOKKOSKERNELS_ENABLE_TPL_CHOLMOD) && defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)

#include "cholmod.h"
#include "KokkosSparse_Utils.hpp"
#include "KokkosSparse_sptrsv_supernode.hpp"

namespace KokkosSparse {
namespace Experimental {

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/* Auxiliary functions for symbolic analysis */

/* =========================================================================================
 */
template <typename cholmod_int_type, typename graph_t, typename KernelHandle>
graph_t read_cholmod_graphL(KernelHandle *kernelHandle, cholmod_factor *L, cholmod_common *cm) {
  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */
  size_t n             = L->n;
  size_t nsuper        = L->nsuper;                    // # of supernodal columns
  cholmod_int_type *mb = (cholmod_int_type *)(L->pi);  // mb[s+1] - mb[s] = total number of rows in the s-th
                                                       // supernodes (diagonal+off-diagonal)
  cholmod_int_type *nb     = (cholmod_int_type *)(L->super);
  cholmod_int_type *colptr = (cholmod_int_type *)(L->px);  // colptr
  cholmod_int_type *rowind = (cholmod_int_type *)(L->s);   // rowind

  bool ptr_by_column = false;
  if (kernelHandle->is_sptrsv_column_major()) {
    int nnzA = colptr[nsuper] - colptr[0];  // overestimated if not block_diag
    return read_supernodal_graphL<graph_t>(kernelHandle, n, nsuper, nnzA, ptr_by_column, mb, nb, rowind);
  } else {
    return read_supernodal_graphLt<graph_t>(kernelHandle, n, nsuper, ptr_by_column, mb, nb, rowind);
  }
}

/* =========================================================================================
 */
template <typename cholmod_int_type>
void compute_etree_cholmod(cholmod_sparse *A, cholmod_common *cm, int **etree) {
  cholmod_factor *L;
  if (std::is_same<cholmod_int_type, long>::value == true) {
    L = cholmod_l_analyze(A, cm);
  } else if (std::is_same<cholmod_int_type, int>::value == true) {
    L = cholmod_analyze(A, cm);
  }

  size_t n                 = L->n;
  size_t nsuper            = L->nsuper;  // # of supernodal columns
  cholmod_int_type *Iwork  = (cholmod_int_type *)(cm->Iwork);
  cholmod_int_type *Parent = Iwork + (2 * n); /* size nfsuper <= n [ */

  *etree = new int[nsuper];
  for (size_t ii = 0; ii < nsuper; ii++) (*etree)[ii] = Parent[ii];
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/* For symbolic analysis */
template <typename cholmod_int_type, typename KernelHandle>
void sptrsv_symbolic(KernelHandle *kernelHandleL, KernelHandle *kernelHandleU, cholmod_factor *L, cholmod_common *cm) {
  // ===================================================================
  // load sptrsv-handles
  auto *handleL = kernelHandleL->get_sptrsv_handle();
  auto *handleU = kernelHandleU->get_sptrsv_handle();

  // ==============================================
  // load supernodes
  int nsuper                  = L->nsuper;
  cholmod_int_type *supercols = (cholmod_int_type *)(L->super);
  // convert supercols into internal-view type
  using integer_view_host_t          = typename KernelHandle::SPTRSVHandleType::integer_view_host_t;
  integer_view_host_t supercols_view = integer_view_host_t("supercols", 1 + nsuper);
  for (int i = 0; i <= nsuper; i++) {
    supercols_view(i) = supercols[i];
  }

  // ==============================================
  // load etree (optional)
  int *etree = handleL->get_etree();

  // ==============================================
  // extract CrsGraph for L from Cholmod
  using host_graph_t = typename KernelHandle::SPTRSVHandleType::host_graph_t;
  auto graphL        = read_cholmod_graphL<cholmod_int_type, host_graph_t>(kernelHandleL, L, cm);

  if (handleU->is_column_major()) {
    // ==============================================
    // extract CrsGraph for U from Cholmod
    handleU->set_column_major(false);
    auto graphU = read_cholmod_graphL<cholmod_int_type, host_graph_t>(kernelHandleU, L, cm);
    handleU->set_column_major(true);

    // ==============================================
    // call supnodal symbolic
    sptrsv_supernodal_symbolic(nsuper, supercols_view.data(), etree, graphL, kernelHandleL, graphU, kernelHandleU);
  } else {
    // ==============================================
    // call supnodal symbolic
    sptrsv_supernodal_symbolic(nsuper, supercols_view.data(), etree, graphL, kernelHandleL, graphL, kernelHandleU);
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/* Auxiliary functions for numeric computation */

/* =========================================================================================
 */
template <typename cholmod_int_type, typename crsmat_t, typename graph_t, typename KernelHandle>
crsmat_t read_cholmod_factor(KernelHandle *kernelHandle, cholmod_factor *L, cholmod_common *cm, graph_t &static_graph) {
  using values_view_t = typename crsmat_t::values_type::non_const_type;
  using scalar_t      = typename values_view_t::value_type;

  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */
  size_t n      = L->n;
  size_t nsuper = L->nsuper;  // # of supernodal columns
  // mb[s+1] - mb[s] = total number of rows in the s-th supernodes
  // (diagonal+off-diagonal)
  cholmod_int_type *mb     = (cholmod_int_type *)(L->pi);
  cholmod_int_type *nb     = (cholmod_int_type *)(L->super);
  cholmod_int_type *colptr = (cholmod_int_type *)(L->px);
  cholmod_int_type *rowind = (cholmod_int_type *)(L->s);
  scalar_t *Lx             = (scalar_t *)(L->x);  // data

  bool ptr_by_column = false;
  if (kernelHandle->is_sptrsv_column_major()) {
    return read_supernodal_values<crsmat_t>(kernelHandle, n, nsuper, ptr_by_column, mb, nb, colptr, rowind, Lx,
                                            static_graph);
  } else {
    return read_supernodal_valuesLt<crsmat_t>(kernelHandle, n, nsuper, ptr_by_column, mb, nb, colptr, rowind, Lx,
                                              static_graph);
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/* For numeric computation */
template <typename cholmod_int_type, typename KernelHandle>
void sptrsv_compute(KernelHandle *kernelHandleL, KernelHandle *kernelHandleU, cholmod_factor *L, cholmod_common *cm) {
  // ==============================================
  // load sptrsv-handles
  auto *handleL = kernelHandleL->get_sptrsv_handle();
  auto *handleU = kernelHandleU->get_sptrsv_handle();

  if (!(handleL->is_symbolic_complete()) || !(handleU->is_symbolic_complete())) {
    std::cout << std::endl
              << " ** needs to call sptrsv_symbolic before calling sptrsv_numeric **" << std::endl
              << std::endl;
    return;
  }

  // ==============================================
  // load options
  bool useSpMV = (handleL->get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
                  handleL->get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG);

  // ==============================================
  // load crsGraph
  auto graph = handleL->get_graph();

  // ==============================================
  // read numerical values of L from Cholmod
  using crsmat_t = typename KernelHandle::SPTRSVHandleType::crsmat_t;
  auto cholmodL  = read_cholmod_factor<cholmod_int_type, crsmat_t>(kernelHandleL, L, cm, graph);

  // ==============================================
  // split the matrix into submatrices for spmv at each level
  if (useSpMV) {
    split_crsmat<crsmat_t>(kernelHandleL, cholmodL);
  }

  // ==============================================
  // save crsmat
  handleL->set_crsmat(cholmodL);
  if (handleU->is_column_major()) {
    auto graphU = handleU->get_graph();

    handleU->set_lower_tri(true);
    handleU->set_column_major(false);
    auto cholmodU = read_cholmod_factor<cholmod_int_type, crsmat_t>(kernelHandleU, L, cm, graphU);

    handleU->set_lower_tri(false);
    handleU->set_column_major(true);
    // ==============================================
    // split the matrix into submatrices for spmv at each level
    if (useSpMV) {
      split_crsmat<crsmat_t>(kernelHandleU, cholmodU);
    }
    handleU->set_crsmat(cholmodU);
  } else {
    handleU->set_crsmat(cholmodL);
    if (useSpMV) {
      if (!handleL->get_invert_offdiagonal()) {
        // copy submatrices to U for SpMV at each level
        auto nlevels = handleL->get_num_levels();
        std::vector<crsmat_t> sub_crsmats(nlevels);
        std::vector<crsmat_t> diag_blocks(nlevels);
        for (int lvl = 0; lvl < nlevels; lvl++) {
          sub_crsmats[lvl] = handleL->get_submatrix(nlevels - lvl - 1);
          diag_blocks[lvl] = handleL->get_diagblock(nlevels - lvl - 1);
        }
        handleU->set_submatrices(sub_crsmats);
        handleU->set_diagblocks(diag_blocks);
      } else {
        // not supported
      }
    }
  }

  // ==============================================
  handleL->set_numeric_complete();
  handleU->set_numeric_complete();
}

}  // namespace Experimental
}  // namespace KokkosSparse

#endif  // KOKKOSKERNELS_ENABLE_TPL_CHOLMOD &&
        // KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
#endif  // KOKKOSSPARSE_SPTRSV_CHOLMOD_HPP_
