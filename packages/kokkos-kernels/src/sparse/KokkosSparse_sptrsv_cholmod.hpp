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

/// \file KokkosSparse_sptrsv.hpp
/// \brief Parallel sparse triangular solve
///
/// This file provides KokkosSparse::sptrsv.  This function performs a
/// local (no MPI) sparse triangular solve on matrices stored in
/// compressed row sparse ("Crs") format.

#ifndef KOKKOSSPARSE_SPTRSV_CHOLMOD_HPP_
#define KOKKOSSPARSE_SPTRSV_CHOLMOD_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_CHOLMOD
#include "cholmod.h"
#include "KokkosSparse_sptrsv_supernode.hpp"

namespace KokkosSparse {
namespace Experimental {


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* Auxiliary functions for symbolic analysis                                                 */

  /* ========================================================================================= */
  template <typename cholmod_int_type, typename graph_t, typename KernelHandle>
  graph_t read_cholmod_graphL(KernelHandle kernelHandle, cholmod_factor *L, cholmod_common *cm) {

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */
    size_t n = L->n;
    size_t nsuper = L->nsuper;     // # of supernodal columns
    cholmod_int_type *mb = (cholmod_int_type*)(L->pi);    // mb[s+1] - mb[s] = total number of rows in the s-th supernodes (diagonal+off-diagonal)
    cholmod_int_type *nb = (cholmod_int_type*)(L->super);
    cholmod_int_type *colptr = (cholmod_int_type*)(L->px);      // colptr
    cholmod_int_type *rowind = (cholmod_int_type*)(L->s);       // rowind

    bool ptr_by_column = false;
    int nnzA = colptr[nsuper] - colptr[0]; // overestimated if not block_diag
    return read_supernodal_graphL<graph_t> (kernelHandle, n, nsuper, nnzA, ptr_by_column, mb, nb, rowind);
  }


  /* ========================================================================================= */
  template <typename cholmod_int_type>
  void compute_etree_cholmod(cholmod_sparse *A, cholmod_common *cm, int **etree) {
    cholmod_factor *L;
    if (std::is_same<cholmod_int_type, long>::value == true) {
      L = cholmod_l_analyze (A, cm);
    } else if (std::is_same<cholmod_int_type, int>::value == true) {
      L = cholmod_analyze (A, cm);
    }

    size_t n = L->n;
    size_t nsuper = L->nsuper;      // # of supernodal columns
    cholmod_int_type *Iwork = (cholmod_int_type*)(cm->Iwork);
    cholmod_int_type *Parent = Iwork + (2*n); /* size nfsuper <= n [ */

    *etree = new int [nsuper];
    for (size_t ii = 0 ; ii < nsuper; ii++) (*etree)[ii] = Parent[ii];
  }


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* For symbolic analysis                                                                     */
  template <typename cholmod_int_type, typename KernelHandle>
  void sptrsv_symbolic(
      KernelHandle *kernelHandleL,
      KernelHandle *kernelHandleU,
      cholmod_factor *L,
      cholmod_common *cm)
  {
    // ===================================================================
    // load sptrsv-handles
    auto *handleL = kernelHandleL->get_sptrsv_handle ();
    auto *handleU = kernelHandleU->get_sptrsv_handle ();

    // load options
    int *etree = handleL->get_etree ();
    bool needEtree = (handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
                      handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_ETREE);
    if (needEtree && etree == nullptr) {
      std::cout << std::endl
                << " ** etree needs to be set before calling sptrsv_symbolic with SuperLU **"
                << std::endl << std::endl;
      return;
    }

    Kokkos::Timer timer;
    // ==============================================
    // extract CrsGraph from Cholmod
    using graph_t  = typename KernelHandle::SPTRSVHandleType::graph_t;
    auto graph = read_cholmod_graphL<cholmod_int_type, graph_t>(kernelHandleL, L, cm);
    auto row_map = graph.row_map;
    auto entries = graph.entries;

    // ==============================================
    // setup supnodal info 
    size_t nsuper = L->nsuper;
    cholmod_int_type *supercols = (cholmod_int_type*)(L->super);
    handleL->set_supernodes (nsuper, supercols, etree);
    handleU->set_supernodes (nsuper, supercols, etree);

    // ==============================================
    // symbolic for L-solve on the host
    timer.reset();
    sptrsv_symbolic (kernelHandleL, row_map, entries);
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    double timeL = timer.seconds ();
    std::cout << " > Lower-TRI: " << std::endl;
    std::cout << "   Symbolic Time: " << timeL << std::endl;
    #endif

    // ==============================================
    // symbolic for L^T-solve on the host
    timer.reset ();
    sptrsv_symbolic (kernelHandleU, row_map, entries);
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    double timeU = timer.seconds ();
    std::cout << " > Upper-TRI: " << std::endl;
    std::cout << "   Symbolic Time: " << timeU << std::endl;
    #endif

    // ==============================================
    // save graphs
    handleL->set_graph (graph);
    handleU->set_graph (graph);

    // ===================================================================
    handleU->set_symbolic_complete ();
    handleU->set_symbolic_complete ();
  }


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* Auxiliary functions for numeric computation                                               */

  /* ========================================================================================= */
  template <typename cholmod_int_type, typename crsmat_t, typename graph_t, typename KernelHandle>
  crsmat_t read_cholmod_factor(KernelHandle kernelHandle, cholmod_factor *L, cholmod_common *cm, graph_t &static_graph) {

    using values_view_t = typename crsmat_t::values_type::non_const_type;
    using scalar_t      = typename values_view_t::value_type;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */
    size_t n = L->n;
    size_t nsuper = L->nsuper;     // # of supernodal columns
    cholmod_int_type *mb = (cholmod_int_type*)(L->pi);    // mb[s+1] - mb[s] = total number of rows in the s-th supernodes (diagonal+off-diagonal)
    cholmod_int_type *nb = (cholmod_int_type*)(L->super);
    cholmod_int_type *colptr = (cholmod_int_type*)(L->px);      // colptr
    cholmod_int_type *rowind = (cholmod_int_type*)(L->s);       // rowind
    scalar_t *Lx = (scalar_t*)(L->x); // data

    bool unit_diag = false;
    bool ptr_by_column = false;
    return read_supernodal_valuesL<crsmat_t> (unit_diag, kernelHandle, n, nsuper,
                                              ptr_by_column, mb, nb, colptr, rowind, Lx, static_graph);
  }


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* For numeric computation                                                                   */
  template <typename cholmod_int_type, typename KernelHandle>
  void sptrsv_compute(
      KernelHandle *kernelHandleL,
      KernelHandle *kernelHandleU,
      cholmod_factor *L,
      cholmod_common *cm)
  {
    // ===================================================================
    // load sptrsv-handles
    auto *handleL = kernelHandleL->get_sptrsv_handle ();
    auto *handleU = kernelHandleU->get_sptrsv_handle ();

    if (!(handleL->is_symbolic_complete()) ||
        !(handleU->is_symbolic_complete())) {
      std::cout << std::endl
                << " ** needs to call sptrsv_symbolic before calling sptrsv_numeric **"
                << std::endl << std::endl;
      return;
    }

    // ==============================================
    // load crsGraph
    auto graph = handleL->get_graph ();

    // ==============================================
    // read numerical values of L from Cholmod
    using crsmat_t = typename KernelHandle::SPTRSVHandleType::crsmat_t;
    auto cholmodL = read_cholmod_factor<cholmod_int_type, crsmat_t> (kernelHandleL, L, cm, graph);

    // ==============================================
    // save crsmat
    handleL->set_crsmat (cholmodL);
    handleU->set_crsmat (cholmodL);

    // ===================================================================
    handleL->set_numeric_complete ();
    handleU->set_numeric_complete ();
  }

} // namespace Experimental
} // namespace KokkosSparse

#endif // KOKKOSKERNELS_ENABLE_TPL_CHOLMOD
#endif // KOKKOSSPARSE_SPTRSV_CHOLMOD_HPP_

