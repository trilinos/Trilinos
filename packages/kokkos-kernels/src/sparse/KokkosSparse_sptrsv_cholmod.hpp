/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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
/* For symbolic analysis                                                                     */

/* ========================================================================================= */
template <typename graph_t>
graph_t read_cholmod_graphL(bool cusparse, cholmod_factor *L, cholmod_common *cm) {

  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */
  int n = L->n;
  int nsuper = L->nsuper;     // # of supernodal columns
  int *mb = (int*)(L->pi);    // mb[s+1] - mb[s] = total number of rows in all the s-th supernodes (diagonal+off-diagonal)
  int *nb = (int*)(L->super);
  int *colptr = (int*)(L->px);      // colptr
  int *rowind = (int*)(L->s);       // rowind

  bool merge = false;
  bool ptr_by_column = false;

  return read_supernodal_graphL<graph_t> (cusparse, merge,
                                          n, nsuper, ptr_by_column, mb, nb, colptr, rowind);
}


/* ========================================================================================= */
void compute_etree_cholmod(cholmod_sparse *A, cholmod_common *cm, int **etree) {
  cholmod_factor *L;
  L = cholmod_analyze (A, cm);

  int n = L->n;
  int nsuper = L->nsuper;      // # of supernodal columns
  int *Iwork = (int*)(cm->Iwork);
  int *Parent = Iwork + (2*((size_t) n)); /* size nfsuper <= n [ */

  *etree = new int [nsuper];
  for (int ii = 0 ; ii < nsuper; ii++) (*etree)[ii] = Parent[ii];
}




/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* For numeric computation                                                                   */

/* ========================================================================================= */
template <typename crsmat_t, typename graph_t>
crsmat_t read_cholmod_factor(bool cusparse, bool invert_diag, cholmod_factor *L, cholmod_common *cm, graph_t &static_graph) {

  using values_view_t = typename crsmat_t::values_type::non_const_type;
  using scalar_t      = typename values_view_t::value_type;

  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */
  int n = L->n;
  int nsuper = L->nsuper;     // # of supernodal columns
  int *mb = (int*)(L->pi);    // mb[s+1] - mb[s] = total number of rows in all the s-th supernodes (diagonal+off-diagonal)
  int *nb = (int*)(L->super);
  int *colptr = (int*)(L->px);      // colptr
  int *rowind = (int*)(L->s);       // rowind
  scalar_t *Lx = (scalar_t*)(L->x); // data

  bool merge = false;
  bool invert_offdiag = false;
  bool unit_diag = false;
  bool ptr_by_column = false;
  return read_supernodal_valuesL<crsmat_t, graph_t, scalar_t> (cusparse, merge, invert_diag, invert_offdiag,
                                                               unit_diag, n, nsuper, ptr_by_column, mb, nb, colptr, rowind, Lx, static_graph);
}


} // namespace Experimental
} // namespace KokkosSparse

#endif // KOKKOSKERNELS_ENABLE_TPL_CHOLMOD
#endif // KOKKOSSPARSE_SPTRSV_CHOLMOD_HPP_

