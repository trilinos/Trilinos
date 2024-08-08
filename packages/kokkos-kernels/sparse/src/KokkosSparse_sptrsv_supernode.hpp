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

#ifndef KOKKOSSPARSE_SPTRSV_SUPERNODE_HPP_
#define KOKKOSSPARSE_SPTRSV_SUPERNODE_HPP_

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)

#include "KokkosBlas3_trmm.hpp"
#include "KokkosLapack_trtri.hpp"

#include "KokkosBatched_Trtri_Decl.hpp"
#include "KokkosBatched_Trtri_Serial_Impl.hpp"

#include "KokkosBatched_Trmm_Decl.hpp"
#include "KokkosBatched_Trmm_Serial_Impl.hpp"

#include "KokkosSparse_SortCrs.hpp"
#include "KokkosSparse_sptrsv.hpp"

namespace KokkosSparse {
namespace Experimental {

template <typename ordinal_type>
class sort_indices {
 public:
  sort_indices(ordinal_type *rowinds) : rowinds_(rowinds) {}
  bool operator()(int i, int j) const { return rowinds_[i] < rowinds_[j]; }

 private:
  ordinal_type *rowinds_;  // rowindices
};

/* =========================================================================================
 */
template <typename host_graph_t, typename graph_t>
graph_t deep_copy_graph(host_graph_t &host_graph) {
  // load graph on host
  auto row_map = host_graph.row_map;
  auto entries = host_graph.entries;
  auto nrows   = host_graph.numRows();
  auto nnz     = row_map(nrows);

  // create graph on device
  using row_map_view_t = typename graph_t::row_map_type::non_const_type;
  using cols_view_t    = typename graph_t::entries_type::non_const_type;
  row_map_view_t rowmap_view("rowmap_view", nrows + 1);
  cols_view_t column_view("colmap_view", nnz);

  // copy graph to device
  Kokkos::deep_copy(rowmap_view, row_map);
  Kokkos::deep_copy(column_view, entries);
  graph_t static_graph(column_view, rowmap_view);
  return static_graph;
}

/* =========================================================================================
 */
template <typename graph_t, typename ptr_type, typename size_type, typename ordinal_type, typename KernelHandle>
graph_t read_supernodal_graphL(KernelHandle *kernelHandle, int n, int nsuper, int nnzA, bool ptr_by_column,
                               ptr_type *mb, size_type *nb, ordinal_type *rowind) {
  using row_map_view_t      = typename graph_t::row_map_type::non_const_type;
  using cols_view_t         = typename graph_t::entries_type::non_const_type;
  using integer_view_host_t = Kokkos::View<ordinal_type *, Kokkos::HostSpace>;

  // load parameters
  auto *handle = kernelHandle->get_sptrsv_handle();
  bool merge   = handle->get_merge_supernodes();

  row_map_view_t rowmap_view("rowmap_view", n + 1);
  cols_view_t column_view("colmap_view", nnzA);
  auto hr = Kokkos::create_mirror_view(rowmap_view);
  auto hc = Kokkos::create_mirror_view(column_view);

  // compute offset for each row
  int j               = 0;
  int max_nnz_per_row = 0;
  hr(j)               = 0;
  for (int s = 0; s < nsuper; s++) {
    int j1 = nb[s];
    int j2 = nb[s + 1];
    // number of columns in the s-th supernode column
    int nscol = j2 - j1;

    int i1, i2;
    if (ptr_by_column) {
      i1 = mb[j1];
      i2 = mb[j1 + 1];
    } else {
      i1 = mb[s];
      i2 = mb[s + 1];
    }
    // "total" number of rows in all the supernodes (diagonal+off-diagonal)
    int nsrow = i2 - i1;

    for (int jj = 0; jj < nscol; jj++) {
      hr(j + 1) = hr(j) + nsrow;
      j++;
    }
    if (nsrow > max_nnz_per_row) {
      max_nnz_per_row = nsrow;
    }
  }

  integer_view_host_t sorted_rowind_view("sorted_rowind", max_nnz_per_row + 1);
  ordinal_type *sorted_rowind = sorted_rowind_view.data();
  // store L in csr
  for (int s = 0; s < nsuper; s++) {
    int j1    = nb[s];
    int j2    = nb[s + 1];
    int nscol = j2 - j1;  // number of columns in the s-th supernode column

    int i1, i2;
    if (ptr_by_column) {
      i1 = mb[j1];
      i2 = mb[j1 + 1];
    } else {
      i1 = mb[s];
      i2 = mb[s + 1];
    }
    int nsrow = i2 - i1;         // "total" number of rows in all the supernodes
                                 // (diagonal+off-diagonal)
    int nsrow2 = nsrow - nscol;  // "total" number of rows in all the off-diagonal supernodes
    int ps2    = i1 + nscol;     // offset into rowind

    /* diagonal block */
    for (int ii = 0; ii < nscol; ii++) {
      // lower-triangular part
      for (int jj = 0; jj < ii; jj++) {
        hc(hr(j1 + jj)) = j1 + ii;
        hr(j1 + jj)++;
      }
      // diagonal
      hc(hr(j1 + ii)) = j1 + ii;
      hr(j1 + ii)++;
      // explicitly store zeros in upper-part
      for (int jj = ii + 1; jj < nscol; jj++) {
        hc(hr(j1 + jj)) = j1 + ii;
        hr(j1 + jj)++;
      }
    }

    /* off-diagonal blocks */
    if (merge) {
      // sort rowind (to merge supernodes)
      for (int ii = 0; ii < nsrow2; ii++) {
        sorted_rowind[ii] = ii;
      }
      std::sort(&(sorted_rowind[0]), &(sorted_rowind[nsrow2]), sort_indices<ordinal_type>(&rowind[ps2]));
    }
    for (int kk = 0; kk < nsrow2; kk++) {
      int ii = (merge ? sorted_rowind[kk] : kk);  // sorted rowind
      int i  = rowind[ps2 + ii];
      for (int jj = 0; jj < nscol; jj++) {
        hc(hr(j1 + jj)) = i;
        hr(j1 + jj)++;
      }
    }
  }

  // fix hr
  for (int i = n; i >= 1; i--) {
    hr(i) = hr(i - 1);
  }
  hr(0) = 0;

#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  std::cout << "    * Matrix size = " << n << std::endl;
  std::cout << "    * Total nnz   = " << hr(n) << std::endl;
  std::cout << "    * nnz / n     = " << hr(n) / n << std::endl;
#endif

  // deepcopy
  Kokkos::deep_copy(rowmap_view, hr);
  Kokkos::deep_copy(column_view, hc);

  // create crs
  graph_t static_graph(column_view, rowmap_view);
  return static_graph;
}

/* =========================================================================================
 */
template <typename graph_t, typename ptr_type, typename size_type, typename ordinal_type, typename KernelHandle>
graph_t read_supernodal_graphLt(KernelHandle *kernelHandle, int n, int nsuper, bool ptr_by_column, ptr_type *mb,
                                size_type *nb, ordinal_type *rowind) {
  using row_map_view_t      = typename graph_t::row_map_type::non_const_type;
  using cols_view_t         = typename graph_t::entries_type::non_const_type;
  using integer_view_host_t = Kokkos::View<int *, Kokkos::HostSpace>;

  // load parameters
  auto *handle = kernelHandle->get_sptrsv_handle();
  bool merge   = handle->get_merge_supernodes();

  /* create a map from row id to supernode id */
  integer_view_host_t map("map", n);
  int supid = 0;
  for (int k = 0; k < nsuper; k++) {
    int j1 = nb[k];
    int j2 = nb[k + 1];
    for (int j = j1; j < j2; j++) {
      map(j) = supid;
    }
    supid++;
  }

  row_map_view_t rowmap_view("rowmap_view", n + 1);
  auto hr = Kokkos::create_mirror_view(rowmap_view);
  Kokkos::deep_copy(hr, 0);

  integer_view_host_t sup("sup", nsuper);
  integer_view_host_t check("check", nsuper);
  Kokkos::deep_copy(check, 0);

  // compute offset for each row
  int nnzA = 0;
  for (int s = 0; s < nsuper; s++) {
    int j1 = nb[s];
    int j2 = nb[s + 1];
    // number of columns in the s-th supernode column
    int nscol = j2 - j1;

    int i1, i2;
    if (ptr_by_column) {
      i1 = mb[j1];
      i2 = mb[j1 + 1];
    } else {
      i1 = mb[s];
      i2 = mb[s + 1];
    }

    // diagonal blocks
    for (int ii = j1; ii < j2; ii++) {
      hr(ii + 1) += nscol;
      nnzA += nscol;
    }

    // offdiagonal blocks
    int nsup   = 0;
    int ps2    = i1 + nscol;  // offset into rowind
    int nsrow  = i2 - i1;
    int nsrow2 = nsrow - nscol;
    for (int kk = 0; kk < nsrow2; kk++) {
      int irow = rowind[ps2 + kk];
      supid    = map(irow);
      if (check(supid) == 0) {
        for (int ii = nb[supid]; ii < nb[supid + 1]; ii++) {
          hr(ii + 1) += nscol;
          nnzA += nscol;
        }
        check(supid) = 1;
        sup(nsup)    = supid;
        nsup++;
      }
    }
    // reset check
    for (int i = 0; i < nsup; i++) {
      check(sup(i)) = 0;
    }
  }
  for (int i = 0; i < n; i++) {
    hr(i + 1) += hr(i);
  }
  cols_view_t column_view("colmap_view", nnzA);
  auto hc = Kokkos::create_mirror_view(column_view);

  // pointer to off-diagonals (diagonal comes first)
  integer_view_host_t off("off", 1 + n);
  for (int s = 0; s < nsuper; s++) {
    int i1 = nb[s];
    int i2 = nb[s + 1];
    // number of columns in the s-th supernode column
    int nscol = i2 - i1;
    for (int ii = i1; ii < i2; ii++) {
      off(ii) = hr(ii) + nscol;
    }
  }

  // store L in csr
  for (int s = 0; s < nsuper; s++) {
    int j1    = nb[s];
    int j2    = nb[s + 1];
    int nscol = j2 - j1;  // number of columns in the s-th supernode column

    int i1, i2;
    if (ptr_by_column) {
      i1 = mb[j1];
      i2 = mb[j1 + 1];
    } else {
      i1 = mb[s];
      i2 = mb[s + 1];
    }
    int nsrow = i2 - i1;         // "total" number of rows in all the supernodes
                                 // (diagonal+off-diagonal)
    int nsrow2 = nsrow - nscol;  // "total" number of rows in all the off-diagonal supernodes
    int ps2    = i1 + nscol;     // offset into rowind

    /* diagonal block */
    for (int ii = 0; ii < nscol; ii++) {
      // explicitly store zeros in upper-part
      for (int jj = 0; jj < nscol; jj++) {
        hc(hr(j1 + ii) + jj) = j1 + jj;
      }
    }

    /* off-diagonal blocks */
    int nsup = 0;
    for (int kk = 0; kk < nsrow2; kk++) {
      int irow = rowind[ps2 + kk];
      supid    = map(irow);
      if (check(supid) == 0) {
        for (int ii = nb[supid]; ii < nb[supid + 1]; ii++) {
          for (int jj = 0; jj < nscol; jj++) {
            hc(off(ii) + jj) = j1 + jj;
          }
          off(ii) += nscol;
        }
        check(supid) = 1;
        sup(nsup)    = supid;
        nsup++;
      }
    }
    // reset check
    for (int i = 0; i < nsup; i++) {
      check(sup(i)) = 0;
    }
  }

  if (merge) {
    // they should be sorted?
  }

#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  std::cout << "    * Matrix size = " << n << std::endl;
  std::cout << "    * Total nnz   = " << hr(n) << std::endl;
  std::cout << "    * nnz / n     = " << hr(n) / n << std::endl;
#endif

  // deepcopy
  Kokkos::deep_copy(rowmap_view, hr);
  Kokkos::deep_copy(column_view, hc);

  // create crs
  graph_t static_graph(column_view, rowmap_view);
  return static_graph;
}

/* =========================================================================================
 */
template <typename input_graph_t, typename input_size_type>
void check_supernode_sizes(const char *title, int n, int nsuper, input_size_type *nb, input_graph_t &graph) {
  auto rowmap_view = graph.row_map;
  auto hr          = Kokkos::create_mirror_view(rowmap_view);
  Kokkos::deep_copy(hr, rowmap_view);

  int min_nsrow = 0, max_nsrow = 0, tot_nsrow = 0;
  int min_nscol = 0, max_nscol = 0, tot_nscol = 0;
  for (int s = 0; s < nsuper; s++) {
    int j1 = nb[s];
    int j2 = nb[s + 1];

    int nscol = j2 - j1;
    int nsrow = hr(j1 + 1) - hr(j1);

    if (s == 0) {
      min_nscol = max_nscol = tot_nscol = nscol;
      min_nsrow = max_nsrow = tot_nsrow = nsrow;
    } else {
      if (min_nsrow > nsrow) {
        min_nsrow = nsrow;
      }
      if (max_nsrow < nsrow) {
        max_nsrow = nsrow;
      }
      tot_nsrow += nsrow;

      if (min_nscol > nscol) {
        min_nscol = nscol;
      }
      if (max_nscol < nscol) {
        max_nscol = nscol;
      }
      tot_nscol += nscol;
    }
  }
  std::cout << std::endl << " ------------------------------------- " << std::endl << std::endl;
  std::cout << " " << title << std::endl;
  std::cout << "  + nsuper = " << nsuper << std::endl;
  std::cout << "  > nsrow: min = " << min_nsrow << ", max = " << max_nsrow << ", avg = " << tot_nsrow / nsuper
            << std::endl;
  std::cout << "  > nscol: min = " << min_nscol << ", max = " << max_nscol << ", avg = " << tot_nscol / nsuper
            << std::endl;
  std::cout << "    + Matrix size = " << n << std::endl;
  std::cout << "    + Total nnz   = " << hr(n) << std::endl;
  std::cout << "    + nnz / n     = " << hr(n) / n << std::endl;
}

/* =========================================================================================
 */
template <typename host_graph_t, typename graph_t, typename input_size_type>
host_graph_t generate_supernodal_graph(bool col_major, graph_t &graph, int nsuper, const input_size_type *nb) {
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  double time_seconds = 0.0;
  Kokkos::Timer timer;
#endif

  using size_type           = typename graph_t::size_type;
  using cols_view_host_t    = typename host_graph_t::entries_type::non_const_type;
  using row_map_view_host_t = typename host_graph_t::row_map_type::non_const_type;
  using integer_view_host_t = Kokkos::View<int *, Kokkos::HostSpace>;

  int n        = graph.numRows();
  auto row_map = graph.row_map;
  auto entries = graph.entries;

  auto row_map_host = Kokkos::create_mirror_view(row_map);
  auto entries_host = Kokkos::create_mirror_view(entries);
  Kokkos::deep_copy(row_map_host, row_map);
  Kokkos::deep_copy(entries_host, entries);

  // map col/row to supernode
  integer_view_host_t map("map", n);
  for (int s = 0; s < nsuper; s++) {
    for (int j = nb[s]; j < nb[s + 1]; j++) {
      map(j) = s;
    }
  }

  // count non-empty supernodal blocks
  row_map_view_host_t hr("rowmap_view", nsuper + 1);
  integer_view_host_t check("check", nsuper);
  integer_view_host_t idxs("idxs", nsuper);
  Kokkos::deep_copy(hr, 0);
  Kokkos::deep_copy(check, 0);

#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  timer.reset();
#endif
  int nblocks = 0;
  for (int s = 0; s < nsuper; s++) {
    int j1 = nb[s];
    int j2 = j1 + 1;  // based on the first row

    size_type nidxs = 0;
    for (size_type i = row_map_host(j1); i < row_map_host(j2); i++) {
      int s2 = map(entries_host(i));
      // supernodal blocks may not be filled with zeros
      // so need to check by each row
      // (also rowids are not sorted)
      if (check(s2) == 0) {
        check(s2) = 1;
        nblocks++;
        // count blocks per row for col_major
        hr(s2 + 1)++;
        // keep track of non-zero block ids
        idxs(nidxs) = s2;
        nidxs++;
      }
    }
    // reset check
    // Kokkos::deep_copy (check, 0);
    for (size_type i = 0; i < nidxs; i++) {
      check(idxs(i)) = 0;
    }
  }

  cols_view_host_t hc("colmap_view", nblocks);
  if (col_major) {
    // convert to offset
    for (int s = 0; s < nsuper; s++) {
      hr(s + 1) += hr(s);
    }
  }
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time_seconds = timer.seconds();
  std::cout << "   > Generate Supernodal Graph: count blocks   : " << time_seconds << std::endl;
  timer.reset();
#endif

  nblocks = 0;
  for (int s = 0; s < nsuper; s++) {
    int j1 = nb[s];
    int j2 = j1 + 1;  // based on the first row

    size_type nidxs = 0;
    for (size_type i = row_map_host(j1); i < row_map_host(j2); i++) {
      int s2 = map(entries_host(i));
      // supernodal blocks may not be filled with zeros
      // so need to check by each row
      // (also rowids are not sorted)
      if (check(s2) == 0) {
        check(s2) = 1;
        if (col_major) {
          hc(hr(s2)) = s;
          hr(s2)++;
        } else {
          hc(nblocks) = s2;
        }
        nblocks++;
        // keep track of non-zero block ids
        idxs(nidxs) = s2;
        nidxs++;
      }
    }
    if (!col_major) {
      hr(s + 1) = nblocks;
    }
    // reset check
    /*if (!col_major) {
      for (size_type s2 = hr(s); s2 < hr(s+1); s2++) {
        check (hc(s2)) = 0;
      }
    } else {
      // NOTE: nonzero supernodes in s-th col are not stored
      Kokkos::deep_copy (check, 0);
    }*/
    for (size_type i = 0; i < nidxs; i++) {
      check(idxs(i)) = 0;
    }
  }
  // fix hr
  if (col_major) {
    for (int s = nsuper; s > 0; s--) {
      hr(s) = hr(s - 1);
    }
    hr(0) = 0;
  }
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time_seconds = timer.seconds();
  std::cout << "   > Generate Supernodal Graph: compress graph : " << time_seconds << " (col_major = " << col_major
            << ")" << std::endl;
  timer.reset();
#endif

  // sort column ids per row
  KokkosSparse::sort_crs_graph<Kokkos::HostSpace::execution_space, row_map_view_host_t, cols_view_host_t>(hr, hc);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time_seconds = timer.seconds();
  std::cout << "   > Generate Supernodal Graph: sort graph     : " << time_seconds << std::endl << std::endl;
#endif

  host_graph_t static_graph(hc, hr);
  return static_graph;
}

template <typename graph_t>
graph_t generate_supernodal_dag(int nsuper, graph_t &supL, graph_t &supU) {
  // graph_t is assumed to be on the host
  auto row_mapL = supL.row_map;
  auto entriesL = supL.entries;
  auto row_mapU = supU.row_map;
  auto entriesU = supU.entries;

  // compute dag
  int totedges         = 0;
  using size_type      = typename graph_t::size_type;
  using row_map_view_t = typename graph_t::row_map_type::non_const_type;
  using cols_view_t    = typename graph_t::entries_type::non_const_type;

  row_map_view_t colptr("rowind", nsuper + 1);
  cols_view_t rowind("colptr", row_mapL(nsuper));  // over-estimate
  cols_view_t edges("edges", nsuper);              // workspace
  cols_view_t check("edges", nsuper);              // workspace
  Kokkos::deep_copy(check, 0);
  colptr(0) = totedges;
  for (int s = 0; s < nsuper; s++) {
    // count # of edges (search for first matching nonzero)
    size_type nedges = 0;
    size_type k1     = 1 + row_mapL(s);  // skip diagonal
    size_type k2     = 1 + row_mapU(s);  // skip diagonal
    for (; k1 < row_mapL(s + 1); k1++) {
      // look for match
      while (k2 + 1 < row_mapU(s + 1) && entriesL(k1) > entriesU(k2)) {
        k2++;
      }
      if (k2 + 1 >= row_mapU(s + 1) || entriesL(k1) <= entriesU(k2)) {
        edges(nedges) = entriesL(k1);
        nedges++;
        if (entriesL(k1) == entriesU(k2)) {
#if 0
           // make sure entriesU(k2) rows include the sparsity structure of s?
           for (int k = row_mapL (entriesU (k2)); k < row_mapL (entriesU (k2) + 1); k++) {
             check (entriesL (k)) = 1;
           }
           for (int k = k1 + 1; k < row_mapL (s+1); k++) {
             if (check (entriesL (k)) != 1) {
               edges (nedges) = entriesL (k1);
               nedges ++;
             }
           }
           for (int k = row_mapL (entriesU (k2)); k < row_mapL (entriesU (k2) + 1); k++) {
             check (entriesL (k)) = 0;
           }
#endif
          break;
        }
      }
    }
    // store the edges
    for (size_type k = 0; k < nedges; k++) {
      rowind(totedges) = edges(k);
      totedges++;
    }
    colptr(s + 1) = totedges;
  }

  // compress dag into graph
  cols_view_t hc("colmap_view", totedges);
  for (size_type k = colptr(0); k < colptr(nsuper); k++) {
    hc(k) = rowind(k);
  }

  graph_t static_graph(hc, colptr);
  return static_graph;
}

/* =========================================================================================
 */
template <typename input_graph_t, typename input_size_type>
void merge_supernodal_graph(int *p_nsuper, input_size_type *nb, bool col_majorL, input_graph_t &graphL, bool col_majorU,
                            input_graph_t &graphU, int *etree) {
  int nsuper = *p_nsuper;

  // ---------------------------------------------------------------
  // looking for supernodes to merge (i.e., dense diagonal blocks)
  int nsuper2 = 0;
  auto supL   = generate_supernodal_graph<input_graph_t, input_graph_t>(!col_majorL, graphL, nsuper, nb);
  auto supU   = generate_supernodal_graph<input_graph_t, input_graph_t>(col_majorU, graphU, nsuper, nb);

  auto row_mapL = supL.row_map;
  auto entriesL = supL.entries;

  auto row_mapU = supU.row_map;
  auto entriesU = supU.entries;

  // map the first supernode
  using integer_view_host_t = Kokkos::View<int *, Kokkos::HostSpace>;
  integer_view_host_t map("map", nsuper);  // map old to new supernodes
  map(0) = 0;
  for (int s = 0; s < nsuper - 1; s++) {
    int s2      = s;
    bool merged = false;
    do {
      // check if L(s2+1:end, s2) and L(s2+1:end, s2+1)are the same
      bool mergedL = false;
      int k1       = row_mapL[s2 + 1] - row_mapL[s2];
      int k2       = row_mapL[s2 + 2] - row_mapL[s2 + 1];
      if (k1 == k2 + 1) {
        mergedL = true;
        for (int k = 0; k < k2 && mergedL; k++) {
          if (entriesL[row_mapL[s2] + k + 1] != entriesL[row_mapL[s2 + 1] + k]) {
            mergedL = false;
          }
        }
      }
      // check if U(s2+1:end, s2) and U(s2+1:end, s2+1) are the same
      bool mergedU = false;
      k1           = row_mapU[s2 + 1] - row_mapU[s2];
      k2           = row_mapU[s2 + 2] - row_mapU[s2 + 1];
      if (k1 == k2 + 1) {
        mergedU = true;
        for (int k = 0; k < k2 && mergedU; k++) {
          if (entriesU[row_mapU[s2] + k + 1] != entriesU[row_mapU[s2 + 1] + k]) {
            mergedU = false;
          }
        }
      }
      merged = (mergedL && mergedU);
      if (merged) {
        // printf( "  >> merge s2+1=%d(%dx%d, row=%d:%d) with s=%d(%dx%d) <\n",
        //              s2+1,nb[s2+2]-nb[s2+1],nb[s2+2]-nb[s2+1],
        //              nb[s2+1],nb[s2+2]-1, s,nb[s+1]-nb[s],nb[s+1]-nb[s]);
        map(s2 + 1) = nsuper2;
        s2++;
      } else {
        // printf( "  -- not merge s2+1=%d(%dx%d, row=%d:%d) with s=%d(%dx%d)
        // --\n",
        //           s2+1,nb[s2+2]-nb[s2+1],nb[s2+2]-nb[s2+1],nb[s2+1],nb[s2+2]-1,
        //           s,nb[s+1]-nb[s],nb[s+1]-nb[s]);
        map(s2 + 1) = nsuper2 + 1;
      }
    } while (merged && s2 < nsuper - 1);
    s = s2;
    nsuper2++;
  }
  nsuper2 = map(nsuper - 1) + 1;
  // printf( " nsuper2 = %d\n",nsuper2 );
  // printf( " map:\n" );
  // for (int s = 0; s < nsuper; s++) printf( "   %d %d\n",s,map (s) );

  // ----------------------------------------------------------
  // make sure each of the merged supernodes has the same parent in the etree
  int nsuper3 = 0;
  integer_view_host_t map2;
  if (etree != nullptr) {
    nsuper3 = 0;
    map2    = integer_view_host_t("map2", nsuper);  // map old to new supernodes
    for (int s2 = 0, s = 0; s2 < nsuper2; s2++) {
      // look for parent of the first supernode
      int s3 = s;
      while (etree[s3] != -1 && map(etree[s3]) == map(s3)) {
        s3++;
      }
      map2(s) = nsuper3;
      int p   = (etree[s3] == -1 ? -1 : map(etree[s3]));

      // go through the rest of the supernode in this merged supernode
      s++;
      while (s < nsuper && map(s) == s2) {
        int q = (etree[s3] == -1 ? -1 : map(etree[s3]));
        while (etree[s3] != -1 && map(etree[s3]) == map(s3)) {
          s3++;
          q = (etree[s3] == -1 ? -1 : map(etree[s3]));
        }

        if (q != p) {
          p = q;
          nsuper3++;
        }
        map2(s) = nsuper3;
        s++;
      }
      nsuper3++;
    }
  } else {
    nsuper3 = nsuper2;
    map2    = map;
  }
  // printf( " nsuper3 = %d\n",nsuper3 );
  // printf( " map:\n" );
  // for (int s = 0; s < nsuper; s++) printf( "   %d %d\n",s,map2 (s) );

  // ----------------------------------------------------------
  // construct new supernodes
  integer_view_host_t nb2("nb2", 1 + nsuper3);
  for (int s2 = 0, s = 0; s2 < nsuper3; s2++) {
    nb2(1 + s2) = 0;
    // merging supernodal rows
    while (s < nsuper && map2(s) == s2) {
      nb2(1 + s2) += (nb[s + 1] - nb[s]);
      s++;
    }
  }

  // copy back the new supernodes "offsets"
  nb2(0) = 0;
  for (int s = 0; s < nsuper3; s++) {
    nb2(s + 1) = nb2(s) + nb2(s + 1);
  }
  // copy nb
  for (int s = 0; s < nsuper3; s++) {
    nb[s + 1] = nb2(s + 1);
  }

  // ----------------------------------------------------------
  // construct new etree
  if (etree != nullptr) {
    integer_view_host_t etree2("etree2", nsuper3);
    for (int s = 0; s < nsuper; s++) {
      // etree
      int s2 = map2(s);
      int p  = (etree[s] == -1 ? -1 : map2(etree[s]));
      if (p != s2) {
        etree2(s2) = p;
      }
    }
    // copy etree
    for (int s = 0; s < nsuper3; s++) {
      etree[s] = etree2(s);
    }
  }

  *p_nsuper = nsuper3;
}

/* =========================================================================================
 */
template <typename output_graph_t, typename input_graph_t, typename input_size_type>
output_graph_t generate_merged_supernodal_graph(bool lower, int nsuper, const input_size_type *nb, int nsuper2,
                                                input_size_type *nb2, input_graph_t &graph, int *nnz) {
  using cols_view_t    = typename output_graph_t::entries_type::non_const_type;
  using row_map_view_t = typename output_graph_t::row_map_type::non_const_type;
  using size_type      = typename input_graph_t::size_type;

  // ----------------------------------------------------------
  // now let me find nsrow for the merged supernode
  auto row_map = graph.row_map;
  auto entries = graph.entries;
  int n        = graph.numRows();

  using integer_view_host_t = Kokkos::View<int *, Kokkos::HostSpace>;
  integer_view_host_t mb2("mb2", nsuper2);
  integer_view_host_t work1("work1", n);
  integer_view_host_t work2("work2", n + 1);
  Kokkos::deep_copy(work1, 0);

  int nnzS = 0;
  int nnzA = 0;
  integer_view_host_t rowind("rowind", row_map(n));  // over-estimate
  integer_view_host_t colptr("colptr", nsuper2 + 1);
  colptr(0) = nnzS;
  for (int s2 = 0, s = 0; s2 < nsuper2; s2++) {
    mb2(s2) = 0;
    // merging supernodal rows
    // NOTE: SuperLU may not fill zeros to fill the supernodes
    //       So, these rows may be just subset of the supernodal rows
    while (s < nsuper && nb[s + 1] <= nb2[s2 + 1]) {
      input_size_type j1 = nb[s];
      for (size_type k = row_map(j1); k < row_map(j1 + 1); k++) {
        // just taking union of rows
        if (work1(entries[k]) == 0) {
          work1(entries[k]) = 1;
          work2(mb2(s2))    = entries[k];
          mb2(s2)++;
        }
      }
      s++;
    }
    // sort such that diagonal come on the top
    std::sort(work2.data(), work2.data() + mb2(s2));

    // save nonzero row ids
    if (lower) {
      // lower in csc, diagonal come on top
      for (int k = 0; k < mb2(s2); k++) {
        rowind(nnzS)    = work2(k);
        work1(work2(k)) = 0;
        nnzS++;
      }
    } else {
      // upper in csc, diagonal block is on bottom right now, but move it to top
      int nd = nb2[s2 + 1] - nb2[s2];  // size of diagonal block
      // > diagonal block
      for (int k = mb2(s2) - nd; k < mb2(s2); k++) {
        rowind(nnzS)    = work2(k);
        work1(work2(k)) = 0;
        nnzS++;
      }
      // > offdiagonal blocks
      for (int k = 0; k < mb2(s2) - nd; k++) {
        rowind(nnzS)    = work2(k);
        work1(work2(k)) = 0;
        nnzS++;
      }
    }
    colptr(s2 + 1) = nnzS;
    nnzA += (nb2[s2 + 1] - nb2[s2]) * mb2(s2);
  }

  // ----------------------------------------------------------
  // now let's create crs graph
  row_map_view_t rowmap_view("rowmap_view", n + 1);
  cols_view_t column_view("colmap_view", nnzA);
  auto hr = Kokkos::create_mirror_view(rowmap_view);
  auto hc = Kokkos::create_mirror_view(column_view);

  nnzA  = 0;
  hr(0) = 0;
  for (int s2 = 0; s2 < nsuper2; s2++) {
    for (int j = nb2[s2]; j < nb2[s2 + 1]; j++) {
      for (int k = colptr(s2); k < colptr(s2 + 1); k++) {
        hc(nnzA) = rowind(k);
        nnzA++;
      }
      hr(j + 1) = nnzA;
    }
  }
  *nnz = nnzA;

  // deepcopy
  Kokkos::deep_copy(rowmap_view, hr);
  Kokkos::deep_copy(column_view, hc);

  // create crs
  output_graph_t static_graph(column_view, rowmap_view);
  return static_graph;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/* For symbolic analysis */
template <typename host_graph_t, typename KernelHandle>
void sptrsv_supernodal_symbolic(int nsuper, int *supercols, int *etree, host_graph_t graphL_host,
                                KernelHandle *kernelHandleL, host_graph_t graphU_host, KernelHandle *kernelHandleU) {
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  int nrows           = graphL_host.numRows();
  double time_seconds = 0.0;
  Kokkos::Timer timer;
  Kokkos::Timer tic;
  timer.reset();
#endif

  // ===================================================================
  // load sptrsv-handles
  auto *handleL = kernelHandleL->get_sptrsv_handle();
  auto *handleU = kernelHandleU->get_sptrsv_handle();

  // store arguments to handle
  handleL->set_graph_host(graphL_host);
  handleU->set_graph_host(graphU_host);
  handleL->set_supernodes(nsuper, supercols, etree);
  handleU->set_supernodes(nsuper, supercols, etree);

  // ===================================================================
  // load options
  bool col_majorL = handleL->is_column_major();
  bool col_majorU = handleU->is_column_major();
  bool merge      = handleL->get_merge_supernodes();
  bool UinCSC     = handleU->is_column_major();
  bool needEtree  = (handleL->get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
                    handleL->get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_ETREE);
  if (needEtree && etree == nullptr) {
    std::cout << std::endl
              << " ** etree needs to be set before calling sptrsv_symbolic "
                 "with SuperLU **"
              << std::endl
              << std::endl;
    return;
  }

  // ===================================================================
  // > make a copy of supercols (merge needs both original and merged supercols)
  using integer_view_host_t = typename KernelHandle::SPTRSVHandleType::integer_view_host_t;
  integer_view_host_t supercols_view("supercols view", 1 + nsuper);
  int *supercols_merged = supercols_view.data();
  for (int i = 0; i <= nsuper; i++) {
    supercols_merged[i] = supercols[i];
  }
  if (merge) {
    // =================================================================
    // merge supernodes
    int nsuper_merged = nsuper;
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    tic.reset();
    check_supernode_sizes("Original L-structure", nrows, nsuper, supercols_merged, graphL_host);
    check_supernode_sizes("Original U-structure", nrows, nsuper, supercols_merged, graphU_host);
#endif
    // etree will be updated
    merge_supernodal_graph(&nsuper_merged, supercols_merged, col_majorL, graphL_host, col_majorU, graphU_host, etree);

// =================================================================
// generate merged graph for L-solve
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    int nnzL = graphL_host.row_map(nrows);
#endif
    int nnzL_merged;
    bool lower = true;
    handleL->set_original_graph_host(graphL_host);  // save graph before merge
    graphL_host = generate_merged_supernodal_graph<host_graph_t>(lower, nsuper, supercols, nsuper_merged,
                                                                 supercols_merged, graphL_host, &nnzL_merged);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    time_seconds = tic.seconds();
    check_supernode_sizes("After Merge", nrows, nsuper_merged, supercols_merged, graphL_host);
    std::cout << " for L factor:" << std::endl;
    std::cout << "   Merge Supernodes Time: " << time_seconds << std::endl;
    std::cout << "   Number of nonzeros   : " << nnzL << " -> " << nnzL_merged << " : "
              << double(nnzL_merged) / double(nnzL) << "x" << std::endl;
#endif

// =================================================================
// generate merged graph for U-solve
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    tic.reset();
    int nnzU = graphU_host.row_map(nrows);
#endif
    int nnzU_merged;
    lower = (UinCSC ? false : true);
    handleU->set_original_graph_host(graphU_host);  // save graph before merge
    graphU_host = generate_merged_supernodal_graph<host_graph_t>(lower, nsuper, supercols, nsuper_merged,
                                                                 supercols_merged, graphU_host, &nnzU_merged);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    time_seconds = tic.seconds();
    check_supernode_sizes("After Merge", nrows, nsuper_merged, supercols_merged, graphU_host);
    std::cout << " for U factor:" << std::endl;
    std::cout << "   Merge Supernodes Time: " << time_seconds << std::endl;
    std::cout << "   Number of nonzeros   : " << nnzU << " -> " << nnzU_merged << " : "
              << double(nnzU_merged) / double(nnzU) << "x" << std::endl;
#endif

    // update the number of supernodes
    nsuper = nsuper_merged;
  }
  // replace the supernodal info with the merged ones
  supercols = supercols_merged;

  // ===================================================================
  // copy graph to device
  using graph_t = typename KernelHandle::SPTRSVHandleType::graph_t;
  auto graphL   = deep_copy_graph<host_graph_t, graph_t>(graphL_host);
  auto graphU   = deep_copy_graph<host_graph_t, graph_t>(graphU_host);

  // ===================================================================
  // save the supernodal info in the handles for L/U solves
  handleL->set_supernodes(nsuper, supercols_view, etree);
  handleU->set_supernodes(nsuper, supercols_view, etree);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time_seconds = tic.seconds();
  std::cout << "   Deep-copy graph Time: " << time_seconds << std::endl;
  tic.reset();
#endif

  if (handleL->get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_DAG ||
      handleL->get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG) {
    // generate supernodal graphs for DAG scheduling
    auto supL = generate_supernodal_graph<host_graph_t>(!col_majorL, graphL_host, nsuper, supercols);
    auto supU = generate_supernodal_graph<host_graph_t>(col_majorU, graphU_host, nsuper, supercols);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    time_seconds = tic.seconds();
    std::cout << "   Compute Supernodal Graph Time: " << time_seconds << std::endl;
    tic.reset();
#endif

    auto dagL = generate_supernodal_dag<host_graph_t>(nsuper, supL, supU);
    auto dagU = generate_supernodal_dag<host_graph_t>(nsuper, supU, supL);
    handleL->set_supernodal_dag(dagL);
    handleU->set_supernodal_dag(dagU);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    time_seconds = tic.seconds();
    std::cout << "   Compute DAG Time: " << time_seconds << std::endl;
    tic.reset();
#endif
  }

  // ===================================================================
  // do symbolic for L solve on the host
  auto row_mapL = graphL.row_map;
  auto entriesL = graphL.entries;
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  tic.reset();
  std::cout << std::endl;
#endif
  sptrsv_symbolic(kernelHandleL, row_mapL, entriesL);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time_seconds = tic.seconds();
  std::cout << " > Lower-TRI: " << std::endl;
  std::cout << "   Symbolic Time: " << time_seconds << std::endl;
  tic.reset();
#endif

  // ===================================================================
  // do symbolic for U solve on the host
  auto row_mapU = graphU.row_map;
  auto entriesU = graphU.entries;
  sptrsv_symbolic(kernelHandleU, row_mapU, entriesU);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time_seconds = tic.seconds();
  std::cout << " > Upper-TRI: " << std::endl;
  std::cout << "   Symbolic Time: " << time_seconds << std::endl;
#endif

  // ===================================================================
  // save options
  handleL->set_merge_supernodes(merge);
  handleU->set_merge_supernodes(merge);

  // ===================================================================
  // save graphs
  handleL->set_graph(graphL);
  handleU->set_graph(graphU);
  // graph on host (merged)
  handleL->set_graph_host(graphL_host);
  handleU->set_graph_host(graphU_host);

  // ===================================================================
  handleL->set_symbolic_complete();
  handleU->set_symbolic_complete();
  handleL->set_etree(etree);
  handleU->set_etree(etree);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time_seconds = timer.seconds();
  std::cout << "   Total Symbolic Time: " << time_seconds << std::endl << std::endl;
  std::cout << "   Total nnz: " << graphL_host.row_map(nrows) << " + " << graphU_host.row_map(nrows) << std::endl;
#endif
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/* Auxiliary functions for numeric computation */

/* =========================================================================================
 */
struct Tag_SupTrtriFunctor {};
struct Tag_SupTrtriTrmmFunctor {};

template <typename UploType, typename DiagType, typename integer_view_host_t, typename input_size_type,
          typename row_map_type, typename index_type, typename values_type>
struct TriSupernodalTrtriFunctor {
  integer_view_host_t supernode_ids;
  const input_size_type *nb;
  row_map_type hr;
  index_type hc;
  values_type hv;

  KOKKOS_INLINE_FUNCTION
  TriSupernodalTrtriFunctor(integer_view_host_t supernode_ids_, const input_size_type *nb_, row_map_type &hr_,
                            index_type &hc_, values_type &hv_)
      : supernode_ids(supernode_ids_), nb(nb_), hr(hr_), hc(hc_), hv(hv_) {}

  // functor: just invert diagonal
  KOKKOS_INLINE_FUNCTION
  void operator()(const Tag_SupTrtriFunctor &, const int i) const {
    using execution_space = typename values_type::execution_space;
    using memory_space    = typename execution_space::memory_space;
    using values_view_t   = typename values_type::non_const_type;
    using scalar_t        = typename values_view_t::value_type;

    using range_type    = Kokkos::pair<int, int>;
    using TrtriAlgoType = KokkosBatched::Algo::Trtri::Unblocked;

    int s     = supernode_ids(i);
    int j1    = nb[s];
    int nsrow = hr(j1 + 1) - hr(j1);
    int nscol = nb[s + 1] - nb[s];

    // invert diagonal
    auto nnzD = hr(j1);
    Kokkos::View<scalar_t **, Kokkos::LayoutLeft, memory_space, Kokkos::MemoryUnmanaged> viewL(&hv(nnzD), nsrow, nscol);
    auto Ljj = Kokkos::subview(viewL, range_type(0, nscol), Kokkos::ALL());
    KokkosBatched::SerialTrtri<UploType, DiagType, TrtriAlgoType>::invoke(Ljj);
  }

  // functor: invert diagonal + apply inverse to off-diagonal
  KOKKOS_INLINE_FUNCTION
  void operator()(const Tag_SupTrtriTrmmFunctor &, const int i) const {
    using execution_space = typename values_type::execution_space;
    using memory_space    = typename execution_space::memory_space;
    using values_view_t   = typename values_type::non_const_type;
    using scalar_t        = typename values_view_t::value_type;

    using range_type    = Kokkos::pair<int, int>;
    using TrtriAlgoType = KokkosBatched::Algo::Trtri::Unblocked;
    using Side          = KokkosBatched::Side;
    using Trans         = KokkosBatched::Trans;

    int s     = supernode_ids(i);
    int j1    = nb[s];
    int nsrow = hr(j1 + 1) - hr(j1);
    int nscol = nb[s + 1] - nb[s];

    // invert diagonal
    auto nnzD = hr(j1);
    Kokkos::View<scalar_t **, Kokkos::LayoutLeft, memory_space, Kokkos::MemoryUnmanaged> viewL(&hv(nnzD), nsrow, nscol);
    auto Ljj = Kokkos::subview(viewL, range_type(0, nscol), Kokkos::ALL());
    KokkosBatched::SerialTrtri<UploType, DiagType, TrtriAlgoType>::invoke(Ljj);

    // apply invse to off-diagonal
    // if (nsrow > nscol && invert_offdiag)
    {
      const scalar_t one(1.0);
      auto Lij = Kokkos::subview(viewL, range_type(nscol, nsrow), Kokkos::ALL());
      KokkosBatched::SerialTrmm<Side::Right, UploType, Trans::NoTranspose, DiagType, TrtriAlgoType>::invoke(one, Ljj,
                                                                                                            Lij);
    }
  }
};

/* =========================================================================================
 */
template <typename KernelHandle, typename input_size_type, typename row_map_type, typename index_type,
          typename values_type, typename integer_view_host_t>
void invert_supernodal_columns_batched(KernelHandle *kernelHandle, bool unit_diag, const input_size_type *nb,
                                       row_map_type &hr, index_type &hc, values_type &hv, int num_batches,
                                       integer_view_host_t supernode_ids) {
  using execution_space = typename values_type::execution_space;

  using Uplo = KokkosBatched::Uplo;
  using Diag = KokkosBatched::Diag;

  // load parameters
  auto *handle        = kernelHandle->get_sptrsv_handle();
  bool invert_diag    = handle->get_invert_diagonal();
  bool invert_offdiag = handle->get_invert_offdiagonal();

  // quick return
  if (!invert_diag && !invert_offdiag) return;

  if (num_batches > 0) {
    // lower is always in CSC, if UinCSC, then lower=false, else lower=true
    bool lower_tri = kernelHandle->is_sptrsv_lower_tri();
    bool lower     = ((lower_tri && handle->is_column_major()) || (!lower_tri && !handle->is_column_major()));

    if (lower) {
      if (unit_diag) {
        if (invert_offdiag) {
          using range_policy = Kokkos::RangePolicy<Tag_SupTrtriTrmmFunctor, execution_space>;
          TriSupernodalTrtriFunctor<Uplo::Lower, Diag::Unit, integer_view_host_t, input_size_type, row_map_type,
                                    index_type, values_type>
              sptrsv_tritri_functor(supernode_ids, nb, hr, hc, hv);
          Kokkos::parallel_for("TriSupernodalTrtriFunctor", range_policy(0, num_batches), sptrsv_tritri_functor);
        } else {
          using range_policy = Kokkos::RangePolicy<Tag_SupTrtriFunctor, execution_space>;
          TriSupernodalTrtriFunctor<Uplo::Lower, Diag::Unit, integer_view_host_t, input_size_type, row_map_type,
                                    index_type, values_type>
              sptrsv_tritri_functor(supernode_ids, nb, hr, hc, hv);
          Kokkos::parallel_for("TriSupernodalTrtriFunctor", range_policy(0, num_batches), sptrsv_tritri_functor);
        }
      } else {
        if (invert_offdiag) {
          using range_policy = Kokkos::RangePolicy<Tag_SupTrtriTrmmFunctor, execution_space>;
          TriSupernodalTrtriFunctor<Uplo::Lower, Diag::NonUnit, integer_view_host_t, input_size_type, row_map_type,
                                    index_type, values_type>
              sptrsv_tritri_functor(supernode_ids, nb, hr, hc, hv);
          Kokkos::parallel_for("TriSupernodalTrtriFunctor", range_policy(0, num_batches), sptrsv_tritri_functor);
        } else {
          using range_policy = Kokkos::RangePolicy<Tag_SupTrtriFunctor, execution_space>;
          TriSupernodalTrtriFunctor<Uplo::Lower, Diag::NonUnit, integer_view_host_t, input_size_type, row_map_type,
                                    index_type, values_type>
              sptrsv_tritri_functor(supernode_ids, nb, hr, hc, hv);
          Kokkos::parallel_for("TriSupernodalTrtriFunctor", range_policy(0, num_batches), sptrsv_tritri_functor);
        }
      }
    } else {
      if (unit_diag) {
        if (invert_offdiag) {
          using range_policy = Kokkos::RangePolicy<Tag_SupTrtriTrmmFunctor, execution_space>;
          TriSupernodalTrtriFunctor<Uplo::Upper, Diag::Unit, integer_view_host_t, input_size_type, row_map_type,
                                    index_type, values_type>
              sptrsv_tritri_functor(supernode_ids, nb, hr, hc, hv);
          Kokkos::parallel_for("TriSupernodalTrtriFunctor", range_policy(0, num_batches), sptrsv_tritri_functor);
        } else {
          using range_policy = Kokkos::RangePolicy<Tag_SupTrtriFunctor, execution_space>;
          TriSupernodalTrtriFunctor<Uplo::Upper, Diag::Unit, integer_view_host_t, input_size_type, row_map_type,
                                    index_type, values_type>
              sptrsv_tritri_functor(supernode_ids, nb, hr, hc, hv);
          Kokkos::parallel_for("TriSupernodalTrtriFunctor", range_policy(0, num_batches), sptrsv_tritri_functor);
        }
      } else {
        if (invert_offdiag) {
          using range_policy = Kokkos::RangePolicy<Tag_SupTrtriTrmmFunctor, execution_space>;
          TriSupernodalTrtriFunctor<Uplo::Upper, Diag::NonUnit, integer_view_host_t, input_size_type, row_map_type,
                                    index_type, values_type>
              sptrsv_tritri_functor(supernode_ids, nb, hr, hc, hv);
          Kokkos::parallel_for("TriSupernodalTrtriFunctor", range_policy(0, num_batches), sptrsv_tritri_functor);
        } else {
          using range_policy = Kokkos::RangePolicy<Tag_SupTrtriFunctor, execution_space>;
          TriSupernodalTrtriFunctor<Uplo::Upper, Diag::NonUnit, integer_view_host_t, input_size_type, row_map_type,
                                    index_type, values_type>
              sptrsv_tritri_functor(supernode_ids, nb, hr, hc, hv);
          Kokkos::parallel_for("TriSupernodalTrtriFunctor", range_policy(0, num_batches), sptrsv_tritri_functor);
        }
      }
    }
  }
}

/* =========================================================================================
 */
template <typename KernelHandle, typename input_size_type, typename row_map_type, typename index_type,
          typename values_type>
void invert_supernodal_columns(KernelHandle *kernelHandle, bool unit_diag, int nsuper, const input_size_type *nb,
                               row_map_type &hr, index_type &hc, values_type &hv) {
  using execution_space     = typename values_type::execution_space;
  using memory_space        = typename execution_space::memory_space;
  using values_view_t       = typename values_type::non_const_type;
  using scalar_t            = typename values_view_t::value_type;
  using range_type          = Kokkos::pair<int, int>;
  using integer_view_host_t = Kokkos::View<int *, Kokkos::HostSpace>;

  const scalar_t one(1.0);

  // load parameters
  auto *handle        = kernelHandle->get_sptrsv_handle();
  bool invert_diag    = handle->get_invert_diagonal();
  bool invert_offdiag = handle->get_invert_offdiagonal();

  // lower is always in CSC, if UinCSC, then lower=false, else lower=true
  bool lower_tri = kernelHandle->is_sptrsv_lower_tri();
  bool lower     = ((lower_tri && handle->is_column_major()) || (!lower_tri && !handle->is_column_major()));

  // quick return
  if (!invert_diag) return;

#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  Kokkos::Timer timer;
  double time1 = 0.0;
  double time2 = 0.0;
  double time3 = 0.0;
#endif

  // ----------------------------------------------------------
  // now let's invert some blocks
  // > first go through all the supernode columns
  // > use KokkosBlas on large blocks, and keep track of small blocks
  // > to call batchedBlas on them
  int num_batches    = 0;
  int size_unblocked = handle->get_supernode_size_unblocked();
  integer_view_host_t supernode_ids("supernode_batch", nsuper);

  // ----------------------------------------------------------
  // If we are running KokkosKernels::trmm on device,
  // then we need to allocate a workspace on device
  using trmm_execution_space = typename KernelHandle::HandleExecSpace;
  using trmm_memory_space    = typename KernelHandle::HandlePersistentMemorySpace;
  using trmm_view_t          = Kokkos::View<scalar_t *, trmm_execution_space>;
#if !defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
  // use KokkosBlas::trmm only with CUBLAS (since deep-copy to host throws an
  // error)
  bool run_trmm_on_device = false;
#else
  bool run_trmm_on_device =
      (handle->get_trmm_on_device() && !std::is_same<trmm_execution_space, execution_space>::value);
#endif

  // figure out largest supernode
  int lwork = 0;
  trmm_view_t trmm_dwork("trmm_dwork", lwork);
  if (run_trmm_on_device) {
    for (int s2 = 0; s2 < nsuper; s2++) {
      int nscol = nb[s2 + 1] - nb[s2];
      if (nscol >= size_unblocked) {
        int j1    = nb[s2];
        int nsrow = hr(j1 + 1) - hr(j1);
        if (lwork < nsrow * nscol) {
          lwork = nsrow * nscol;
        }
      }
    }
    try {
      Kokkos::resize(trmm_dwork, lwork);
    } catch (...) {
      // something went wrong allocating device memory
      // so we'll just do trmm on host
      run_trmm_on_device = false;
    }
  }

  // ----------------------------------------------------------
  // now go through the supernode columns and invert "large" supernodes
  // using KokkosKernels::trtri (host) and KokkosKernels::trmm (host or device)
  for (int s2 = 0; s2 < nsuper; s2++) {
    int nscol = nb[s2 + 1] - nb[s2];

    if (nscol >= size_unblocked) {
      int j1    = nb[s2];
      int nsrow = hr(j1 + 1) - hr(j1);

      auto nnzD      = hr(j1);
      char uplo_char = (lower ? 'L' : 'U');
      char diag_char = (unit_diag ? 'U' : 'N');

      // NOTE: we currently supports only default_layout = LayoutLeft
      Kokkos::View<scalar_t **, default_layout, memory_space, Kokkos::MemoryUnmanaged> viewL(&hv(nnzD), nsrow, nscol);
      auto Ljj = Kokkos::subview(viewL, range_type(0, nscol), Kokkos::ALL());

#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
      timer.reset();
#endif
#if defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
      if (run_trmm_on_device) {
        Kokkos::View<scalar_t **, Kokkos::LayoutLeft, trmm_memory_space, Kokkos::MemoryUnmanaged> dViewL(
            trmm_dwork.data(), nsrow, nscol);

        // deep-copy the whole supernode column to device
        Kokkos::deep_copy(dViewL, viewL);

        // call trtri on device
        auto dViewLjj = Kokkos::subview(dViewL, range_type(0, nscol), Kokkos::ALL());
        KokkosLapack::trtri(&uplo_char, &diag_char, dViewLjj);
      } else
#endif
      {
        // call trtri on host
        KokkosLapack::trtri(&uplo_char, &diag_char, Ljj);
      }
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
      time1 += timer.seconds();
#endif

      if (nsrow > nscol && invert_offdiag) {
        char side_char = 'R';
        char tran_char = 'N';
        auto Lij       = Kokkos::subview(viewL, range_type(nscol, nsrow), Kokkos::ALL());

#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
        timer.reset();
#endif
        if (run_trmm_on_device) {
          Kokkos::View<scalar_t **, Kokkos::LayoutLeft, trmm_memory_space, Kokkos::MemoryUnmanaged> dViewL(
              trmm_dwork.data(), nsrow, nscol);

#if !defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
          // deep-copy the whole supernode column to device
          Kokkos::deep_copy(dViewL, viewL);
#endif

          // NOTE: we currently supports only default_layout = LayoutLeft
          auto dViewLjj = Kokkos::subview(dViewL, range_type(0, nscol), Kokkos::ALL());
          auto dViewLij = Kokkos::subview(dViewL, range_type(nscol, nsrow), Kokkos::ALL());

          KokkosBlas::trmm(&side_char, &uplo_char, &tran_char, &diag_char, one, dViewLjj, dViewLij);

#if !defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
          // deep-copy the whole panel back to host (since I cannot just
          // deep-copy Lij)
          Kokkos::deep_copy(viewL, dViewL);
#endif
        } else {
          KokkosBlas::trmm(&side_char, &uplo_char, &tran_char, &diag_char, one, Ljj, Lij);
        }
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
        time2 += timer.seconds();
#endif
      }

#if defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
      if (run_trmm_on_device) {
        // deep-copy the whole supernode column back to host
        Kokkos::View<scalar_t **, Kokkos::LayoutLeft, trmm_memory_space, Kokkos::MemoryUnmanaged> dViewL(
            trmm_dwork.data(), nsrow, nscol);
        Kokkos::deep_copy(viewL, dViewL);
      }
#endif
    } else {
      supernode_ids(num_batches) = s2;
      num_batches++;
    }
  }

// ----------------------------------------------------------
// now call batchedBLAS on "small" supernodes
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  timer.reset();
#endif
  invert_supernodal_columns_batched(kernelHandle, unit_diag, nb, hr, hc, hv, num_batches, supernode_ids);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time3 = timer.seconds();
#endif

  if (run_trmm_on_device) {
    // to make sure the data is deep-copied to host..
    Kokkos::fence();
  }
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  std::cout << "   invert_supernodes" << std::endl;
  std::cout << "   + num supernodes = " << nsuper << " num batchs = " << num_batches << std::endl;
  std::cout << "   > Time for inversion::trtri : " << time1 << std::endl;
  std::cout << "   > Time for inversion::trmm  : " << time2 << std::endl;
  std::cout << "   > Time for batchs           : " << time3 << std::endl;
#endif
}

/* =========================================================================================
 */
template <typename crsmat_t, typename input_crsmat_t, typename input_ptr_type, typename graph_t, typename KernelHandle>
crsmat_t read_merged_supernodes(KernelHandle *kernelHandle, int nsuper, const input_ptr_type *mb, bool unit_diag,
                                input_crsmat_t &L, graph_t &static_graph) {
  using values_view_t      = typename crsmat_t::values_type::non_const_type;
  using scalar_t           = typename values_view_t::value_type;
  using scalar_view_host_t = Kokkos::View<scalar_t *, Kokkos::HostSpace>;

  const scalar_t zero(0.0);

  // original matrix
  auto graphL   = L.graph;  // in_graph
  auto row_mapL = graphL.row_map;
  auto entriesL = graphL.entries;
  auto valuesL  = L.values;
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  Kokkos::Timer timer;
  Kokkos::Timer timer2;
  timer.reset();
  timer2.reset();
#endif

  // merged graph
  auto rowmap_view = static_graph.row_map;
  auto column_view = static_graph.entries;

  auto hr = Kokkos::create_mirror_view(rowmap_view);
  auto hc = Kokkos::create_mirror_view(column_view);
  Kokkos::deep_copy(hr, rowmap_view);
  Kokkos::deep_copy(hc, column_view);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  double time_copy = timer2.seconds();
  timer2.reset();
#endif

  // ----------------------------------------------------------
  // now let's merge supernodes
  int n = graphL.numRows();
  scalar_view_host_t dwork("dwork", n);
  Kokkos::deep_copy(dwork, zero);

  auto nnzA = hr(n);
  values_view_t values_view("values_view", nnzA);
  auto hv = Kokkos::create_mirror_view(values_view);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  double time_mirror = timer2.seconds();
  timer2.reset();
#endif

  for (int s2 = 0; s2 < nsuper; s2++) {
    for (int j = mb[s2]; j < mb[s2 + 1]; j++) {
      for (int k = row_mapL[j]; k < row_mapL[j + 1]; k++) {
        dwork(entriesL[k]) = valuesL[k];
      }
      for (int k = hr(j); k < hr(j + 1); k++) {
        hv(k) = dwork(hc(k));
      }
      for (int k = row_mapL[j]; k < row_mapL[j + 1]; k++) {
        dwork(entriesL[k]) = zero;
      }
    }
  }
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  double time_merge = timer2.seconds();
  timer2.reset();
#endif

  // invert blocks (TODO done on host for now)
  invert_supernodal_columns(kernelHandle, unit_diag, nsuper, mb, hr, hc, hv);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  double time_invert = timer2.seconds();
  timer2.reset();
#endif
  // deepcopy
  Kokkos::deep_copy(values_view, hv);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time_copy += timer2.seconds();
#endif

#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  double time = timer.seconds();
  std::cout << "   read_merged_supernodes" << std::endl;
  std::cout << "   > Time       : " << time << std::endl;
  std::cout << "    + copy   time : " << time_copy << std::endl;
  std::cout << "    + mirror time : " << time_mirror << std::endl;
  std::cout << "    + merge  time : " << time_merge << std::endl;
  std::cout << "    + invert time : " << time_invert << std::endl;
#endif

  // create crs
  crsmat_t crsmat("CrsMatrix", n, values_view, static_graph);

  return crsmat;
}

/* =========================================================================================
 */
template <typename crsmat_t, typename graph_t, typename scalar_t, typename input_size_type, typename input_ptr_type,
          typename size_type, typename ordinal_type, typename KernelHandle>
crsmat_t read_supernodal_values(KernelHandle *kernelHandle, int n, int nsuper, bool ptr_by_column,
                                const input_size_type *mb, const input_ptr_type *nb, const size_type *colptr,
                                ordinal_type *rowind, scalar_t *Lx, graph_t &static_graph) {
  using values_view_t       = typename crsmat_t::values_type::non_const_type;
  using integer_view_host_t = Kokkos::View<ordinal_type *, Kokkos::HostSpace>;

  const scalar_t zero(0.0);
  const scalar_t one(1.0);

  Kokkos::Timer timer;
  timer.reset();

  // load parameters
  auto *handle   = kernelHandle->get_sptrsv_handle();
  bool unit_diag = handle->is_unit_diagonal();
  bool merge     = handle->get_merge_supernodes();

  // lower is always in CSC, if UinCSC, then lower=false, else lower=true
  bool lower_tri = kernelHandle->is_sptrsv_lower_tri();
  bool lower     = ((lower_tri && handle->is_column_major()) || (!lower_tri && !handle->is_column_major()));

  // load graph
  auto rowmap_view = static_graph.row_map;
  auto column_view = static_graph.entries;
  auto hr          = Kokkos::create_mirror_view(rowmap_view);
  auto hc          = Kokkos::create_mirror_view(column_view);
  Kokkos::deep_copy(hr, rowmap_view);
  Kokkos::deep_copy(hc, column_view);

  // total nnz
  int nnzL = hr(n);
  values_view_t values_view("values_view", nnzL);
  auto hv = Kokkos::create_mirror_view(values_view);
  Kokkos::deep_copy(hv,
                    zero);  // seems to be needed (instead of zeroing out upper)

  // compute max nnz per row
  int max_nnz_per_row = 0;
  for (int s = 0; s < nsuper; s++) {
    int i1, i2;
    if (ptr_by_column) {
      int j1 = nb[s];

      i1 = mb[j1];
      i2 = mb[j1 + 1];
    } else {
      i1 = mb[s];
      i2 = mb[s + 1];
    }
    // "total" number of rows in all the supernodes (diagonal+off-diagonal)
    int nsrow = i2 - i1;
    if (nsrow > max_nnz_per_row) {
      max_nnz_per_row = nsrow;
    }
  }

  integer_view_host_t sorted_rowind("sorted_rowind", max_nnz_per_row + 1);
  // store L in csr
  for (int s = 0; s < nsuper; s++) {
    int j1    = nb[s];
    int j2    = nb[s + 1];
    int nscol = j2 - j1;  // number of columns in the s-th supernode column

    int i1, i2;
    if (ptr_by_column) {
      i1 = mb[j1];
      i2 = mb[j1 + 1];
    } else {
      i1 = mb[s];
      i2 = mb[s + 1];
    }
    int nsrow = i2 - i1;         // "total" number of rows in all the supernodes
                                 // (diagonal+off-diagonal)
    int nsrow2 = nsrow - nscol;  // "total" number of rows in all the off-diagonal supernodes
    int ps2    = i1 + nscol;     // offset into rowind

    int psx;  // offset into data,   Lx[s][s]
    if (ptr_by_column) {
      psx = colptr[j1];
    } else {
      psx = colptr[s];
    }

    /* diagonal block */
    // for each column (or row due to symmetry), the diagonal supernodal block
    // is stored (in ascending order of row indexes) first so that we can do
    // TRSM on the diagonal block
    for (int jj = 0; jj < nscol; jj++) {
      if (lower) {
        // shift for explicitly store zeros in upper-part
        hr(j1 + jj) += jj;
        // diagonal
        if (unit_diag) {
          hv(hr(j1 + jj)) = one;
        } else {
          hv(hr(j1 + jj)) = Lx[psx + (jj + jj * nsrow)];
        }
        hr(j1 + jj)++;
        // lower-triangular part
        for (int ii = jj + 1; ii < nscol; ii++) {
          hv(hr(j1 + jj)) = Lx[psx + (ii + jj * nsrow)];
          hr(j1 + jj)++;
        }
      } else {
        // upper-triangular part
        for (int ii = 0; ii < jj; ii++) {
          hv(hr(j1 + jj)) = Lx[psx + (ii + jj * nsrow)];
          hr(j1 + jj)++;
        }
        // diagonal
        if (unit_diag) {
          hv(hr(j1 + jj)) = one;
        } else {
          hv(hr(j1 + jj)) = Lx[psx + (jj + jj * nsrow)];
        }
        hr(j1 + jj)++;
        // shift for explicitly store zeros in lower-part
        hr(j1 + jj) += (nscol - jj - 1);
      }
    }
    /* off-diagonal blocks */
    if (merge) {
      // sort rowind (to merge supernodes)
      for (int ii = 0; ii < nsrow2; ii++) {
        sorted_rowind(ii) = ii;
      }
      std::sort(sorted_rowind.data(), sorted_rowind.data() + nsrow2, sort_indices<ordinal_type>(&rowind[ps2]));
    }
    for (int jj = 0; jj < nscol; jj++) {
      for (int kk = 0; kk < nsrow2; kk++) {
        int ii          = (merge ? sorted_rowind(kk) : kk);  // sorted rowind
        hv(hr(j1 + jj)) = Lx[psx + (nscol + ii + jj * nsrow)];
        hr(j1 + jj)++;
      }
    }
  }

  // fix hr
  for (int i = n; i >= 1; i--) {
    hr(i) = hr(i - 1);
  }
  hr(0) = 0;

#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  std::cout << "    read_supernodal_values(" << (lower ? "lower)" : "upper)") << std::endl;
  std::cout << "    * Matrix size = " << n << std::endl;
  std::cout << "    * Total nnz   = " << hr(n) << std::endl;
  std::cout << "    * nnz / n     = " << hr(n) / n << std::endl;

  double time = timer.seconds();
  std::cout << "    > Time : " << time << std::endl;
#endif

  // invert blocks (TODO done on host for now)
  invert_supernodal_columns(kernelHandle, unit_diag, nsuper, nb, hr, hc, hv);
  // deepcopy
  Kokkos::deep_copy(values_view, hv);

  // create crs
  crsmat_t crsmat("CrsMatrix", n, values_view, static_graph);

  return crsmat;
}

/* =========================================================================================
 */
template <typename crsmat_t, typename graph_t, typename scalar_t, typename input_size_type, typename input_ptr_type,
          typename size_type, typename ordinal_type, typename KernelHandle>
crsmat_t read_supernodal_valuesLt(KernelHandle *kernelHandle, int n, int nsuper, bool ptr_by_column,
                                  const input_size_type *mb, const input_ptr_type *nb, const size_type *colptr,
                                  ordinal_type *rowind, scalar_t *Lx, graph_t &static_graph) {
  using values_view_t       = typename crsmat_t::values_type::non_const_type;
  using integer_view_host_t = Kokkos::View<int *, Kokkos::HostSpace>;

  const scalar_t zero(0.0);
  const scalar_t one(1.0);

  Kokkos::Timer timer;
  timer.reset();

  // load parameters
  auto *handle   = kernelHandle->get_sptrsv_handle();
  bool unit_diag = handle->is_unit_diagonal();

  // load graph
  auto rowmap_view = static_graph.row_map;
  auto column_view = static_graph.entries;
  auto hr          = Kokkos::create_mirror_view(rowmap_view);
  auto hc          = Kokkos::create_mirror_view(column_view);
  Kokkos::deep_copy(hr, rowmap_view);
  Kokkos::deep_copy(hc, column_view);

  // total nnz
  int nnzL = hr(n);
  values_view_t values_view("values_view", nnzL);
  auto hv = Kokkos::create_mirror_view(values_view);
  Kokkos::deep_copy(hv,
                    zero);  // seems to be needed (instead of zeroing out upper)

  /* create a map from row id to supernode id */
  integer_view_host_t map("map", n);
  int supid = 0;
  for (int k = 0; k < nsuper; k++) {
    int j1 = nb[k];
    int j2 = nb[k + 1];
    for (int j = j1; j < j2; j++) {
      map(j) = supid;
    }
    supid++;
  }

  // pointer to off-diagonals (diagonal comes first)
  integer_view_host_t off("off", n + 1);
  for (int s = 0; s < nsuper; s++) {
    int i1 = nb[s];
    int i2 = nb[s + 1];
    // number of columns in the s-th supernode column
    int nscol = i2 - i1;
    for (int ii = i1; ii < i2; ii++) {
      off(ii) = hr(ii) + nscol;
    }
  }

  // store L in csr
  integer_view_host_t sup("sup", nsuper);
  integer_view_host_t check("check", nsuper);
  Kokkos::deep_copy(check, 0);
  for (int s = 0; s < nsuper; s++) {
    int j1    = nb[s];
    int j2    = nb[s + 1];
    int nscol = j2 - j1;  // number of columns in the s-th supernode column

    int i1, i2;
    if (ptr_by_column) {
      i1 = mb[j1];
      i2 = mb[j1 + 1];
    } else {
      i1 = mb[s];
      i2 = mb[s + 1];
    }
    int nsrow = i2 - i1;         // "total" number of rows in all the supernodes
                                 // (diagonal+off-diagonal)
    int nsrow2 = nsrow - nscol;  // "total" number of rows in all the off-diagonal supernodes
    int ps2    = i1 + nscol;     // offset into rowind

    int psx;  // offset into data,   Lx[s][s]
    if (ptr_by_column) {
      psx = colptr[j1];
    } else {
      psx = colptr[s];
    }

    /* diagonal block */
    // for each column (or row due to symmetry), the diagonal supernodal block
    // is stored (in ascending order of row indexes) first so that we can do
    // TRSM on the diagonal block
    for (int ii = 0; ii < nscol; ii++) {
      // lower-triangular part
      for (int jj = 0; jj < ii; jj++) {
        hv(hr(j1 + ii) + jj) = Lx[psx + (ii + jj * nsrow)];
      }
      // diagonal
      if (unit_diag) {
        hv(hr(j1 + ii) + ii) = one;
      } else {
        hv(hr(j1 + ii) + ii) = Lx[psx + (ii + ii * nsrow)];
      }
    }
    /* off-diagonal blocks */
    int nsup = 0;
    for (int jj = 0; jj < nscol; jj++) {
      for (int kk = 0; kk < nsrow2; kk++) {
        int irow = rowind[ps2 + kk];

        hv(off(irow) + jj) = Lx[psx + (nscol + kk + jj * nsrow)];

        supid = map(irow);
        if (check(supid) == 0) {
          check(supid) = 1;
          sup(nsup)    = supid;
          nsup++;
        }
      }
    }
    // shift pointers, and reset check
    for (int i = 0; i < nsup; i++) {
      supid = sup(i);
      for (int ii = nb[supid]; ii < nb[supid + 1]; ii++) {
        off(ii) += nscol;
      }
      check(supid) = 0;
    }
  }

#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  std::cout << "    read_supernodal_valuesLt" << std::endl;
  std::cout << "    * Matrix size = " << n << std::endl;
  std::cout << "    * Total nnz   = " << hr(n) << std::endl;
  std::cout << "    * nnz / n     = " << hr(n) / n << std::endl;

  double time = timer.seconds();
  std::cout << "    > Time : " << time << std::endl;
#endif

  // invert blocks (TODO done on host for now)
  invert_supernodal_columns(kernelHandle, unit_diag, nsuper, nb, hr, hc, hv);
  // deepcopy
  Kokkos::deep_copy(values_view, hv);

  // create crs
  crsmat_t crsmat("CrsMatrix", n, values_view, static_graph);

  return crsmat;
}

/* =========================================================================================
 */
template <typename crsmat_t, typename KernelHandle, typename host_crsmat_t>
void split_crsmat(KernelHandle *kernelHandleL, host_crsmat_t superluL) {
  using graph_t        = typename crsmat_t::StaticCrsGraphType;
  using row_map_view_t = typename graph_t::row_map_type::non_const_type;
  using cols_view_t    = typename graph_t::entries_type::non_const_type;
  using values_view_t  = typename crsmat_t::values_type::non_const_type;

  using row_map_view_host_t = typename row_map_view_t::HostMirror;
  using cols_view_host_t    = typename cols_view_t::HostMirror;
  using values_view_host_t  = typename values_view_t::HostMirror;

  using scalar_t  = typename KernelHandle::nnz_scalar_t;
  using size_type = typename KernelHandle::size_type;

  const scalar_t zero(0.0);

  // get sparse-triangular solve handle
  auto *handleL = kernelHandleL->get_sptrsv_handle();

#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  Kokkos::Timer timer;
  Kokkos::Timer timer2;
  double time1 = 0.0;
  double time2 = 0.0;
  double time3 = 0.0;
  double time4 = 0.0;
#endif
  // ===================================================================
  // number of supernodes per level
  auto nodes_per_level  = handleL->get_nodes_per_level();
  auto hnodes_per_level = Kokkos::create_mirror_view(nodes_per_level);
  Kokkos::deep_copy(hnodes_per_level, nodes_per_level);

  // id of supernodes at each level
  auto nodes_grouped_by_level      = handleL->get_nodes_grouped_by_level();
  auto nodes_grouped_by_level_host = Kokkos::create_mirror_view(nodes_grouped_by_level);
  Kokkos::deep_copy(nodes_grouped_by_level_host, nodes_grouped_by_level);

  // load graphs
  auto graphL = handleL->get_graph_host();

  // crsgraph for L
  int nrows     = graphL.numRows();
  auto row_mapL = graphL.row_map;
  auto entriesL = graphL.entries;

  auto values  = superluL.values;
  auto valuesL = Kokkos::create_mirror_view(values);
  Kokkos::deep_copy(valuesL, values);

  int node_count = 0;  // number of supernodes processed
  int nlevels    = handleL->get_num_levels();

  bool invert_offdiag       = handleL->get_invert_offdiagonal();
  const int *supercols_host = handleL->get_supercols_host();

#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  timer.reset();
#endif
  // count total nnz
  int newNnz = 0;
  for (int j = 0; j < nrows; j++) {
    for (size_type k = row_mapL(j); k < row_mapL(j + 1); k++) {
      if (valuesL(k) != zero) {
        newNnz++;
      }
    }
  }
  // allocate for all the subgraphs
  row_map_view_t total_rowmap_view(Kokkos::view_alloc(Kokkos::WithoutInitializing, "rowmap_view"),
                                   2 * nlevels * (nrows + 1));
  cols_view_t total_column_view(Kokkos::view_alloc(Kokkos::WithoutInitializing, "colmap_view"), newNnz);
  values_view_t total_values_view(Kokkos::view_alloc(Kokkos::WithoutInitializing, "values_view"), newNnz);
  // create host-mirrors
  row_map_view_host_t total_hr = Kokkos::create_mirror_view(total_rowmap_view);
  cols_view_host_t total_hc    = Kokkos::create_mirror_view(total_column_view);
  values_view_host_t total_hv  = Kokkos::create_mirror_view(total_values_view);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time4 = timer.seconds();
#endif

// form crsgraph for each submatrix at each level
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  int oldNnz = row_mapL(nrows);
#endif
  newNnz = 0;
  std::vector<crsmat_t> sub_crsmats(nlevels);
  std::vector<crsmat_t> diag_blocks(nlevels);
  int offset_view = 0;
  for (int lvl = 0; lvl < nlevels; ++lvl) {
// > count nnz
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    timer.reset();
#endif
    int nnzL      = 0;
    int nnzD      = 0;
    int lvl_nodes = hnodes_per_level(lvl);  // number of supernodes at this level
    for (int league_rank = 0; league_rank < lvl_nodes; league_rank++) {
      auto s = nodes_grouped_by_level_host(node_count + league_rank);

      // supernodal column size
      int j1    = supercols_host[s];
      int j2    = supercols_host[s + 1];
      int nscol = j2 - j1;  // number of columns in the s-th supernode column
      for (int j = j1; j < j2; j++) {
        for (size_type k = row_mapL(j); k < row_mapL(j + 1); k++) {
          if (valuesL(k) != zero) {
            if (invert_offdiag || k >= row_mapL(j) + nscol) {
              nnzL++;
            } else {
              nnzD++;
            }
          }
        }
      }
    }

#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    timer2.reset();
#endif
    // create subviews for the subgraph
    using range_type  = Kokkos::pair<int, int>;
    int offset_rowmap = lvl * 2 * (nrows + 1);
    row_map_view_t rowmap_view =
        Kokkos::subview(total_rowmap_view, range_type(offset_rowmap, offset_rowmap + (nrows + 1)));
    cols_view_t column_view   = Kokkos::subview(total_column_view, range_type(offset_view, offset_view + nnzL));
    values_view_t values_view = Kokkos::subview(total_values_view, range_type(offset_view, offset_view + nnzL));

    row_map_view_host_t hr = Kokkos::subview(total_hr, range_type(offset_rowmap, offset_rowmap + (nrows + 1)));
    cols_view_host_t hc    = Kokkos::subview(total_hc, range_type(offset_view, offset_view + nnzL));
    values_view_host_t hv  = Kokkos::subview(total_hv, range_type(offset_view, offset_view + nnzL));
    offset_view += nnzL;

    // create subviews for the subgraph, just for diagonal blocks
    offset_rowmap += nrows + 1;
    row_map_view_t rowmapD_view =
        Kokkos::subview(total_rowmap_view, range_type(offset_rowmap, offset_rowmap + (nrows + 1)));
    cols_view_t columnD_view   = Kokkos::subview(total_column_view, range_type(offset_view, offset_view + nnzD));
    values_view_t valuesD_view = Kokkos::subview(total_values_view, range_type(offset_view, offset_view + nnzD));

    row_map_view_host_t hrD = Kokkos::subview(total_hr, range_type(offset_rowmap, offset_rowmap + (nrows + 1)));
    cols_view_host_t hcD    = Kokkos::subview(total_hc, range_type(offset_view, offset_view + nnzD));
    values_view_host_t hvD  = Kokkos::subview(total_hv, range_type(offset_view, offset_view + nnzD));
    offset_view += nnzD;
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    time3 += timer2.seconds();
#endif

    // create subgraph
    hr(0)  = 0;
    hrD(0) = 0;
    if (handleL->transpose_spmv()) {  // NOTE: it always transpose, for now
      // >> store submatrix with transpose (CSR -> CSC or CSC -> CSR) <<
      // count nnz / row
      for (int j = 0; j <= nrows; j++) {
        hr(j)  = 0;
        hrD(j) = 0;
      }
      nnzL = 0;
      nnzD = 0;
      for (int league_rank = 0; league_rank < lvl_nodes; league_rank++) {
        auto s    = nodes_grouped_by_level_host(node_count + league_rank);
        int j1    = supercols_host[s];      // start of this supernode
        int j2    = supercols_host[s + 1];  // start of next supernode
        int nscol = j2 - j1;                // number of columns in the s-th supernode column
        for (int j = j1; j < j2; j++) {
          for (size_type k = row_mapL(j); k < row_mapL(j + 1); k++) {
            if (valuesL(k) != zero) {
              if (invert_offdiag || k >= row_mapL(j) + nscol) {
                hr(1 + entriesL(k))++;
                nnzL++;
              } else {
                hrD(1 + entriesL(k))++;
                nnzD++;
              }
            }
          }
        }
      }
      for (int j = 0; j < nrows; j++) {
        hr(j + 1) += hr(j);
        hrD(j + 1) += hrD(j);
      }
      // insert nzs
      for (int league_rank = 0; league_rank < lvl_nodes; league_rank++) {
        auto s = nodes_grouped_by_level_host(node_count + league_rank);

        // start/end column id for this supernodal column at this level
        // (these column ids are sorted in ascending order at each level)
        int j1    = supercols_host[s];      // start of this supernode
        int j2    = supercols_host[s + 1];  // start of next supernode
        int nscol = j2 - j1;                // number of columns in the s-th supernode column
        for (int j = j1; j < j2; j++) {
          // diagonals
          for (size_type k = row_mapL(j); k < row_mapL(j) + nscol; k++) {
            // remove zeros
            if (valuesL(k) != zero) {
              if (invert_offdiag) {
                hc(hr(entriesL(k))) = j;
                hv(hr(entriesL(k))) = valuesL(k);
                hr(entriesL(k))++;
              } else {
                hcD(hrD(entriesL(k))) = j;
                hvD(hrD(entriesL(k))) = valuesL(k);
                hrD(entriesL(k))++;
              }
            }
          }

          // off-diagonals
          for (size_type k = row_mapL(j) + nscol; k < row_mapL(j + 1); k++) {
            // remove zeros, and minus for updating off-diagonal elements with
            // Spmv
            if (valuesL(k) != zero) {
              hc(hr(entriesL(k))) = j;
              hv(hr(entriesL(k))) = -valuesL(k);
              hr(entriesL(k))++;
            }
          }
        }
      }
      // fix pointers
      for (int j = nrows; j > 0; j--) {
        hr(j)  = hr(j - 1);
        hrD(j) = hrD(j - 1);
      }
      hr(0)  = 0;
      hrD(0) = 0;
    } else {
      // >> store submatrix without transpose (CSC -> CSC or CSR -> CSR) <<
      nnzL   = 0;
      nnzD   = 0;
      int j0 = 0;  // end of previous supernode at this level (not including
                   // this column)
      for (int league_rank = 0; league_rank < lvl_nodes; league_rank++) {
        auto s = nodes_grouped_by_level_host(node_count + league_rank);

        // start/end column id for this supernodal column at this level
        // (these column ids are sorted in ascending order at each level)
        int j1    = supercols_host[s];      // start of this supernode
        int j2    = supercols_host[s + 1];  // start of next supernode
        int nscol = j2 - j1;                // number of columns in the s-th supernode column

        // insert empty columns for the columns skipped (between the previous
        // and this supernodes)
        for (int j = j0 + 1; j <= j1; j++) {
          hr(j)  = hr(j - 1);
          hrD(j) = hrD(j - 1);
        }

        // insert the columns in this supernode
        for (int j = j1; j < j2; j++) {
          // diagonals
          for (size_type k = row_mapL(j); k < row_mapL(j) + nscol; k++) {
            // remove zeros
            if (valuesL(k) != zero) {
              if (invert_offdiag) {
                hc(nnzL) = entriesL(k);
                hv(nnzL) = valuesL(k);
                nnzL++;
              } else {
                hcD(nnzD) = entriesL(k);
                hvD(nnzD) = valuesL(k);
                nnzD++;
              }
            }
          }

          // off-diagonals
          for (size_type k = row_mapL(j) + nscol; k < row_mapL(j + 1); k++) {
            // remove zeros, and minus for updating off-diagonal elements with
            // Spmv
            if (valuesL(k) != zero) {
              hc(nnzL) = entriesL(k);
              hv(nnzL) = -valuesL(k);
              nnzL++;
            }
          }
          hr(j + 1)  = nnzL;
          hrD(j + 1) = nnzD;
        }
        j0 = j2;  // update the last column of the processed supernode (not
                  // including this column)
      }

      // insert empty columns at the end
      for (int j = j0 + 1; j <= nrows; j++) {
        hr(j)  = hr(j - 1);
        hrD(j) = hrD(j - 1);
      }
    }
    newNnz += nnzL + nnzD;
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    time1 += timer.seconds();
    timer.reset();
#endif

    // create crs-graph
    graph_t sub_graph(column_view, rowmap_view);
    sub_crsmats[lvl] = crsmat_t("CrsMatrix", nrows, values_view, sub_graph);
    if (!invert_offdiag) {
      graph_t diag_graph(columnD_view, rowmapD_view);
      diag_blocks[lvl] = crsmat_t("DiagMatrix", nrows, valuesD_view, diag_graph);
    }
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    // std::cout << "   > split nnz(" << lvl << ") = " << nnzL+nnzD <<
    // std::endl;
    time2 += timer.seconds();
#endif

    // update the number of supernodes processed
    node_count += lvl_nodes;
  }
// deep-copy all the subviews
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  timer.reset();
#endif
  Kokkos::deep_copy(total_rowmap_view, total_hr);
  Kokkos::deep_copy(total_column_view, total_hc);
  Kokkos::deep_copy(total_values_view, total_hv);
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time2 += timer.seconds();
#endif
  handleL->set_submatrices(sub_crsmats);
  if (!invert_offdiag) {
    handleL->set_diagblocks(diag_blocks);
  }
#ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  std::cout << "   split_crsmat" << std::endl;
  std::cout << "   > Time to split to submatrices       : " << time1 << std::endl;
  std::cout << "      + allocate submatrices            : " << time4 << std::endl;
  std::cout << "      + create subviews                 : " << time3 << std::endl;
  std::cout << "   > Time to copy submatrices to device : " << time2 << std::endl;
  std::cout << "   > Total NNZ                          : " << oldNnz << " -> " << newNnz << std::endl << std::endl;
#endif
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/* For numeric computation */
template <typename crsmat_input_t, typename KernelHandle>
void sptrsv_compute(KernelHandle *kernelHandleL, crsmat_input_t L) {
  // ===================================================================
  // load sptrsv-handles
  auto *handleL = kernelHandleL->get_sptrsv_handle();
  if (!(handleL->is_symbolic_complete())) {
    std::cout << std::endl
              << " ** needs to call sptrsv_symbolic before calling sptrsv_numeric **" << std::endl
              << std::endl;
    return;
  }
  bool merged = handleL->get_merge_supernodes();
  if (merged) {
    // TODO: follow what's done in sptrsv_compute in superlu
    std::cout << std::endl << " ** merge is not supported through this interface, yet **" << std::endl << std::endl;
    return;
  }

  // ===================================================================
  // load supernodes
  int nsuper           = handleL->get_num_supernodes();
  const int *supercols = handleL->get_supercols_host();

  // ==============================================
  // load crsGraph
  // auto graph = handleL->get_original_graph_host (); // graph stored in handle
  // (before merge)
  auto graph      = handleL->get_graph();       // graph stored in handle (before merge)
  auto graph_host = handleL->get_graph_host();  // graph stored in handle (before merge)
  auto row_map    = graph_host.row_map;
  auto entries    = graph_host.entries;
  auto nrows      = graph_host.numRows();

  // from input CrsMatrix
  auto values = L.values;  // numerical values from input (host), output will be
                           // stored in handle

  // ==============================================
  // read numerical values of L from Cholmod
  using crsmat_t     = typename KernelHandle::SPTRSVHandleType::crsmat_t;
  bool ptr_by_column = true;
  auto crsmatL       = read_supernodal_values<crsmat_t>(kernelHandleL, nrows, nsuper, ptr_by_column, row_map.data(),
                                                  supercols, row_map.data(), entries.data(), values.data(), graph);

  // ===================================================================
  bool useSpMV = (handleL->get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
                  handleL->get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG);
  if (useSpMV) {
    // ----------------------------------------------------
    // split the matrix into submatrices for spmv at each level
    split_crsmat<crsmat_t>(kernelHandleL, crsmatL);
  }

  // ==============================================
  // save crsmat
  handleL->set_crsmat(crsmatL);

  // ===================================================================
  handleL->set_numeric_complete();
}

}  // namespace Experimental
}  // namespace KokkosSparse

#endif  // KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
#endif  // KOKKOSSPARSE_SPTRSV_SUPERNODE_HPP_
