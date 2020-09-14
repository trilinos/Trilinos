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

#ifndef KOKKOSSPARSE_SPTRSV_SUPERNODE_HPP_
#define KOKKOSSPARSE_SPTRSV_SUPERNODE_HPP_

// trmm & trtri are called on Host
#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS)

#include "KokkosBlas3_trmm.hpp"
#include "KokkosBlas_trtri.hpp"
#include "KokkosSparse_sptrsv.hpp"


namespace KokkosSparse {
namespace Experimental {

template<typename ordinal_type>
class sort_indices {
   public:
     sort_indices(ordinal_type* rowinds) : rowinds_(rowinds){}
     bool operator()(int i, int j) const { return rowinds_[i] < rowinds_[j]; }
   private:
     ordinal_type* rowinds_; // rowindices
};

/* ========================================================================================= */
template <typename host_graph_t, typename graph_t>
graph_t deep_copy_graph (host_graph_t &host_graph) {
  // load graph on host
  auto row_map = host_graph.row_map;
  auto entries = host_graph.entries;
  auto nrows = host_graph.numRows ();
  auto nnz = row_map (nrows);

  // create graph on device
  using row_map_view_t = typename graph_t::row_map_type::non_const_type;
  using cols_view_t    = typename graph_t::entries_type::non_const_type;
  row_map_view_t rowmap_view ("rowmap_view", nrows+1);
  cols_view_t    column_view ("colmap_view", nnz);

  // copy graph to device
  Kokkos::deep_copy (rowmap_view, row_map);
  Kokkos::deep_copy (column_view, entries);
  graph_t static_graph (column_view, rowmap_view);
  return static_graph;
}

/* ========================================================================================= */
template <typename graph_t, typename ptr_type, typename size_type, typename ordinal_type, typename KernelHandle>
graph_t
read_supernodal_graphL(KernelHandle kernelHandle, int n, int nsuper, int nnzA, bool ptr_by_column,
                       ptr_type *mb, size_type *nb, ordinal_type *rowind) {

  using row_map_view_t = typename graph_t::row_map_type::non_const_type;
  using cols_view_t    = typename graph_t::entries_type::non_const_type;
  using integer_view_host_t = Kokkos::View<ordinal_type*, Kokkos::HostSpace>;

  // load parameters
  auto *handle = kernelHandle->get_sptrsv_handle ();
  bool merge = handle->get_merge_supernodes ();

  row_map_view_t rowmap_view ("rowmap_view", n+1);
  cols_view_t    column_view ("colmap_view", nnzA);
  auto hr = Kokkos::create_mirror_view (rowmap_view);
  auto hc = Kokkos::create_mirror_view (column_view);

  // compute offset for each row
  int j = 0;
  int max_nnz_per_row = 0;
  hr(j) = 0;
  for (int s = 0 ; s < nsuper ; s++) {
    int j1 = nb[s];
    int j2 = nb[s+1];
    // number of columns in the s-th supernode column
    int nscol = j2 - j1;

    int i1, i2;
    if (ptr_by_column) {
      i1 = mb[j1];
      i2 = mb[j1+1];
    } else {
      i1 = mb[s];
      i2 = mb[s+1];
    }
    // "total" number of rows in all the supernodes (diagonal+off-diagonal)
    int nsrow = i2 - i1;

    for (int jj = 0; jj < nscol; jj++) {
      hr(j+1) = hr(j) + nsrow;
      j++;
    }
    if (nsrow > max_nnz_per_row) {
      max_nnz_per_row = nsrow;
    }
  }

  integer_view_host_t sorted_rowind_view ("sorted_rowind", max_nnz_per_row);
  ordinal_type *sorted_rowind = sorted_rowind_view.data ();
  // store L in csr
  for (int s = 0 ; s < nsuper ; s++) {
    int j1 = nb[s];
    int j2 = nb[s+1];
    int nscol = j2 - j1;      // number of columns in the s-th supernode column

    int i1, i2;
    if (ptr_by_column) {
      i1 = mb[j1];
      i2 = mb[j1+1];
    } else {
      i1 = mb[s];
      i2 = mb[s+1];
    }
    int nsrow  = i2 - i1;    // "total" number of rows in all the supernodes (diagonal+off-diagonal)
    int nsrow2 = nsrow - nscol;  // "total" number of rows in all the off-diagonal supernodes
    int ps2    = i1 + nscol;     // offset into rowind

    /* diagonal block */
    for (int ii = 0; ii < nscol; ii++) {
      // lower-triangular part
      for (int jj = 0; jj < ii; jj++) {
        hc(hr(j1+jj)) = j1+ii;
        hr(j1+jj) ++;
      }
      // diagonal
      hc(hr(j1+ii)) = j1+ii;
      hr(j1+ii) ++;
      // explicitly store zeros in upper-part
      for (int jj = ii+1; jj < nscol; jj++) {
        hc(hr(j1+jj)) = j1+ii;
        hr(j1+jj) ++;
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
      int ii = (merge ? sorted_rowind[kk] : kk); // sorted rowind
      int i = rowind[ps2 + ii];
      for (int jj = 0; jj < nscol; jj++) {
        hc(hr(j1+jj)) = i;
        hr(j1+jj) ++;
      }
    }
  }

  // fix hr
  for (int i = n; i >= 1; i--) {
    hr(i) = hr(i-1);
  }
  hr(0) = 0;

  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  std::cout << "    * Matrix size = " << n << std::endl;
  std::cout << "    * Total nnz   = " << hr (n) << std::endl;
  std::cout << "    * nnz / n     = " << hr (n)/n << std::endl;
  #endif

  // deepcopy
  Kokkos::deep_copy (rowmap_view, hr);
  Kokkos::deep_copy (column_view, hc);

  // create crs
  graph_t static_graph (column_view, rowmap_view);
  return static_graph;
}


/* ========================================================================================= */
template <typename graph_t, typename ptr_type, typename size_type, typename ordinal_type, typename KernelHandle>
graph_t
read_supernodal_graphLt(KernelHandle kernelHandle, int n, int nsuper, bool ptr_by_column,
                        ptr_type *mb, size_type *nb, ordinal_type *rowind) {

  using row_map_view_t = typename graph_t::row_map_type::non_const_type;
  using cols_view_t    = typename graph_t::entries_type::non_const_type;
  using integer_view_host_t = Kokkos::View<int*, Kokkos::HostSpace>;

  // load parameters
  auto *handle = kernelHandle->get_sptrsv_handle ();
  bool merge = handle->get_merge_supernodes ();

  /* create a map from row id to supernode id */
  integer_view_host_t map ("map", n);
  int supid = 0;
  for (int k = 0; k < nsuper; k++) {
    int j1 = nb[k];
    int j2 = nb[k+1];
    for (int j = j1; j < j2; j++) {
      map (j) = supid;
    }
    supid ++;
  }

  row_map_view_t rowmap_view ("rowmap_view", n+1);
  auto hr = Kokkos::create_mirror_view (rowmap_view);
  Kokkos::deep_copy (hr, 0);

  integer_view_host_t sup ("sup", nsuper);
  integer_view_host_t check ("check", nsuper);
  Kokkos::deep_copy (check, 0);

  // compute offset for each row
  int nnzA = 0;
  for (int s = 0; s < nsuper; s++) {
    int j1 = nb[s];
    int j2 = nb[s+1];
    // number of columns in the s-th supernode column
    int nscol = j2 - j1;

    int i1, i2;
    if (ptr_by_column) {
      i1 = mb[j1];
      i2 = mb[j1+1];
    } else {
      i1 = mb[s];
      i2 = mb[s+1];
    }

    // diagonal blocks
    for (int ii = j1; ii < j2; ii++) {
      hr(ii+1) += nscol;
      nnzA += nscol;
    }

    // offdiagonal blocks
    int nsup = 0;
    int ps2   = i1 + nscol; // offset into rowind
    int nsrow = i2 - i1;
    int nsrow2 = nsrow-nscol;
    for (int kk = 0; kk < nsrow2; kk++) {
      int irow = rowind[ps2 + kk];
      supid = map (irow);
      if (check (supid) == 0) {
        for (int ii = nb[supid]; ii < nb[supid+1]; ii++) {
          hr (ii + 1) += nscol;
          nnzA += nscol;
        }
        check (supid) = 1;
        sup (nsup) = supid;
        nsup ++;
      }
    }
    // reset check
    for (int i = 0; i < nsup; i++ ) {
      check (sup (i)) = 0;
    }
  }
  for (int i = 0; i < n; i++) {
    hr(i+1) += hr(i);
  }
  cols_view_t    column_view ("colmap_view", nnzA);
  auto hc = Kokkos::create_mirror_view (column_view);

  // pointer to off-diagonals (diagonal comes first)
  integer_view_host_t off ("off", 1+n);
  for (int s = 0; s < nsuper; s++) {
    int i1 = nb[s];
    int i2 = nb[s+1];
    // number of columns in the s-th supernode column
    int nscol = i2 - i1;
    for (int ii = i1; ii < i2; ii++) {
      off (ii) = hr (ii) + nscol;
    }
  }

  // store L in csr
  for (int s = 0; s < nsuper; s++) {
    int j1 = nb[s];
    int j2 = nb[s+1];
    int nscol = j2 - j1;      // number of columns in the s-th supernode column

    int i1, i2;
    if (ptr_by_column) {
      i1 = mb[j1];
      i2 = mb[j1+1];
    } else {
      i1 = mb[s];
      i2 = mb[s+1];
    }
    int nsrow  = i2 - i1;    // "total" number of rows in all the supernodes (diagonal+off-diagonal)
    int nsrow2 = nsrow - nscol;  // "total" number of rows in all the off-diagonal supernodes
    int ps2    = i1 + nscol;     // offset into rowind

    /* diagonal block */
    for (int ii = 0; ii < nscol; ii++) {
      // explicitly store zeros in upper-part
      for (int jj = 0; jj < nscol; jj++) {
        hc(hr(j1+ii)+jj) = j1+jj;
      }
    }

    /* off-diagonal blocks */
    int nsup = 0;
    for (int kk = 0; kk < nsrow2; kk++) {
      int irow = rowind[ps2 + kk];
      supid = map (irow);
      if (check (supid) == 0) {
        for (int ii = nb[supid]; ii < nb[supid+1]; ii++) {
          for (int jj = 0; jj < nscol; jj++) {
            hc (off(ii)+jj) = j1+jj;
          }
          off(ii) += nscol;
        }
       check (supid) = 1;
       sup (nsup) = supid;
       nsup ++;
      }
    }
    // reset check
    for (int i = 0; i < nsup; i++ ) {
      check (sup (i)) = 0;
    }
  }

  if (merge) {
    // they should be sorted?
  }

  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  std::cout << "    * Matrix size = " << n << std::endl;
  std::cout << "    * Total nnz   = " << hr (n) << std::endl;
  std::cout << "    * nnz / n     = " << hr (n)/n << std::endl;
  #endif

  // deepcopy
  Kokkos::deep_copy (rowmap_view, hr);
  Kokkos::deep_copy (column_view, hc);

  // create crs
  graph_t static_graph (column_view, rowmap_view);
  return static_graph;
}


/* ========================================================================================= */
template <typename input_graph_t, typename input_size_type>
void check_supernode_sizes(const char *title, int n, int nsuper, input_size_type *nb, input_graph_t &graph) {

  auto rowmap_view = graph.row_map;
  auto hr = Kokkos::create_mirror_view (rowmap_view);
  Kokkos::deep_copy (hr, rowmap_view);

  int min_nsrow = 0, max_nsrow = 0, tot_nsrow = 0;
  int min_nscol = 0, max_nscol = 0, tot_nscol = 0;
  for (int s = 0; s <nsuper; s++) {
    int j1 = nb[s];
    int j2 = nb[s+1];

    int nscol = j2 - j1;
    int nsrow = hr(j1+1) - hr(j1);

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
  std::cout << std::endl
            << " ------------------------------------- "
            << std::endl << std::endl;
  std::cout << " " << title << std::endl;
  std::cout << "  + nsuper = " << nsuper << std::endl;
  std::cout << "  > nsrow: min = " << min_nsrow << ", max = " << max_nsrow
            << ", avg = " << tot_nsrow/nsuper << std::endl;
  std::cout << "  > nscol: min = " <<  min_nscol << ", max = " << max_nscol
            << ", avg = " << tot_nscol/nsuper << std::endl;
  std::cout << "    + Matrix size = " << n << std::endl;
  std::cout << "    + Total nnz   = " << hr(n) << std::endl;
  std::cout << "    + nnz / n     = " << hr(n)/n << std::endl;
}


/* ========================================================================================= */
template <typename host_graph_t, typename graph_t, typename input_size_type>
host_graph_t
generate_supernodal_graph(bool col_major, graph_t &graph, int nsuper, const input_size_type *nb) {

  using size_type = typename graph_t::size_type;
  using cols_view_host_t    = typename host_graph_t::entries_type::non_const_type;
  using row_map_view_host_t = typename host_graph_t::row_map_type::non_const_type;
  using integer_view_host_t = Kokkos::View<int*, Kokkos::HostSpace>;

  int n = graph.numRows ();
  auto row_map = graph.row_map;
  auto entries = graph.entries;

  auto row_map_host = Kokkos::create_mirror_view (row_map);
  auto entries_host = Kokkos::create_mirror_view (entries);
  Kokkos::deep_copy (row_map_host, row_map);
  Kokkos::deep_copy (entries_host, entries);

  // map col/row to supernode
  integer_view_host_t map ("map", n);
  for (int s = 0; s < nsuper; s++) {
    for (int j = nb[s]; j < nb[s+1]; j++) {
      map (j) = s;
    }
  }

  // count non-empty supernodal blocks
  row_map_view_host_t hr ("rowmap_view", nsuper+1);
  integer_view_host_t check ("check", nsuper);
  Kokkos::deep_copy (hr, 0);
  Kokkos::deep_copy (check, 0);

  int nblocks = 0;
  for (int s = 0; s < nsuper; s++) {
    int j1 = nb[s];
    int j2 = j1+1;  // based on the first row
    for (size_type i = row_map_host (j1); i < row_map_host (j2); i++) {
      int s2 = map (entries_host (i));
      // supernodal blocks may not be filled with zeros
      // so need to check by each row
      // (also rowids are not sorted)
      if (check (s2) == 0) {
        check (s2) = 1;
        nblocks ++;
        // count blocks per row for col_major
        hr (s2+1) ++;
      }
    }
    // reset check
    Kokkos::deep_copy (check, 0);
  }

  cols_view_host_t hc ("colmap_view", nblocks);
  if (col_major) {
    // convert to offset
    for (int s = 0; s < nsuper; s++) {
      hr (s+1) += hr (s);
    }
  }

  nblocks = 0;
  for (int s = 0; s < nsuper; s++) {
    int j1 = nb[s];
    int j2 = j1+1;  // based on the first row
    for (size_type i = row_map_host (j1); i < row_map_host (j2); i++) {
      int s2 = map (entries_host (i));
      // supernodal blocks may not be filled with zeros
      // so need to check by each row
      // (also rowids are not sorted)
      if (check (s2) == 0) {
        check (s2) = 1;
        if (col_major) {
          hc (hr (s2)) = s;
          hr (s2) ++;
        } else {
          hc (nblocks) = s2;
        }
        nblocks ++;
      }
    }
    if (!col_major) {
      hr (s+1) = nblocks;
    }
    // reset check
    if (!col_major) {
      for (size_type s2 = hr(s); s2 < hr(s+1); s2++) {
        check (hc(s2)) = 0;
      }
    } else {
      // NOTE: nonzero supernodes in s-th col are not stored
      Kokkos::deep_copy (check, 0);
    }
  }
  // fix hr
  if (col_major) {
    for (int s = nsuper; s > 0; s--) {
      hr (s) = hr (s-1);
    }
    hr (0) = 0;
  }
  // sort column ids per row
  for (int s = 0; s < nsuper; s++) {
    std::sort(&(hc (hr (s))), &(hc (hr (s+1))));
  }

  host_graph_t static_graph (hc, hr);
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
  int totedges = 0;
  using size_type      = typename graph_t::size_type;
  using row_map_view_t = typename graph_t::row_map_type::non_const_type;
  using cols_view_t    = typename graph_t::entries_type::non_const_type;

  row_map_view_t colptr ("rowind", nsuper+1);
  cols_view_t rowind ("colptr", row_mapL (nsuper)); // over-estimate
  cols_view_t edges ("edges", nsuper); // workspace
  cols_view_t check ("edges", nsuper); // workspace
  Kokkos::deep_copy (check, 0);
  colptr (0) = totedges;
  for (int s = 0; s < nsuper; s ++) {
    // count # of edges (search for first matching nonzero)
    size_type nedges = 0;
    size_type k1 = 1 + row_mapL (s); // skip diagonal
    size_type k2 = 1 + row_mapU (s); // skip diagonal
    for (; k1 < row_mapL (s+1); k1++) {
       // look for match
       while (k2+1 < row_mapU (s+1) && entriesL (k1) > entriesU (k2)) {
         k2 ++;
       }
       if (k2+1 >= row_mapU (s+1) || entriesL (k1) <= entriesU (k2)) {
         edges (nedges) = entriesL (k1);
         nedges ++;
         if (entriesL (k1) == entriesU (k2)) {
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
      rowind (totedges) = edges (k);
      totedges ++;
    }
    colptr (s+1) = totedges;
  }

  // compress dag into graph
  cols_view_t hc ("colmap_view", totedges);
  for (size_type k = colptr (0); k < colptr (nsuper); k++) {
    hc (k) = rowind (k);
  }

  graph_t static_graph (hc, colptr);
  return static_graph;
}



/* ========================================================================================= */
template <typename input_graph_t, typename input_size_type>
void merge_supernodal_graph(int *p_nsuper, input_size_type *nb,
                            bool col_majorL, input_graph_t &graphL,
                            bool col_majorU, input_graph_t &graphU,
                            int *etree) {

  int nsuper = *p_nsuper;

  // ---------------------------------------------------------------
  // looking for supernodes to merge (i.e., dense diagonal blocks)
  int nsuper2 = 0;
  auto supL = generate_supernodal_graph<input_graph_t, input_graph_t> (!col_majorL, graphL, nsuper, nb);
  auto supU = generate_supernodal_graph<input_graph_t, input_graph_t> ( col_majorU, graphU, nsuper, nb);

  auto row_mapL = supL.row_map;
  auto entriesL = supL.entries;

  auto row_mapU = supU.row_map;
  auto entriesU = supU.entries;

  // map the first supernode
  using integer_view_host_t = Kokkos::View<int*, Kokkos::HostSpace>;
  integer_view_host_t map ("map", nsuper); // map old to new supernodes
  map (0) = 0;
  for (int s = 0; s < nsuper-1; s++) {
    int s2 = s;
    bool merged = false;
    do {
      // check if L(s2+1:end, s2) and L(s2+1:end, s2+1)are the same
      bool mergedL = false;
      int k1 = row_mapL[s2+1] - row_mapL[s2];
      int k2 = row_mapL[s2+2] - row_mapL[s2+1];
      if (k1 == k2+1) {
        mergedL = true;
        for (int k = 0; k < k2 && mergedL; k++) {
          if (entriesL[row_mapL[s2]+k+1] != entriesL[row_mapL[s2+1]+k]) {
            mergedL = false;
          }
        }
      }
      // check if U(s2+1:end, s2) and U(s2+1:end, s2+1) are the same
      bool mergedU = false;
      k1 = row_mapU[s2+1] - row_mapU[s2];
      k2 = row_mapU[s2+2] - row_mapU[s2+1];
      if (k1 == k2+1) {
        mergedU = true;
        for (int k = 0; k < k2 && mergedU; k++) {
          if (entriesU[row_mapU[s2]+k+1] != entriesU[row_mapU[s2+1]+k]) {
            mergedU = false;
          }
        }
      }
      merged = (mergedL && mergedU);
      if (merged) {
        //printf( "  >> merge s2+1=%d(%dx%d, row=%d:%d) with s=%d(%dx%d) <\n",
        //              s2+1,nb[s2+2]-nb[s2+1],nb[s2+2]-nb[s2+1], nb[s2+1],nb[s2+2]-1,
        //              s,nb[s+1]-nb[s],nb[s+1]-nb[s]);
        map (s2+1) = nsuper2;
        s2 ++;
      } else {
        //printf( "  -- not merge s2+1=%d(%dx%d, row=%d:%d) with s=%d(%dx%d) --\n",
        //           s2+1,nb[s2+2]-nb[s2+1],nb[s2+2]-nb[s2+1],nb[s2+1],nb[s2+2]-1,
        //           s,nb[s+1]-nb[s],nb[s+1]-nb[s]);
        map (s2+1) = nsuper2+1;
      }
    } while (merged && s2 < nsuper-1);
    s = s2;
    nsuper2 ++;
  }
  nsuper2 = map (nsuper-1)+1;
  //printf( " nsuper2 = %d\n",nsuper2 );
  //printf( " map:\n" );
  //for (int s = 0; s < nsuper; s++) printf( "   %d %d\n",s,map (s) );

  // ----------------------------------------------------------
  // make sure each of the merged supernodes has the same parent in the etree
  int nsuper3 = 0;
  integer_view_host_t map2;
  if (etree != nullptr) {
    nsuper3 = 0;
    map2 = integer_view_host_t ("map2", nsuper); // map old to new supernodes
    for (int s2 = 0, s = 0; s2 < nsuper2; s2++) {
      // look for parent of the first supernode
      int s3 = s;
      while (etree[s3] != -1 && map (etree[s3]) == map (s3)) {
        s3 ++;
      }
      map2 (s) = nsuper3;
      int p = (etree[s3] == -1 ? -1 : map (etree[s3]));

      // go through the rest of the supernode in this merged supernode
      s++;
      while (s < nsuper && map (s) == s2) {
        int q = (etree[s3] == -1 ? -1 : map (etree[s3]));
        while (etree[s3] != -1 && map (etree[s3]) == map (s3)) {
          s3 ++;
          q = (etree[s3] == -1 ? -1 : map (etree[s3]));
        }

        if (q != p) {
          p = q;
          nsuper3 ++;
        }
        map2 (s) = nsuper3;
        s ++;
      }
      nsuper3 ++;
    }
  } else {
    nsuper3 = nsuper2;
    map2 = map;
  }
  //printf( " nsuper3 = %d\n",nsuper3 );
  //printf( " map:\n" );
  //for (int s = 0; s < nsuper; s++) printf( "   %d %d\n",s,map2 (s) );

  // ----------------------------------------------------------
  // construct new supernodes
  integer_view_host_t nb2 ("nb2", 1+nsuper3);
  for (int s2 = 0, s = 0; s2 < nsuper3; s2++) {
    nb2 (1+s2) = 0;
    // merging supernodal rows
    while(s < nsuper && map2 (s) == s2) {
      nb2 (1+s2) += (nb[s+1]-nb[s]);
      s ++;
    }
  }

  // copy back the new supernodes "offsets"
  nb2 (0) = 0;
  for (int s = 0; s < nsuper3; s++) {
    nb2 (s+1) = nb2 (s) + nb2 (s+1);
  }
  // copy nb
  for (int s = 0; s <nsuper3; s++) {
    nb[s+1] = nb2 (s+1);
  }

  // ----------------------------------------------------------
  // construct new etree
  if (etree != nullptr) {
    integer_view_host_t etree2 ("etree2", nsuper3);
    for (int s = 0; s < nsuper; s++) {
      // etree
      int s2 = map2 (s);
      int p = (etree[s] == -1 ? -1 : map2 (etree[s]));
      if (p != s2) {
        etree2 (s2) = p;
      }
    }
    // copy etree
    for (int s = 0; s <nsuper3; s++) {
      etree[s] = etree2 (s);
    }
  }

  *p_nsuper = nsuper3;
}


/* ========================================================================================= */
template <typename output_graph_t, typename input_graph_t, typename input_size_type>
output_graph_t
generate_merged_supernodal_graph(bool lower, 
                                 int nsuper, const input_size_type *nb,
                                 int nsuper2,      input_size_type *nb2,
                                 input_graph_t &graph, int *nnz) {

  using cols_view_t    = typename output_graph_t::entries_type::non_const_type;
  using row_map_view_t = typename output_graph_t::row_map_type::non_const_type;
  using size_type      = typename input_graph_t::size_type;

  // ----------------------------------------------------------
  // now let me find nsrow for the merged supernode
  auto row_map = graph.row_map;
  auto entries = graph.entries;
  int n = graph.numRows ();

  using integer_view_host_t = Kokkos::View<int*, Kokkos::HostSpace>;
  integer_view_host_t mb2 ("mb2", nsuper2);
  integer_view_host_t work1 ("work1", n);
  integer_view_host_t work2 ("work2", n);
  Kokkos::deep_copy (work1, 0);

  int nnzS = 0;
  int nnzA = 0;
  integer_view_host_t rowind ("rowind", row_map(n)); // over-estimate
  integer_view_host_t colptr ("colptr", nsuper2+1);
  colptr (0) = nnzS;
  for (int s2 = 0, s = 0; s2 < nsuper2; s2++) {
    mb2 (s2) = 0;
    // merging supernodal rows
    // NOTE: SuperLU may not fill zeros to fill the supernodes
    //       So, these rows may be just subset of the supernodal rows
    while (s < nsuper && nb[s+1] <= nb2[s2+1]) {
      input_size_type j1 = nb[s];
      for (size_type k = row_map (j1); k < row_map (j1+1); k++) {
        // just taking union of rows
        if (work1 (entries[k]) == 0) {
          work1 (entries[k]) = 1;
          work2 (mb2 (s2)) = entries[k];
          mb2 (s2) ++;
        }
      }
      s++;
    }
    // sort such that diagonal come on the top
    std::sort(work2.data (), work2.data () + mb2 (s2));

    // save nonzero row ids
    if (lower) {
      // lower in csc, diagonal come on top
      for (int k = 0; k < mb2 (s2); k++) {
        rowind (nnzS) = work2 (k);
        work1 (work2 (k)) = 0;
        nnzS ++;
      }
    } else {
      // upper in csc, diagonal block is on bottom right now, but move it to top
      int nd = nb2[s2+1]-nb2[s2]; // size of diagonal block
      // > diagonal block
      for (int k = mb2 (s2) - nd; k < mb2 (s2); k++) {
        rowind (nnzS) = work2 (k);
        work1 (work2 (k)) = 0;
        nnzS ++;
      }
      // > offdiagonal blocks
      for (int k = 0; k < mb2 (s2) - nd; k++) {
        rowind (nnzS) = work2 (k);
        work1 (work2 (k)) = 0;
        nnzS ++;
      }
    }
    colptr (s2+1) = nnzS;
    nnzA += (nb2[s2+1]-nb2[s2]) * mb2 (s2);
  }

  // ----------------------------------------------------------
  // now let's create crs graph
  row_map_view_t rowmap_view ("rowmap_view", n+1);
  cols_view_t    column_view ("colmap_view", nnzA);
  auto hr = Kokkos::create_mirror_view (rowmap_view);
  auto hc = Kokkos::create_mirror_view (column_view);

  nnzA = 0;
  hr(0) = 0;
  for (int s2 = 0; s2 < nsuper2; s2++) {
    for (int j = nb2[s2]; j < nb2[s2+1]; j++) {
      for (int k = colptr (s2); k < colptr (s2+1); k++) {
        hc(nnzA) = rowind (k);
        nnzA ++;
      }
      hr(j+1) = nnzA;
    }
  }
  *nnz = nnzA;

  // deepcopy
  Kokkos::deep_copy (rowmap_view, hr);
  Kokkos::deep_copy (column_view, hc);

  // create crs
  output_graph_t static_graph (column_view, rowmap_view);
  return static_graph;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* For symbolic analysis                                                                     */
template <typename host_graph_t, typename KernelHandle>
void sptrsv_supernodal_symbolic(
    int nsuper, int *supercols, int *etree,
    host_graph_t graphL_host, KernelHandle *kernelHandleL,
    host_graph_t graphU_host, KernelHandle *kernelHandleU) {

  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  double time_seconds = 0.0;
  Kokkos::Timer timer;
  Kokkos::Timer tic;
  timer.reset ();
  #endif

  // ===================================================================
  // load sptrsv-handles
  auto *handleL = kernelHandleL->get_sptrsv_handle ();
  auto *handleU = kernelHandleU->get_sptrsv_handle ();

  // store arguments to handle
  handleL->set_graph_host (graphL_host);
  handleU->set_graph_host (graphU_host);
  handleL->set_supernodes (nsuper, supercols, etree);
  handleU->set_supernodes (nsuper, supercols, etree);

  // ===================================================================
  // load options
  bool col_majorL = handleL->is_column_major ();
  bool col_majorU = handleU->is_column_major ();
  bool merge = handleL->get_merge_supernodes ();
  bool UinCSC = handleU->is_column_major ();
  bool needEtree = (handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
                    handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_ETREE);
  if (needEtree && etree == nullptr) {
    std::cout << std::endl
              << " ** etree needs to be set before calling sptrsv_symbolic with SuperLU **"
              << std::endl << std::endl;
    return;
  }

  // ===================================================================
  // > make a copy of supercols (merge needs both original and merged supercols)
  using integer_view_host_t = typename KernelHandle::SPTRSVHandleType::integer_view_host_t;
  integer_view_host_t supercols_view ("supercols view", 1+nsuper);
  int *supercols_merged = supercols_view.data ();
  for (int i = 0; i <= nsuper; i++) {
    supercols_merged[i] = supercols[i];
  }
  if (merge) {
    // =================================================================
    // merge supernodes
    int nsuper_merged = nsuper;
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    tic.reset ();
    int nrows = graphL_host.numRows ();
    check_supernode_sizes("Original L-structure", nrows, nsuper, supercols_merged, graphL_host);
    check_supernode_sizes("Original U-structure", nrows, nsuper, supercols_merged, graphU_host);
    #endif
    // etree will be updated
    merge_supernodal_graph (&nsuper_merged, supercols_merged,
                            col_majorL, graphL_host, col_majorU, graphU_host,
                            etree);

    // =================================================================
    // generate merged graph for L-solve
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    int nnzL = graphL_host.row_map (nrows);
    #endif
    int nnzL_merged;
    bool lower = true;
    handleL->set_original_graph_host (graphL_host); // save graph before merge
    graphL_host = generate_merged_supernodal_graph<host_graph_t> (lower, nsuper, supercols,
                                                                  nsuper_merged, supercols_merged,
                                                                  graphL_host, &nnzL_merged);
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    time_seconds = tic.seconds ();
    check_supernode_sizes ("After Merge", nrows, nsuper_merged, supercols_merged, graphL_host);
    std::cout << " for L factor:" << std::endl;
    std::cout << "   Merge Supernodes Time: " << time_seconds << std::endl;
    std::cout << "   Number of nonzeros   : " << nnzL << " -> " << nnzL_merged
              << " : " << double(nnzL_merged) / double(nnzL) << "x" << std::endl;
    #endif

    // =================================================================
    // generate merged graph for U-solve
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    tic.reset ();
    int nnzU = graphU_host.row_map (nrows);
    #endif
    int nnzU_merged;
    lower = (UinCSC ? false : true);
    handleU->set_original_graph_host (graphU_host); // save graph before merge
    graphU_host = generate_merged_supernodal_graph<host_graph_t> (lower, nsuper, supercols,
                                                                  nsuper_merged, supercols_merged,
                                                                  graphU_host, &nnzU_merged);
    #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
    time_seconds = tic.seconds ();
    check_supernode_sizes("After Merge", nrows, nsuper_merged, supercols_merged, graphU_host);
    std::cout << " for U factor:" << std::endl;
    std::cout << "   Merge Supernodes Time: " << time_seconds << std::endl;
    std::cout << "   Number of nonzeros   : " << nnzU << " -> " << nnzU_merged
              << " : " << double(nnzU_merged) / double(nnzU) << "x" << std::endl;
    #endif

    // update the number of supernodes
    nsuper = nsuper_merged;
  }
  // replace the supernodal info with the merged ones
  supercols = supercols_merged;

  // ===================================================================
  // copy graph to device
  using graph_t = typename KernelHandle::SPTRSVHandleType::graph_t;
  auto graphL = deep_copy_graph<host_graph_t, graph_t> (graphL_host);
  auto graphU = deep_copy_graph<host_graph_t, graph_t> (graphU_host);

  // ===================================================================
  // save the supernodal info in the handles for L/U solves
  handleL->set_supernodes (nsuper, supercols_view, etree);
  handleU->set_supernodes (nsuper, supercols_view, etree);

  if (handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_DAG ||
      handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG) {
    // generate supernodal graphs for DAG scheduling
    auto supL = generate_supernodal_graph<host_graph_t> (!col_majorL, graphL_host, nsuper, supercols);
    auto supU = generate_supernodal_graph<host_graph_t> ( col_majorU, graphU_host, nsuper, supercols);

    auto dagL = generate_supernodal_dag<host_graph_t> (nsuper, supL, supU);
    auto dagU = generate_supernodal_dag<host_graph_t> (nsuper, supU, supL);
    handleL->set_supernodal_dag (dagL);
    handleU->set_supernodal_dag (dagU);
  }

  // ===================================================================
  // do symbolic for L solve on the host
  auto row_mapL = graphL.row_map;
  auto entriesL = graphL.entries;
  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  tic.reset();
  std::cout << std::endl;
  #endif
  sptrsv_symbolic (kernelHandleL, row_mapL, entriesL);
  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time_seconds = tic.seconds ();
  std::cout << " > Lower-TRI: " << std::endl;
  std::cout << "   Symbolic Time: " << time_seconds << std::endl;
  tic.reset ();
  #endif

  // ===================================================================
  // do symbolic for U solve on the host
  auto row_mapU = graphU.row_map;
  auto entriesU = graphU.entries;
  sptrsv_symbolic (kernelHandleU, row_mapU, entriesU);
  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time_seconds = tic.seconds ();
  std::cout << " > Upper-TRI: " << std::endl;
  std::cout << "   Symbolic Time: " << time_seconds << std::endl;
  #endif

  // ===================================================================
  // save options
  handleL->set_merge_supernodes (merge);
  handleU->set_merge_supernodes (merge);

  // ===================================================================
  // save graphs
  handleL->set_graph (graphL);
  handleU->set_graph (graphU);
  // graph on host (merged)
  handleL->set_graph_host (graphL_host);
  handleU->set_graph_host (graphU_host);

  // ===================================================================
  handleL->set_symbolic_complete ();
  handleU->set_symbolic_complete ();
  handleL->set_etree (etree);
  handleU->set_etree (etree);
  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  time_seconds = timer.seconds ();
  std::cout << "   Total Symbolic Time: " << time_seconds << std::endl << std::endl;
  std::cout << "   Total nnz: " << graphL_host.row_map (nrows) << " + " << graphU_host.row_map (nrows) << std::endl;
  #endif
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* Auxiliary functions for numeric computation                                               */

/* ========================================================================================= */
template <typename KernelHandle, typename input_size_type,
          typename row_map_type, typename index_type, typename values_type>
void
invert_supernodal_columns(KernelHandle kernelHandle, bool unit_diag, int nsuper, const input_size_type *nb, 
                          row_map_type& hr, index_type& hc, values_type& hv) {

  using execution_space = typename values_type::execution_space;
  using memory_space = typename execution_space::memory_space;
  using values_view_t  = typename values_type::non_const_type;
  using scalar_t = typename values_view_t::value_type;
  using range_type = Kokkos::pair<int, int>;

  const scalar_t one (1.0);

  // load parameters
  auto *handle = kernelHandle->get_sptrsv_handle ();
  bool invert_diag = handle->get_invert_diagonal ();
  bool invert_offdiag = handle->get_invert_offdiagonal ();

  // lower is always in CSC, if UinCSC, then lower=false, else lower=true
  bool lower_tri = kernelHandle->is_sptrsv_lower_tri ();
  bool lower = ((lower_tri && handle->is_column_major ()) || (!lower_tri && !handle->is_column_major ()));

  // quick return
  if (!invert_diag) return;

  Kokkos::Timer timer;
  double time1 = 0.0;
  double time2 = 0.0;

  // ----------------------------------------------------------
  // now let's invert some blocks
  for (int s2 = 0; s2 < nsuper; s2++) {
    int j1 = nb[s2];
    int nsrow = hr(j1+1) - hr(j1);
    int nscol = nb[s2+1] - nb[s2];

    auto nnzD = hr (j1);
    char uplo_char = (lower ? 'L' : 'U');
    char diag_char = (unit_diag ? 'U' : 'N');

    Kokkos::View<scalar_t**, Kokkos::LayoutLeft, memory_space, Kokkos::MemoryUnmanaged>
      viewL (&hv(nnzD), nsrow, nscol);
    auto Ljj = Kokkos::subview (viewL, range_type (0, nscol), Kokkos::ALL ());

    timer.reset ();
    KokkosBlas::trtri(&uplo_char, &diag_char, Ljj);
    time1 += timer.seconds ();

    if (nsrow > nscol && invert_offdiag) {
      char side_char = 'R';
      char tran_char = 'N';
      auto Lij = Kokkos::subview (viewL, range_type (nscol, nsrow), Kokkos::ALL ());

      timer.reset ();
      KokkosBlas::trmm (&side_char, &uplo_char,
                        &tran_char, &diag_char,
                        one, Ljj, Lij);
      time2 += timer.seconds ();
    }
  }

  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  std::cout << "   invert_supernodes" << std::endl;
  std::cout << "   > Time for inversion::trtri : " << time1 << std::endl;
  std::cout << "   > Time for inversion::trmm  : " << time2 << std::endl;
  #endif
}


/* ========================================================================================= */
template <typename crsmat_t, typename input_crsmat_t, typename input_ptr_type, typename graph_t,
          typename KernelHandle>
crsmat_t
read_merged_supernodes(KernelHandle kernelHandle, int nsuper, const input_ptr_type *mb,
                       bool unit_diag, input_crsmat_t &L, graph_t &static_graph) {

  using values_view_t  = typename crsmat_t::values_type::non_const_type;
  using scalar_t = typename values_view_t::value_type;
  using scalar_view_host_t = Kokkos::View<scalar_t*, Kokkos::HostSpace>;

  const scalar_t zero (0.0);

  // original matrix
  auto graphL = L.graph; // in_graph
  auto row_mapL = graphL.row_map;
  auto entriesL = graphL.entries;
  auto valuesL  = L.values;

  Kokkos::Timer timer;
  timer.reset ();

  // merged graph
  auto rowmap_view = static_graph.row_map;
  auto column_view = static_graph.entries;

  auto hr = Kokkos::create_mirror_view (rowmap_view);
  auto hc = Kokkos::create_mirror_view (column_view);
  Kokkos::deep_copy (hr, rowmap_view);
  Kokkos::deep_copy (hc, column_view);

  // ----------------------------------------------------------
  // now let's merge supernodes
  int n = graphL.numRows ();
  scalar_view_host_t dwork ("dwork", n);
  Kokkos::deep_copy (dwork, zero);

  auto nnzA = hr (n);
  values_view_t values_view ("values_view", nnzA);
  auto hv = Kokkos::create_mirror_view (values_view);

  for (int s2 = 0; s2 < nsuper; s2++) {
    for (int j = mb[s2]; j < mb[s2+1]; j++) {
      for (int k = row_mapL[j]; k < row_mapL[j+1]; k++) {
        dwork (entriesL[k]) = valuesL[k];
      }
      for (int k = hr (j); k < hr (j+1); k++) {
        hv(k) = dwork (hc(k));
      }
      for (int k = row_mapL[j]; k < row_mapL[j+1]; k++) {
        dwork (entriesL[k]) = zero;
      }
    }
  }

  // invert blocks (TODO done on host for now)
  invert_supernodal_columns (kernelHandle, unit_diag, nsuper, mb, hr, hc, hv);
  // deepcopy
  Kokkos::deep_copy (values_view, hv);

  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  double time = timer.seconds ();
  std::cout << "   read_merged_supernodes" << std::endl;
  std::cout << "   > Time : " << time << std::endl;
  #endif

  // create crs
  crsmat_t crsmat("CrsMatrix", n, values_view, static_graph);

  return crsmat;
}


/* ========================================================================================= */
template <typename crsmat_t, typename graph_t, typename scalar_t,
          typename input_size_type, typename input_ptr_type,
          typename size_type, typename ordinal_type, typename KernelHandle>
crsmat_t
read_supernodal_valuesL(KernelHandle kernelHandle,
                        int n, int nsuper, bool ptr_by_column, const input_size_type *mb, const input_ptr_type *nb,
                        const size_type *colptr, ordinal_type *rowind, scalar_t *Lx, graph_t &static_graph) {

  using  values_view_t = typename crsmat_t::values_type::non_const_type;
  using integer_view_host_t = Kokkos::View<ordinal_type*, Kokkos::HostSpace>;

  const scalar_t zero (0.0);
  const scalar_t one (1.0);

  Kokkos::Timer timer;
  timer.reset ();

  // load parameters
  auto *handle = kernelHandle->get_sptrsv_handle ();
  bool unit_diag = handle->is_unit_diagonal ();
  bool merge = handle->get_merge_supernodes ();

  // load graph
  auto rowmap_view = static_graph.row_map;
  auto column_view = static_graph.entries;
  auto hr = Kokkos::create_mirror_view (rowmap_view);
  auto hc = Kokkos::create_mirror_view (column_view);
  Kokkos::deep_copy (hr, rowmap_view);
  Kokkos::deep_copy (hc, column_view);

  // total nnz
  int nnzL = hr (n);
  values_view_t values_view ("values_view", nnzL);
  auto hv = Kokkos::create_mirror_view (values_view);
  Kokkos::deep_copy (hv, zero); // seems to be needed (instead of zeroing out upper)

  // compute max nnz per row
  int max_nnz_per_row = 0;
  for (int s = 0; s < nsuper; s++) {
    int i1, i2;
    if (ptr_by_column) {
      int j1 = nb[s];

      i1 = mb[j1];
      i2 = mb[j1+1];
    } else {
      i1 = mb[s];
      i2 = mb[s+1];
    }
    // "total" number of rows in all the supernodes (diagonal+off-diagonal)
    int nsrow = i2 - i1;
    if (nsrow > max_nnz_per_row) {
      max_nnz_per_row = nsrow;
    }
  }

  integer_view_host_t sorted_rowind ("sorted_rowind", max_nnz_per_row);
  // store L in csr
  for (int s = 0; s < nsuper; s++) {
    int j1 = nb[s];
    int j2 = nb[s+1];
    int nscol = j2 - j1;      // number of columns in the s-th supernode column

    int i1, i2;
    if (ptr_by_column) {
      i1 = mb[j1];
      i2 = mb[j1+1];
    } else {
      i1 = mb[s];
      i2 = mb[s+1];
    }
    int nsrow  = i2 - i1;    // "total" number of rows in all the supernodes (diagonal+off-diagonal)
    int nsrow2 = nsrow - nscol;  // "total" number of rows in all the off-diagonal supernodes
    int ps2    = i1 + nscol;     // offset into rowind

    int psx;                 // offset into data,   Lx[s][s]
    if (ptr_by_column) {
      psx = colptr[j1];
    } else {
      psx = colptr[s];
    }

    /* diagonal block */
    // for each column (or row due to symmetry), the diagonal supernodal block is stored (in ascending order of row indexes) first
    // so that we can do TRSM on the diagonal block
    for (int jj = 0; jj < nscol; jj++) {
      // shift for explicitly store zeros in upper-part
      hr(j1+jj) += jj;
      // diagonal
      if (unit_diag) {
        hv(hr(j1+jj)) = one;
      } else {
        hv(hr(j1+jj)) = Lx[psx + (jj + jj*nsrow)];
      }
      hr(j1+jj) ++;
      // lower-triangular part
      for (int ii = jj+1; ii < nscol; ii++) {
        hv(hr(j1+jj)) = Lx[psx + (ii + jj*nsrow)];
        hr(j1+jj) ++;
      }
    }
    /* off-diagonal blocks */
    if (merge) {
      // sort rowind (to merge supernodes)
      for (int ii = 0; ii < nsrow2; ii++) {
        sorted_rowind (ii) = ii;
      }
      std::sort(sorted_rowind.data (), sorted_rowind.data () + nsrow2, sort_indices<ordinal_type>(&rowind[ps2]));
    }
    for (int jj = 0; jj < nscol; jj++) {
      for (int kk = 0; kk < nsrow2; kk++) {
      int ii = (merge ? sorted_rowind (kk) : kk); // sorted rowind
        hv(hr(j1+jj)) = Lx[psx + (nscol+ii + jj*nsrow)];
        hr(j1+jj) ++;
      }
    }
  }

  // fix hr
  for (int i = n; i >= 1; i--) {
    hr(i) = hr(i-1);
  }
  hr(0) = 0;

  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  std::cout << "    read_supernodal_valuesL" << std::endl;
  std::cout << "    * Matrix size = " << n << std::endl;
  std::cout << "    * Total nnz   = " << hr (n) << std::endl;
  std::cout << "    * nnz / n     = " << hr (n)/n << std::endl;

  double time = timer.seconds ();
  std::cout << "    > Time : " << time << std::endl;
  #endif

  // invert blocks (TODO done on host for now)
  invert_supernodal_columns (kernelHandle, unit_diag, nsuper, nb, hr, hc, hv);
  // deepcopy
  Kokkos::deep_copy (values_view, hv);

  // create crs
  crsmat_t crsmat ("CrsMatrix", n, values_view, static_graph);

  return crsmat;
}


/* ========================================================================================= */
template <typename crsmat_t, typename graph_t, typename scalar_t,
          typename input_size_type, typename input_ptr_type,
          typename size_type, typename ordinal_type, typename KernelHandle>
crsmat_t
read_supernodal_valuesLt(KernelHandle kernelHandle,
                         int n, int nsuper, bool ptr_by_column, const input_size_type *mb, const input_ptr_type *nb,
                         const size_type *colptr, ordinal_type *rowind, scalar_t *Lx, graph_t &static_graph) {

  using  values_view_t = typename crsmat_t::values_type::non_const_type;
  using integer_view_host_t = Kokkos::View<int*, Kokkos::HostSpace>;

  const scalar_t zero (0.0);
  const scalar_t one (1.0);

  Kokkos::Timer timer;
  timer.reset ();

  // load parameters
  auto *handle = kernelHandle->get_sptrsv_handle ();
  bool unit_diag = handle->is_unit_diagonal ();

  // load graph
  auto rowmap_view = static_graph.row_map;
  auto column_view = static_graph.entries;
  auto hr = Kokkos::create_mirror_view (rowmap_view);
  auto hc = Kokkos::create_mirror_view (column_view);
  Kokkos::deep_copy (hr, rowmap_view);
  Kokkos::deep_copy (hc, column_view);

  // total nnz
  int nnzL = hr (n);
  values_view_t values_view ("values_view", nnzL);
  auto hv = Kokkos::create_mirror_view (values_view);
  Kokkos::deep_copy (hv, zero); // seems to be needed (instead of zeroing out upper)

  /* create a map from row id to supernode id */
  integer_view_host_t map ("map", n);
  int supid = 0;
  for (int k = 0; k < nsuper; k++) {
    int j1 = nb[k];
    int j2 = nb[k+1];
    for (int j = j1; j < j2; j++) {
      map (j) = supid;
    }
    supid ++;
  }

  // pointer to off-diagonals (diagonal comes first)
  integer_view_host_t off ("off", n+1);
  for (int s = 0; s < nsuper; s++) {
    int i1 = nb[s];
    int i2 = nb[s+1];
    // number of columns in the s-th supernode column
    int nscol = i2 - i1;
    for (int ii = i1; ii < i2; ii++) {
      off (ii) = hr (ii) + nscol;
    }
  }

  // store L in csr
  integer_view_host_t sup ("sup", nsuper);
  integer_view_host_t check ("check", nsuper);
  Kokkos::deep_copy (check, 0);
  for (int s = 0; s < nsuper; s++) {

    int j1 = nb[s];
    int j2 = nb[s+1];
    int nscol = j2 - j1;      // number of columns in the s-th supernode column

    int i1, i2;
    if (ptr_by_column) {
      i1 = mb[j1];
      i2 = mb[j1+1];
    } else {
      i1 = mb[s];
      i2 = mb[s+1];
    }
    int nsrow  = i2 - i1;    // "total" number of rows in all the supernodes (diagonal+off-diagonal)
    int nsrow2 = nsrow - nscol;  // "total" number of rows in all the off-diagonal supernodes
    int ps2    = i1 + nscol;     // offset into rowind

    int psx;                 // offset into data,   Lx[s][s]
    if (ptr_by_column) {
      psx = colptr[j1];
    } else {
      psx = colptr[s];
    }

    /* diagonal block */
    // for each column (or row due to symmetry), the diagonal supernodal block is stored (in ascending order of row indexes) first
    // so that we can do TRSM on the diagonal block
    for (int ii = 0; ii < nscol; ii++) {
      // lower-triangular part
      for (int jj = 0; jj < ii; jj++) {
        hv(hr(j1+ii)+jj) = Lx[psx + (ii + jj*nsrow)];
      }
      // diagonal
      if (unit_diag) {
        hv(hr(j1+ii)+ii) = one;
      } else {
        hv(hr(j1+ii)+ii) = Lx[psx + (ii + ii*nsrow)];
      }
    }
    /* off-diagonal blocks */
    int nsup = 0;
    for (int jj = 0; jj < nscol; jj++) {
      for (int kk = 0; kk < nsrow2; kk++) {
        int irow = rowind[ps2 + kk];

        hv(off(irow)+jj) = Lx[psx + (nscol+kk + jj*nsrow)];

        supid = map (irow);
        if (check (supid) == 0) {
          check (supid) = 1;
          sup (nsup) = supid;
          nsup ++;
        }
      }
    }
    // shift pointers, and reset check
    for (int i = 0; i < nsup; i++ ) {
      supid = sup (i);
      for (int ii = nb[supid]; ii < nb[supid+1]; ii++) {
        off (ii) += nscol;
      }
      check (supid) = 0;
    }
  }

  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  std::cout << "    read_supernodal_valuesL" << std::endl;
  std::cout << "    * Matrix size = " << n << std::endl;
  std::cout << "    * Total nnz   = " << hr (n) << std::endl;
  std::cout << "    * nnz / n     = " << hr (n)/n << std::endl;

  double time = timer.seconds ();
  std::cout << "    > Time : " << time << std::endl;
  #endif

  // invert blocks (TODO done on host for now)
  invert_supernodal_columns (kernelHandle, unit_diag, nsuper, nb, hr, hc, hv);
  // deepcopy
  Kokkos::deep_copy (values_view, hv);

  // create crs
  crsmat_t crsmat ("CrsMatrix", n, values_view, static_graph);

  return crsmat;
}


/* ========================================================================================= */
template <typename crsmat_t, typename KernelHandle, typename host_crsmat_t>
void split_crsmat(KernelHandle *kernelHandleL, host_crsmat_t superluL) {

  using        graph_t = typename crsmat_t::StaticCrsGraphType;
  using row_map_view_t = typename graph_t::row_map_type::non_const_type;
  using    cols_view_t = typename graph_t::entries_type::non_const_type;
  using  values_view_t = typename crsmat_t::values_type::non_const_type;

  using row_map_view_host_t = typename row_map_view_t::HostMirror;
  using    cols_view_host_t = typename cols_view_t::HostMirror;
  using  values_view_host_t = typename values_view_t::HostMirror;

  using scalar_t = typename KernelHandle::nnz_scalar_t;
  using size_type = typename KernelHandle::size_type;

  const scalar_t zero (0.0);

  // get sparse-triangular solve handle
  auto *handleL = kernelHandleL->get_sptrsv_handle ();

  Kokkos::Timer timer;
  double time1 = 0.0;
  double time2 = 0.0;
  // ===================================================================
  // number of supernodes per level
  auto nodes_per_level = handleL->get_nodes_per_level ();
  auto hnodes_per_level = Kokkos::create_mirror_view (nodes_per_level);
  Kokkos::deep_copy (hnodes_per_level, nodes_per_level);

  // id of supernodes at each level
  auto nodes_grouped_by_level = handleL->get_nodes_grouped_by_level ();
  auto nodes_grouped_by_level_host = Kokkos::create_mirror_view (nodes_grouped_by_level);
  Kokkos::deep_copy (nodes_grouped_by_level_host, nodes_grouped_by_level);

  // load graphs
  auto graphL = handleL->get_graph_host ();

  // crsgraph for L
  int nrows = graphL.numRows ();
  auto row_mapL = graphL.row_map;
  auto entriesL = graphL.entries;

  auto values = superluL.values;
  auto valuesL = Kokkos::create_mirror_view (values);
  Kokkos::deep_copy (valuesL, values);

  int node_count = 0; // number of supernodes processed
  int nlevels = handleL->get_num_levels ();

  bool invert_offdiag = handleL->get_invert_offdiagonal ();
  const int* supercols_host = handleL->get_supercols_host ();

  // form crsgraph for each submatrix at each level
  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  int oldNnz = row_mapL (nrows);
  #endif
  int newNnz = 0;
  std::vector <crsmat_t> sub_crsmats (nlevels);
  std::vector <crsmat_t> diag_blocks (nlevels);
  for (int lvl = 0; lvl < nlevels; ++lvl) {
    timer.reset ();
    // > count nnz
    int nnzL = 0;
    int nnzD = 0;
    int lvl_nodes = hnodes_per_level (lvl); // number of supernodes at this level
    for (int league_rank = 0; league_rank < lvl_nodes; league_rank++) {
      auto s = nodes_grouped_by_level_host (node_count + league_rank);

      // supernodal column size
      int j1 = supercols_host[s];
      int j2 = supercols_host[s+1];
      int nscol = j2 - j1 ;        // number of columns in the s-th supernode column
      for (int j = j1; j < j2; j++) {
        for (size_type k = row_mapL (j); k < row_mapL (j+1); k++) {
          if (valuesL (k) != zero) {
            if (invert_offdiag || k >= row_mapL (j) + nscol) {
              nnzL ++;
            } else {
              nnzD ++;
            }
          }
        }
      }
    }

    // allocate subgraph
    row_map_view_t rowmap_view ("rowmap_view", nrows+1);
    cols_view_t    column_view ("colmap_view", nnzL);
    values_view_t  values_view ("values_view", nnzL);
    row_map_view_host_t hr = Kokkos::create_mirror_view (rowmap_view);
    cols_view_host_t    hc = Kokkos::create_mirror_view (column_view);
    values_view_host_t  hv = Kokkos::create_mirror_view (values_view);

    // allocate subgraph, just for diagonal blocks
    row_map_view_t rowmapD_view ("rowmapD_view", nrows+1);
    cols_view_t    columnD_view ("colmapD_view", nnzD);
    values_view_t  valuesD_view ("valuesD_view", nnzD);
    row_map_view_host_t hrD = Kokkos::create_mirror_view (rowmapD_view);
    cols_view_host_t    hcD = Kokkos::create_mirror_view (columnD_view);
    values_view_host_t  hvD = Kokkos::create_mirror_view (valuesD_view);

    // create subgraph
    hr (0) = 0;
    hrD (0) = 0;
    if (handleL->transpose_spmv()) { // NOTE: it always transpose, for now
      // >> store submatrix with transpose (CSR -> CSC or CSC -> CSR) <<
      // count nnz / row
      for (int j = 0; j <= nrows; j++) {
        hr (j) = 0;
        hrD (j) = 0;
      }
      nnzL = 0;
      nnzD = 0;
      for (int league_rank = 0; league_rank < lvl_nodes; league_rank++) {
        auto s = nodes_grouped_by_level_host (node_count + league_rank);
        int j1 = supercols_host[s];                  // start of this supernode
        int j2 = supercols_host[s+1];                // start of next supernode
        int nscol = j2 - j1 ;        // number of columns in the s-th supernode column
        for (int j = j1; j < j2; j++) {
          for (size_type k = row_mapL (j); k < row_mapL (j+1); k++) {
            if (valuesL (k) != zero) {
              if (invert_offdiag || k >= row_mapL (j) + nscol) {
                hr (1 + entriesL (k)) ++;
                nnzL ++;
              } else {
                hrD (1 + entriesL (k)) ++;
                nnzD ++;
              }
            }
          }
        }
      }
      for (int j = 0; j < nrows; j++) {
        hr (j+1) += hr (j);
        hrD (j+1) += hrD (j);
      }
      // insert nzs
      for (int league_rank = 0; league_rank < lvl_nodes; league_rank++) {
        auto s = nodes_grouped_by_level_host (node_count + league_rank);

        // start/end column id for this supernodal column at this level
        // (these column ids are sorted in ascending order at each level)
        int j1 = supercols_host[s];                  // start of this supernode
        int j2 = supercols_host[s+1];                // start of next supernode
        int nscol = j2 - j1 ;        // number of columns in the s-th supernode column
        for (int j = j1; j < j2; j++) {
          // diagonals
          for (size_type k = row_mapL (j); k < row_mapL (j) + nscol; k++) {
            // remove zeros
            if (valuesL (k) != zero) {
              if (invert_offdiag) {
                hc (hr (entriesL (k))) = j;
                hv (hr (entriesL (k))) = valuesL (k);
                hr (entriesL (k)) ++;
              } else {
                hcD (hrD (entriesL (k))) = j;
                hvD (hrD (entriesL (k))) = valuesL (k);
                hrD (entriesL (k)) ++;
              }
            }
          }

          // off-diagonals
          for (size_type k = row_mapL (j) + nscol; k < row_mapL (j+1); k++) {
            // remove zeros, and minus for updating off-diagonal elements with Spmv
            if (valuesL (k) != zero) {
              hc (hr (entriesL (k))) = j;
              hv (hr (entriesL (k))) = -valuesL (k);
              hr (entriesL (k)) ++;
            }
          }
        }
      }
      // fix pointers
      for (int j = nrows; j > 0; j--) {
        hr (j) = hr (j-1);
        hrD (j) = hrD (j-1);
      }
      hr (0) = 0;
      hrD (0) = 0;
    } else {
      // >> store submatrix without transpose (CSC -> CSC or CSR -> CSR) <<
      nnzL = 0;
      nnzD = 0;
      int j0 = 0; // end of previous supernode at this level (not including this column)
      for (int league_rank = 0; league_rank < lvl_nodes; league_rank++) {
        auto s = nodes_grouped_by_level_host (node_count + league_rank);

        // start/end column id for this supernodal column at this level
        // (these column ids are sorted in ascending order at each level)
        int j1 = supercols_host[s];                  // start of this supernode
        int j2 = supercols_host[s+1];                // start of next supernode
        int nscol = j2 - j1 ;        // number of columns in the s-th supernode column

        // insert empty columns for the columns skipped (between the previous and this supernodes)
        for (int j = j0+1; j <= j1; j++) {
          hr (j) = hr (j-1);
          hrD (j) = hrD (j-1);
        }

        // insert the columns in this supernode
        for (int j = j1; j < j2; j++) {
          // diagonals
          for (size_type k = row_mapL (j); k < row_mapL (j) + nscol; k++) {
            // remove zeros
            if (valuesL (k) != zero) {
              if (invert_offdiag) {
                hc (nnzL) = entriesL (k);
                hv (nnzL) = valuesL (k);
                nnzL ++;
              } else {
                hcD (nnzD) = entriesL (k);
                hvD (nnzD) = valuesL (k);
                nnzD ++;
              }
            }
          }

          // off-diagonals
          for (size_type k = row_mapL (j) + nscol; k < row_mapL (j+1); k++) {
            // remove zeros, and minus for updating off-diagonal elements with Spmv
            if (valuesL (k) != zero) {
              hc (nnzL) =  entriesL (k);
              hv (nnzL) = -valuesL (k);
              nnzL ++;
            }
          }
          hr (j+1) = nnzL;
          hrD (j+1) = nnzD;
        }
        j0 = j2; // update the last column of the processed supernode (not including this column)
      }

      // insert empty columns at the end
      for (int j = j0+1; j <= nrows; j++) {
        hr (j) = hr (j-1);
        hrD (j) = hrD (j-1);
      }
    }
    newNnz += nnzL+nnzD;
    time1 += timer.seconds ();

    // create crs-graph
    timer.reset ();
    Kokkos::deep_copy (rowmap_view, hr);
    Kokkos::deep_copy (column_view, hc);
    Kokkos::deep_copy (values_view, hv);
    graph_t sub_graph(column_view, rowmap_view);
    sub_crsmats[lvl] = crsmat_t("CrsMatrix", nrows, values_view, sub_graph);
    if (!invert_offdiag) {
      Kokkos::deep_copy (rowmapD_view, hrD);
      Kokkos::deep_copy (columnD_view, hcD);
      Kokkos::deep_copy (valuesD_view, hvD);
      graph_t diag_graph(columnD_view, rowmapD_view);
      diag_blocks[lvl] = crsmat_t("DiagMatrix", nrows, valuesD_view, diag_graph);
    }
    //std::cout << "   > split nnz(" << lvl << ") = " << nnzL+nnzD << std::endl; 
    time2 += timer.seconds ();

    // update the number of supernodes processed
    node_count += lvl_nodes;
  }
  handleL->set_submatrices (sub_crsmats);
  if (!invert_offdiag) {
    handleL->set_diagblocks (diag_blocks);
  }
  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  std::cout << "   split_crsmat" << std::endl;
  std::cout << "   > Time to split to submatrices      : " << time1 << std::endl;
  std::cout << "   > Time to copy submatrices to device: " << time2 << std::endl;
  std::cout << "   > Total NNZ                         : " << oldNnz << " -> " << newNnz
            << std::endl << std::endl;
  #endif
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* For numeric computation                                                                   */
  template <typename crsmat_input_t,
            typename KernelHandle>
  void sptrsv_compute(
      KernelHandle *kernelHandleL,
      crsmat_input_t L)
  {
    // ===================================================================
    // load sptrsv-handles
    auto *handleL = kernelHandleL->get_sptrsv_handle ();
    if (!(handleL->is_symbolic_complete())) {
      std::cout << std::endl
                << " ** needs to call sptrsv_symbolic before calling sptrsv_numeric **"
                << std::endl << std::endl;
      return;
    }
    bool merged = handleL->get_merge_supernodes ();
    if (merged) {
      // TODO: follow what's done in sptrsv_compute in superlu
      std::cout << std::endl
                << " ** merge is not supported through this interface, yet **"
                << std::endl << std::endl;
      return;
    }

    // ===================================================================
    // load supernodes
    int nsuper = handleL->get_num_supernodes ();
    const int *supercols = handleL->get_supercols_host ();

    // ==============================================
    // load crsGraph
    //auto graph = handleL->get_original_graph_host (); // graph stored in handle (before merge)
    auto graph = handleL->get_graph (); // graph stored in handle (before merge)
    auto graph_host = handleL->get_graph_host (); // graph stored in handle (before merge)
    auto row_map = graph_host.row_map;
    auto entries = graph_host.entries;
    auto nrows = graph_host.numRows ();

    // from input CrsMatrix
    auto values = L.values;            // numerical values from input (host), output will be stored in handle

    // ==============================================
    // read numerical values of L from Cholmod
    using crsmat_t = typename KernelHandle::SPTRSVHandleType::crsmat_t;
    bool ptr_by_column = true;
    auto crsmatL = read_supernodal_valuesL<crsmat_t> (kernelHandleL, nrows, nsuper, ptr_by_column, row_map.data (), supercols,
                                                      row_map.data (), entries.data (), values.data (), graph);

    // ===================================================================
    bool useSpMV = (handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
                    handleL->get_algorithm () == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG);
    if (useSpMV) {
      // ----------------------------------------------------
      // split the matrix into submatrices for spmv at each level
      split_crsmat<crsmat_t> (kernelHandleL, crsmatL);
    }

    // ==============================================
    // save crsmat
    handleL->set_crsmat (crsmatL);

    // ===================================================================
    handleL->set_numeric_complete ();
  }

} // namespace Experimental
} // namespace KokkosSparse

#endif // KOKKOSKERNELS_ENABLE_TPL_BLAS
#endif // KOKKOSSPARSE_SPTRSV_SUPERNODE_HPP_

