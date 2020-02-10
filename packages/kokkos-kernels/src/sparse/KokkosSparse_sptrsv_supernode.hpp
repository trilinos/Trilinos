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

#ifndef KOKKOSSPARSE_SPTRSV_SUPERNODE_HPP_
#define KOKKOSSPARSE_SPTRSV_SUPERNODE_HPP_

#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACKE) && defined(KOKKOSKERNELS_ENABLE_TPL_CBLAS)
#include "cblas.h"
#include "lapacke.h"

namespace KokkosSparse {
namespace Experimental {

class sort_indices {
   public:
     sort_indices(int* rowinds) : rowinds_(rowinds){}
     bool operator()(int i, int j) const { return rowinds_[i] < rowinds_[j]; }
   private:
     int* rowinds_; // rowindices
};

/* ========================================================================================= */
template <typename graph_t>
graph_t
read_supernodal_graphL(bool cusparse, bool merge,
                       int n, int nsuper, bool ptr_by_column, int *mb,
                       int *nb, int *colptr, int *rowind) {

  using row_map_view_t = typename graph_t::row_map_type::non_const_type;
  using cols_view_t    = typename graph_t::entries_type::non_const_type;
  using integer_view_host_t = Kokkos::View<int*, Kokkos::HostSpace>;

  int nnzA;
  if (ptr_by_column) {
    nnzA = colptr[n] - colptr[0];
  } else {
    nnzA = colptr[nsuper] - colptr[0];
  }

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
      if (cusparse) {
        hr(j+1) = hr(j) + nsrow - jj;
      } else {
        hr(j+1) = hr(j) + nsrow;
      }
      j++;
    }
    if (nsrow > max_nnz_per_row) {
      max_nnz_per_row = nsrow;
    }
  }

  integer_view_host_t sorted_rowind_view ("sorted_rowind", max_nnz_per_row);
  int *sorted_rowind = sorted_rowind_view.data ();
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
      if (!cusparse) {
        // explicitly store zeros in upper-part
        for (int jj = ii+1; jj < nscol; jj++) {
          hc(hr(j1+jj)) = j1+ii;
          hr(j1+jj) ++;
        }
      }
    }

    /* off-diagonal blocks */
    if (merge) {
      // sort rowind (to merge supernodes)
      for (int ii = 0; ii < nsrow2; ii++) {
        sorted_rowind[ii] = ii;
      }
      std::sort(&(sorted_rowind[0]), &(sorted_rowind[nsrow2]), sort_indices(&rowind[ps2]));
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
template <typename input_graph_t>
void check_supernode_sizes(const char *title, int n, int nsuper, int *nb, input_graph_t &graph) {

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
template <typename host_graph_t, typename graph_t>
host_graph_t
generate_supernodal_graph(bool col_major, graph_t &graph, int nsuper, int *nb) {

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
  //int *map = new int[n];
  integer_view_host_t map ("map", n);
  for (int s = 0; s < nsuper; s++) {
    for (int j = nb[s]; j < nb[s+1]; j++) {
      map (j) = s;
    }
  }

  // count non-empty supernodal blocks
  row_map_view_host_t hr ("rowmap_view", nsuper+1);
  for (int s = 0; s < nsuper; s++ ) {
    hr (s) = 0;
  }

  integer_view_host_t check ("check", nsuper);
  Kokkos::deep_copy (check, 0);

  int nblocks = 0;
  for (int s = 0; s < nsuper; s++) {
    int j1 = nb[s];
    for (int i = row_map_host (j1); i < row_map_host (j1+1);) {
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
      i ++;
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
    for (int i = row_map_host (j1); i < row_map_host (j1+1);) {
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
      i ++;
    }
    if (!col_major) {
      hr (s+1) = nblocks;
    }
    // reset check
    if (!col_major) {
      for (int s2 = hr(s); s2 < hr(s+1); s2++) {
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
  using row_map_view_t = typename graph_t::row_map_type::non_const_type;
  using cols_view_t    = typename graph_t::entries_type::non_const_type;
  row_map_view_t colptr ("rowind", nsuper+1);
  cols_view_t rowind ("colptr", row_mapL (nsuper)); // over-estimate
  cols_view_t edges ("edges", nsuper); // workspace
  colptr (0) = totedges;
  for (int s = 0; s < nsuper; s ++) {
    // count # of edges (search for first matching nonzero)
    int nedges = 0;
    int k1 = 1 + row_mapL (s); // skip diagonal
    int k2 = 1 + row_mapU (s); // skip diagonal
    for (; k1 < row_mapL (s+1); k1++) {
       // look for match
       while (entriesL (k1) > entriesU (k2) && k2 < row_mapU (s+1)) {
         k2 ++;
       }
       if (entriesL (k1) <= entriesU (k2)) {
         edges (nedges) = entriesL (k1);
         nedges ++;
         if (entriesL (k1) == entriesU (k2)) {
           break;
         }
      }
    }
    // store the edges
    for (int k = 0; k < nedges; k++) {
      rowind (totedges) = edges (k);
      totedges ++;
    }
    colptr (s+1) = totedges;
  }

  // compress dag into graph
  cols_view_t hc ("colmap_view", totedges);
  for (int k = colptr (0); k < colptr (nsuper); k++) {
    hc (k) = rowind (k);
  }

  graph_t static_graph (hc, colptr);
  return static_graph;
}



/* ========================================================================================= */
template <typename input_graph_t>
void merge_supernodal_graph(int *p_nsuper, int *nb,
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
template <typename output_graph_t, typename input_graph_t>
output_graph_t
generate_merged_supernodal_graph(bool lower, 
                                 int nsuper, int *nb,
                                 int nsuper2, int *nb2,
                                 input_graph_t &graph, int *nnz) {


  //auto graphL = L.graph; // in_graph
  auto row_map = graph.row_map;
  auto entries = graph.entries;

  // ----------------------------------------------------------
  // now let me find nsrow for the merged supernode
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
      int j1 = nb[s];
      for (int k = row_map[j1]; k < row_map[j1+1]; k++) {
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
  using row_map_view_t = typename output_graph_t::row_map_type::non_const_type;
  using cols_view_t    = typename output_graph_t::entries_type::non_const_type;

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
/* For numeric computation                                                                   */

/* ========================================================================================= */
template <typename crsmat_t, typename input_crsmat_t, typename graph_t>
crsmat_t
read_merged_supernodes(int nsuper, const int *nb,
                       bool lower, bool unit_diag, bool invert_diag, bool invert_offdiag,
                       input_crsmat_t &L, graph_t &static_graph) {

  using values_view_t  = typename crsmat_t::values_type::non_const_type;
  using scalar_t = typename values_view_t::value_type;
  using scalar_view_host_t = Kokkos::View<scalar_t*, Kokkos::HostSpace>;

  const scalar_t zero (0.0);
  const scalar_t one (1.0);

  // original matrix
  auto graphL = L.graph; // in_graph
  auto row_mapL = graphL.row_map;
  auto entriesL = graphL.entries;
  auto valuesL  = L.values;

  Kokkos::Timer tic;
  double time1 = 0.0;
  double time2 = 0.0;

  // merged graph
  auto rowmap_view = static_graph.row_map;
  auto column_view = static_graph.entries;

  auto hr = Kokkos::create_mirror_view (rowmap_view);
  auto hc = Kokkos::create_mirror_view (column_view);
  Kokkos::deep_copy (hr, rowmap_view);
  Kokkos::deep_copy (hc, column_view);

  // ----------------------------------------------------------
  // now let's copy numerical values
  int n = graphL.numRows ();
  scalar_view_host_t dwork ("dwork", n);
  Kokkos::deep_copy (dwork, zero);

  auto nnzA = hr (n);
  values_view_t values_view ("values_view", nnzA);
  auto hv = Kokkos::create_mirror_view (values_view);

  for (int s2 = 0; s2 < nsuper; s2++) {
    for (int j = nb[s2]; j < nb[s2+1]; j++) {
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

    int j1 = nb[s2];
    int nsrow = hr(j1+1) - hr(j1);
    int nscol = nb[s2+1]-nb[s2];
    if (invert_diag) {
      auto nnzD = hr (j1);
      char uplo_char = (lower ? 'L' : 'U');
      char diag_char = (unit_diag ? 'U' : 'N');

      tic.reset ();
      if (std::is_same<scalar_t, double>::value) {
        LAPACKE_dtrtri (LAPACK_COL_MAJOR,
                        uplo_char, diag_char, nscol,
                        reinterpret_cast <double*> (&hv(nnzD)), nsrow);
      }
      else if (std::is_same<scalar_t, std::complex<double>>::value ||
               std::is_same<scalar_t, Kokkos::complex<double>>::value) {
        LAPACKE_ztrtri (LAPACK_COL_MAJOR,
                        uplo_char, diag_char, nscol, 
                        reinterpret_cast <lapack_complex_double*> (&hv(nnzD)), nsrow);
      }
      else {
        throw std::runtime_error( "Unsupported scalar type for calling trtri");
      }

      time1 += tic.seconds ();
      if (nsrow > nscol && invert_offdiag) {
        CBLAS_UPLO uplo_cblas = (lower ? CblasLower : CblasUpper);
        CBLAS_DIAG diag_cblas = (unit_diag ? CblasUnit : CblasNonUnit);

        tic.reset ();
        if (std::is_same<scalar_t, double>::value) {
          cblas_dtrmm (CblasColMajor,
                CblasRight, uplo_cblas, CblasNoTrans, diag_cblas,
                nsrow-nscol, nscol,
                1.0, reinterpret_cast <double*> (&hv(nnzD)), nsrow,
                     reinterpret_cast <double*> (&hv(nnzD+nscol)), nsrow);
        } else {
          // NOTE: use double pointers
          scalar_t alpha = one;
          cblas_ztrmm (CblasColMajor,
                CblasRight, uplo_cblas, CblasNoTrans, diag_cblas,
                nsrow-nscol, nscol,
                reinterpret_cast <double*> (&alpha),
                reinterpret_cast <double*> (&hv(nnzD)), nsrow,
                reinterpret_cast <double*> (&hv(nnzD+nscol)), nsrow);
        }
        time2 += tic.seconds ();
      }
    }
  }

  // deepcopy
  Kokkos::deep_copy (values_view, hv);
  #ifdef KOKKOS_SPTRSV_SUPERNODE_PROFILE
  std::cout << "   read_merged_supernodes" << std::endl;
  std::cout << "   > Time for inversion::trtri : " << time1 << std::endl;
  std::cout << "   > Time for inversion::trmm  : " << time2 << std::endl;
  #endif

  // create crs
  crsmat_t crsmat("CrsMatrix", n, values_view, static_graph);
  return crsmat;
}


/* ========================================================================================= */
template <typename crsmat_t, typename graph_t, typename scalar_t>
crsmat_t
read_supernodal_valuesL(bool cusparse, bool merge, bool invert_diag, bool invert_offdiag,
                        bool unit_diag, int n, int nsuper, bool ptr_by_column, int *mb,
                        int *nb, int *colptr, int *rowind, scalar_t *Lx,
                        graph_t &static_graph) {

  using  values_view_t = typename crsmat_t::values_type::non_const_type;
  using integer_view_host_t = Kokkos::View<int*, Kokkos::HostSpace>;

  const scalar_t zero (0.0);
  const scalar_t one (1.0);

  Kokkos::Timer tic;
  double time1 = 0.0; // time for trtri
  double time2 = 0.0; // time for trmm to offdiagonal

  // total nnz
  int nnzL;
  if (ptr_by_column) {
    nnzL = colptr[n] - colptr[0];
  } else {
    nnzL = colptr[nsuper] - colptr[0];
  }
  auto rowmap_view = static_graph.row_map;
  auto column_view = static_graph.entries;
  values_view_t values_view ("values_view", nnzL);

  auto hr = Kokkos::create_mirror_view (rowmap_view);
  auto hc = Kokkos::create_mirror_view (column_view);
  auto hv = Kokkos::create_mirror_view (values_view);
  Kokkos::deep_copy (hr, rowmap_view);
  Kokkos::deep_copy (hc, column_view);
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
    if (invert_diag) {
      tic.reset ();
      char diag_char = (unit_diag ? 'U' : 'N');
      if (std::is_same<scalar_t, double>::value) {
        LAPACKE_dtrtri (LAPACK_COL_MAJOR,
                        'L', diag_char, nscol,
                        reinterpret_cast <double*> (&Lx[psx]), nsrow);
      } else {
        LAPACKE_ztrtri (LAPACK_COL_MAJOR,
                        'L', diag_char, nscol,
                        reinterpret_cast <lapack_complex_double*> (&Lx[psx]), nsrow);
      }
      time1 += tic.seconds ();

      if (nsrow2 > 0 && invert_offdiag) {
        tic.reset ();
        CBLAS_DIAG diag_int = (unit_diag ? CblasUnit : CblasNonUnit);
        if (std::is_same<scalar_t, double>::value) {
          cblas_dtrmm (CblasColMajor,
                CblasRight, CblasLower, CblasNoTrans, diag_int,
                nsrow2, nscol,
                1.0, reinterpret_cast <double*> (&Lx[psx]), nsrow,
                     reinterpret_cast <double*> (&Lx[psx+nscol]), nsrow);
        } else {
          // NOTE: use double pointers
          scalar_t alpha = one;
          cblas_ztrmm (CblasColMajor,
                CblasRight, CblasLower, CblasNoTrans, diag_int,
                nsrow2, nscol,
                reinterpret_cast <double*> (&alpha),
                reinterpret_cast <double*> (&Lx[psx]), nsrow,
                reinterpret_cast <double*> (&Lx[psx+nscol]), nsrow);
        }
        time2 += tic.seconds ();
      }
    }
    for (int jj = 0; jj < nscol; jj++) {
      if (!cusparse) {
        // explicitly store zeros in upper-part
#if 0
        for (int ii = 0; ii < jj; ii++) {
          hv(hr(j1+jj) + ii) = zero;
        }
#endif
        hr(j1+jj) += jj;
      }
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
      std::sort(sorted_rowind.data (), sorted_rowind.data () + nsrow2, sort_indices(&rowind[ps2]));
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

  std::cout << "    > Time for trtri on diagonal    : " << time1 << std::endl;
  std::cout << "    > Time for trmm  on off-diagonal: " << time2 << std::endl;
  #endif

  // deepcopy
  Kokkos::deep_copy (values_view, hv);
  // create crs
  crsmat_t crsmat ("CrsMatrix", n, values_view, static_graph);
  return crsmat;
}


/* ========================================================================================= */
template <typename crsmat_t, typename KernelHandle, typename host_crsmat_t>
void split_crsmat(KernelHandle *kernelHandleL, host_crsmat_t superluL) {

  using graph_t        = typename crsmat_t::StaticCrsGraphType;
  using row_map_view_t = typename graph_t::row_map_type::non_const_type;
  using    cols_view_t = typename graph_t::entries_type::non_const_type;
  using  values_view_t = typename crsmat_t::values_type::non_const_type;

  using row_map_view_host_t = typename row_map_view_t::HostMirror;
  using    cols_view_host_t = typename cols_view_t::HostMirror;
  using  values_view_host_t = typename values_view_t::HostMirror;

  using scalar_t = typename KernelHandle::nnz_scalar_t;

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
  values_view_host_t valuesL = Kokkos::create_mirror_view (values);
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
        for (int k = row_mapL (j); k < row_mapL (j+1); k++) {
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
          for (int k = row_mapL (j); k < row_mapL (j+1); k++) {
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
          for (int k = row_mapL (j); k < row_mapL (j) + nscol; k++) {
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
          for (int k = row_mapL (j) + nscol; k < row_mapL (j+1); k++) {
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
          for (int k = row_mapL (j); k < row_mapL (j) + nscol; k++) {
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
          for (int k = row_mapL (j) + nscol; k < row_mapL (j+1); k++) {
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

} // namespace Experimental
} // namespace KokkosSparse

#endif // KOKKOSKERNELS_ENABLE_TPL_LAPACKE & KOKKOSKERNELS_ENABLE_TPL_CBLAS
#endif // KOKKOSSPARSE_SPTRSV_SUPERNODE_HPP_

