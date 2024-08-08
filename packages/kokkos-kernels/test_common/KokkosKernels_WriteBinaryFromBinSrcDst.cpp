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
#include <cstdlib>
#include <iostream>
#include "KokkosKernels_IOUtils.hpp"
#include <string.h>

// typedef long long size_type;
// typedef int lno_t;
// typedef double wt;
typedef size_t size_type;
typedef size_t lno_t;
typedef double wt;

typedef size_t input_lno_t;
int main(int argc, char **argv) {
  char *in_src = NULL, *in_dst = NULL;
  for (int i = 1; i < argc; ++i) {
    if (0 == Test::string_compare_no_case(argv[i], "in_src")) {
      in_src = argv[++i];
    } else if (0 == Test::string_compare_no_case(argv[i], "in_dst")) {
      in_dst = argv[++i];
    } else {
      std::cerr << "Usage:" << argv[0] << " in_src srcs.bin in_dst dsts.bin" << std::endl;
      exit(1);
    }
  }
  if (in_src == NULL || in_dst == NULL) {
    std::cerr << "Usage:" << argv[0] << " in_src srcs.bin in_dst dsts.bin" << std::endl;
    exit(1);
  }

  size_t numEdges = 0;
  size_t *srcs, *dst;  // this type is hard coded
  KokkosKernels::Impl::buildEdgeListFromBinSrcTarg_undirected<size_t>(in_src, in_dst, numEdges, &srcs, &dst);
  std::cout << "read numEdges:" << numEdges << std::endl;
  size_t num_vertex = 0;
  for (size_t i = 0; i < numEdges; ++i) {
    if (srcs[i] == 0 || dst[i] == 0) std::cout << "i:" << i << " src:" << srcs[i] << " dst:" << dst[i] << std::endl;
    if (num_vertex < srcs[i]) num_vertex = srcs[i];
    if (num_vertex < dst[i]) num_vertex = dst[i];
  }
  num_vertex += 1;
  std::cout << "num_vertex:" << num_vertex << std::endl;

  lno_t nv     = num_vertex;
  size_type ne = numEdges * 2;
  // std::vector<wt> ew1(ne);
  // wt *ew = &(ew1[0]);
  wt *ew;
  KokkosKernels::Impl::md_malloc<wt>(&ew, ne);

  size_type *xadj;
  lno_t *adj;

  KokkosKernels::Impl::md_malloc<size_type>(&xadj, nv + 1);
  KokkosKernels::Impl::md_malloc<lno_t>(&adj, ne);

  std::cout << "converting" << std::endl;
  KokkosKernels::Impl::convert_undirected_edge_list_to_csr<size_t, size_type, lno_t>(
      nv, numEdges,  // numEdges should be num undirected edges.
      srcs, dst, xadj, adj);
  delete[] srcs;
  delete[] dst;
  size_t num_cols = 0;
  for (size_t i = 0; i < numEdges * 2; ++i) {
    if (num_cols < adj[i]) num_cols = adj[i];
  }
  std::cout << "num_cols:" << num_cols << std::endl;
  // std::vector<size_type> i_xadj(ne / 2 + 1);
  size_type *i_xadj;
  KokkosKernels::Impl::md_malloc<size_type>(&i_xadj, ne / 2 + 1);
  lno_t *i_adj;
  KokkosKernels::Impl::md_malloc<lno_t>(&i_adj, ne);
  // std::vector<lno_t> i_adj(ne);

  std::cout << "writing original" << std::endl;
  KokkosKernels::Impl::write_graph_bin(nv, ne, xadj, adj, ew, "actual.bin");

  std::cout << "calculating incidence transpose" << std::endl;
  KokkosKernels::Impl::kk_sequential_create_incidence_matrix_transpose(nv, ne, xadj, adj,
                                                                       &(i_xadj[0]),  // output. preallocated
                                                                       &(i_adj[0])    // output. preallocated
  );
  std::cout << "writing bin incidence transpose" << std::endl;
  KokkosKernels::Impl::write_graph_bin(lno_t(ne / 2), ne, &(i_xadj[0]), &(i_adj[0]), ew, "incidence-transpose.bin");
  size_type *i_adj2;
  KokkosKernels::Impl::md_malloc<size_type>(&i_adj2, ne);

  std::cout << "calculating incidence " << std::endl;
  KokkosKernels::Impl::kk_sequential_create_incidence_matrix(nv, xadj, adj,
                                                             &(i_adj2[0])  // output. preallocated
  );
  std::cout << "writing bin incidence" << std::endl;
  KokkosKernels::Impl::write_graph_bin(nv, ne, xadj, i_adj2, ew, "incidence.bin");

  lno_t average_degree = ne / nv;
  std::vector<lno_t> row_sizes(nv, 0);
  for (lno_t i = 0; i < nv; ++i) {
    size_type row_s = xadj[i + 1] - xadj[i];
    if (row_s > 1000)
      std::cout << "row:" << i << " size:" << row_s << " average_degree:" << average_degree << std::endl;
    row_sizes[row_s] += 1;
  }

  for (lno_t i = 0; i < nv; ++i) {
    if (row_sizes[i] != 0) {
      std::cout << row_sizes[i] << " rows has " << i << " nonzeroes" << std::endl;
    }
  }
  delete[] i_xadj;

  delete[] i_adj;

  delete[] xadj;

  delete[] adj;

  delete[] ew;
}
