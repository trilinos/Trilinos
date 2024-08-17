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
#include "KokkosKernels_Utils.hpp"
#include "KokkosKernels_Sorting.hpp"

#include <string.h>
#include "KokkosSparse_CrsMatrix.hpp"

int main(int argc, char *argv[]) {
  typedef int size_type;
  typedef int idx;
  typedef double wt;

  bool symmetrize = false, remove_diagonal = false, transpose = false;
  char *in_mtx = NULL, *out_bin = NULL;
  // bool create_incidence = false;
  for (int i = 1; i < argc; ++i) {
    if (0 == Test::string_compare_no_case(argv[i], "--symmetrize")) {
      symmetrize = true;
    } else if (0 == Test::string_compare_no_case(argv[i], "--remove_diagonal")) {
      remove_diagonal = true;
    } else if (0 == Test::string_compare_no_case(argv[i], "--transpose")) {
      transpose = true;
    } else if (0 == Test::string_compare_no_case(argv[i], "--in_mtx")) {
      in_mtx = argv[++i];
    } else if (0 == Test::string_compare_no_case(argv[i], "--out_mtx")) {
      out_bin = argv[++i];
    } else {
      std::cerr << "Usage:" << argv[0]
                << " --in_mtx matrixfile --out_mtx output_file [--symmetrize] "
                   "[--remove_diagonal] [--transpose]"
                << std::endl;
      std::cerr << "Input format .mtx for matrix market, .bin for binary, .crs "
                   "for crs format"
                << std::endl;
      std::cerr << "Output format .mtx for matrix market, .bin for binary, "
                   ".crs for crs format, .ligra for ligra output format"
                << std::endl;

      exit(1);
    }
  }
  if (in_mtx == NULL || out_bin == NULL) {
    std::cerr << "Usage:" << argv[0]
              << " --in_mtx matrixfile --out_mtx output_file [--symmetrize] "
                 "[--remove_diagonal] [--transpose]"
              << std::endl;
    std::cerr << "Input format .mtx for matrix market, .bin for binary, .crs "
                 "for crs format"
              << std::endl;
    std::cerr << "Output format .mtx for matrix market, .bin for binary, .crs "
                 "for crs format, .ligra for ligra output format"
              << std::endl;

    exit(1);
  }

  Kokkos::initialize(argc, argv);
  {
    typedef Kokkos::DefaultHostExecutionSpace MyExecSpace;

    typedef typename KokkosSparse::CrsMatrix<wt, idx, MyExecSpace, void, size_type> crstmat_t;
    typedef typename crstmat_t::StaticCrsGraphType graph_t;
    typedef typename graph_t::row_map_type::non_const_type row_map_view_t;
    typedef typename graph_t::entries_type::non_const_type cols_view_t;
    typedef typename crstmat_t::values_type::non_const_type values_view_t;

    typedef typename graph_t::row_map_type::const_type c_row_map_view_t;
    typedef typename graph_t::entries_type::const_type c_cols_view_t;
    typedef typename crstmat_t::values_type::const_type c_values_view_t;

    crstmat_t a_crsmat = KokkosKernels::Impl::read_kokkos_crst_matrix<crstmat_t>(in_mtx);

    c_row_map_view_t orm    = a_crsmat.graph.row_map;
    c_cols_view_t oentries  = a_crsmat.graph.entries;
    c_values_view_t ovalues = a_crsmat.values;

    const size_type *prm = orm.data();
    const idx *pentries  = oentries.data();
    const wt *pvals      = ovalues.data();

    idx numrows = a_crsmat.numRows();
    // idx numcols = a_crsmat.numCols();
    idx nnz = ovalues.extent(0);
    std::cout << "numrows :" << numrows << " nnz:" << nnz << std::endl;
    // Kokkos::deep_copy(new_rowmap, a_crsmat.graph.row_map);

    if (remove_diagonal) {
      std::vector<size_type> nrm(numrows + 1, 0);
      std::vector<idx> nentries(nnz + 1);
      std::vector<wt> nvals(nnz + 1);

      for (idx i = 0; i < numrows; ++i) {
        size_type begin = prm[i];
        size_type end   = prm[i + 1];
        for (size_type j = begin; j < end; ++j) {
          idx col = pentries[j];
          // wt val = pvals[j];

          if (i == col) {
            nrm[i] = 1;
            break;
          }
        }
      }

      size_type prefix = 0;
      for (idx i = 0; i <= numrows; ++i) {
        size_type current = nrm[i];
        nrm[i]            = prefix;
        prefix += current;
      }

      for (idx i = 0; i <= numrows; ++i) {
        nrm[i] = prm[i] - nrm[i];
      }

      for (idx i = 0; i < numrows; ++i) {
        size_type begin = prm[i];
        size_type end   = prm[i + 1];

        size_type obegin = nrm[i];

        for (size_type j = begin; j < end; ++j) {
          idx col = pentries[j];
          wt val  = pvals[j];
          if (i != col) {
            nentries[obegin] = col;
            nvals[obegin++]  = val;
          }
        }
        if (obegin != nrm[i + 1]) {
          std::cout << "i:" << i << " nrm[i+1]:" << nrm[i + 1] << " obegin:" << obegin << std::endl;
          exit(1);
        }
      }

      row_map_view_t new_rowmap("new rowmap", numrows + 1);

      cols_view_t new_entries("new colmap", nrm[numrows]);
      values_view_t new_values("new values", nrm[numrows]);

      for (idx i = 0; i <= numrows; ++i) {
        new_rowmap(i) = nrm[i];
      }

      for (size_type i = 0; i < nrm[numrows]; ++i) {
        new_entries(i) = nentries[i];
        new_values(i)  = nvals[i];
      }

      graph_t transpose_graph(new_entries, new_rowmap);
      crstmat_t transpose_matrix("transpose", numrows, new_values, transpose_graph);
      a_crsmat = transpose_matrix;

      orm      = a_crsmat.graph.row_map;
      oentries = a_crsmat.graph.entries;
      ovalues  = a_crsmat.values;

      prm      = orm.data();
      pentries = oentries.data();
      pvals    = ovalues.data();

      numrows = a_crsmat.numRows();
      // numcols = a_crsmat.numCols();
      nnz = ovalues.extent(0);
    }

    if (symmetrize) {
      row_map_view_t new_rowmap;
      cols_view_t new_entries;

      KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap<c_row_map_view_t, c_cols_view_t, row_map_view_t,
                                                             cols_view_t, MyExecSpace>(numrows, orm, oentries,
                                                                                       new_rowmap, new_entries);
      values_view_t new_values("new_values", new_entries.extent(0));

      KokkosKernels::sort_crs_matrix<MyExecSpace, row_map_view_t, cols_view_t, values_view_t>(new_rowmap, new_entries,
                                                                                              new_values);

      graph_t symmetric_graph(new_entries, new_rowmap);
      crstmat_t symmetric_marix("transpose", numrows, new_values, symmetric_graph);
      a_crsmat = symmetric_marix;

      orm      = a_crsmat.graph.row_map;
      oentries = a_crsmat.graph.entries;
      ovalues  = a_crsmat.values;

      prm      = orm.data();
      pentries = oentries.data();
      pvals    = ovalues.data();

      numrows = a_crsmat.numRows();
      // numcols = a_crsmat.numCols();
      nnz = ovalues.extent(0);
    }
    if (transpose) {
      row_map_view_t new_rowmap("new_rowmap", a_crsmat.numCols() + 1);
      cols_view_t new_entries("new_rowmap", a_crsmat.nnz());
      values_view_t new_values("new_rowmap", a_crsmat.nnz());

      KokkosSparse::Impl::transpose_matrix<c_row_map_view_t, c_cols_view_t, c_values_view_t, row_map_view_t,
                                           cols_view_t, values_view_t, row_map_view_t, MyExecSpace>(
          a_crsmat.numRows(), a_crsmat.numCols(), a_crsmat.graph.row_map, a_crsmat.graph.entries, a_crsmat.values,
          new_rowmap, new_entries, new_values);

      std::cout << 1 << std::endl;
      std::cout << 2 << std::endl;

      KokkosKernels::sort_crs_matrix<MyExecSpace, row_map_view_t, cols_view_t, values_view_t>(new_rowmap, new_entries,
                                                                                              new_values);

      std::cout << 3 << std::endl;
      MyExecSpace().fence();
      KokkosKernels::Impl::kk_print_1Dview(new_entries);
      KokkosKernels::Impl::kk_print_1Dview(new_values);

      graph_t transpose_graph(new_entries, new_rowmap);
      crstmat_t transpose_matrix("transpose", a_crsmat.numRows(), new_values, transpose_graph);
      a_crsmat = transpose_matrix;

      orm      = a_crsmat.graph.row_map;
      oentries = a_crsmat.graph.entries;
      ovalues  = a_crsmat.values;

      prm      = orm.data();
      pentries = oentries.data();
      pvals    = ovalues.data();

      numrows = a_crsmat.numRows();
      // numcols = a_crsmat.numCols();
      nnz = ovalues.extent(0);
    }

    KokkosKernels::Impl::write_kokkos_crst_matrix(a_crsmat, out_bin);
  }

  Kokkos::finalize();
}
