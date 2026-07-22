// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"

//#include "Tpetra_RowMatrixTransposer.hpp"

#include "KokkosKernels_Utils.hpp"
#include "KokkosSparse_Utils.hpp"

namespace Zoltan2
{

// Some utility functions to help with coloring
namespace Impl
{

// Check coloring is valid for a given graph
template <typename LO, typename GO, typename NO, typename list_of_colors_t>
bool
check_coloring(
  const Tpetra::CrsGraph<LO, GO, NO> &graph, 
  const list_of_colors_t &list_of_colors)
{
  typedef typename list_of_colors_t::execution_space execution_space;

  Teuchos::RCP<const Teuchos::Comm<int>> comm = graph.getRowMap()->getComm();
  const int rank = comm->getRank();
  auto local_graph = graph.getLocalGraphDevice();
  const size_t num_rows = graph.getLocalNumRows();
  size_t num_conflict = 0;

  Kokkos::parallel_reduce(
      "check_coloring()", Kokkos::RangePolicy<execution_space>(0, num_rows),
      KOKKOS_LAMBDA(const size_t row, size_t &lcl_conflict) {
        const size_t entry_begin = local_graph.row_map(row);
        const size_t entry_end   = local_graph.row_map(row + 1);
        for (size_t entry = entry_begin; entry < entry_end; entry++)
        {
          const size_t col = local_graph.entries(entry);
          for (size_t entry2 = entry_begin; entry2 < entry_end; ++entry2)
          {
            const size_t col2 = local_graph.entries(entry2);
            if (col != col2 && list_of_colors[col] == list_of_colors[col2])
            {
              ++lcl_conflict;
              Kokkos::printf(
                "proc = %i : Invalid coloring!  Local row %zu"
                " and columns %zu, %zu have the same color %i\n",
                rank, row, col, col2, list_of_colors[col]);
            }
          }
        }
      },
      num_conflict);

  size_t global_num_conflict = 0;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &num_conflict, 
                     &global_num_conflict);

  return global_num_conflict == 0;
}

//////////////////////////////////////////////////////////////////////////////
template <typename LocalCrsGraphType>
LocalCrsGraphType
compute_local_transpose_graph(
  const LocalCrsGraphType &local_graph, 
  const size_t num_cols)
{
  using KokkosSparse::Impl::transpose_graph;

  typedef LocalCrsGraphType graph_t;
  typedef typename graph_t::row_map_type lno_view_t;
  typedef typename graph_t::entries_type lno_nnz_view_t;
  typedef typename lno_view_t::non_const_type trans_row_view_t;
  typedef typename lno_nnz_view_t::non_const_type trans_nnz_view_t;
  typedef typename lno_view_t::execution_space exec_t;

  // Build tranpose graph
  const size_t num_rows = local_graph.row_map.extent(0) - 1;
  const size_t num_nnz  = local_graph.entries.extent(0);

  trans_row_view_t trans_row_map("trans_row_map", num_cols + 1);
  trans_nnz_view_t trans_entries("trans_entries", num_nnz);

  transpose_graph<lno_view_t, lno_nnz_view_t, trans_row_view_t,
                  trans_nnz_view_t, trans_row_view_t, exec_t>(
      num_rows, num_cols, local_graph.row_map, local_graph.entries,
      trans_row_map, trans_entries);

  graph_t local_trans_graph(trans_entries, trans_row_map);

  return local_trans_graph;
}

//////////////////////////////////////////////////////////////////////////////
// template <typename LO, typename GO, typename NO>
// Teuchos::RCP<const Tpetra::CrsGraph<LO,GO,NO> >
// compute_transpose_graph(const Tpetra::CrsGraph<LO,GO,NO>& graph)
// {
//   typedef double SC;
//   Teuchos::RCP< Tpetra::CrsMatrix<SC,LO,GO,NO> > matrix =
//     Teuchos::rcp(new Tpetra::CrsMatrix<double,LO,GO,NO>(Teuchos::rcp(&graph,false)));
//   Tpetra::RowMatrixTransposer<SC,LO,GO,NO> transposer(matrix);
//   Teuchos::RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> > trans_matrix =
//     transposer.createTranspose();
//   return trans_matrix->getCrsGraph();
// }

//////////////////////////////////////////////////////////////////////////////
template <typename LO, typename GO, typename NO>
Teuchos::RCP<Tpetra::CrsGraph<LO, GO, NO>>
compute_transpose_graph(const Tpetra::CrsGraph<LO, GO, NO> &graph)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::CrsGraph<LO, GO, NO> graph_t;
  typedef typename graph_t::local_graph_device_type local_graph_t;

  // Transpose local graph
  local_graph_t local_graph = graph.getLocalGraphDevice();
  local_graph_t local_trans_graph = 
                   compute_local_transpose_graph(local_graph,
                                                 graph.getLocalNumCols());

  // Build (possibly overlapped) transpose graph using original graph's
  // column map as the new row map, and vice versa
  //
  // Normally, for a transpose graph, we'd use
  //     original matrix's RangeMap as the DomainMap, and
  //     original matrix's DomainMap as the RangeMap.
  // Moreover, we're not doing SpMV with the transpose graph, so
  // you might think we wouldn't care about RangeMap and DomainMap.
  // But we DO care, because exportAndFillCompleteCrsGraph uses the
  // shared (overlapped) transpose matrix's RangeMap as the RowMap of the
  // transpose matrix, while the Zoltan callbacks are assuming the
  // transpose matrix's RowMap is the same as the original matrix's.
  // So we'll use the original matrix's RowMap as the RangeMap.
  RCP<graph_t> trans_graph_shared = rcp(new graph_t(
                  local_trans_graph, graph.getColMap(), graph.getRowMap(),
                  graph.getRangeMap(), graph.getRowMap()));

  RCP<graph_t> trans_graph;

  // Export graph to non-overlapped distribution if necessary.
  // If the exporter is null, we don't need to export
  RCP<const Tpetra::Export<LO, GO, NO>> exporter = 
                    trans_graph_shared->getExporter();
  if (exporter == Teuchos::null)
    trans_graph = trans_graph_shared;
  else
  {
    RCP<const graph_t> trans_graph_shared_const = trans_graph_shared;
    trans_graph = Tpetra::exportAndFillCompleteCrsGraph(
                                trans_graph_shared_const, *exporter,
                                Teuchos::null, Teuchos::null, Teuchos::null);
  }

  return trans_graph;
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
Teuchos::RCP<Tpetra::CrsGraph<LO, GO, NO>>
compute_transpose_graph(const Tpetra::CrsMatrix<SC, LO, GO, NO> &matrix)
{
  return compute_transpose_graph(*(matrix.getCrsGraph()));
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
Teuchos::RCP<Tpetra::CrsGraph<LO, GO, NO>>
compute_transpose_graph(const Tpetra::BlockCrsMatrix<SC, LO, GO, NO> &matrix)
{
  return compute_transpose_graph(matrix.getCrsGraph());
}
} // namespace Impl
} // namespace Zoltan2
