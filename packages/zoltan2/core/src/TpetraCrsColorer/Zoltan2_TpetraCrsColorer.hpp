// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once

#include <stdexcept>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_BlockMultiVector.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"

#include "Zoltan2_TpetraCrsColorerUtils.hpp"
#include "Zoltan2_TpetraCrsColorer_Zoltan.hpp"
#include "Zoltan2_TpetraCrsColorer_Zoltan2.hpp"

namespace Zoltan2
{

// Base class for coloring Tpetra::CrsMatrix for use in column compression
template <typename CrsMatrixType>
class TpetraCrsColorer
{
public:
  typedef CrsMatrixType matrix_t;
  typedef typename matrix_t::crs_graph_type graph_t;
  typedef typename matrix_t::scalar_type scalar_t;
  typedef typename matrix_t::local_ordinal_type lno_t;
  typedef typename matrix_t::global_ordinal_type gno_t;
  typedef typename matrix_t::node_type node_t;
  typedef typename node_t::device_type device_t;
  typedef typename device_t::execution_space execution_space;
  typedef Kokkos::View<int *, device_t> list_of_colors_t;
  typedef typename list_of_colors_t::HostMirror list_of_colors_host_t;

  // Constructor
  TpetraCrsColorer(const Teuchos::RCP<matrix_t> &matrix_);

  // Destructor
  ~TpetraCrsColorer() {}

  // Compute coloring data
  void
  computeColoring(Teuchos::ParameterList &coloring_params);

  // Compute seed matrix
  template <typename MultiVectorType>
  void
  computeSeedMatrix(MultiVectorType &V, const int color_beg = 0) const;

  // Compute seed matrix with distribution fitted to the graph's column map
  template <typename MultiVectorType>
  void
  computeSeedMatrixFitted(MultiVectorType &V, const int color_beg = 0) const;

  // Reconstruct matrix from supplied compressed vector
  template <typename MultiVectorType>
  void
  reconstructMatrix(MultiVectorType &W, const int color_beg = 0) const;
  template <typename MultiVectorType>
  void
  reconstructMatrix(MultiVectorType &W, matrix_t &mat, const int color_beg = 0) const;

  // Reconstruct matrix from supplied compressed vector fitted to the graph's
  // row map
  template <typename MultiVectorType>
  void
  reconstructMatrixFitted(MultiVectorType &W, const int color_beg = 0) const;
  template <typename MultiVectorType>
  void
  reconstructMatrixFitted(MultiVectorType &W, matrix_t &mat, const int color_beg = 0) const;

  // Return number of colors
  // KDD should num_colors be int or size_t?
  int
  getNumColors() const
  {
    return num_colors;
  }

  // Get color given a column index
  int
  getColor(const size_t col) const
  {
    return list_of_colors_host(col) - 1;
  }

  // Check coloring is valid, i.e., no row has two nonzero columns with
  // same color
  bool
  checkColoring() const
  {
    return Impl::check_coloring(*graph, list_of_colors);
  }

protected:
  Teuchos::RCP<matrix_t> matrix;
  Teuchos::RCP<const graph_t> graph;
  list_of_colors_t list_of_colors;
  list_of_colors_host_t list_of_colors_host;
  int num_colors;
};

//////////////////////////////////////////////////////////////////////////////
// Specialization of TpetraCrsColorer for BlockCrsMatrix
template <typename SC, typename LO, typename GO, typename NO>
class TpetraCrsColorer<Tpetra::BlockCrsMatrix<SC,LO,GO,NO> >
{
public:
  typedef Tpetra::BlockCrsMatrix<SC, LO, GO, NO> matrix_t;
  typedef Tpetra::BlockMultiVector<SC, LO, GO, NO> multivector_t;
  typedef typename matrix_t::crs_graph_type graph_t;
  typedef typename matrix_t::scalar_type scalar_t;
  typedef typename matrix_t::local_ordinal_type lno_t;
  typedef typename matrix_t::global_ordinal_type gno_t;
  typedef typename matrix_t::node_type node_t;
  typedef typename node_t::device_type device_t;
  typedef typename device_t::execution_space execution_space;
  typedef Kokkos::View<int *, device_t> list_of_colors_t;
  typedef typename list_of_colors_t::HostMirror list_of_colors_host_t;

  // Constructor
  TpetraCrsColorer(const Teuchos::RCP<matrix_t> &matrix_);

  // Destructor
  ~TpetraCrsColorer() {}

  // Compute coloring data
  void
  computeColoring(Teuchos::ParameterList &coloring_params);

  // Compute seed matrix
  void
  computeSeedMatrix(multivector_t &V, const int color_beg = 0) const;

  // Compute seed matrix with distribution fitted to the graph's column map
  void
  computeSeedMatrixFitted(multivector_t &V, const int color_beg = 0) const;

  // Reconstruct matrix from supplied compressed vector
  void
  reconstructMatrix(multivector_t &W, const int color_beg = 0) const;
  void
  reconstructMatrix(multivector_t &W, matrix_t &mat, const int color_beg = 0) const;

  // Reconstruct matrix from supplied compressed vector fitted to the graph's
  // row map
  void
  reconstructMatrixFitted(multivector_t &W, const int color_beg = 0) const;
  void
  reconstructMatrixFitted(multivector_t &W, matrix_t &mat, const int color_beg = 0) const;

  // Return number of colors
  int
  getNumColors() const
  {
    return num_colors * matrix->getBlockSize();
  }

  // Get color given a column index
  int
  getColor(const size_t col) const
  {
    const int block_size   = matrix->getBlockSize();
    const size_t block_col = col / block_size;
    const int offset       = col - block_col * block_size;
    return (list_of_colors_host(block_col) - 1) * block_size + offset;
  }

  // Check coloring is valid, i.e., no row has two nonzero columns with
  // same color
  bool
  checkColoring() const
  {
    return Impl::check_coloring(*graph, list_of_colors);
  }

protected:
  Teuchos::RCP<matrix_t> matrix;
  Teuchos::RCP<const graph_t> graph;
  list_of_colors_t list_of_colors;
  list_of_colors_host_t list_of_colors_host;
  int num_colors;
};

//////////////////////////////////////////////////////////////////////////////

template <typename CrsMatrixType>
TpetraCrsColorer<CrsMatrixType>::TpetraCrsColorer(
  const Teuchos::RCP<matrix_t> &matrix_
)
  : matrix(matrix_),
    graph(matrix->getCrsGraph()),
    list_of_colors(),
    list_of_colors_host(),
    num_colors(0)
{

}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
void
TpetraCrsColorer<CrsMatrixType>::computeColoring(
  Teuchos::ParameterList &coloring_params
)
{
  const std::string library = coloring_params.get("library", "zoltan");

  if (library == "zoltan") {
    // Use Zoltan's coloring
    ZoltanCrsColorer<matrix_t> zz(matrix);
    zz.computeColoring(coloring_params,
                       num_colors, list_of_colors_host, list_of_colors);
  }
  else {
    // Use Zoltan2's coloring
    Zoltan2CrsColorer<matrix_t> zz2(matrix);
    zz2.computeColoring(coloring_params,
                        num_colors, list_of_colors_host, list_of_colors);
  }
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
template <typename MultiVectorType>
void
TpetraCrsColorer<CrsMatrixType>::computeSeedMatrix(
  MultiVectorType &V, const int color_beg) const
{
  MultiVectorType V_fitted(graph->getColMap(), V.getNumVectors());

  computeSeedMatrixFitted(V_fitted, color_beg);

  Tpetra::Import<lno_t, gno_t, node_t>
          importer(graph->getColMap(), V.getMap());

  V.doImport(V_fitted, importer, Tpetra::INSERT);
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
template <typename MultiVectorType>
void
TpetraCrsColorer<CrsMatrixType>::computeSeedMatrixFitted(
  MultiVectorType &V, const int color_beg) const
{
  // Check V's map is locally fitted to the graph's column map
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(graph->getColMap()->isLocallyFitted(*(V.getMap()))), std::domain_error,
      "Map for supplied vector is not locally fitted to the column map of the"
      " Jacobian graph.  "
      "You must call the non-fitted version of this function.");

  auto V_view_dev = V.getLocalViewDevice(Tpetra::Access::OverwriteAll);
  const size_t num_local_cols = graph->getLocalNumCols();
  const int chunk_size = V.getNumVectors();
  list_of_colors_t my_list_of_colors = list_of_colors;

  Kokkos::parallel_for(
      "TpetraCrsColorer::computeSeedMatrixFitted()",
      Kokkos::RangePolicy<execution_space>(0, num_local_cols),
      KOKKOS_LAMBDA(const size_t i) {
        const int color = my_list_of_colors[i] - 1 - color_beg;
        if (color >= 0 && color < chunk_size)
          V_view_dev(i, color) = scalar_t(1.0);
      });

}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
template <typename MultiVectorType>
void
TpetraCrsColorer<CrsMatrixType>::reconstructMatrix(
  MultiVectorType &W, const int color_beg) const
{
  reconstructMatrix(W, *matrix, color_beg);
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
template <typename MultiVectorType>
void
TpetraCrsColorer<CrsMatrixType>::reconstructMatrix(
  MultiVectorType &W,
  matrix_t &mat,
  const int color_beg) const
{
  MultiVectorType W_fitted(graph->getRowMap(), W.getNumVectors());

  Tpetra::Import<lno_t, gno_t, node_t>
          importer(W.getMap(), graph->getRowMap());

  W_fitted.doImport(W, importer, Tpetra::INSERT);

  reconstructMatrixFitted(W_fitted, mat, color_beg);
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
template <typename MultiVectorType>
void
TpetraCrsColorer<CrsMatrixType>::reconstructMatrixFitted(
  MultiVectorType &W, const int color_beg) const
{
  reconstructMatrixFitted(W, *matrix, color_beg);
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
template <typename MultiVectorType>
void
TpetraCrsColorer<CrsMatrixType>::reconstructMatrixFitted(
  MultiVectorType &W,
  matrix_t &mat,
  const int color_beg) const
{
  // Check the graph's row map is locally fitted to W's map
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(W.getMap()->isLocallyFitted(*(graph->getRowMap()))), std::domain_error,
      "Row map of the Jacobian graph is not locally fitted to the vector's map."
      "  You must call the non-fitted version of this function.");

  auto W_view_dev = W.getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto local_matrix = mat.getLocalMatrixDevice();
  auto local_graph = local_matrix.graph;
  const size_t num_local_rows = graph->getLocalNumRows();
  const int chunk_size = W.getNumVectors();
  list_of_colors_t my_list_of_colors = list_of_colors;

  Kokkos::parallel_for(
      "TpetraCrsColorer::reconstructMatrixFitted()",
      Kokkos::RangePolicy<execution_space>(0, num_local_rows),
      KOKKOS_LAMBDA(const size_t row) {
        const size_t entry_begin = local_graph.row_map(row);
        const size_t entry_end   = local_graph.row_map(row + 1);
        for (size_t entry = entry_begin; entry < entry_end; entry++)
        {
          const size_t col = local_graph.entries(entry);
          const int color  = my_list_of_colors[col] - 1 - color_beg;
          if (color >= 0 && color < chunk_size)
            local_matrix.values(entry) = W_view_dev(row, color);
        }
      });
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
TpetraCrsColorer<Tpetra::BlockCrsMatrix<SC,LO,GO,NO> >::TpetraCrsColorer(
  const Teuchos::RCP<matrix_t> &matrix_
)
    : matrix(matrix_)
    , graph(Teuchos::rcp(&(matrix->getCrsGraph()), false))
    , list_of_colors()
    , list_of_colors_host()
    , num_colors(0)
{
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
void
TpetraCrsColorer<Tpetra::BlockCrsMatrix<SC,LO,GO,NO> >::computeColoring(
  Teuchos::ParameterList &coloring_params
)
{
  const std::string library = coloring_params.get("library", "zoltan");

  if (library == "zoltan") {
    // Use Zoltan's coloring
    ZoltanCrsColorer<matrix_t> zz(matrix);
    zz.computeColoring(coloring_params,
                       num_colors, list_of_colors_host, list_of_colors);
  }
  else {
    // Use Zoltan2's coloring
    Zoltan2CrsColorer<matrix_t> zz2(matrix);
    zz2.computeColoring(coloring_params,
                        num_colors, list_of_colors_host, list_of_colors);
  }
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
void
TpetraCrsColorer<Tpetra::BlockCrsMatrix<SC,LO,GO,NO> >::computeSeedMatrix(
  multivector_t &block_V, const int color_beg) const
{
  multivector_t block_V_fitted(*(graph->getColMap()), matrix->getBlockSize(),
                                    block_V.getNumVectors());

  computeSeedMatrixFitted(block_V_fitted, color_beg);

  const lno_t block_size = block_V.getBlockSize();
  auto col_point_map = multivector_t::makePointMap(*graph->getColMap(),
                                                      block_size);
  auto blockV_point_map = multivector_t::makePointMap(*block_V.getMap(),
                                                         block_size);
  const auto col_point_map_rcp =
             Teuchos::rcp(new Tpetra::Map<LO, GO, NO>(col_point_map));
  const auto blockV_point_map_rcp =
             Teuchos::rcp(new Tpetra::Map<LO, GO, NO>(blockV_point_map));

  Tpetra::Import<LO, GO, NO> importer_point(col_point_map_rcp,
                                            blockV_point_map_rcp);

  block_V.getMultiVectorView().doImport(block_V_fitted.getMultiVectorView(),
                                        importer_point, Tpetra::INSERT);

  // Tpetra::Import<LO,GO,NO> importer(graph->getColMap(), block_V.getMap());
  // block_V.doImport(block_V_fitted, importer, Tpetra::INSERT);
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
void
TpetraCrsColorer<Tpetra::BlockCrsMatrix<SC,LO,GO,NO> >::computeSeedMatrixFitted(
  multivector_t &blockV, const int color_beg) const
{
  // Check blockV's map is locally fitted to the graph's column map
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(graph->getColMap()->isLocallyFitted(*(blockV.getMap()))),
       std::domain_error,
      "Map for supplied vector is not locally fitted to the column map of the"
      " Jacobian graph.  "
      "You must call the non-fitted version of this function.");

  const lno_t block_size = blockV.getBlockSize();
  auto V = blockV.getMultiVectorView();
  auto V_view_dev = V.getLocalViewDevice(Tpetra::Access::ReadWrite);
  const size_t num_local_cols = V_view_dev.extent(0) / block_size;
  const int chunk_size = V.getNumVectors() / block_size;
  const int block_color_beg = color_beg / block_size;
  list_of_colors_t my_list_of_colors = list_of_colors;

  Kokkos::parallel_for(
      "TpetraCrsColorer::computeSeedMatrixFitted()",
      Kokkos::RangePolicy<execution_space>(0, num_local_cols),
      KOKKOS_LAMBDA(const size_t i) {
        const int color = my_list_of_colors[i] - 1 - block_color_beg;
        if (color >= 0 && color < chunk_size)
          for (lno_t j = 0; j < block_size; ++j)
            V_view_dev(i*block_size+j, color*block_size+j) = scalar_t(1.0);
      });

}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
void
TpetraCrsColorer<Tpetra::BlockCrsMatrix<SC,LO,GO,NO> >::reconstructMatrix(
  multivector_t &W, const int color_beg) const
{
  reconstructMatrix(W, *matrix, color_beg);
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
void
TpetraCrsColorer<Tpetra::BlockCrsMatrix<SC,LO,GO,NO> >::reconstructMatrix(
  multivector_t &block_W,
  matrix_t &mat,
  const int color_beg) const
{
  multivector_t block_W_fitted(*(graph->getRowMap()), matrix->getBlockSize(),
                                  block_W.getNumVectors());

  const lno_t block_size = block_W.getBlockSize();
  auto row_point_map = multivector_t::makePointMap(*graph->getRowMap(),
                                                      block_size);
  auto blockW_point_map = multivector_t::makePointMap(*block_W.getMap(),
                                                         block_size);
  const auto row_point_map_rcp =
             Teuchos::rcp(new Tpetra::Map<LO, GO, NO>(row_point_map));
  const auto blockW_point_map_rcp =
             Teuchos::rcp(new Tpetra::Map<LO, GO, NO>(blockW_point_map));

  Tpetra::Import<lno_t, gno_t, node_t>
          importer_point(blockW_point_map_rcp, row_point_map_rcp);

  block_W_fitted.getMultiVectorView().doImport(block_W.getMultiVectorView(),
                                               importer_point, Tpetra::INSERT);
  reconstructMatrixFitted(block_W_fitted, mat, color_beg);
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
void
TpetraCrsColorer<Tpetra::BlockCrsMatrix<SC,LO,GO,NO> >::reconstructMatrixFitted(
  multivector_t &W, const int color_beg) const
{
  reconstructMatrixFitted(W, *matrix, color_beg);
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
void
TpetraCrsColorer<Tpetra::BlockCrsMatrix<SC,LO,GO,NO> >::reconstructMatrixFitted(
  multivector_t &block_W,
  matrix_t &mat,
  const int color_beg) const
{
  // Check the graph's row map is locally fitted to W's map
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(block_W.getMap()->isLocallyFitted(*(graph->getRowMap()))),
      std::domain_error,
      "Row map of the Jacobian graph is not locally fitted to the vector's map."
      "  You must call the non-fitted version of this function.");

  // Blocks are block_size x block_size stored with LayoutRight
  const lno_t block_size       = block_W.getBlockSize();
  const lno_t block_stride     = block_size * block_size;
  const lno_t block_row_stride = block_size;

  auto W = block_W.getMultiVectorView();
  auto W_view_dev                       = W.getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto matrix_vals                      = mat.getValuesDeviceNonConst();
  auto local_graph                      = graph->getLocalGraphDevice();
  const size_t num_local_rows           = graph->getLocalNumRows();
  const int chunk_size                  = W.getNumVectors() / block_size;
  const int block_color_beg             = color_beg / block_size;
  list_of_colors_t my_list_of_colors = list_of_colors;

  Kokkos::parallel_for(
      "TpetraCrsColorer::reconstructMatrix()",
      Kokkos::RangePolicy<execution_space>(0, num_local_rows),
      KOKKOS_LAMBDA(const size_t block_row) {
        const size_t entry_begin = local_graph.row_map(block_row);
        const size_t entry_end   = local_graph.row_map(block_row + 1);
        for (size_t block_entry = entry_begin; block_entry < entry_end;
                    block_entry++)
        {
          const size_t block_col    = local_graph.entries(block_entry);
          const size_t block_offset = block_stride * block_entry;
          const int block_color     = my_list_of_colors[block_col] - 1 - block_color_beg;
          if (block_color >= 0 && block_color < chunk_size)
          {
            for (lno_t i = 0; i < block_size; ++i)
            {
              const size_t row = block_row * block_size + i;
              for (lno_t j = 0; j < block_size; ++j)
              {
                const size_t entry = block_offset + block_row_stride * i + j;
                matrix_vals(entry) = W_view_dev(row, block_color*block_size+j);
              }
            }
          }
        }
      });
}
} // namespace Tpetra
