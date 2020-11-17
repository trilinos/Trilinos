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

#include "Zoltan2_CrsColorerUtils.hpp"

namespace Zoltan2
{

// Base class for coloring Tpetra::CrsMatrix for use in column compression
template <typename CrsMatrixType>
class CrsColorer
{
public:
  typedef CrsMatrixType matrix_type;
  typedef typename matrix_type::crs_graph_type graph_type;
  typedef typename matrix_type::scalar_type scalar_type;
  typedef typename matrix_type::local_ordinal_type local_ordinal_type;
  typedef typename matrix_type::global_ordinal_type global_ordinal_type;
  typedef typename matrix_type::node_type node_type;
  typedef typename node_type::device_type device_type;
  typedef typename device_type::execution_space execution_space;
  typedef Kokkos::View<int *, device_type> list_of_colors_type;
  typedef typename list_of_colors_type::HostMirror list_of_colors_host_type;

  // Constructor
  CrsColorer(const Teuchos::RCP<matrix_type> &matrix_);

  // Destructor
  virtual ~CrsColorer() {}

  // Compute coloring data
  virtual void
  computeColoring(Teuchos::ParameterList &coloring_params) = 0;

  // Compute seed matrix
  template <typename MultiVectorType>
  void
  computeSeedMatrix(MultiVectorType &V) const;

  // Compute seed matrix with distribution fitted to the graph's column map
  template <typename MultiVectorType>
  void
  computeSeedMatrixFitted(MultiVectorType &V) const;

  // Reconstruct matrix from supplied compressed vector
  template <typename MultiVectorType>
  void
  reconstructMatrix(MultiVectorType &W) const;
  template <typename MultiVectorType>
  void
  reconstructMatrix(MultiVectorType &W, matrix_type &mat) const;

  // Reconstruct matrix from supplied compressed vector fitted to the graph's
  // row map
  template <typename MultiVectorType>
  void
  reconstructMatrixFitted(MultiVectorType &W) const;
  template <typename MultiVectorType>
  void
  reconstructMatrixFitted(MultiVectorType &W, matrix_type &mat) const;

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
  Teuchos::RCP<matrix_type> matrix;
  Teuchos::RCP<const graph_type> graph;
  list_of_colors_type list_of_colors;
  list_of_colors_host_type list_of_colors_host;
  int num_colors;
};

//////////////////////////////////////////////////////////////////////////////
// Specialization of Tpetra::CrsColorer for BlockCrsMatrix
template <typename SC, typename LO, typename GO, typename NO>
class CrsColorer<BlockCrsMatrix<SC, LO, GO, NO>>
{
public:
  typedef BlockCrsMatrix<SC, LO, GO, NO> matrix_type;
  typedef BlockMultiVector<SC, LO, GO, NO> multivector_type;
  typedef typename matrix_type::crs_graph_type graph_type;
  typedef typename matrix_type::scalar_type scalar_type;
  typedef typename matrix_type::local_ordinal_type local_ordinal_type;
  typedef typename matrix_type::global_ordinal_type global_ordinal_type;
  typedef typename matrix_type::node_type node_type;
  typedef typename node_type::device_type device_type;
  typedef typename device_type::execution_space execution_space;
  typedef Kokkos::View<int *, device_type> list_of_colors_type;
  typedef typename list_of_colors_type::HostMirror list_of_colors_host_type;

  // Constructor
  CrsColorer(const Teuchos::RCP<matrix_type> &matrix_);

  // Destructor
  virtual ~CrsColorer() {}

  // Compute coloring data
  virtual void
  computeColoring(Teuchos::ParameterList &coloring_params) = 0;

  // Compute seed matrix
  void
  computeSeedMatrix(multivector_type &V) const;

  // Compute seed matrix with distribution fitted to the graph's column map
  void
  computeSeedMatrixFitted(multivector_type &V) const;

  // Reconstruct matrix from supplied compressed vector
  void
  reconstructMatrix(multivector_type &W) const;
  void
  reconstructMatrix(multivector_type &W, matrix_type &mat) const;

  // Reconstruct matrix from supplied compressed vector fitted to the graph's
  // row map
  void
  reconstructMatrixFitted(multivector_type &W) const;
  void
  reconstructMatrixFitted(multivector_type &W, matrix_type &mat) const;

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
  Teuchos::RCP<matrix_type> matrix;
  Teuchos::RCP<const graph_type> graph;
  list_of_colors_type list_of_colors;
  list_of_colors_host_type list_of_colors_host;
  int num_colors;
};

//////////////////////////////////////////////////////////////////////////////

template <typename CrsMatrixType>
CrsColorer<CrsMatrixType>::CrsColorer(const Teuchos::RCP<matrix_type> &matrix_)
  : matrix(matrix_), 
    graph(matrix->getCrsGraph()), 
    list_of_colors(), 
    list_of_colors_host(), 
    num_colors(0)
{

}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
CrsColorer<CrsMatrixType>::computeColoring(
  Teuchos::ParameterList &coloring_params
)
{
  const std::string library = coloring_params.get("library", "zoltan");
  const std::string matrixType = coloring_params.get("matrixType", "Jacobian");

  // TODO:  Check the logic here

  if (matrixType == "Jacobian") {
    if (library == "zoltan") {
      // Use Zoltan's partial distance 2 coloring
      ZoltanCrsColorer zz(J);
      zz.computeColoring(coloring_params, "PARTIAL-DISTANCE-2", 
                         num_colors, list_of_colors_host, list_of_colors);
    }
    else {
      // Use Zoltan2's partial distance 2 coloring when it is ready
      Zoltan2CrsColorer zz2(J);
      zz2.computeColoring(coloring_params, "PARTIAL-DISTANCE-2", 
                          num_colors, list_of_colors_host, list_of_colors);
    }
  else if (matrixType == "Hessian") {
    // TODO: Figure out whether this is the right thing to do and whether
    // TODO: the code supports it
    // Hessian is already symmetric; 
    // code currently still uses partial distance 2.
    // Should use Distance2 instead
    if (library == "zoltan") {
      ZoltanCrsColorer zz(J);
      zz.computeColoring(coloring_params, "DISTANCE-2", 
                         num_colors, list_of_colors_host, list_of_colors);
    }
    else {
      Zoltan2CrsColorer zz2(J);
      zz2.computeColoring(coloring_params, "DISTANCE-2", 
                          num_colors, list_of_colors_host, list_of_colors);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
template <typename MultiVectorType>
void
CrsColorer<CrsMatrixType>::computeSeedMatrix(MultiVectorType &V) const
{
  MultiVectorType V_fitted(graph->getColMap(), num_colors);

  computeSeedMatrixFitted(V_fitted);

  Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> 
          importer(graph->getColMap(), V.getMap());

  V.doImport(V_fitted, importer, Tpetra::INSERT);
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
template <typename MultiVectorType>
void
CrsColorer<CrsMatrixType>::computeSeedMatrixFitted(MultiVectorType &V) const
{
  // Check V's map is locally fitted to the graph's column map
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(graph->getColMap()->isLocallyFitted(*(V.getMap()))), std::domain_error,
      "Map for supplied vector is not locally fitted to the column map of the"
      " Jacobian graph.  "
      "You must call the non-fitted version of this function.");

  V.sync_device();

  auto V_view_dev = V.getLocalViewDevice();
  const size_t num_local_cols = graph->getNodeNumCols();
  list_of_colors_type my_list_of_colors = list_of_colors;

  Kokkos::parallel_for(
      Kokkos::RangePolicy<execution_space>(0, num_local_cols),
      KOKKOS_LAMBDA(const size_t i) { 
        V_view_dev(i, my_list_of_colors[i] - 1) = scalar_type(1.0); },
        "CrsColorer::computeSeedMatrixFitted()");

  V.modify_device();
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
template <typename MultiVectorType>
void
CrsColorer<CrsMatrixType>::reconstructMatrix(MultiVectorType &W) const
{
  reconstructMatrix(W, *matrix);
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
template <typename MultiVectorType>
void
CrsColorer<CrsMatrixType>::reconstructMatrix(
  MultiVectorType &W, 
  matrix_type &mat) const
{
  MultiVectorType W_fitted(graph->getRowMap(), num_colors);

  Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> 
          importer(W.getMap(), graph->getRowMap());

  W_fitted.doImport(W, importer, Tpetra::INSERT);

  reconstructMatrixFitted(W_fitted, mat);
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
template <typename MultiVectorType>
void
CrsColorer<CrsMatrixType>::reconstructMatrixFitted(MultiVectorType &W) const
{
  reconstructMatrixFitted(W, *matrix);
}

//////////////////////////////////////////////////////////////////////////////
template <typename CrsMatrixType>
template <typename MultiVectorType>
void
CrsColorer<CrsMatrixType>::reconstructMatrixFitted(
  MultiVectorType &W, 
  matrix_type &mat) const
{
  // Check the graph's row map is locally fitted to W's map
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(W.getMap()->isLocallyFitted(*(graph->getRowMap()))), std::domain_error,
      "Row map of the Jacobian graph is not locally fitted to the vector's map."
      "  You must call the non-fitted version of this function.");

  W.sync_device();

  auto W_view_dev = W.getLocalViewDevice();
  auto local_matrix = mat.getLocalMatrix();
  auto local_graph = graph->getLocalGraph();
  const size_t num_local_rows = graph->getNodeNumRows();
  list_of_colors_type my_list_of_colors = list_of_colors;

  Kokkos::parallel_for(
      Kokkos::RangePolicy<execution_space>(0, num_local_rows),
      KOKKOS_LAMBDA(const size_t row) {
        const size_t entry_begin = local_graph.row_map(row);
        const size_t entry_end   = local_graph.row_map(row + 1);
        for (size_t entry = entry_begin; entry < entry_end; entry++)
        {
          const size_t col           = local_graph.entries(entry);
          local_matrix.values(entry) = W_view_dev(row,my_list_of_colors[col]-1);
        }
      },
      "CrsColorer::reconstructMatrixFitted()");
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
CrsColorer<BlockCrsMatrix<SC, LO, GO, NO>>::CrsColorer(
  const Teuchos::RCP<matrix_type> &matrix_
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
CrsColorer<BlockCrsMatrix<SC, LO, GO, NO>>::computeSeedMatrix(
  multivector_type &block_V) const
{
  multivector_type block_V_fitted(*(graph->getColMap()), matrix->getBlockSize(),
                                  num_colors * matrix->getBlockSize());

  computeSeedMatrixFitted(block_V_fitted);

  const local_ordinal_type block_size = block_V.getBlockSize();
  auto col_point_map = multivector_type::makePointMap(*graph->getColMap(),
                                                      block_size);
  auto blockV_point_map = multivector_type::makePointMap(*block_V.getMap(),
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
CrsColorer<BlockCrsMatrix<SC, LO, GO, NO> >::computeSeedMatrixFitted(
  multivector_type &blockV) const
{
  // Check blockV's map is locally fitted to the graph's column map
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(graph->getColMap()->isLocallyFitted(*(blockV.getMap()))),
       std::domain_error,
      "Map for supplied vector is not locally fitted to the column map of the"
      " Jacobian graph.  "
      "You must call the non-fitted version of this function.");

  const local_ordinal_type block_size = blockV.getBlockSize();
  auto V = blockV.getMultiVectorView();

  V.sync_device();
  auto V_view_dev = V.getLocalViewDevice();
  const size_t num_local_cols = V_view_dev.extent(0) / block_size;
  list_of_colors_type my_list_of_colors = list_of_colors;

  Kokkos::parallel_for(
      Kokkos::RangePolicy<execution_space>(0, num_local_cols),
      KOKKOS_LAMBDA(const size_t i) {
        for (local_ordinal_type j = 0; j < block_size; ++j)
          V_view_dev(i*block_size+j, (my_list_of_colors[i]-1)*block_size+j) = 
                                     scalar_type(1.0);
      },
      "CrsColorer::computeSeedMatrixOverlapped()");

  V.modify_device();
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
void
CrsColorer<BlockCrsMatrix<SC, LO, GO, NO>>::reconstructMatrix(
  multivector_type &W) const
{
  reconstructMatrix(W, *matrix);
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
void
CrsColorer<BlockCrsMatrix<SC, LO, GO, NO> >::reconstructMatrix(
  multivector_type &block_W, 
  matrix_type &mat) const
{
  multivector_type block_W_fitted(*(graph->getRowMap()), matrix->getBlockSize(),
                                  num_colors * matrix->getBlockSize());

  const local_ordinal_type block_size = block_W.getBlockSize();
  auto row_point_map = multivector_type::makePointMap(*graph->getRowMap(),
                                                      block_size);
  auto blockW_point_map = multivector_type::makePointMap(*block_W.getMap(),
                                                         block_size);
  const auto row_point_map_rcp = 
             Teuchos::rcp(new Tpetra::Map<LO, GO, NO>(row_point_map));
  const auto blockW_point_map_rcp = 
             Teuchos::rcp(new Tpetra::Map<LO, GO, NO>(blockW_point_map));

  Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> 
          importer_point(blockW_point_map_rcp, row_point_map_rcp);

  block_W_fitted.getMultiVectorView().doImport(block_W.getMultiVectorView(),
                                               importer_point, Tpetra::INSERT);
  reconstructMatrixFitted(block_W_fitted, mat);
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
void
CrsColorer<BlockCrsMatrix<SC, LO, GO, NO> >::reconstructMatrixFitted(
  multivector_type &W) const
{
  reconstructMatrixFitted(W, *matrix);
}

//////////////////////////////////////////////////////////////////////////////
template <typename SC, typename LO, typename GO, typename NO>
void
CrsColorer<BlockCrsMatrix<SC, LO, GO, NO>>::reconstructMatrixFitted(
  multivector_type &block_W,
  matrix_type &mat) const
{
  // Check the graph's row map is locally fitted to W's map
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(block_W.getMap()->isLocallyFitted(*(graph->getRowMap()))),
      std::domain_error,
      "Row map of the Jacobian graph is not locally fitted to the vector's map."
      "  You must call the non-fitted version of this function.");

  // Blocks are block_size x block_size stored with LayoutRight
  const local_ordinal_type block_size       = block_W.getBlockSize();
  const local_ordinal_type block_stride     = block_size * block_size;
  const local_ordinal_type block_row_stride = block_size;

  auto W = block_W.getMultiVectorView();
  W.sync_device();
  auto W_view_dev                       = W.getLocalViewDevice();
  auto matrix_vals                      = matrix->getValuesDevice();
  auto local_graph                      = graph->getLocalGraph();
  const size_t num_local_rows           = graph->getNodeNumRows();
  list_of_colors_type my_list_of_colors = list_of_colors;

  Kokkos::parallel_for(
      Kokkos::RangePolicy<execution_space>(0, num_local_rows),
      KOKKOS_LAMBDA(const size_t block_row) {
        const size_t entry_begin = local_graph.row_map(block_row);
        const size_t entry_end   = local_graph.row_map(block_row + 1);
        for (size_t block_entry = entry_begin; block_entry < entry_end;
                    block_entry++)
        {
          const size_t block_col    = local_graph.entries(block_entry);
          const size_t block_offset = block_stride * block_entry;
          const int block_color     = my_list_of_colors[block_col] - 1;
          for (local_ordinal_type i = 0; i < block_size; ++i)
          {
            const size_t row = block_row * block_size + i;
            for (local_ordinal_type j = 0; j < block_size; ++j)
            {
              const size_t entry = block_offset + block_row_stride * i + j;
              matrix_vals(entry) = W_view_dev(row, block_color*block_size+j);
            }
          }
        }
      },
      "CrsColorer::reconstructMatrix()");
}
} // namespace Tpetra
