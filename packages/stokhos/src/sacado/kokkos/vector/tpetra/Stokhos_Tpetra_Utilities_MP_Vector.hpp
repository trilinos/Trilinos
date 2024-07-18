// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_TPETRA_UTILITIES_MP_VECTOR_HPP
#define STOKHOS_TPETRA_UTILITIES_MP_VECTOR_HPP

#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace Stokhos {

  // Create a flattened map for a map representing a distribution for an
  // embedded scalar type
  template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
  Teuchos::RCP< Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  create_flat_map(const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& map,
                  const LocalOrdinal block_size) {
    using Tpetra::global_size_t;
    using Teuchos::ArrayView;
    using Teuchos::Array;
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Map;

    // Get map info
    const global_size_t num_global_entries = map.getGlobalNumElements();
    const size_t num_local_entries = map.getLocalNumElements();
    const GlobalOrdinal index_base = map.getIndexBase();
    ArrayView<const GlobalOrdinal> element_list =
      map.getLocalElementList();

    // Create new elements
    const global_size_t flat_num_global_entries = num_global_entries*block_size;
    const size_t flat_num_local_entries = num_local_entries * block_size;
    const GlobalOrdinal flat_index_base = index_base;
    Array<GlobalOrdinal> flat_element_list(flat_num_local_entries);
    for (size_t i=0; i<num_local_entries; ++i)
      for (LocalOrdinal j=0; j<block_size; ++j)
        flat_element_list[i*block_size+j] = element_list[i]*block_size+j;

    // Create new map
    RCP<Map> flat_map =
      rcp(new Map(flat_num_global_entries, flat_element_list(),
                  flat_index_base, map.getComm()));

    return flat_map;
  }

  // Create a flattened graph for a graph from a matrix with the
  // MP::Vector scalar type (each block is an identity matrix)
  // If flat_domain_map and/or flat_range_map are null, they will be computed,
  // otherwise they will be used as-is.
  template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
  Teuchos::RCP< Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node > >
  create_flat_mp_graph(
    const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node >& graph,
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node > >& flat_domain_map,
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node > >& flat_range_map,
    const LocalOrdinal block_size) {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Map;
    typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Graph;
    typedef typename Graph::local_graph_device_type::row_map_type::non_const_type RowPtrs;
    typedef typename Graph::local_graph_device_type::entries_type::non_const_type LocalIndices;

    // Build domain map if necessary
    if (flat_domain_map == Teuchos::null)
      flat_domain_map = create_flat_map(*(graph.getDomainMap()), block_size);

    // Build range map if necessary
    if (flat_range_map == Teuchos::null)
      flat_range_map = create_flat_map(*(graph.getRangeMap()), block_size);

    // Build column map
    RCP<const Map> flat_col_map =
      create_flat_map(*(graph.getColMap()), block_size);

    // Build row map if necessary
    // Check if range_map == row_map, then we can use flat_range_map
    // as the flattened row map
    RCP<const Map> flat_row_map;
    if (graph.getRangeMap() == graph.getRowMap())
      flat_row_map = flat_range_map;
    else
      flat_row_map = create_flat_map(*(graph.getRowMap()), block_size);

    // Build flattened row offsets and column indices
    auto row_offsets = graph.getLocalRowPtrsHost();
    auto col_indices = graph.getLocalIndicesHost();
    const size_t num_row = graph.getLocalNumRows();
    const size_t num_col_indices = col_indices.size();
    RowPtrs flat_row_offsets("row_ptrs", num_row*block_size+1);
    LocalIndices flat_col_indices("col_indices", num_col_indices * block_size);
    auto flat_row_offsets_host = Kokkos::create_mirror_view(flat_row_offsets);
    auto flat_col_indices_host = Kokkos::create_mirror_view(flat_col_indices);
    for (size_t row=0; row<num_row; ++row) {
      const size_t row_beg = row_offsets[row];
      const size_t row_end = row_offsets[row+1];
      const size_t num_col = row_end - row_beg;
      for (LocalOrdinal j=0; j<block_size; ++j) {
        const size_t flat_row = row*block_size + j;
        const size_t flat_row_beg = row_beg*block_size + j*num_col;
        flat_row_offsets_host[flat_row] = flat_row_beg;
        for (size_t entry=0; entry<num_col; ++entry) {
          const LocalOrdinal col = col_indices[row_beg+entry];
          const LocalOrdinal flat_col = col*block_size + j;
          flat_col_indices_host[flat_row_beg+entry] = flat_col;
        }
      }
    }
    flat_row_offsets_host[num_row*block_size] = num_col_indices*block_size;
    Kokkos::deep_copy(flat_row_offsets, flat_row_offsets_host);
    Kokkos::deep_copy(flat_col_indices, flat_col_indices_host);

    // Build flattened graph
    RCP<Graph> flat_graph =
      rcp(new Graph(flat_row_map, flat_col_map,
                    flat_row_offsets, flat_col_indices));
    flat_graph->fillComplete(flat_domain_map, flat_range_map);

    return flat_graph;
  }

  // Create a flattened vector by unrolling the MP::Vector scalar type.  The
  // returned vector is a view of the original
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
  Teuchos::RCP<Tpetra::MultiVector<typename Storage::value_type, LocalOrdinal, GlobalOrdinal, Node>>
  create_flat_vector_view(
    const Tpetra::MultiVector<Sacado::MP::Vector<Storage>, LocalOrdinal, GlobalOrdinal, Node>& vec,
    const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& flat_map
  ) {
    using BaseScalar = typename Storage::value_type;
    using FlatVector = Tpetra::MultiVector<BaseScalar, LocalOrdinal, GlobalOrdinal, Node>;

    // Create flattenend view using reshaping conversion copy constructor
    typename FlatVector::wrapped_dual_view_type flat_vals(vec.getWrappedDualView());

    // Create flat vector
    return Teuchos::make_rcp<FlatVector>(flat_map, flat_vals);
  }

  // This version creates the map if necessary
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
  Teuchos::RCP<Tpetra::MultiVector<typename Storage::value_type, LocalOrdinal, GlobalOrdinal, Node>>
  create_flat_vector_view(
    const Tpetra::MultiVector<Sacado::MP::Vector<Storage>, LocalOrdinal, GlobalOrdinal, Node>& vec,
          Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& flat_map
  ) {
    if (flat_map == Teuchos::null) {
      const LocalOrdinal mp_size = Storage::static_size;
      flat_map = create_flat_map(*(vec.getMap()), mp_size);
    }
    const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> const_flat_map = flat_map;
    return create_flat_vector_view(vec, const_flat_map);
  }

  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
  Teuchos::RCP<Tpetra::Vector<typename Storage::value_type, LocalOrdinal, GlobalOrdinal, Node>>
  create_flat_vector_view(
    const Tpetra::Vector<Sacado::MP::Vector<Storage>, LocalOrdinal, GlobalOrdinal, Node>& vec,
    const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& flat_map
  ) {
    return create_flat_vector_view(
      static_cast<const Tpetra::MultiVector<Sacado::MP::Vector<Storage>, LocalOrdinal, GlobalOrdinal, Node>&>(vec),
      flat_map
    )->getVectorNonConst(0);
  }

  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
  Teuchos::RCP<Tpetra::Vector<typename Storage::value_type, LocalOrdinal, GlobalOrdinal, Node>>
  create_flat_vector_view(
    const Tpetra::Vector<Sacado::MP::Vector<Storage>, LocalOrdinal, GlobalOrdinal, Node>& vec,
          Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& flat_map
  ) {
    if (flat_map == Teuchos::null) {
      const LocalOrdinal mp_size = Storage::static_size;
      flat_map = create_flat_map(*(vec.getMap()), mp_size);
    }
    const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> const_flat_map = flat_map;
    return create_flat_vector_view(vec, const_flat_map);
  }

  // Create a flattened matrix by unrolling the MP::Vector scalar type.  The
  // returned matrix is NOT a view of the original (and can't be)
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Node>
  Teuchos::RCP< Tpetra::CrsMatrix<typename Storage::value_type,
                                  LocalOrdinal,GlobalOrdinal,Node> >
  create_flat_matrix(
    const Tpetra::CrsMatrix<Sacado::MP::Vector<Storage>,
                            LocalOrdinal,GlobalOrdinal,Node>& mat,
    const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& flat_graph,
    const LocalOrdinal block_size) {
    using Teuchos::ArrayView;
    using Teuchos::Array;
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Sacado::MP::Vector<Storage> Scalar;
    typedef typename Storage::value_type BaseScalar;
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Matrix;
    typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FlatMatrix;

    // Create flat matrix
    RCP<FlatMatrix> flat_mat = rcp(new FlatMatrix(flat_graph));

    // Set values
    const size_t num_rows = mat.getLocalNumRows();
    const size_t max_cols = mat.getLocalMaxNumRowEntries();
    typename Matrix::local_inds_host_view_type indices, flat_indices;
    typename Matrix::values_host_view_type values;
    Array<BaseScalar> flat_values(max_cols);
    for (size_t row=0; row<num_rows; ++row) {
      mat.getLocalRowView(row, indices, values);
      const size_t num_col = mat.getNumEntriesInLocalRow(row);
      for (LocalOrdinal i=0; i<block_size; ++i) {
        const LocalOrdinal flat_row = row*block_size + i;
        for (size_t j=0; j<num_col; ++j)
          flat_values[j] = values[j].coeff(i);
        flat_graph->getLocalRowView(flat_row, flat_indices);
        flat_mat->replaceLocalValues(flat_row, Kokkos::Compat::getConstArrayView(flat_indices),
                                     flat_values(0, num_col));
      }
    }
    flat_mat->fillComplete(flat_graph->getDomainMap(),
                           flat_graph->getRangeMap());

    return flat_mat;
  }

} // namespace Stokhos

#endif // STOKHOS_TPETRA_UTILITIES_MP_VECTOR_HPP
