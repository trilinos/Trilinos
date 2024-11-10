// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_TPETRA_UTILITIES_UQ_PCE_HPP
#define STOKHOS_TPETRA_UTILITIES_UQ_PCE_HPP

#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#include "Stokhos_Tpetra_Utilities_MP_Vector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace Stokhos {


  // Build a CRS graph from a sparse Cijk tensor
  template <typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            typename CijkType>
  Teuchos::RCP< Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node > >
  create_cijk_crs_graph(const CijkType& cijk_dev,
                        const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                        const size_t matrix_pce_size) {
    using Teuchos::RCP;
    using Teuchos::arrayView;

    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Map;
    typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Graph;

    // Code below accesses cijk entries on the host, so make sure it is
    // accessible there
    auto cijk = create_mirror_view(cijk_dev);
    deep_copy(cijk, cijk_dev);

    const size_t pce_sz = cijk.dimension();
    RCP<const Map> map =
      Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(pce_sz, comm);
    RCP<Graph> graph;
    if (matrix_pce_size == 1) {
      graph =  Tpetra::createCrsGraph(map, 1);
      // Mean-based case -- graph is diagonal
      for (size_t i=0; i<pce_sz; ++i) {
        const GlobalOrdinal row = i;
        graph->insertGlobalIndices(row, arrayView(&row, 1));
      }
    }
    else {
      // General case

      // Get max num entries
      size_t max_num_entry = 0;
      for (size_t i=0; i<pce_sz; ++i) {
        const size_t num_entry = cijk.num_entry(i);
        max_num_entry = (num_entry > max_num_entry) ? num_entry : max_num_entry;
      }
      max_num_entry *= 2; // 1 entry each for j, k coord
      graph =  Tpetra::createCrsGraph(map, max_num_entry);

      for (size_t i=0; i<pce_sz; ++i) {
        const GlobalOrdinal row = i;
        const size_t num_entry = cijk.num_entry(i);
        const size_t entry_beg = cijk.entry_begin(i);
        const size_t entry_end = entry_beg + num_entry;
        for (size_t entry = entry_beg; entry < entry_end; ++entry) {
          const GlobalOrdinal j = cijk.coord(entry,0);
          const GlobalOrdinal k = cijk.coord(entry,1);
          graph->insertGlobalIndices(row, arrayView(&j, 1));
          graph->insertGlobalIndices(row, arrayView(&k, 1));
        }
      }
    }
    graph->fillComplete();
    return graph;
  }

  // Create a flattened graph for a graph from a matrix with the
  // UQ::PCE scalar type
  // If flat_domain_map and/or flat_range_map are null, they will be computed,
  // otherwise they will be used as-is.
  template <typename LocalOrdinal, typename GlobalOrdinal, typename Node,
            typename CijkType>
  Teuchos::RCP< Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node > >
  create_flat_pce_graph(
    const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node >& graph,
    const CijkType& cijk,
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node > >& flat_domain_map,
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node > >& flat_range_map,
    Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node > >& cijk_graph,
    const size_t matrix_pce_size) {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::Array;
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Map;
    typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Graph;

    const LocalOrdinal block_size = cijk.dimension();

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

    // Build Cijk graph if necessary
    if (cijk_graph == Teuchos::null)
      cijk_graph = create_cijk_crs_graph<LocalOrdinal,GlobalOrdinal,Node>(
        cijk,
        flat_domain_map->getComm(),
        matrix_pce_size);

    // Build flattened graph that is the Kronecker product of the given
    // graph and cijk_graph

    // Loop over outer rows
    typename Graph::local_inds_host_view_type outer_cols;
    typename Graph::local_inds_host_view_type inner_cols;
    size_t max_num_row_entries = graph.getLocalMaxNumRowEntries()*block_size;
    Array<LocalOrdinal> flat_col_indices;
    flat_col_indices.reserve(max_num_row_entries);
    RCP<Graph> flat_graph = rcp(new Graph(flat_row_map, flat_col_map, max_num_row_entries));
    const LocalOrdinal num_outer_rows = graph.getLocalNumRows();
    for (LocalOrdinal outer_row=0; outer_row < num_outer_rows; outer_row++) {

      // Get outer columns for this outer row
      Kokkos::fence();
      graph.getLocalRowView(outer_row, outer_cols);
      const LocalOrdinal num_outer_cols = outer_cols.size();

      // Loop over inner rows
      for (LocalOrdinal inner_row=0; inner_row < block_size; inner_row++) {

        // Compute flat row index
        const LocalOrdinal flat_row = outer_row*block_size + inner_row;

        // Get inner columns for this inner row
        Kokkos::fence();
        cijk_graph->getLocalRowView(inner_row, inner_cols);
        const LocalOrdinal num_inner_cols = inner_cols.size();

        // Create space to store all column indices for this flat row
        flat_col_indices.resize(0);

        // Loop over outer cols
        for (LocalOrdinal outer_entry=0; outer_entry<num_outer_cols;
             ++outer_entry) {
          const LocalOrdinal outer_col = outer_cols[outer_entry];

          // Loop over inner cols
          for (LocalOrdinal inner_entry=0; inner_entry<num_inner_cols;
             ++inner_entry) {
            const LocalOrdinal inner_col = inner_cols[inner_entry];

            // Compute and store flat column index
            const LocalOrdinal flat_col = outer_col*block_size + inner_col;
            flat_col_indices.push_back(flat_col);
          }

        }

        // Insert all indices for this flat row
        flat_graph->insertLocalIndices(flat_row, flat_col_indices());

      }

    }
    flat_graph->fillComplete(flat_domain_map, flat_range_map);

    return flat_graph;
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Node>
  Teuchos::RCP< const Tpetra::MultiVector<typename Storage::value_type,
                                          LocalOrdinal,GlobalOrdinal,Node > >
  create_flat_vector_view(
    const Tpetra::MultiVector<Sacado::UQ::PCE<Storage>,
                              LocalOrdinal,GlobalOrdinal,Node >& vec,
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node > >& flat_map) {
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef typename Storage::value_type BaseScalar;
    typedef Tpetra::MultiVector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FlatVector;
    typedef typename FlatVector::dual_view_type::t_dev flat_view_type;

    // Have to do a nasty const-cast because getLocalViewDevice(ReadWrite) is a
    // non-const method, yet getLocalViewDevice(ReadOnly) returns a const-view
    // (i.e., with a constant scalar type), and there is no way to make a
    // MultiVector out of it!
    typedef Tpetra::MultiVector<Sacado::UQ::PCE<Storage>, LocalOrdinal,GlobalOrdinal, Node > mv_type;
    mv_type& vec_nc = const_cast<mv_type&>(vec);

    // Create flattenend view using special reshaping view assignment operator
    flat_view_type flat_vals = vec_nc.getLocalViewDevice(Tpetra::Access::ReadWrite);

    // Create flat vector
    RCP<FlatVector> flat_vec = rcp(new FlatVector(flat_map, flat_vals));

    return flat_vec;
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Node>
  Teuchos::RCP< Tpetra::MultiVector<typename Storage::value_type,
                                    LocalOrdinal,GlobalOrdinal,Node > >
  create_flat_vector_view(
    Tpetra::MultiVector<Sacado::UQ::PCE<Storage>,
                        LocalOrdinal,GlobalOrdinal,Node >& vec,
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node > >& flat_map) {
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef typename Storage::value_type BaseScalar;
    typedef Tpetra::MultiVector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FlatVector;
    typedef typename FlatVector::dual_view_type::t_dev flat_view_type;

    // Create flattenend view using special reshaping view assignment operator
    flat_view_type flat_vals = vec.getLocalViewDevice(Tpetra::Access::ReadWrite);

    // Create flat vector
    RCP<FlatVector> flat_vec = rcp(new FlatVector(flat_map, flat_vals));

    return flat_vec;
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original.  This version creates the
  // map if necessary
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Node>
  Teuchos::RCP< const Tpetra::MultiVector<typename Storage::value_type,
                                          LocalOrdinal,GlobalOrdinal,Node > >
  create_flat_vector_view(
    const Tpetra::MultiVector<Sacado::UQ::PCE<Storage>,
                              LocalOrdinal,GlobalOrdinal,Node >& vec,
    Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node > >& flat_map) {
    typedef typename Node::device_type Device;
    if (flat_map == Teuchos::null) {
      const LocalOrdinal pce_size =
        Kokkos::dimension_scalar(vec.template getLocalView<Device>(Tpetra::Access::ReadOnly));
      flat_map = create_flat_map(*(vec.getMap()), pce_size);
    }
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > const_flat_map = flat_map;
    return create_flat_vector_view(vec, const_flat_map);
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original.  This version creates the
  // map if necessary
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Node>
  Teuchos::RCP< Tpetra::MultiVector<typename Storage::value_type,
                                    LocalOrdinal,GlobalOrdinal,Node > >
  create_flat_vector_view(
    Tpetra::MultiVector<Sacado::UQ::PCE<Storage>,
                        LocalOrdinal,GlobalOrdinal,Node >& vec,
    Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node > >& flat_map) {
    typedef typename Node::device_type Device;
    if (flat_map == Teuchos::null) {
      const LocalOrdinal pce_size =
        Kokkos::dimension_scalar(vec.template getLocalView<Device>(Tpetra::Access::ReadOnly));
      flat_map = create_flat_map(*(vec.getMap()), pce_size);
    }
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > const_flat_map = flat_map;
    return create_flat_vector_view(vec, const_flat_map);
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Node>
  Teuchos::RCP< const Tpetra::Vector<typename Storage::value_type,
                                     LocalOrdinal,GlobalOrdinal,Node > >
  create_flat_vector_view(
    const Tpetra::Vector<Sacado::UQ::PCE<Storage>,
                         LocalOrdinal,GlobalOrdinal,Node >& vec_const,
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node > >& flat_map) {
    const Tpetra::MultiVector<Sacado::UQ::PCE<Storage>,LocalOrdinal,GlobalOrdinal,Node>& mv = vec_const;
    Teuchos::RCP< Tpetra::MultiVector<typename Storage::value_type,LocalOrdinal,GlobalOrdinal,Node> > flat_mv = create_flat_vector_view(mv, flat_map);
    return flat_mv->getVector(0);
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original.  This version creates the
  // map if necessary
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Node>
  Teuchos::RCP< const Tpetra::Vector<typename Storage::value_type,
                                     LocalOrdinal,GlobalOrdinal,Node > >
  create_flat_vector_view(
    const Tpetra::Vector<Sacado::UQ::PCE<Storage>,
                         LocalOrdinal,GlobalOrdinal,Node >& vec,
    Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node > >& flat_map) {
    typedef typename Node::device_type Device;
    if (flat_map == Teuchos::null) {
      const LocalOrdinal pce_size =
        Kokkos::dimension_scalar(vec.template getLocalView<Device>(Tpetra::Access::ReadOnly));
      flat_map = create_flat_map(*(vec.getMap()), pce_size);
    }
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > const_flat_map = flat_map;
    return create_flat_vector_view(vec, const_flat_map);
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Node>
  Teuchos::RCP< Tpetra::Vector<typename Storage::value_type,
                               LocalOrdinal,GlobalOrdinal,Node > >
  create_flat_vector_view(
    Tpetra::Vector<Sacado::UQ::PCE<Storage>,
                   LocalOrdinal,GlobalOrdinal,Node >& vec,
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node > >& flat_map) {
    Tpetra::MultiVector<Sacado::UQ::PCE<Storage>,LocalOrdinal,GlobalOrdinal,Node>& mv = vec;
    Teuchos::RCP< Tpetra::MultiVector<typename Storage::value_type,LocalOrdinal,GlobalOrdinal,Node> > flat_mv = create_flat_vector_view(mv, flat_map);
    return flat_mv->getVectorNonConst(0);
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original.  This version creates the
  // map if necessary
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Node>
  Teuchos::RCP< Tpetra::Vector<typename Storage::value_type,
                               LocalOrdinal,GlobalOrdinal,Node > >
  create_flat_vector_view(
    Tpetra::Vector<Sacado::UQ::PCE<Storage>,
                   LocalOrdinal,GlobalOrdinal,Node >& vec,
    Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node > >& flat_map) {
    typedef typename Node::device_type Device;
    if (flat_map == Teuchos::null) {
      const LocalOrdinal pce_size =
        Kokkos::dimension_scalar(vec.template getLocalView<Device>(Tpetra::Access::ReadOnly));
      flat_map = create_flat_map(*(vec.getMap()), pce_size);
    }
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > const_flat_map = flat_map;
    return create_flat_vector_view(vec, const_flat_map);
  }

  // Create a flattened matrix by unrolling the UQ::PCE scalar type.  The
  // returned matrix is NOT a view of the original (and can't be)
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Node, typename CijkType>
  Teuchos::RCP< Tpetra::CrsMatrix<typename Storage::value_type,
                                  LocalOrdinal,GlobalOrdinal,Node > >
  create_flat_matrix(
    const Tpetra::CrsMatrix<Sacado::UQ::PCE<Storage>,
                            LocalOrdinal,GlobalOrdinal,Node >& mat,
    const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node > >& flat_graph,
    const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node > >& cijk_graph,
    const CijkType& cijk_dev) {
    using Teuchos::ArrayView;
    using Teuchos::Array;
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Sacado::UQ::PCE<Storage> Scalar;
    typedef typename Storage::value_type BaseScalar;
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Matrix;
    typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FlatMatrix;

    // Code below accesses cijk entries on the host, so make sure it is
    // accessible there
    auto cijk = create_mirror_view(cijk_dev);
    deep_copy(cijk, cijk_dev);

    const LocalOrdinal block_size = cijk.dimension();
    const LocalOrdinal matrix_pce_size =
      Kokkos::dimension_scalar(mat.getLocalMatrixDevice().values);

    // Create flat matrix
    RCP<FlatMatrix> flat_mat = rcp(new FlatMatrix(flat_graph));

    // Fill flat matrix
    typename Matrix::values_host_view_type outer_values;
    typename Matrix::local_inds_host_view_type outer_cols;
    typename Matrix::local_inds_host_view_type inner_cols;
    typename Matrix::local_inds_host_view_type flat_cols;
    Array<BaseScalar> flat_values;
    flat_values.reserve(flat_graph->getLocalMaxNumRowEntries());
    const LocalOrdinal num_outer_rows = mat.getLocalNumRows();
    for (LocalOrdinal outer_row=0; outer_row < num_outer_rows; outer_row++) {

      // Get outer columns and values for this outer row
      mat.getLocalRowView(outer_row, outer_cols, outer_values);
      const LocalOrdinal num_outer_cols = outer_cols.size();

      // Loop over inner rows
      for (LocalOrdinal inner_row=0; inner_row < block_size; inner_row++) {

        // Compute flat row index
        const LocalOrdinal flat_row = outer_row*block_size + inner_row;

        // Get cijk column indices for this row
        cijk_graph->getLocalRowView(inner_row, inner_cols);
        const LocalOrdinal num_inner_cols = inner_cols.size();
        ArrayView<const LocalOrdinal> inner_cols_av =
          Kokkos::Compat::getConstArrayView(inner_cols);

        // Create space to store all values for this flat row
        const LocalOrdinal num_flat_indices = num_outer_cols*num_inner_cols;
        //flat_values.resize(num_flat_indices);
        flat_values.assign(num_flat_indices, BaseScalar(0));

        if (matrix_pce_size == 1) {
          // Mean-based case

          // Loop over outer cols
          for (LocalOrdinal outer_entry=0; outer_entry<num_outer_cols;
               ++outer_entry) {

            // Extract mean PCE entry for each outer column
            flat_values[outer_entry] = outer_values[outer_entry].coeff(0);

          }

        }
        else {

          // Loop over cijk non-zeros for this inner row
          const size_t num_entry = cijk.num_entry(inner_row);
          const size_t entry_beg = cijk.entry_begin(inner_row);
          const size_t entry_end = entry_beg + num_entry;
          for (size_t entry = entry_beg; entry < entry_end; ++entry) {
            const LocalOrdinal j = cijk.coord(entry,0);
            const LocalOrdinal k = cijk.coord(entry,1);
            const BaseScalar   c = cijk.value(entry);

            // Find column offset for j
            typedef typename ArrayView<const LocalOrdinal>::iterator iterator;
            iterator ptr_j =
              std::find(inner_cols_av.begin(), inner_cols_av.end(), j);
            iterator ptr_k =
              std::find(inner_cols_av.begin(), inner_cols_av.end(), k);
            const LocalOrdinal j_offset = ptr_j - inner_cols_av.begin();
            const LocalOrdinal k_offset = ptr_k - inner_cols_av.begin();

            // Loop over outer cols
            for (LocalOrdinal outer_entry=0; outer_entry<num_outer_cols;
                 ++outer_entry) {

              // Add contributions for each outer column
              flat_values[outer_entry*num_inner_cols + j_offset] +=
                c*outer_values[outer_entry].coeff(k);
              flat_values[outer_entry*num_inner_cols + k_offset] +=
                c*outer_values[outer_entry].coeff(j);

            }

          }

        }

        // Set values in flat matrix
        flat_graph->getLocalRowView(flat_row, flat_cols);
        flat_mat->replaceLocalValues(flat_row, Kokkos::Compat::getConstArrayView(flat_cols), flat_values());

      }

    }
    flat_mat->fillComplete(flat_graph->getDomainMap(),
                           flat_graph->getRangeMap());

    return flat_mat;
  }


} // namespace Stokhos

#endif // STOKHOS_TPETRA_UTILITIES_MP_VECTOR_HPP
