// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
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

#if defined(TPETRA_HAVE_KOKKOS_REFACTOR)

  // Build a CRS graph from a sparse Cijk tensor
  template <typename LocalOrdinal, typename GlobalOrdinal, typename Device,
            typename CijkType>
  Teuchos::RCP< Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,
                                 Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >
  create_cijk_crs_graph(const CijkType& cijk,
                        const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                        const Teuchos::RCP<Kokkos::Compat::KokkosDeviceWrapperNode<Device> >& node,
                        const size_t matrix_pce_size) {
    using Teuchos::RCP;
    using Teuchos::arrayView;

    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Device> Node;
    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Map;
    typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Graph;

    const size_t pce_sz = cijk.dimension();
    RCP<const Map> map =
      Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal>(pce_sz, comm, node);
    RCP<Graph> graph = Tpetra::createCrsGraph(map);
    if (matrix_pce_size == 1) {
      // Mean-based case -- graph is diagonal
      for (size_t i=0; i<pce_sz; ++i) {
        const GlobalOrdinal row = i;
        graph->insertGlobalIndices(row, arrayView(&row, 1));
      }
    }
    else {
      // General case
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
  template <typename LocalOrdinal, typename GlobalOrdinal, typename Device,
            typename CijkType>
  Teuchos::RCP< Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,
                                 Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >
  create_flat_pce_graph(
    const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<Device> >& graph,
    const CijkType& cijk,
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >& flat_domain_map,
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >& flat_range_map,
    Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >& cijk_graph,
    const size_t matrix_pce_size) {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::Array;
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Device> Node;
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
      cijk_graph = create_cijk_crs_graph<LocalOrdinal,GlobalOrdinal>(
        cijk,
        flat_domain_map->getComm(),
        flat_domain_map->getNode(),
        matrix_pce_size);

    // Build flattened graph that is the Kronecker product of the given
    // graph and cijk_graph
    RCP<Graph> flat_graph = rcp(new Graph(flat_row_map, flat_col_map, 0));

    // Loop over outer rows
    ArrayView<const LocalOrdinal> outer_cols;
    ArrayView<const LocalOrdinal> inner_cols;
    Array<LocalOrdinal> flat_col_indices;
    flat_col_indices.reserve(graph.getNodeMaxNumRowEntries()*block_size);
    const LocalOrdinal num_outer_rows = graph.getNodeNumRows();
    for (LocalOrdinal outer_row=0; outer_row < num_outer_rows; outer_row++) {

      // Get outer columns for this outer row
      graph.getLocalRowView(outer_row, outer_cols);
      const LocalOrdinal num_outer_cols = outer_cols.size();

      // Loop over inner rows
      for (LocalOrdinal inner_row=0; inner_row < block_size; inner_row++) {

        // Compute flat row index
        const LocalOrdinal flat_row = outer_row*block_size + inner_row;

        // Get inner columns for this inner row
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
            typename Device>
  Teuchos::RCP< const Tpetra::MultiVector<typename Storage::value_type,
                                          LocalOrdinal,GlobalOrdinal,
                                          Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >
  create_flat_vector_view(
    const Tpetra::MultiVector<Sacado::UQ::PCE<Storage>,
                              LocalOrdinal,GlobalOrdinal,
                              Kokkos::Compat::KokkosDeviceWrapperNode<Device> >& vec,
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,
                                          Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >& flat_map) {
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef typename Storage::value_type BaseScalar;
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Device> Node;
    typedef Tpetra::MultiVector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FlatVector;
    typedef typename FlatVector::dual_view_type flat_view_type;

    // Create flattenend view using special reshaping view assignment operator
    flat_view_type flat_vals = vec.getDualView();

    // Create flat vector
    RCP<FlatVector> flat_vec = rcp(new FlatVector(flat_map, flat_vals));

    return flat_vec;
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Device>
  Teuchos::RCP< Tpetra::MultiVector<typename Storage::value_type,
                                    LocalOrdinal,GlobalOrdinal,
                                    Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >
  create_flat_vector_view(
    Tpetra::MultiVector<Sacado::UQ::PCE<Storage>,
                        LocalOrdinal,GlobalOrdinal,
                        Kokkos::Compat::KokkosDeviceWrapperNode<Device> >& vec,
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,
                                          Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >& flat_map) {
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef typename Storage::value_type BaseScalar;
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Device> Node;
    typedef Tpetra::MultiVector<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FlatVector;
    typedef typename FlatVector::dual_view_type flat_view_type;

    // Create flattenend view using special reshaping view assignment operator
    flat_view_type flat_vals = vec.getDualView();

    // Create flat vector
    RCP<FlatVector> flat_vec = rcp(new FlatVector(flat_map, flat_vals));

    return flat_vec;
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original.  This version creates the
  // map if necessary
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Device>
  Teuchos::RCP< const Tpetra::MultiVector<typename Storage::value_type,
                                          LocalOrdinal,GlobalOrdinal,
                                          Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >
  create_flat_vector_view(
    const Tpetra::MultiVector<Sacado::UQ::PCE<Storage>,
                              LocalOrdinal,GlobalOrdinal,
                              Kokkos::Compat::KokkosDeviceWrapperNode<Device> >& vec,
    Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,
                                    Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >& flat_map) {
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Device> Node;
    if (flat_map == Teuchos::null) {
      const LocalOrdinal pce_size =
        vec.template getLocalView<Device>().sacado_size();
      flat_map = create_flat_map(*(vec.getMap()), pce_size);
    }
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > const_flat_map = flat_map;
    return create_flat_vector_view(vec, const_flat_map);
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original.  This version creates the
  // map if necessary
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Device>
  Teuchos::RCP< Tpetra::MultiVector<typename Storage::value_type,
                                    LocalOrdinal,GlobalOrdinal,
                                    Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >
  create_flat_vector_view(
    Tpetra::MultiVector<Sacado::UQ::PCE<Storage>,
                        LocalOrdinal,GlobalOrdinal,
                        Kokkos::Compat::KokkosDeviceWrapperNode<Device> >& vec,
    Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,
                                    Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >& flat_map) {
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Device> Node;
    if (flat_map == Teuchos::null) {
      const LocalOrdinal pce_size =
        vec.template getLocalView<Device>().sacado_size();
      flat_map = create_flat_map(*(vec.getMap()), pce_size);
    }
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > const_flat_map = flat_map;
    return create_flat_vector_view(vec, const_flat_map);
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Device>
  Teuchos::RCP< const Tpetra::Vector<typename Storage::value_type,
                                     LocalOrdinal,GlobalOrdinal,
                                     Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >
  create_flat_vector_view(
    const Tpetra::Vector<Sacado::UQ::PCE<Storage>,
                         LocalOrdinal,GlobalOrdinal,
                         Kokkos::Compat::KokkosDeviceWrapperNode<Device> >& vec_const,
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,
                                          Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >& flat_map) {
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Device> Node;
    const Tpetra::MultiVector<Sacado::UQ::PCE<Storage>,LocalOrdinal,GlobalOrdinal,Node>& mv = vec_const;
    Teuchos::RCP< Tpetra::MultiVector<typename Storage::value_type,LocalOrdinal,GlobalOrdinal,Node> > flat_mv = create_flat_vector_view(mv, flat_map);
    return flat_mv->getVector(0);
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original.  This version creates the
  // map if necessary
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Device>
  Teuchos::RCP< const Tpetra::Vector<typename Storage::value_type,
                                     LocalOrdinal,GlobalOrdinal,
                                     Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >
  create_flat_vector_view(
    const Tpetra::Vector<Sacado::UQ::PCE<Storage>,
                         LocalOrdinal,GlobalOrdinal,
                         Kokkos::Compat::KokkosDeviceWrapperNode<Device> >& vec,
    Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,
                                    Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >& flat_map) {
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Device> Node;
    if (flat_map == Teuchos::null) {
      const LocalOrdinal pce_size =
        vec.template getLocalView<Device>().sacado_size();
      flat_map = create_flat_map(*(vec.getMap()), pce_size);
    }
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > const_flat_map = flat_map;
    return create_flat_vector_view(vec, const_flat_map);
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Device>
  Teuchos::RCP< Tpetra::Vector<typename Storage::value_type,
                               LocalOrdinal,GlobalOrdinal,
                               Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >
  create_flat_vector_view(
    Tpetra::Vector<Sacado::UQ::PCE<Storage>,
                   LocalOrdinal,GlobalOrdinal,
                   Kokkos::Compat::KokkosDeviceWrapperNode<Device> >& vec,
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,
                                          Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >& flat_map) {
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Device> Node;
    Tpetra::MultiVector<Sacado::UQ::PCE<Storage>,LocalOrdinal,GlobalOrdinal,Node>& mv = vec;
    Teuchos::RCP< Tpetra::MultiVector<typename Storage::value_type,LocalOrdinal,GlobalOrdinal,Node> > flat_mv = create_flat_vector_view(mv, flat_map);
    return flat_mv->getVectorNonConst(0);
  }

  // Create a flattened vector by unrolling the UQ::PCE scalar type.  The
  // returned vector is a view of the original.  This version creates the
  // map if necessary
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Device>
  Teuchos::RCP< Tpetra::Vector<typename Storage::value_type,
                               LocalOrdinal,GlobalOrdinal,
                               Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >
  create_flat_vector_view(
    Tpetra::Vector<Sacado::UQ::PCE<Storage>,
                   LocalOrdinal,GlobalOrdinal,
                   Kokkos::Compat::KokkosDeviceWrapperNode<Device> >& vec,
    Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,
                                    Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >& flat_map) {
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Device> Node;
    if (flat_map == Teuchos::null) {
      const LocalOrdinal pce_size =
        vec.template getLocalView<Device>().sacado_size();
      flat_map = create_flat_map(*(vec.getMap()), pce_size);
    }
    const Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > const_flat_map = flat_map;
    return create_flat_vector_view(vec, const_flat_map);
  }

  // Create a flattened matrix by unrolling the UQ::PCE scalar type.  The
  // returned matrix is NOT a view of the original (and can't be)
  template <typename Storage, typename LocalOrdinal, typename GlobalOrdinal,
            typename Device, typename CijkType>
  Teuchos::RCP< Tpetra::CrsMatrix<typename Storage::value_type,
                                  LocalOrdinal,GlobalOrdinal,
                                  Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >
  create_flat_matrix(
    const Tpetra::CrsMatrix<Sacado::UQ::PCE<Storage>,
                            LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<Device> >& mat,
    const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >& flat_graph,
    const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<Device> > >& cijk_graph,
    const CijkType& cijk) {
    using Teuchos::ArrayView;
    using Teuchos::Array;
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Kokkos::Compat::KokkosDeviceWrapperNode<Device> Node;
    typedef Sacado::UQ::PCE<Storage> Scalar;
    typedef typename Storage::value_type BaseScalar;
    typedef Tpetra::CrsMatrix<BaseScalar,LocalOrdinal,GlobalOrdinal,Node> FlatMatrix;

    const LocalOrdinal block_size = cijk.dimension();
    const LocalOrdinal matrix_pce_size =
      mat.getLocalMatrix().values.sacado_size();

    // Create flat matrix
    RCP<FlatMatrix> flat_mat = rcp(new FlatMatrix(flat_graph));

    // Fill flat matrix
    ArrayView<const Scalar> outer_values;
    ArrayView<const LocalOrdinal> outer_cols;
    ArrayView<const LocalOrdinal> inner_cols;
    ArrayView<const LocalOrdinal> flat_cols;
    Array<BaseScalar> flat_values;
    flat_values.reserve(flat_graph->getNodeMaxNumRowEntries());
    const LocalOrdinal num_outer_rows = mat.getNodeNumRows();
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
              std::find(inner_cols.begin(), inner_cols.end(), j);
            iterator ptr_k =
              std::find(inner_cols.begin(), inner_cols.end(), k);
            const LocalOrdinal j_offset = ptr_j - inner_cols.begin();
            const LocalOrdinal k_offset = ptr_k - inner_cols.begin();

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
        flat_mat->replaceLocalValues(flat_row, flat_cols, flat_values());

      }

    }
    flat_mat->fillComplete(flat_graph->getDomainMap(),
                           flat_graph->getRangeMap());

    return flat_mat;
  }

#endif

} // namespace Stokhos

#endif // STOKHOS_TPETRA_UTILITIES_MP_VECTOR_HPP
