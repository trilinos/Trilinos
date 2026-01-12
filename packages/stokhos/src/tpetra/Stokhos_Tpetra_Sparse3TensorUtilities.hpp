// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_TPETRA_SPARSE_3_TENSOR_UTILITIES_HPP
#define STOKHOS_TPETRA_SPARSE_3_TENSOR_UTILITIES_HPP

#include "Stokhos_Sparse3Tensor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "MatrixMarket_Tpetra.hpp"

namespace Stokhos {

  //! Build an Tpetra::CrsGraph from a sparse 3 tensor
  /*!
   * Builds a sparse graph from a sparse 3 tensor by summing over the third
   * index.  This graph then represents the sparsity pattern of the stochastic
   * part of the block stochastic Galerkin operator.  Redistributing the graph
   * should then provide a suitable parallel distribution for block
   * stochastic Galerkin linear solves.
   */
  template <typename TpetraLocalOrdinal, typename TpetraGlobalOrdinal, typename TpetraNode,
            typename StokhosOrdinal, typename StokhosValue>
  Teuchos::RCP< Tpetra::CrsGraph<TpetraLocalOrdinal,TpetraGlobalOrdinal,TpetraNode > >
  sparse3Tensor2TpetraCrsGraph(
    const Stokhos::OrthogPolyBasis<StokhosOrdinal,StokhosValue>& basis,
    const Stokhos::Sparse3Tensor<StokhosOrdinal,StokhosValue>& Cijk,
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm) 
  {
    typedef Stokhos::Sparse3Tensor<StokhosOrdinal,StokhosValue> Cijk_type;
    typedef Tpetra::Map<TpetraLocalOrdinal,TpetraGlobalOrdinal,TpetraNode> Map;
    typedef Tpetra::CrsGraph<TpetraLocalOrdinal,TpetraGlobalOrdinal,TpetraNode> Graph;
    using Teuchos::RCP;
    using Teuchos::arrayView;

    // Number of stochastic rows
    StokhosOrdinal num_rows = basis.size();

    // Compute max num entries -- note this is a substantial overestimate in most cases because it doesn't
    // take into account overlap in j entries for each k. Also, while Cijk is symmetric with respect to i and
    // j, it is not necessarily symmetric with respect to k.  So this is the best we can do without adding
    // and i -> j -> k traversal to the Cijk tensor.
    size_t max_num_entry = 0;
    for (typename Cijk_type::i_iterator i_it=Cijk.i_begin(); i_it!=Cijk.i_end(); ++i_it) {
      size_t num_entry = 0;
      for (typename Cijk_type::ik_iterator k_it = Cijk.k_begin(i_it); k_it != Cijk.k_end(i_it); ++k_it) {
        for (typename Cijk_type::ikj_iterator j_it = Cijk.j_begin(k_it); j_it != Cijk.j_end(k_it); ++j_it) {
          ++num_entry;
        }
      max_num_entry = (num_entry > max_num_entry) ? num_entry : max_num_entry;
      }
    }

    // Replicated local map
    RCP<const Map> map =
      Tpetra::createLocalMapWithNode<TpetraLocalOrdinal,TpetraGlobalOrdinal,TpetraNode>(num_rows, comm);

    // Graph to be created
    RCP<Graph> graph = Tpetra::createCrsGraph(map, max_num_entry);
    
    // Loop over Cijk entries including a non-zero in the graph at
    // indices (i,j) if there is any k for which Cijk is non-zero
    for (typename Cijk_type::i_iterator i_it=Cijk.i_begin(); i_it!=Cijk.i_end(); ++i_it) {
      TpetraGlobalOrdinal i = index(i_it);
      for (typename Cijk_type::ik_iterator k_it = Cijk.k_begin(i_it); k_it != Cijk.k_end(i_it); ++k_it) {
        for (typename Cijk_type::ikj_iterator j_it = Cijk.j_begin(k_it); j_it != Cijk.j_end(k_it); ++j_it) {
          TpetraGlobalOrdinal j = index(j_it);
          graph->insertGlobalIndices(i, arrayView(&j, 1));
	      }
      }
    }

    // Sort, remove redundencies, transform to local, ...
    graph->fillComplete();

    return graph;
  }

  template <typename TpetraLocalOrdinal, typename TpetraGlobalOrdinal, typename TpetraNode,
            typename StokhosOrdinal, typename StokhosValue>
  void
  sparse3Tensor2TpetraMatrixMarket(
    const Stokhos::OrthogPolyBasis<StokhosOrdinal,StokhosValue>& basis,
    const Stokhos::Sparse3Tensor<StokhosOrdinal,StokhosValue>& Cijk,
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const std::string& file)
  {
    typedef Tpetra::CrsGraph<TpetraLocalOrdinal,TpetraGlobalOrdinal,TpetraNode> Graph;
    typedef Tpetra::CrsMatrix<StokhosValue,TpetraLocalOrdinal,TpetraGlobalOrdinal,TpetraNode> Matrix;
    typedef Tpetra::MatrixMarket::Writer<Matrix> Writer;

    Teuchos::RCP<Graph> graph =
      Stokhos::sparse3Tensor2TpetraCrsGraph<TpetraLocalOrdinal,TpetraGlobalOrdinal,TpetraNode>(
        basis, Cijk, comm);
    Writer::writeSparseGraphFile(file, *graph);
  }
} // namespace Stokhos

#endif // STOKHOS_TPETRA_SPARSE_3_TENSOR_UTILITIES_HPP
