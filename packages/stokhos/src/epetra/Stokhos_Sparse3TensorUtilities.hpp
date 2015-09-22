// $Id$ 
// $Source$ 
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

#ifndef STOKHOS_SPARSE_3_TENSOR_UTILITIES_HPP
#define STOKHOS_SPARSE_3_TENSOR_UTILITIES_HPP

#include "Stokhos_Sparse3Tensor.hpp"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_BlockMap.h"
#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

namespace Stokhos {

  //! Build an Epetra_CrsGraph from a sparse 3 tensor
  /*!
   * Builds a sparse graph from a sparse 3 tensor by summing over the third
   * index.  This graph then represents the sparsity pattern of the stochastic
   * part of the block stochastic Galerkin operator.  Redistributing the graph
   * should then provide a suitable parallel distribution for block
   * stochastic Galerkin linear solves.
   */
  template <typename ordinal_type, typename value_type>
  Teuchos::RCP<Epetra_CrsGraph>
  sparse3Tensor2CrsGraph(
    const Stokhos::OrthogPolyBasis<ordinal_type,value_type>& basis,
    const Stokhos::Sparse3Tensor<ordinal_type,value_type>& Cijk,
    const Epetra_Comm& comm) 
  {
    typedef Stokhos::Sparse3Tensor<ordinal_type,value_type> Cijk_type;

    // Number of stochastic rows
    ordinal_type num_rows = basis.size();

    // Replicated local map
    Epetra_LocalMap map(num_rows, 0, comm);

    // Graph to be created
    Teuchos::RCP<Epetra_CrsGraph> graph = 
      Teuchos::rcp(new Epetra_CrsGraph(Copy, map, 0));
    
    // Loop over Cijk entries including a non-zero in the graph at
    // indices (i,j) if there is any k for which Cijk is non-zero
    for (typename Cijk_type::k_iterator k_it=Cijk.k_begin(); 
	 k_it!=Cijk.k_end(); ++k_it) {
      for (typename Cijk_type::kj_iterator j_it = Cijk.j_begin(k_it); 
	   j_it != Cijk.j_end(k_it); ++j_it) {
	ordinal_type j = index(j_it);
	for (typename Cijk_type::kji_iterator i_it = Cijk.i_begin(j_it);
	     i_it != Cijk.i_end(j_it); ++i_it) {
	  ordinal_type i = index(i_it);
	  graph->InsertGlobalIndices(i, 1, &j);
	}
      }
    }

    // Sort, remove redundencies, transform to local, ...
    graph->FillComplete();

    return graph;
  }

  //! Build an Epetra_CrsGraph from a sparse 3 tensor
  /*!
   * Builds a sparse graph from a sparse 3 tensor by summing over the third
   * index.  This graph then represents the sparsity pattern of the stochastic
   * part of the block stochastic Galerkin operator.  Redistributing the graph
   * should then provide a suitable parallel distribution for block
   * stochastic Galerkin linear solves.
   */
  template <typename ordinal_type, typename value_type>
  Teuchos::RCP<Epetra_CrsGraph>
  sparse3Tensor2CrsGraph(
    const Stokhos::Sparse3Tensor<ordinal_type,value_type>& Cijk,
    const Epetra_BlockMap& map) 
  {
    typedef Stokhos::Sparse3Tensor<ordinal_type,value_type> Cijk_type;

    // Graph to be created
    Teuchos::RCP<Epetra_CrsGraph> graph = 
      Teuchos::rcp(new Epetra_CrsGraph(Copy, map, 0));
    
    // Loop over Cijk entries including a non-zero in the graph at
    // indices (i,j) if there is any k for which Cijk is non-zero
    for (typename Cijk_type::k_iterator k_it=Cijk.k_begin(); 
	 k_it!=Cijk.k_end(); ++k_it) {
      for (typename Cijk_type::kj_iterator j_it = Cijk.j_begin(k_it); 
	   j_it != Cijk.j_end(k_it); ++j_it) {
	ordinal_type j = index(j_it);
	for (typename Cijk_type::kji_iterator i_it = Cijk.i_begin(j_it);
	     i_it != Cijk.i_end(j_it); ++i_it) {
	  ordinal_type i = index(i_it);
	  graph->InsertGlobalIndices(i, 1, &j);
	}
      }
    }

    // Sort, remove redundencies, transform to local, ...
    graph->FillComplete();

    return graph;
  }

  template <typename ordinal_type, typename value_type>
  void
  sparse3Tensor2MatrixMarket(
    const Stokhos::OrthogPolyBasis<ordinal_type,value_type>& basis,
    const Stokhos::Sparse3Tensor<ordinal_type,value_type>& Cijk,
    const Epetra_Comm& comm,
    const std::string& file)
  {
    Teuchos::RCP<Epetra_CrsGraph> graph = 
      Stokhos::sparse3Tensor2CrsGraph(basis, Cijk, comm);
    Epetra_CrsMatrix mat(Copy, *graph);
    mat.FillComplete();
    mat.PutScalar(1.0);
    EpetraExt::RowMatrixToMatrixMarketFile(file.c_str(), mat);
  }

  template <typename ordinal_type, typename value_type>
  void
  sparse3Tensor2MatrixMarket(
    const Stokhos::Sparse3Tensor<ordinal_type,value_type>& Cijk,
    const Epetra_BlockMap& map,
    const std::string& file)
  {
    Teuchos::RCP<Epetra_CrsGraph> graph = 
      Stokhos::sparse3Tensor2CrsGraph(Cijk, map);
    Epetra_CrsMatrix mat(Copy, *graph);
    mat.FillComplete();
    mat.PutScalar(1.0);
    EpetraExt::RowMatrixToMatrixMarketFile(file.c_str(), mat);
  }

} // namespace Stokhos

#endif // SPARSE_3_TENSOR_UTILITIES_HPP
