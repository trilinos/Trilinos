//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
//@HEADER

#ifndef _Isorropia_EpetraLibrary_hpp_
#define _Isorropia_EpetraLibrary_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_EpetraCostDescriber.hpp>

#ifdef HAVE_EPETRA
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_LinearProblem;

namespace Isorropia {

namespace Epetra {
  class CostDescriber;

/** An implementation of the Partitioner interface that operates on
    Epetra matrices and linear systems.

*/

class Library {
public:

  Library(Teuchos::RCP<const Epetra_CrsGraph> input_graph, int itype = unspecified_input_);
  Library(Teuchos::RCP<const Epetra_CrsGraph> input_graph, Teuchos::RCP<const Epetra_MultiVector> input_coords,
          int itype = unspecified_input_);
  Library(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	  Teuchos::RCP<CostDescriber> costs, int itype = unspecified_input_);
  Library(Teuchos::RCP<const Epetra_CrsGraph> input_graph, Teuchos::RCP<CostDescriber> costs, 
	  Teuchos::RCP<const Epetra_MultiVector> input_coords, Teuchos::RCP<const Epetra_MultiVector> weights,
          int itype = unspecified_input_);
  Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, int itype = unspecified_input_);
  Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, Teuchos::RCP<const Epetra_MultiVector> input_coords,
          int itype = unspecified_input_);
  Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
 	  Teuchos::RCP<CostDescriber> costs, int itype = unspecified_input_);
  Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, Teuchos::RCP<CostDescriber> costs, 
	  Teuchos::RCP<const Epetra_MultiVector> input_coords, Teuchos::RCP<const Epetra_MultiVector> weights,
          int itype = unspecified_input_);
  Library(Teuchos::RCP<const Epetra_MultiVector> input_coords, int itype = unspecified_input_);
  Library(Teuchos::RCP<const Epetra_MultiVector> input_coords,
          Teuchos::RCP<const Epetra_MultiVector> weights, int itype = unspecified_input_);
  Library(Teuchos::RCP<const Epetra_BlockMap> input_map, int itype = unspecified_input_);

  virtual ~Library();

  virtual int
  repartition(Teuchos::ParameterList& paramlist,
	      std::vector<int>& myNewElements,
	      int& exportsSize,
	      std::vector<int>& imports) = 0;

  virtual int
  color(Teuchos::ParameterList& paramlist,
	std::vector<int>& colorAssignment) = 0 ;

  virtual int
  order(Teuchos::ParameterList& paramlist,
	std::vector<int>& orderAssignment) = 0 ;

  /** input_type_ == hgraph_input_
      This indicates that the matrix or graph represents a hypergraph.  Columns
      represent hyperedges, and row (vertex) partitioning is to be performed.
    */

  static const int hgraph_input_ = 1;

  /** input_type_ == hgraph2d_finegrain_input_
      This indicates that the matrix or graph represents a hypergraph.  Columns
      represent hyperedges, and non-zeroes are to be partitioned.
    */
  static const int hgraph2d_finegrain_input_ = 2;

  /** input_type_ == graph_input_ 
      This indicates that the square symmetric matrix or graph represents a graph
      in the sense that row/column IDs are vertices and non-zeroes represent
      edges.  The vertices are to be partitioned.
    */
  static const int graph_input_ = 3;

  /** input_type_ == geometric_input_
      This indicates that the Epetra_MultiVector represents geometric
      coordinates.  The MultiVector should have 1, 2 or 3 vectors,
      representing 1, 2 or 3 dimensional coordinates.  The coordinates
      are to be partitioned.
    */
  static const int geometric_input_ = 4;

  /** input_type_ == hgraph_graph_input_
      This indicates that the Epetra_MultiVector represents a hypergraph
      and graph (see above).  This is necessary for hierarchical partitioning
      with both hypergraph and graph methods.
    */
  static const int hgraph_graph_input_ = 5;

  /** input_type_ == hgraph_geom_input_
      This indicates that the Epetra_MultiVector represents a hypergraph
      and graph (see above).  This is necessary for hierarchical partitioning
      with both hypergraph and geometric methods.
    */
  static const int hgraph_geometric_input_ = 6;

  /** input_type_ == graph_geom_input_
      This indicates that the Epetra_MultiVector represents a hypergraph
      and graph (see above).  This is necessary for hierarchical partitioning
      with both graph and geometric methods.
    */
  static const int graph_geometric_input_ = 7;

  /** input_type_ == hgraph_graph_geom_input_
      This indicates that the Epetra_MultiVector represents a hypergraph
      and graph (see above).  This is necessary for hierarchical partitioning
      using hypergraph, graph, and geometric methods.
    */
  static const int hgraph_graph_geometric_input_ = 8;


  /** input_type_ == simple_input_
      This is used to indicate that a simple partitiong method
      (block, cyclic, or random) will be used.
    */
  static const int simple_input_ = 9;


  /** input_type_ == unspecified_input_ 
      This value is the "unset" state for the input_type_ instance variable.
    */

  static const int unspecified_input_ = 10;

  int input_type_;

  int numPartSizes;
  int *partGIDs;
  float *partSizes;

protected:

  Teuchos::RCP<const Epetra_BlockMap> input_map_;
  Teuchos::RCP<const Epetra_CrsGraph> input_graph_;
  Teuchos::RCP<const Epetra_RowMatrix> input_matrix_;
  Teuchos::RCP<const Epetra_MultiVector> input_coords_;
  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs_;
  Teuchos::RCP<const Epetra_MultiVector> weights_;

  virtual int precompute();

  virtual int postcompute() = 0;

};//class Library

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

