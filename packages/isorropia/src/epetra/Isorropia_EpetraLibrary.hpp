//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

************************************************************************
*/
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
  Library(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	  Teuchos::RCP<CostDescriber> costs, int itype = unspecified_input_);

  Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, int itype = unspecified_input_);
  Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
 	  Teuchos::RCP<CostDescriber> costs, int itype = unspecified_input_);

  Library(Teuchos::RCP<const Epetra_MultiVector> input_coords, int itype = unspecified_input_);
  Library(Teuchos::RCP<const Epetra_MultiVector> input_coords,
          Teuchos::RCP<const Epetra_MultiVector> weights, int itype = unspecified_input_);

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

  /** input_type_ == unspecified_input_ 
      This value is the "unset" state for the input_type_ instance variable.
    */
  static const int unspecified_input_ = 5;

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

