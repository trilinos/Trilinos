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

  Library(Teuchos::RCP<const Epetra_CrsGraph> input_graph);
  Library(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	  Teuchos::RCP<CostDescriber> costs);

  Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix);
  Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
 	  Teuchos::RCP<CostDescriber> costs);

  Library(Teuchos::RCP<const Epetra_MultiVector> input_coords);
  Library(Teuchos::RCP<const Epetra_MultiVector> input_coords,
          Teuchos::RCP<const Epetra_MultiVector> weights);

  virtual ~Library();

  virtual int
  repartition(Teuchos::ParameterList& paramlist,
	      std::vector<int>& myNewElements,
	      int& exportsSize,
	      std::vector<int>& imports) = 0;
// 	      std::map<int,int>& exports,
// 	      std::map<int,int>& imports) = 0;

  virtual int
  color(Teuchos::ParameterList& paramlist,
	std::vector<int>& myNewElements) = 0 ;

  virtual int
  order(Teuchos::ParameterList& paramlist,
	std::vector<int>& myNewElements) = 0 ;

  /** Object on which library acts - does it represent a graph,
      a hypergraph, or a collection of geometric coordinates?
    */

  static const int hypergraph_input_ = 1;
  static const int graph_input_ = 2;
  static const int geometric_input_ = 3;
  static const int unspecified_input_ = 4;

  int input_type_;

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

