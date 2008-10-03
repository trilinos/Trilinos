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

#include <Isorropia_EpetraLibrary.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

#endif

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>

namespace Isorropia {

#ifdef HAVE_EPETRA

namespace Epetra {

Library::
Library(Teuchos::RCP<const Epetra_CrsGraph> input_graph)
  : input_map_(),
    input_graph_(input_graph),
    input_matrix_(),
    costs_()
{
  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
}

Library::
Library(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	 Teuchos::RCP<CostDescriber> costs)
  : input_map_(),
    input_graph_(input_graph),
    input_matrix_(),
    costs_(costs)
{

  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
}

Library::
Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix)
  : input_map_(),
    input_graph_(),
    input_matrix_(input_matrix),
    costs_()
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
}

Library::
Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	Teuchos::RCP<CostDescriber> costs)
  : input_map_(),
    input_graph_(),
    input_matrix_(input_matrix),
    costs_(costs)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
}

void Library::
setInput(Teuchos::RCP<const Epetra_CrsGraph> input_graph)
{
  input_graph_ = input_graph;
  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
}

void Library::
setInput(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	 Teuchos::RCP<CostDescriber> costs)
{
  input_graph_ = input_graph;
  costs_ = costs;
  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
}

void Library::
setInput(Teuchos::RCP<const Epetra_RowMatrix> input_matrix)
{
  input_matrix_ = input_matrix;
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
}

void Library::
setInput(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	 Teuchos::RCP<CostDescriber> costs)
{
  input_matrix_ = input_matrix;
  costs_ = costs;
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
}


Library::~Library()
{
}


int Library::precompute()
{
  std::string str1("Isorropia::Epetra::Operator::precompute ");
  std::string str2;

  if ((input_graph_.get() == 0 && input_matrix_.get() == 0)
      || (input_graph_.get() != 0 && input_matrix_.get() != 0)) {
    str2 = "ERROR: not holding valid input graph (x)OR matrix.";
    throw Isorropia::Exception(str1+str2);
    return (-1);
  }

  if (inputType_ == "GRAPH") {
    bool square = false;
    bool symmetric = false;
    if (input_graph_.get() != 0){
      if (input_graph_->NumGlobalRows() == input_graph_->NumGlobalCols()){
	square = true;
	// TODO - is there a quick way to figure out if the graph is
	// symmetric?  I can't see a way to do it.  For now we let
	// Zoltan figure this out.
	symmetric = true;
      }
    }
    else{
      if (input_matrix_->NumGlobalRows() == input_matrix_->NumGlobalCols()){
	square = true;
	// TODO - is there a quick way to figure out if the matrix is
	// symmetric?  I can't see a way to do it.  For now we let
	// Zoltan figure this out.
	symmetric = true;
      }
    }
    if (!square){
      str2 = "LB_METHOD=GRAPH, matrix or graph must be square";
      throw Isorropia::Exception(str1+str2);
      return (-1);
    }
    if (!symmetric){
      str2 = "LB_METHOD=GRAPH, matrix or graph must be symmetric";
      throw Isorropia::Exception(str1+str2);
      return (-1);
    }
  }

  return (0);
}


} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

