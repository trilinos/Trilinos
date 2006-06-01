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
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

#include <Isorropia_Redistributor.hpp>
#include <Isorropia_Rebalance.hpp>
#include <Isorropia_Zoltan_Rebalance.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra_utils.hpp>
#include <Isorropia_Partitioner.hpp>

#include <Teuchos_RefCountPtr.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#endif

/** Isorropia is the namespace that contains isorropia's declarations
  for classes and functions.
*/
namespace Isorropia {

#ifdef HAVE_EPETRA

Redistributor::Redistributor(Teuchos::RefCountPtr<const Partitioner> partitioner)
  : partitioner_(partitioner),
  importer_(NULL),
  target_map_(NULL),
  created_importer_(false)
{
}

Redistributor::~Redistributor()
{
  delete importer_;
  delete target_map_;
}

Teuchos::RefCountPtr<Epetra_CrsGraph>
Redistributor::redistribute(const Epetra_CrsGraph& input_graph)
{
  if (!created_importer_) {
    create_importer(input_graph.RowMap());
  }

  Teuchos::RefCountPtr<Epetra_CrsGraph> new_graph =
    Teuchos::rcp(new Epetra_CrsGraph(Copy, *target_map_, 0));

  new_graph->Import(input_graph, *importer_, Insert);

  return( new_graph );
}

Teuchos::RefCountPtr<Epetra_CrsMatrix>
Redistributor::redistribute(const Epetra_CrsMatrix& input_matrix)
{
  if (!created_importer_) {
    create_importer(input_matrix.RowMap());
  }

  Teuchos::RefCountPtr<Epetra_CrsMatrix> new_matrix =
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, *target_map_, 0));

  new_matrix->Import(input_matrix, *importer_, Insert);

  return( new_matrix );
}

Teuchos::RefCountPtr<Epetra_CrsMatrix>
Redistributor::redistribute(const Epetra_RowMatrix& input_matrix)
{
  if (!created_importer_) {
    create_importer(input_matrix.RowMatrixRowMap());
  }

  Teuchos::RefCountPtr<Epetra_CrsMatrix> new_matrix =
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, *target_map_, 0));

  new_matrix->Import(input_matrix, *importer_, Insert);

  return( new_matrix );
}

Teuchos::RefCountPtr<Epetra_MultiVector>
Redistributor::redistribute(const Epetra_MultiVector& input_vector)
{
  if (!created_importer_) {
    create_importer(input_vector.Map());
  }

  Teuchos::RefCountPtr<Epetra_MultiVector> new_vec =
    Teuchos::rcp(new Epetra_MultiVector(*target_map_, input_vector.NumVectors()));

  new_vec->Import(input_vector, *importer_, Insert);

  return( new_vec );
}

void Redistributor::create_importer(const Epetra_BlockMap& src_map)
{
  if (created_importer_) return;

  int numMyNewElements = partitioner_->newNumMyElements();
  std::vector<int> myElements(numMyNewElements);
  partitioner_->myNewElements(&myElements[0], numMyNewElements);

  target_map_ = new Epetra_Map(-1, numMyNewElements, &myElements[0],
			       0, src_map.Comm());

  importer_ = new Epetra_Import(*target_map_, src_map);

  created_importer_ = true;
}

#endif //HAVE_EPETRA

}//namespace Isorropia

