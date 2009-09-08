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

#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraPartitioner.hpp>

#include <Teuchos_RCP.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Comm.h>
#endif

namespace Isorropia {

namespace Epetra {

#ifdef HAVE_EPETRA

Redistributor::Redistributor(Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner)
  : partitioner_(partitioner),
  importer_(),
  target_map_(),
  created_importer_(false)
{
  if (!partitioner_->alreadyComputed()) {
    partitioner_->partition();
  }
}

Redistributor::~Redistributor()
{
}

void
Redistributor::redistribute(const Epetra_SrcDistObject& src,
				 Epetra_DistObject& target)
{
  if (!created_importer_) {
    create_importer(src.Map());
  }
  else {
    if (!src.Map().PointSameAs(importer_->SourceMap()) ||
	!target.Map().PointSameAs(importer_->TargetMap())) {
      created_importer_ = false;
      create_importer(src.Map());
    }
  }

  target.Import(src, *importer_, Insert);
}

Teuchos::RCP<Epetra_CrsGraph>
Redistributor::redistribute(const Epetra_CrsGraph& input_graph, bool callFillComplete)
{
  if (!created_importer_) {
    create_importer(input_graph.RowMap());
  }

  // First obtain the length of each of my new rows

  int myOldRows = input_graph.NumMyRows();
  int myNewRows = target_map_->NumMyElements();

  double *nnz = new double [myOldRows];
  for (int i=0; i < myOldRows; i++){
    nnz[i] = input_graph.NumMyIndices(i);
  }

  Epetra_Vector oldRowSizes(Copy, input_graph.RowMap(), nnz);

  if (myOldRows)
    delete [] nnz;

  Epetra_Vector newRowSizes(*target_map_);

  newRowSizes.Import(oldRowSizes, *importer_, Insert);

  int *rowSize = new int [myNewRows];
  for (int i=0; i< myNewRows; i++){
    rowSize[i] = static_cast<int>(newRowSizes[i]);
  }

  // Receive new rows, send old rows

  Teuchos::RCP<Epetra_CrsGraph> new_graph =
    Teuchos::rcp(new Epetra_CrsGraph(Copy, *target_map_, rowSize, true));

  if (myNewRows)
    delete [] rowSize;

  new_graph->Import(input_graph, *importer_, Insert);

  // Set the new domain map such that
  // (a) if old DomainMap == old RangeMap, preserve this property,
  // (b) otherwise, let the new DomainMap be the old DomainMap 
  const Epetra_BlockMap *newDomainMap;
  if (input_graph.DomainMap().SameAs(input_graph.RangeMap()))
     newDomainMap = &(new_graph->RangeMap());
  else
     newDomainMap = &(input_graph.DomainMap());

  if (callFillComplete && (!new_graph->Filled()))
    new_graph->FillComplete(*newDomainMap, *target_map_);

  return( new_graph );
}

Teuchos::RCP<Epetra_CrsMatrix>
Redistributor::redistribute(const Epetra_CrsMatrix& input_matrix, bool callFillComplete)
{
  if (!created_importer_) {
    create_importer(input_matrix.RowMap());
  }

  // First obtain the length of each of my new rows

  int myOldRows = input_matrix.NumMyRows();
  int myNewRows = target_map_->NumMyElements();

  double *nnz = new double [myOldRows];
  for (int i=0; i < myOldRows; i++){
    nnz[i] = input_matrix.NumMyEntries(i);
  }

  Epetra_Vector oldRowSizes(Copy, input_matrix.RowMap(), nnz);

  if (myOldRows)
    delete [] nnz;

  Epetra_Vector newRowSizes(*target_map_);

  newRowSizes.Import(oldRowSizes, *importer_, Insert);

  int *rowSize=0;
  if(myNewRows){
    rowSize = new int [myNewRows];
    for (int i=0; i< myNewRows; i++){
      rowSize[i] = static_cast<int>(newRowSizes[i]);
    }
  }

  // Receive new rows, send old rows

  Teuchos::RCP<Epetra_CrsMatrix> new_matrix =
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, *target_map_, rowSize, true));

  if (myNewRows)
    delete [] rowSize;

  new_matrix->Import(input_matrix, *importer_, Insert);

  // Set the new domain map such that
  // (a) if old DomainMap == old RangeMap, preserve this property,
  // (b) otherwise, let the new DomainMap be the old DomainMap 
  const Epetra_Map *newDomainMap;
  if (input_matrix.DomainMap().SameAs(input_matrix.RangeMap()))
     newDomainMap = &(new_matrix->RangeMap());
  else
     newDomainMap = &(input_matrix.DomainMap());

  if (callFillComplete && (!new_matrix->Filled()))
    new_matrix->FillComplete(*newDomainMap,  *target_map_);

  return( new_matrix );
}

Teuchos::RCP<Epetra_CrsMatrix>
Redistributor::redistribute(const Epetra_RowMatrix& input_matrix, bool callFillComplete)
{
  if (!created_importer_) {
    create_importer(input_matrix.RowMatrixRowMap());
  }
 // First obtain the length of each of my new rows

  int myOldRows = input_matrix.NumMyRows();
  int myNewRows = target_map_->NumMyElements();


  double *nnz = new double [myOldRows];
  int val;
  for (int i=0; i < myOldRows; i++){
    input_matrix.NumMyRowEntries(i, val);
    nnz[i] = static_cast<double>(val);
  }

  Epetra_Vector oldRowSizes(Copy, input_matrix.RowMatrixRowMap(), nnz);

  if (myOldRows)
    delete [] nnz;

  Epetra_Vector newRowSizes(*target_map_);

  newRowSizes.Import(oldRowSizes, *importer_, Insert);

  int *rowSize = new int [myNewRows];
  for (int i=0; i< myNewRows; i++){
    rowSize[i] = static_cast<int>(newRowSizes[i]);
  }

  // Receive new rows, send old rows

  Teuchos::RCP<Epetra_CrsMatrix> new_matrix =
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, *target_map_, rowSize, true));

  if (myNewRows)
    delete [] rowSize;

  new_matrix->Import(input_matrix, *importer_, Insert);

  // Set the new domain map such that
  // (a) if old DomainMap == old RangeMap, preserve this property,
  // (b) otherwise, use the original OperatorDomainMap
  //if (input_matrix.NumGlobalRows() == input_matrix.NumGlobalCols()){
  if (input_matrix.OperatorDomainMap().SameAs(input_matrix.OperatorRangeMap())){
    if (callFillComplete && (!new_matrix->Filled()))
      new_matrix->FillComplete();
  }
  else {
    if (callFillComplete && (!new_matrix->Filled()))
      new_matrix->FillComplete(input_matrix.OperatorDomainMap(), *target_map_);
  }

  return( new_matrix );
}

Teuchos::RCP<Epetra_Vector>
Redistributor::redistribute(const Epetra_Vector& input_vector)
{
  if (!created_importer_) {
    create_importer(input_vector.Map());
  }

  Teuchos::RCP<Epetra_Vector> new_vec = Teuchos::rcp(new Epetra_Vector(*target_map_));

  new_vec->Import(input_vector, *importer_, Insert);

  return( new_vec );
}

Teuchos::RCP<Epetra_MultiVector>
Redistributor::redistribute(const Epetra_MultiVector& input_vector)
{
  if (!created_importer_) {
    create_importer(input_vector.Map());
  }

  Teuchos::RCP<Epetra_MultiVector> new_vec =
    Teuchos::rcp(new Epetra_MultiVector(*target_map_, input_vector.NumVectors()));

  new_vec->Import(input_vector, *importer_, Insert);

  return( new_vec );
}

// Reverse redistribute methods (for vectors). 

void
Redistributor::redistribute_reverse(const Epetra_Vector& input_vector, Epetra_Vector& output_vector)
{
  if (!created_importer_) {
    create_importer(input_vector.Map());
  }

  // Export using the importer
  output_vector.Export(input_vector, *importer_, Insert);

}

void
Redistributor::redistribute_reverse(const Epetra_MultiVector& input_vector, Epetra_MultiVector& output_vector)
{
  if (!created_importer_) {
    create_importer(input_vector.Map());
  }

  // Export using the importer
  output_vector.Export(input_vector, *importer_, Insert);

}

void Redistributor::create_importer(const Epetra_BlockMap& src_map)
{
  if (created_importer_) return;

  if (partitioner_->numProperties() > src_map.Comm().NumProc()) {
    throw Isorropia::Exception("Cannot redistribute: Too many parts for too few processors.");
  }

  target_map_ = partitioner_->createNewMap();

  importer_ = Teuchos::rcp(new Epetra_Import(*target_map_, src_map));

  created_importer_ = true;
}

#endif //HAVE_EPETRA

}//namespace Epetra

}//namespace Isorropia

