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

#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_Exception.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_BlockMap.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

namespace Isorropia {
namespace Epetra {

CostDescriber::CostDescriber()
  : vertex_weights_(),
    graph_edge_weights_(),
    hypergraph_edge_weights_(),
    paramlist_()
{
}

CostDescriber::~CostDescriber()
{
}

void CostDescriber::setParameters(const Teuchos::ParameterList& paramlist)
{
  paramlist_ = paramlist;
}

void CostDescriber::setVertexWeights(Teuchos::RefCountPtr<Epetra_Vector> vwts)
{
  vertex_weights_ = vwts;
}

void
CostDescriber::setGraphEdgeWeights(Teuchos::RefCountPtr<Epetra_CrsMatrix> gewts)
{
  graph_edge_weights_ = gewts;
}

void
CostDescriber::setHypergraphEdgeWeights(Teuchos::RefCountPtr<Epetra_Vector> hgewts)
{
  hypergraph_edge_weights_ = hgewts;
}

bool CostDescriber::haveVertexWeights() const
{
  return( vertex_weights_.get() != 0 );
}

int CostDescriber::getNumVertices() const
{
  return( vertex_weights_.get()==0 ? 0 : vertex_weights_->MyLength() );
}

void CostDescriber::getVertexWeights(int numVertices,
				     int* global_ids,
				     float* weights) const
{
  if (vertex_weights_.get() == 0) {
    throw Isorropia::Exception("CostDescriber::getVertexWeights: no vertex weights");
  }

  const Epetra_BlockMap& map = vertex_weights_->Map();

  if (numVertices != map.NumMyElements()) {
    throw Isorropia::Exception("CostDescriber::getVertexWeights: wrong numVertices");
  }

  map.MyGlobalElements(global_ids);

  double* vals = vertex_weights_->Values();
  for(int i=0; i<numVertices; ++i) {
    weights[i] = vals[i];
  }
}

bool CostDescriber::haveGraphEdgeWeights() const
{
  return( graph_edge_weights_.get() != 0);
}

int CostDescriber::getNumGraphEdges(int vertex_global_id) const
{
  if (graph_edge_weights_.get() == 0) {
    throw Isorropia::Exception("CostDescriber::getNumGraphEdges: no graph edge weights");
  }

  return( graph_edge_weights_->NumGlobalEntries(vertex_global_id) );
}

void CostDescriber::getGraphEdgeWeights(int vertex_global_id,
				       int num_neighbors,
				       int* neighbor_global_ids,
				       float* weights) const
{
  if (graph_edge_weights_.get() == 0) {
    throw Isorropia::Exception("CostDescriber::getGraphEdgeWeights: no graph edge weights");
  }

  int rowlen = graph_edge_weights_->NumGlobalEntries(vertex_global_id);
  if (rowlen != num_neighbors) {
    throw Isorropia::Exception("CostDescriber::getGraphEdgeWeights: wrong num_neighbors");
  }

  int err = graph_edge_weights_->Graph().ExtractGlobalRowCopy(vertex_global_id,
							      num_neighbors, rowlen,
							      neighbor_global_ids);

  double* coefptr = 0;
  err = graph_edge_weights_->ExtractGlobalRowView(vertex_global_id, rowlen,
						  coefptr);

  for(int i=0; i<rowlen; ++i) {
    weights[i] = coefptr[i];
  }
}

bool CostDescriber::haveHypergraphEdgeWeights() const
{
  return( hypergraph_edge_weights_.get() != 0);
}

int CostDescriber::getNumHypergraphEdges() const
{
  return( hypergraph_edge_weights_.get()==0 ?
	  0 : hypergraph_edge_weights_->MyLength());
}

void CostDescriber::getHypergraphEdgeWeights(int numEdges,
					     int* global_ids,
					     float* weights) const
{
  if (hypergraph_edge_weights_.get() == 0) {
    throw Isorropia::Exception("CostDescriber::getHypergraphEdgeWeights: no hypergraph edge weights");
  }

  const Epetra_BlockMap& map = hypergraph_edge_weights_->Map();

  if (numEdges != map.NumMyElements()) {
    throw Isorropia::Exception("CostDescriber::getHypergraphEdgeWeights: wrong numEdges");
  }

  map.MyGlobalElements(global_ids);

  double* vals = hypergraph_edge_weights_->Values();
  for(int i=0; i<numEdges; ++i) {
    weights[i] = vals[i];
  }
}

}//namespace Epetra
}//namespace Isorropia
#endif

