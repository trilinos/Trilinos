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

#ifndef _Isorropia_EpetraCostDescriber_hpp_
#define _Isorropia_EpetraCostDescriber_hpp_

#include <Isorropia_configdefs.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_CostDescriber.hpp>

#ifdef HAVE_EPETRA
class Epetra_Vector;
class Epetra_CrsMatrix;

/** Isorropia is the namespace that contains isorropia's declarations
  for classes and functions.
*/
namespace Isorropia {

/** The Epetra namespace contains Isorropia's Epetra-specific
  classes and functions.
*/
namespace Epetra {

/** An Epetra-specific implementation of Isorropia::CostDescriber.
 */
class CostDescriber : public Isorropia::CostDescriber {
public:
  /** Constructor */
  CostDescriber();

  /** Destructor */
  virtual ~CostDescriber();

  /** Set parameters for the CostDescriber instance. The contents of the
      input paramlist object are copied into an internal ParameterList
      attribute. This class does not retain a reference
      to the input ParameterList after this method returns.
  */
  void setParameters(const Teuchos::ParameterList& paramlist);

  /** Vertex Weights; these typically correspond to matrix rows.
   */
  void setVertexWeights(Teuchos::RefCountPtr<Epetra_Vector> vwts);

  /** Graph edge weights; these typically correspond to matrix nonzeros.
   */
  void setGraphEdgeWeights(Teuchos::RefCountPtr<Epetra_CrsMatrix> gewts);

  /** Hypergraph edge weights; these typically correspond to matrix columns.
   */
  void setHypergraphEdgeWeights(Teuchos::RefCountPtr<Epetra_Vector> hgewts);

  /** Query whether non-default vertex weights are present. If this
    function returns false, the caller can assume that vertex weights
    are all 1.0.
  */
  bool haveVertexWeights() const;

  /** Get the number of vertices. Vertices typically correspond to
    matrix rows.
  */
  int getNumVertices() const;

  /** Get the lists of vertex ids and weights.
  */
  void getVertexWeights(int numVertices,
                                int* global_ids,
                                float* weights) const;

  /** Query whether non-default graph edge weights are present.
     If this function returns false, the caller can assume that
     graph edge weights are all 1.0.
  */
  bool haveGraphEdgeWeights() const;

  /** Get the number of graph edges for a specified vertex.
     Graph edges typically correspond to matrix nonzeros.
  */
  int getNumGraphEdges(int vertex_global_id) const;

  /** Get the graph edge weights for a specified vertex.
  */
  void getGraphEdgeWeights(int vertex_global_id,
                                   int num_neighbors,
                                   int* neighbor_global_ids,
                                   float* weights) const;

  /** Query whether non-default hypergraph edge weights are present.
     If this function returns false, the caller can assume that
     hypergraph edge weights are all 1.0.
  */
  bool haveHypergraphEdgeWeights() const;

  /** Get the number of Hypergraph edges. Hypergraph edges typically
     correspond to matrix columns.
  */
  int getNumHypergraphEdges() const;

  /** Get the hypergraph edge weights.
  */
  void getHypergraphEdgeWeights(int numEdges,
                                        int* global_ids,
                                        float* weights) const;

private:
  Teuchos::RefCountPtr<Epetra_Vector> vertex_weights_;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> graph_edge_weights_;
  Teuchos::RefCountPtr<Epetra_Vector> hypergraph_edge_weights_;
  Teuchos::ParameterList paramlist_;
};//class CostDescriber

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif

