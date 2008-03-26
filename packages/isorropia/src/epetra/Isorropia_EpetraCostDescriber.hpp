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

#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_CostDescriber.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <map>

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

class CostDescriber : public Isorropia::CostDescriber {
public:
  /** Constructor */
  CostDescriber();

  /** Destructor */
  ~CostDescriber();

  /** Set parameters for the CostDescriber instance. The contents of the
      input paramlist object are copied into an internal ParameterList
      attribute. This class does not retain a reference
      to the input ParameterList after this method returns.
  */
  void setParameters(const Teuchos::ParameterList& paramlist);

  /** Vertex Weights; these typically correspond to matrix rows.
   ** Each vertex (row) weight should be supplied by exactly one process
   ** (the process that "owns" that row).
   */
  void setVertexWeights(Teuchos::RefCountPtr<Epetra_Vector> vwts);

  /** Graph edge weights; these typically correspond to nonzeros of a
   ** square symmetric matrix.
   **
   ** Edges belonging to a vertex are represented by non-zeroes in the
   ** row that corresponds to the vertex.  All of the edges for a given
   ** vertex should be supplied by the one process that "owns" that row.
   */
  void setGraphEdgeWeights(Teuchos::RefCountPtr<Epetra_CrsMatrix> gewts);

  /** Hypergraph edge weights; these typically correspond to matrix columns.
   **
   ** Matrices that represent hypergraphs are not in general square.
   ** There may be more or fewer hyperedges (columns) than vertices (rows).
   **
   ** More than one process can supply a weight for the same column.  (So
   ** there is no concept of a process owning a hyperedge.)  Zoltan
   ** combines these weights according to the setting of the
   ** PHG_EDGE_WEIGHT_OPERATION parameter.
   */
  void setHypergraphEdgeWeights(Teuchos::RefCountPtr<Epetra_Vector> hgewts);

  /** Supplly a list of hypergraph edge weights and corresponding hypergraph
   ** global IDs.
   */
  void setHypergraphEdgeWeights(int numHGedges, int *hgGIDs, float *hgEwgts);

  /** Supplly a list of hypergraph edge weights and corresponding hypergraph
   ** global IDs.
   */
  void setHypergraphEdgeWeights(int numHGedges, int *hgGIDs, double *hgEwgts);

  /** Query whether non-default vertex weights are present. If this
    function returns false, the caller can assume that vertex weights
    are all 1.0, or this process has no vertices.
  */
  bool haveVertexWeights() const;

  /** Get the number of vertices. Vertices typically correspond to
    matrix rows.
  */
  int getNumVertices() const;

  /** Create a map from each vertex global ID to its weight.  Return the
      number of vertices in the map.  If vertex weights are defined, there
      is one weight for each vertex owned by this process. 
  */
  int getVertexWeights(std::map<int, float> &wgtMap) const;

 /** Get the lists of vertex ids and weights.
  */
  void getVertexWeights(int numVertices,
                      int* global_ids, float* weights) const; 

  /** Query whether non-default graph edge weights are present.
     If this function returns false, the caller can assume that
     graph edge weights are all 1.0, or that this process has
     no rows.
  */
  bool haveGraphEdgeWeights() const;

  /** Get the number of graph edges for a specified vertex.
     Graph edges typically correspond to matrix nonzeros.
  */
  int getNumGraphEdges(int vertex_global_id) const;

  /** Create a map from neighbor global ID to edge weight, and return
   ** the number of neighbors.
  */
  int getGraphEdgeWeights(int vertex_global_id, std::map<int, float> &wgtMap) const;

  /** Get the graph edge weights for a specified vertex.
  */
  void getGraphEdgeWeights(int vertex_global_id,
                                   int num_neighbors,
                                   int* neighbor_global_ids,
                                   float* weights) const; 

  /** Query whether non-default hypergraph edge weights are present.
     If this function returns false, the caller can assume that
     hypergraph edge weights are all 1.0, or that all hyperedge
     weights were supplied by other processes.
  */
  bool haveHypergraphEdgeWeights() const;

  /** Get the number of Hypergraph edge weights. Hypergraph edges typically
     correspond to matrix columns.
  */
  int getNumHypergraphEdgeWeights() const;

  /** Get the hypergraph edge weights.
  */
  void getHypergraphEdgeWeights(int numEdges,
                                        int* global_ids,
                                        float* weights) const;

  /*
   * It's possible that (haveVertexWeights() == false) because this
   * process has no rows of the matrix, but that
   * the application overall has supplied vertex weights.
   */
  bool haveGlobalVertexWeights() const;
  void setNumGlobalVertexWeights(int num);

  /*
   * It's possible that (haveGraphEdgeWeights() == false) because this
   * process has no rows of the matrix, but that
   * the application overall has supplied graph edge weights.
   */
  bool haveGlobalGraphEdgeWeights() const;
  void setNumGlobalGraphEdgeWeights(int num);

  /*
   * It's possible that (haveHypergraphEdgeWeights() == false) because 
   * the hypergraph edge weights were all supplied by other processes,
   * but that the application overall has supplied hypergraph edge weights.
   */
  bool haveGlobalHypergraphEdgeWeights() const;
  void setNumGlobalHypergraphEdgeWeights(int num);

private:
  void allocate_hg_edge_weights_(int n);
  void free_hg_edge_weights_();

  Teuchos::RefCountPtr<Epetra_Vector> vertex_weights_;

  Teuchos::RefCountPtr<Epetra_CrsMatrix> graph_edge_weights_;
  std::set<int> graph_self_edges_;

  int *hg_edge_gids_;
  float *hg_edge_weights_;
  int num_hg_edge_weights_;
  Teuchos::ParameterList paramlist_;

  int numGlobalVertexWeights_;
  int numGlobalGraphEdgeWeights_;
  int numGlobalHypergraphEdgeWeights_;

  // Create an array of the global IDs for the neighbors of the
  // given vertex, and also array of the edge weight for each edge.
  // Return the number of edges.  Self edges are not included.

  int getEdges(int vertexGID, int len, int *nborGID, float *weights) const;


};//class CostDescriber

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif

