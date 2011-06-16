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

#ifndef _Isorropia_EpetraCostDescriber_hpp_
#define _Isorropia_EpetraCostDescriber_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_CostDescriber.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <map>
#include <list>
#include <set>
#include <iostream>

#ifdef HAVE_EPETRA
class Epetra_Vector;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_CrsGraph;

namespace Isorropia {

/** The Isorropia::Epetra namespace implements a parallel partitioner
    that operates on Epetra objects.

    The objects to be partitioned may be the rows of an Epetra matrix
    or linear system.  There are two modes of operation.

    In one, the user supplies
    an Epetra object to Isorropia::Epetra::create_balanced_copy() and a
    new, rebalanced Epetra object is returned.

    In the other mode, the user supplies
    an Epetra object to Isorropia::Epetra::create_partitioner() and a
    Isorropia::Epetra::Partitioner object is returned.  This Partitioner
    object may be used to create an Isorropia::Epetra::Redistributor object.
    Finally the "redistribute" method of the Redistributor object
    may be used multiple times to redistribute Epetra objects with the
    same distribution as the Epetra object given to create_partitioner().

    In both cases, the application may supply row and edge weights with
    an Isorropia::Epetra::CostDescriber and may supply parameters with
    an Teuchos::ParameterList.

    If Trilinos has been built with the Zoltan parallel dynamic load
    balancing library (http://www.cs.sandia.gov/Zoltan),
    then Isorropia will use Zoltan to partition the Epetra object.
    The application can set Zoltan parameters with the ParameterList object.

    If Zoltan is not available, or the application has set the parameter
    "PARTITIONING_METHOD" to "SIMPLE_LINEAR", Isorropia will
    perform a simple linear partitioning
    of the object.  If the application supplied vertex weights with a
    CostDescriber, the linear partitioner will observe the vertex (row) weights.
    Otherwise the linear partitioner will use the number of non-zeroes in
    the row as the row weight.  The linear partitioner does not
    consider edge weights in the partitioning.
*/

class Operator;

namespace Epetra {

/** The CostDescriber class describes the vertex, edge and/or hyperedge weights.

    It is instantiated by the application to define
    weights, and then supplied to Isorropia with the
    Isorropia::Epetra::create_balanced_copy method or the
    Isorropia::Epetra::create_partitioner method.

    The CostDescriber can hold vertex (row) weights.
    For graph partitioning with Zoltan it may
    also define graph edge (nonzero) weights.  For hypergraph partitioning
    with Zoltan it may describe hyperedge (column) weights.

    If the application does not provide weights, reasonable defaults will
    be used.  Zoltan parameters (supplied with an Isorropia::Epetra::ParameterList)
    are available to override those defaults.
*/

// Forward declarations of friends

namespace ZoltanLib{
  class QueryObject;
}
class ZoltanLibClass;

class CostDescriber : public Isorropia::CostDescriber {

  // public methods are part of API, private methods are used by different
  // classes in isorropia

  friend class Isorropia::Operator;
  friend class Isorropia::Epetra::ZoltanLib::QueryObject;
  friend class Isorropia::Epetra::ZoltanLibClass;

public:
  /** Constructor */
  CostDescriber();

  /** Destructor */
  ~CostDescriber();

  /** Copy constructor */
  CostDescriber(const CostDescriber &costs);

  /**  Overloaded << operator for CostDescriber object
   */
  friend std::ostream& operator <<(std::ostream &, const Isorropia::Epetra::CostDescriber &cd);

  /** setVertexWeights is called by a process to supply the
      weight of each vertex (row) in the graph or hypergraph.
      If the object to be partitioned is instead an Epetra_MultiVector 
      representing real
      coordinates, then the weights represent the weight assigned to each coordinate.

      \param[in] vwgts  vector of weights, one for each vertex
   */
  void setVertexWeights(Teuchos::RCP<const Epetra_Vector> vwgts);

  /** setGraphEdgeWeights is called by a process to supply the weights for
      each of the edges of its vertices.  An edge corresponds to a non-zero
      in the row representing the vertex.

      This method is called only when performing graph partitioning with a
      square symmetric matrix.  For hypergraph partitioning call the equivalent
      hypergraph method.

      \param gewts an Epetra_CrsMatrix supplied by the application, each non-zero
             represents a weight for an edge
   */
  void setGraphEdgeWeights(Teuchos::RCP<const Epetra_CrsMatrix> gewts);

  /** setHypergraphEdgeWeights is called by processes in an application to
      supply weights for the hyperedges, which are represented by the columns
      of the matrix.  (A hyperedge can in general link more than one vertex.)
     
      Matrices that represent hypergraphs are not in general square.
      There may be more or fewer hyperedges (columns) than vertices (rows).

      \param hgewts  an Epetra_Vector containing the weights for each hyperedge.
   */
  void setHypergraphEdgeWeights(Teuchos::RCP<const Epetra_Vector> hgewts);

  /** setHypergraphEdgeWeights is called by processes in an application to
      supply weights for the hyperedges, which are represented by the columns
      of the matrix.  (A hyperedge can in general link more than one vertex.)

      Matrices that represent hypergraphs are not in general square.
      There may be more or fewer hyperedges (columns) than vertices (rows).

      More than one process can supply a weight for the same column.  (So
      there is no concept of a process owning a hyperedge.)  Zoltan
      combines these weights according to the setting of the
      PHG_EDGE_WEIGHT_OPERATION parameter.
     
       \param numHGedges  the length of the hgGIDs and heEwgts arrays
       \param hgGIDs      the global ID for each hyperedge this process will supply a weight for
       \param hgEwgts   the hyperedge weight corresponding to each hyperedge listed in hgGIDs
   */
  void setHypergraphEdgeWeights(int numHGedges, const int *hgGIDs, const float *hgEwgts);

  /** \copydoc Isorropia::Epetra::CostDescriber::setHypergraphEdgeWeights(int numHGedges, const int *hgGIDs, const float *hgEwgts)
   */

  void setHypergraphEdgeWeights(int numHGedges, const int *hgGIDs, const double *hgEwgts);

  /** Get the contents of this CostDescriber

      \param vertexWeights is set to a mapping from vertex global IDs to their weights
      \param graphEdgeWeights is a mapping from vertex global IDs to a map from neighboring
                 IDs to edge weights
      \param hypergraphEdgeWeights is a mapping from hyperedge (column) global IDs to hyperedge weights

   */
  void getCosts(std::map<int, float > &vertexWeights,
                std::map<int, std::map<int, float > > &graphEdgeWeights,
                std::map<int, float > &hypergraphEdgeWeights) const;

  /** Print out the contents of this CostDescriber
   */
  void show_cd(std::ostream &) const;

  // Add documentation

#ifdef USE_UTILS
  int compareBeforeAndAfterHypergraph( const Epetra_RowMatrix &in_m, const Epetra_RowMatrix &out_m,
            const Epetra_Import &importer,
            std::vector<double> &balance, std::vector<double> &cutn, std::vector<double> &cutl) const;

  int compareBeforeAndAfterGraph( const Epetra_RowMatrix &in_m, const Epetra_RowMatrix &out_m,
            const Epetra_Import &importer,
            std::vector<double> &balance, std::vector<int> &numCuts, std::vector<double> &cutWgt,
            std::vector<double> &cutn, std::vector<double> &cutl) const;

  int compareBeforeAndAfterGraph( const Epetra_CrsGraph &in_g, const Epetra_CrsGraph &out_g,
            const Epetra_Import &importer,
            std::vector<double> &balance, std::vector<int> &numCuts, std::vector<double> &cutWgt,
            std::vector<double> &cutn, std::vector<double> &cutl) const;

  int compareBeforeAndAfterImbalance(const Epetra_MultiVector &mv, const  Epetra_Import &importer,
                  std::vector<double> &min, std::vector<double> &max, std::vector<double> &avg) const;
#endif

private:

  void _transformWeights(const Epetra_Import &importer);

  int _compareBeforeAndAfterGraph( const Epetra_RowMatrix *in_m, const Epetra_RowMatrix *out_m,
            const Epetra_CrsGraph *in_g, const Epetra_CrsGraph *out_g,
            const Epetra_Import &importer,
            std::vector<double> &balance, std::vector<int> &numCuts, std::vector<double> &cutWgt,
            std::vector<double> &cutn, std::vector<double> &cutl) const;

  /** \copydoc Isorropia::CostDescriber::setParameters
   */
  void setParameters(const Teuchos::ParameterList& paramlist);

  /** \copydoc Isorropia::CostDescriber::haveVertexWeights
   */
  bool haveVertexWeights() const;
  /** \copydoc Isorropia::CostDescriber::getNumVertices
   */
  int getNumVertices() const;
  /** \copydoc Isorropia::CostDescriber::getVertexWeights
   */
  void getVertexWeights(int numVertices,
                      int* global_ids, float* weights) const; 
  /** \copydoc Isorropia::CostDescriber::haveGraphEdgeWeights
   */
  bool haveGraphEdgeWeights() const;
  /** \copydoc Isorropia::CostDescriber::getNumGraphEdges
   */
  int getNumGraphEdges(int vertex_global_id) const;

  /** Get the set of global IDs for the vertices that we have edge information for.

      \param gids will be set to the global IDs of the vertices for which neighbor and edge weight information have been provided
   */
  int getGraphEdgeVertices(std::set<int> &gids) const;


  /** \copydoc Isorropia::CostDescriber::getGraphEdgeWeights
   */
  void getGraphEdgeWeights(int vertex_global_id,
                                   int num_neighbors,
                                   int* neighbor_global_ids,
                                   float* weights) const; 
  /** \copydoc Isorropia::CostDescriber::haveHypergraphEdgeWeights
   */
  bool haveHypergraphEdgeWeights() const;
  /** \copydoc Isorropia::CostDescriber::getNumHypergraphEdgeWeights
   */
  int getNumHypergraphEdgeWeights() const;
  /** \copydoc Isorropia::CostDescriber::getHypergraphEdgeWeights
   */
  void getHypergraphEdgeWeights(int numEdges,
                                        int* global_ids,
                                        float* weights) const;

   /** Return the CostDescribers hypergraph edge weights as a map from hyperedge (column)
       global ID to weight.
  
       \param wgtMap will be set to a map from hyperedge global ID to hyperedge weight
   */
   int getHypergraphEdgeWeights(std::map<int, float> &wgtMap) const;


  /** Get vertex weights in the form of a map from vertex global ID to vertex weight.

      \param wgtMap a map supplied by the caller, the vertex weights will be added to this map
      \return the size of wgtMap
  */
  int getVertexWeights(std::map<int, float> &wgtMap) const;

  // TODO documentation
  const Epetra_Vector &getVertexWeights() { return *vertex_weights_;}


  /** getGraphEdgeWeights is called to obtain
      the graph edge weights for a given vertex (row).
 
      \param vertex_global_id the global ID of the vertex the caller wants edge weights for
      \param wgtMap a map from the global ID of each vertex neighbor to the weight of the edge formed by the vertex and this neighbor
      \return  the count of the neighbors of vertex_global_id
  */
  int getGraphEdgeWeights(int vertex_global_id, std::map<int, float> &wgtMap) const;



  /** haveGlobalVertexWeights returns true if any process in the application has
        supplied vertex weights, it returns false otherwise.
   */
  bool haveGlobalVertexWeights() const;

  /** setNumGlobalVertexWeights may be used to set the count of the
        global number of vertex weights supplied to the CostDescriber
   */
  void setNumGlobalVertexWeights(int num);

  /** haveGlobalGraphEdgeWeights returns true if any process in the application has
        supplied graph edge weights, it returns false otherwise.
   */
  bool haveGlobalGraphEdgeWeights() const;

  /** setNumGlobalGraphEdgeWeights may be used to set the count of the
        global number of graph edge weights supplied to the CostDescriber
   */
  void setNumGlobalGraphEdgeWeights(int num);

  /** haveGlobalHypergraphEdgeWeights returns true if any process in the application has
        supplied hyperedge weights, it returns false otherwise.
   */
  bool haveGlobalHypergraphEdgeWeights() const;

  /** setNumGlobalHypergraphEdgeWeights may be used to set the count of the
        global number of hyperedge weights supplied to the CostDescriber
   */
  void setNumGlobalHypergraphEdgeWeights(int num);

  /** Dynamically allocate storage for hypergraph edge weights.
   */
  void allocate_hg_edge_weights_(int n);

  /** Free storage used by hypergraph edge weights.
   */
  void free_hg_edge_weights_();

  Teuchos::RCP<const Epetra_Vector> vertex_weights_;
  Teuchos::RCP<const Epetra_CrsMatrix> graph_edge_weights_;
  std::set<int> graph_self_edges_;

  Teuchos::ParameterList paramlist_;

  int *hg_edge_gids_;
  float *hg_edge_weights_;
  int num_hg_edge_weights_;

  int numGlobalVertexWeights_;
  int numGlobalGraphEdgeWeights_;
  int numGlobalHypergraphEdgeWeights_;

  /** getEdges creates an array of the neighbors and edge weights for given vertex.
      Self edges are not included.

      \param vertexGID the global ID of the vertex (must be one owned by calling process)
      \param len of preallocated nborGID and weights arrays
      \param nborGID on return contains the global ID of each vertex neighboring vertexGID,
                         allocated by caller      
      \param weights on return contains the weight for each edge formed by the vertices in nborGID

      \return the number of neighbors in nborGID is returned

   */
  int getEdges(int vertexGID, int len, int *nborGID, float *weights) const;

};//class CostDescriber

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif

