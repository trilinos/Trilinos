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

#ifndef _Isorropia_TpetraCostDescriber_hpp_
#define _Isorropia_TpetraCostDescriber_hpp_

//#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_CostDescriber.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <EEP_Isorropia_TpetraZoltanLib.hpp>

#include <map>
#include <list>
#include <set>
#include <iostream>

#ifdef USE_UTILS
#include "ispatest_lbeval_utils.hpp"
#endif

#include <Tpetra_Vector_decl.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_CrsGraph_decl.hpp>

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

namespace Tpetra {

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

//namespace ZoltanLib{
//  class QueryObject;
//}

template<class LocalOrdinal,
	 class GlobalOrdinal,
	 class Node> class ZoltanLibClass;

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class CostDescriber : public Isorropia::CostDescriber {

  // public methods are part of API, private methods are used by different
  // classes in isorropia

  friend class Isorropia::Operator;
  friend class Isorropia::Tpetra::ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>; // EEP
  friend class Isorropia::Tpetra::ZoltanLibClass<LocalOrdinal, GlobalOrdinal, Node>;

public:
  /** Constructor */
  CostDescriber();

  /** Destructor */
  ~CostDescriber();

  /** Copy constructor */
  CostDescriber(const CostDescriber &costs);

  /**  Overloaded << operator for CostDescriber object
   */
  friend std::ostream& operator <<(std::ostream &, const Isorropia::Tpetra::CostDescriber<LocalOrdinal, GlobalOrdinal, Node> &cd);

  /** setVertexWeights is called by a process to supply the
      weight of each vertex (row) in the graph or hypergraph.
      If the object to be partitioned is instead an Epetra_MultiVector 
      representing real
      coordinates, then the weights represent the weight assigned to each coordinate.

      \param vwgts [in]  vector of weights, one for each vertex
   */
  void setVertexWeights(Teuchos::RCP< const ::Tpetra::Vector<double, LocalOrdinal, GlobalOrdinal, Node> > vwgts); // EEP__

  /** setGraphEdgeWeights is called by a process to supply the weights for
      each of the edges of its vertices.  An edge corresponds to a non-zero
      in the row representing the vertex.

      This method is called only when performing graph partitioning with a
      square symmetric matrix.  For hypergraph partitioning call the equivalent
      hypergraph method.

      \param gewts an Epetra_CrsMatrix supplied by the application, each non-zero
             represents a weight for an edge
   */
  void setGraphEdgeWeights(Teuchos::RCP< const ::Tpetra::CrsMatrix<double, LocalOrdinal, GlobalOrdinal, Node> > gewts); // EEP__

  /** setHypergraphEdgeWeights is called by processes in an application to
      supply weights for the hyperedges, which are represented by the columns
      of the matrix.  (A hyperedge can in general link more than one vertex.)
     
      Matrices that represent hypergraphs are not in general square.
      There may be more or fewer hyperedges (columns) than vertices (rows).

      \param hgewts  an Epetra_Vector containing the weights for each hyperedge.
   */
  void setHypergraphEdgeWeights(Teuchos::RCP< const ::Tpetra::Vector<double, LocalOrdinal, GlobalOrdinal, Node> > hgewts); // EEP__

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

  /** \copydoc Isorropia::Tpetra::CostDescriber::setHypergraphEdgeWeights(int numHGedges, const int *hgGIDs, const float *hgEwgts)
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

private:

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


  // TODO documentation
  const ::Tpetra::Vector<double, LocalOrdinal, GlobalOrdinal, Node> &getVertexWeights() { return *vertex_weights_;}


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

  Teuchos::RCP< const ::Tpetra::Vector<double, LocalOrdinal, GlobalOrdinal, Node> > vertex_weights_;
  Teuchos::RCP< const ::Tpetra::CrsMatrix<double, LocalOrdinal, GlobalOrdinal, Node> > graph_edge_weights_;
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

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::CostDescriber()
  : vertex_weights_(),
    graph_edge_weights_(),
    graph_self_edges_(),
    paramlist_(),
    hg_edge_gids_(NULL),
    hg_edge_weights_(NULL),
    num_hg_edge_weights_(0),
    numGlobalVertexWeights_(0),
    numGlobalGraphEdgeWeights_(0),
    numGlobalHypergraphEdgeWeights_(0)
{
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::~CostDescriber()
{
  free_hg_edge_weights_();
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::setParameters(const Teuchos::ParameterList& paramlist)
{
  paramlist_ = paramlist;
}

#if 0 // EEP
/** Supply a vector of vertex (row) weights.  If rows are distributed, then
    each process must supply a weight for each of its rows.  (Alternatively
    the application can supply no vertex weights at all.)  The weights should
    be in the same order as the rows in the Epetra object being partitioned.
*/
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void CostDescriber::setVertexWeights(Teuchos::RCP<const Epetra_Vector> vwts)
{
  if (vertex_weights_.get() != 0){
    vertex_weights_.release();
  }
  vertex_weights_ = vwts;
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
CostDescriber::setGraphEdgeWeights(Teuchos::RCP<const Epetra_CrsMatrix> gewts)
{
  if (graph_edge_weights_.get() != 0){
    graph_edge_weights_.release();
    graph_self_edges_.clear();
  }
  graph_edge_weights_ = gewts;

  if (gewts->NumMyDiagonals() > 0){

    // Save list of self edges - we omit them in the Zoltan query functions

    const Epetra_Map &rowmap = gewts->RowMap();

    Epetra_Vector diag(rowmap);
    
    gewts->ExtractDiagonalCopy(diag);

    int nvals = gewts->NumMyRows();
    double *entry;
    diag.ExtractView(&entry);
    for (int i=0; i<nvals; i++){
      if (entry[i] != 0){
        graph_self_edges_.insert( rowmap.GID(i));
      }
    }
  }
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
CostDescriber::setHypergraphEdgeWeights(Teuchos::RCP<const Epetra_Vector> hgewts)
{
  free_hg_edge_weights_();
  const Epetra_BlockMap& map = hgewts->Map();

  int numWeights = map.NumMyElements();

  if (numWeights > 0){
    allocate_hg_edge_weights_(numWeights);
    map.MyGlobalElements(hg_edge_gids_);
    double *v;
    int stride;
    hgewts->ExtractView(&v, &stride);
    for (int i=0; i<numWeights; i++){
      hg_edge_weights_[i] = (float)v[i];
    }
  }
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
CostDescriber::setHypergraphEdgeWeights(int numHGedges, const int *hgGIDs, const float *hgEwgts)
{
  free_hg_edge_weights_();
  if (numHGedges > 0){
    allocate_hg_edge_weights_(numHGedges);
    for (int i=0; i<numHGedges; i++){
      hg_edge_weights_[i] = hgEwgts[i];
      hg_edge_gids_[i] = hgGIDs[i];
    }
  }
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
CostDescriber::setHypergraphEdgeWeights(int numHGedges, const int *hgGIDs, const double *hgEwgts)
{
  free_hg_edge_weights_();
  if (numHGedges > 0){
    allocate_hg_edge_weights_(numHGedges);
    for (int i=0; i<numHGedges; i++){
      hg_edge_weights_[i] = (float)hgEwgts[i];
      hg_edge_gids_[i] = hgGIDs[i];
    }
  }
}
#endif // EEP

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
bool CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::haveVertexWeights() const
{
  const int n = getNumVertices();
  return( n > 0);
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::getNumVertices() const
{
  return( vertex_weights_.get()==0 ? 0 : vertex_weights_->getLocalLength() ); // EEP MyLength()
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::getVertexWeights(int numVertices,
                                     int* global_ids,
                                     float* weights) const
{
  if (getNumVertices() == 0){
    return;
  }

  if (numVertices != vertex_weights_->getMap()->getLocalNumElements()) {
    throw std::runtime_error/*Isorropia::Exception*/("CostDescriber::getVertexWeights: wrong numVertices");
  }

  //vertex_weights_->getMap()->getGlobalElements(global_ids); // EEP___

  double* vals = nullptr; // vertex_weights_->Values(); // EEP___
  for(int i=0; i<numVertices; ++i) {
    weights[i] = vals[i];
  }
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
bool CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::haveGraphEdgeWeights() const
{
  int n = 0;
  if (graph_edge_weights_.get()){
    n = graph_edge_weights_->getLocalNumEntries(); // NumMyNonzeros(); // EEP
  }
  return( n > 0);
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::getNumGraphEdges(int vertex_global_id) const
{
  int n = 0;
  if (graph_edge_weights_.get() != 0) {
    int lrid = graph_edge_weights_->getRowMap()->getLocalElement(vertex_global_id); // LRID(vertex_global_id); // EEP
    if (lrid >= 0){   
      n = graph_edge_weights_->getNumEntriesInLocalRow(lrid); // NumMyEntries(lrid); // EEP

      if (graph_self_edges_.size() > 0){
        std::set<int>::const_iterator it = graph_self_edges_.find(vertex_global_id);
        if (it != graph_self_edges_.end()){
          n--;     // don't count self edges
        }
      }

    }
  }

  return n;
}

// ONLY called if graph_self_edge_ is not null

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::getEdges(int vertexGID, int len, int *nborGID, float *weights) const
{
  //int vertexLID = 0; // graph_edge_weights_->getRowMap()->LID(vertexGID); // EEP___
  int numRealEdges = getNumGraphEdges(vertexGID);  //excluding self edges

  if (numRealEdges < 1){
    return 0;
  }

  if (len < numRealEdges ){
    throw std::runtime_error/*Isorropia::Exception*/("CostDescriber::getEdges: length of allocated arrays");
  }

  int self_edge = 0;
  std::set<int>::const_iterator it = graph_self_edges_.find(vertexGID);
  if (it != graph_self_edges_.end()){
    self_edge = 1;
  }

  //int *viewIds=NULL;
  double *viewWgts=NULL;
  int numedges;         // including self edges

  int rc = 0; // graph_edge_weights_->ExtractMyRowView(vertexLID, numedges, viewWgts, viewIds); // EEP___

  if (rc){
    throw std::runtime_error/*Isorropia::Exception*/("CostDescriber::getEdges: Extract matrix row view");
  }

  if (numedges != (numRealEdges + self_edge)){
    throw std::runtime_error/*Isorropia::Exception*/("CostDescriber::getEdges: Extract matrix count");
  }

  int nextID = 0;

  for (int j=0; j < numedges; j++){
    int gid = 0; // graph_edge_weights_->ColMap()->GID(viewIds[j]); // EEP___
    if (gid == vertexGID) continue;   // skip the self edges

    nborGID[nextID] = gid;
    weights[nextID] = (float)viewWgts[j];

    nextID++;
  }

  return nextID;
}

#if 0 // EEP
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int CostDescriber::getGraphEdgeVertices(std::set<int> &gids) const
{
  gids.clear();
  int ngids = 0;

  if (haveGraphEdgeWeights()){
    const Epetra_Map &rowmap = graph_edge_weights_->RowMap();
    ngids = rowmap.NumMyElements();
    for (int i=0; i<ngids; i++){
      gids.insert(rowmap.GID(i));
    }
  }
  return ngids;
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int CostDescriber::getGraphEdgeWeights(int vertex_global_id, std::map<int, float> &wgtMap) const
{
  int rowlen = getNumGraphEdges(vertex_global_id);

  if (rowlen < 1){
    return 0;
  }

  float *wgt = new float [rowlen];
  int *nborGID = new int [rowlen];

  int numEdges = getEdges(vertex_global_id, rowlen, nborGID, wgt);

  for (int i=0; i<numEdges; i++){
    wgtMap[nborGID[i]] = wgt[i];
  }

  if (rowlen > 0){
    delete [] nborGID;
    delete [] wgt;
  }
  return numEdges;
}
#endif // EEP

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::getGraphEdgeWeights(int vertex_global_id,
                                       int num_neighbors,
                                       int* neighbor_global_ids,
                                       float* weights) const
{
  int rowlen = getNumGraphEdges(vertex_global_id);

  if (rowlen < 1){
    return;
  }

  if (rowlen > num_neighbors) {
    throw std::runtime_error/*Isorropia::Exception*/("CostDescriber::getGraphEdgeWeights: wrong num_neighbors");
  }

  getEdges(vertex_global_id, num_neighbors, neighbor_global_ids, weights);
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
bool CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::haveHypergraphEdgeWeights() const
{
  return(num_hg_edge_weights_ > 0);
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::getNumHypergraphEdgeWeights() const
{
  return num_hg_edge_weights_;
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::getHypergraphEdgeWeights(int numEdges,
					     int* global_ids,
					     float* weights) const
{
  if (numEdges != num_hg_edge_weights_) {
    throw std::runtime_error/*Isorropia::Exception*/("CostDescriber::getHypergraphEdgeWeights: wrong numEdges");
  }

  for(int i=0; i<numEdges; ++i) {
    weights[i] = hg_edge_weights_[i];
    global_ids[i] = hg_edge_gids_[i];
  }
}

#if 0 // EEP
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int CostDescriber::getHypergraphEdgeWeights(std::map<int, float> &wgtMap) const
{
  int nEdges = num_hg_edge_weights_;
  if (nEdges < 1) return 0;

  for(int i=0; i<nEdges; ++i) {
    wgtMap[hg_edge_gids_[i]] = hg_edge_weights_[i];
  }
  return nEdges;
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void CostDescriber::getCosts(std::map<int, float> &vertexWeights,
                           std::map<int, std::map<int, float > > &graphEdgeWeights, 
                           std::map<int, float> &hypergraphEdgeWeights) const
{
  if (haveVertexWeights()){
    getVertexWeights(vertexWeights);
  }

  if (haveHypergraphEdgeWeights()){
    getHypergraphEdgeWeights(hypergraphEdgeWeights);
  }

  if (haveGraphEdgeWeights()){
    std::set<int> vgids;
//     int ngids = getGraphEdgeVertices(vgids);
    std::set<int>::iterator curr;
    std::set<int>::iterator end = vgids.end();
    curr = vgids.begin();
    while (curr != end){
      std::map<int, float> nborMap;
      getGraphEdgeWeights(*curr, nborMap);
      graphEdgeWeights[*curr] = nborMap;
      curr++;
    }
  }
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
bool CostDescriber::haveGlobalVertexWeights() const
{
  return (numGlobalVertexWeights_ > 0);
}
#endif // EEP

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::setNumGlobalVertexWeights(int num)
{
  this->numGlobalVertexWeights_ = num;
}

#if 0 // EEP
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
bool CostDescriber::haveGlobalGraphEdgeWeights() const
{
  return (numGlobalGraphEdgeWeights_ > 0);
}
#endif // EEP

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::setNumGlobalGraphEdgeWeights(int num)
{
  this->numGlobalGraphEdgeWeights_ = num;
}

#if 0 // EEP
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
bool CostDescriber::haveGlobalHypergraphEdgeWeights() const
{
  return (numGlobalHypergraphEdgeWeights_ > 0);
}
#endif // EEP

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::setNumGlobalHypergraphEdgeWeights(int num)
{
  this->numGlobalHypergraphEdgeWeights_ = num;
}

#if 0 // EEP
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void CostDescriber::allocate_hg_edge_weights_(int n)
{
  free_hg_edge_weights_();
  if (n > 0){
    hg_edge_gids_ = new int [n];
    hg_edge_weights_ = new float [n];
    num_hg_edge_weights_ = n;
  }
}
#endif // EEP

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void CostDescriber<LocalOrdinal, GlobalOrdinal, Node>::free_hg_edge_weights_()
{
  if (hg_edge_gids_){
    delete [] hg_edge_gids_;
    delete [] hg_edge_weights_;
    hg_edge_gids_ = NULL;
    hg_edge_weights_ = NULL;
    num_hg_edge_weights_ = 0;
  }
}

#if 0 // EEP
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void CostDescriber::show_cd(std::ostream &os) const
{
  int nv = getNumVertices();
  int nhge = getNumHypergraphEdgeWeights();

  int *gids = NULL;
  if (nv){
    os << "Vertices and weights" << std::endl << "  ";
    gids = new int [nv];
    float *w = new float [nv];

    getVertexWeights(nv, gids, w);
    for (int j=0; j<nv; j++){
      os << gids[j] << " (" << w[j] << ") ";
    }
    os << std::endl;
    delete [] w;
  }
  else{
    os << "No vertex weights" << std::endl;
  }
  if (gids && haveGraphEdgeWeights()){
    os << "Graph edge (non zero) weights for each vertex (row)" << std::endl;
    for (int i=0; i < nv; i++){
      int vid = gids[i];
      std::map<int, float> wgts;

      getGraphEdgeWeights(vid, wgts);

      os << "  Vertex (row) GID " << vid << std::endl << "    ";
      std::map<int, float>::iterator curr;

      for(curr = wgts.begin(); curr != wgts.end(); curr++){
        os << curr->first << " (" << curr->second << ") ";
      } 
      os << std::endl;
    }
  }
  else{
    os << "No graph edge weights" << std::endl;
  }
  if (nhge){
    int *colgids = new int [nhge];
    float *wgts = new float [nhge];

    getHypergraphEdgeWeights(nhge, colgids, wgts);

    os << "Hypergraph Edge (column) weights" << std::endl << "  ";

    for (int j=0; j < nhge; j++){
      os << colgids[j] << " (" << wgts[j] << ") ";
    }
    os << std::endl;

    delete [] colgids;
    delete [] wgts;
  }
  else{
    os << "No hypergraph edge weights" << std::endl;
  }
  
  if (gids) delete [] gids;

  nv = numGlobalVertexWeights_;
  int nge = numGlobalGraphEdgeWeights_;
  nhge = numGlobalHypergraphEdgeWeights_;

  if (paramlist_.begin() == paramlist_.end()){
    os << "No parameters set" << std::endl;
  }
  else{
    os << "Have some parameters set" << std::endl;
  }

  if (haveGlobalVertexWeights()){
    os << "Number of global vertices " << nv << std::endl;
  }
  else{
    os << "Don't know number of global vertices " << std::endl;
  }

  if (haveGlobalGraphEdgeWeights()){
    os << "Number of global graph edge weights " << nge << std::endl;
  }
  else{
    os << "Don't know number of global graph edge weights " << std::endl;
  }

  if (haveGlobalHypergraphEdgeWeights()){
    os << "Number of global hypergraph edge weights " << nhge << std::endl;
  }
  else{
    os << "Don't know number of global hypergraph edge weights " << std::endl;
  }
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
std::ostream& operator <<(std::ostream& os, const Isorropia::Epetra::CostDescriber &cd)
{
  cd.show_cd(os);
  return os;
}

#endif // EEP

}//namespace Epetra
}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

