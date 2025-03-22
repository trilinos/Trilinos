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

#ifndef ISORROPIA_TPETRA_ZOLTAN_QUERYOBJECT_H
#define ISORROPIA_TPETRA_ZOLTAN_QUERYOBJECT_H

//#include "Isorropia_ConfigDefs.hpp"

#include <Teuchos_RCP.hpp>

#include <zoltan_cpp.h>

#include <set>
#include <map>

#include <EEP_Isorropia_TpetraCostDescriber.hpp>

#include <Tpetra_CrsGraph_decl.hpp>
//#include <Epetra_Comm_decl.hpp>
#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Vector_decl.hpp>
#include <Tpetra_RowMatrix_decl.hpp>

#include <map>

// if non-zero, print out query function responses
//#define DEBUG_QUERIES 1
// process to print out responses, or -1 for all
//#define DEBUG_PROC -1

namespace Isorropia {

namespace Tpetra {

  //class CostDescriber;

/** The ZoltanLib namespace within the Tpetra namespace contains the
    classes and functions that use the Zoltan library to partition an
    Tpetra object.
*/

namespace ZoltanLib {

/** QueryObject is a class that contains the query functions required
    by the Zoltan library.

    These methods are not part of the Isorropia API (except to Zoltan).
    They are called by Isorropia itself and by Zoltan.

    For a better understanding of Zoltan's query functions, see the
    Zoltan User's Guide at http://www.cs.sandia.gov/zoltan
 */
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class QueryObject
{

  /** haveGraph is true if we have CrsGraph, and not a CrsMatrix or
      a MultiVector.
    */
  const bool haveGraph_;

  /** The CrsGraph.  The QueryObject must be constructed with one of
      an Tpetra::CrsGraph, an Tpetra::RowMatrix or an Tpetra::MultiVector.
    */
  Teuchos::RCP<const ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> graph_;

  /** The CrsMatrix. 
      The QueryObject must be constructed with one of
      an Tpetra::CrsGraph, an Tpetra::RowMatrix or an Tpetra::MultiVector.
    */
  //Teuchos::RCP<const ::Tpetra::RowMatrix<LocalOrdinal, GlobalOrdinal, Node>> matrix_; // EEP

  /** The MultiVector containing 1, 2 or 3 dimensional coordinates.  If
      supplied, we will perform geometric partitioning.
      The QueryObject must be constructed with one of
      an Tpetra::CrsGraph, an Tpetra::RowMatrix or an Tpetra::MultiVector.
    */
  //Teuchos::RCP<const ::Tpetra::MultiVector<LocalOrdinal, GlobalOrdinal, Node>> coords_; // EEP

  /** The graph or matrix row map, or the MultiVector map
    */
  const ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> *rowMap_;

  /** The graph or matrix column map
    */
  const ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> *colMap_;

  /** The CostDescriber contains optional vertex and/or edge weights for
      graph and hypergraph partitioning.
    */

  Teuchos::RCP< const ::Isorropia::Tpetra::CostDescriber<LocalOrdinal, GlobalOrdinal, Node> > costs_;

  /** The MultiVector contains optional object (point) weights for
      geometric partitioning.  Zoltan currently will use only the
      weights in the first vector (1 dimensional point weights).
    */

  //Teuchos::RCP<const ::Tpetra::MultiVector<LocalOrdinal, GlobalOrdinal, Node>> weights_; // EEP

  std::map<int,int> procmap_;
  std::set<int> graph_self_edges_;

  unsigned int myProc_;
  unsigned int base_;

  void fill_procmap();

  /** My_Number_Objects() returns the number of objects currently
      assigned to this process.  (The objects are interpreted as
      graph vertices for Graph partitioning, as hypergraph
      vertices for hypergraph partitioning, or as coordinates for
      geometric partitioning.)
   */
  int My_Number_Objects(int *ierr);

  /** My_ObjectList() returns to Zoltan the global ID and weight of the
      objects currently assigned to this process.
   */
  void My_Object_List  (int num_gid_entries, int num_lid_entries,
		     ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
		     int weight_dim, float * object_weights, int * ierr );

  /** My_Number_Edges_Multi() is a query function used for graph partitioning
      only.  It returns to Zoltan the number of edges (non-zeroes) that each
      vertex (row) has.
   */
  void My_Number_Edges_Multi  (int num_gid_entries, int num_lid_entries,
	       int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
	       int *num_edges, int * ierr );

  /** My_Edge_List_Multi() is a query function used for graph partitioning
      only.  For each vertex (row), it returns a list of the global ID of
      each neighbor (non-zero) and the process owning that neighbor (that row).
   */
  void My_Edge_List_Multi(int num_gid_entries, int num_lid_entries,
	       int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
	       int *num_edges, ZOLTAN_ID_PTR neighbor_global_ids, int * neighbor_procs,
	       int weight_dim, float * edge_weights, int * ierr );

  /** My_HG_Size_CS() is a query function used for hypergraph partitioning
      only. Zoltan calls this query to get size of the non-zeros lists from the QueryObject. 
   */
  void My_HG_Size_CS (int* num_lists, int* num_pins, int* format,
			  int * ierr );

  /** My_HG_CS() is a query function used for hypergraph partitioning
      only. Zoltan calls this query to get the non-zeros lists from the QueryObject. 
   */
  void My_HG_CS (int num_gid_entries, int num_row_or_col, int num_pins,
	   int format, ZOLTAN_ID_PTR vtxedge_GID, int* vtxedge_ptr, ZOLTAN_ID_PTR pin_GID,
				       int * ierr );


  /** My_FGHG_CS() is a query function used for fine-grain hypergraph partitioning
      only. Zoltan calls this query to get the non-zeros lists from the QueryObject. 
   */
  void My_FGHG_CS (int num_gid_entries, int num_row_or_col, int num_pins,
	   int format, ZOLTAN_ID_PTR vtxedge_GID, int* vtxedge_ptr, ZOLTAN_ID_PTR pin_GID,
				       int * ierr );

  /** My_HG_Size_Edge_Weights() is a query function used for hypergraph partitioning
      only. Zoltan calls this query to get number of hyperedge weights that this
      QueryObject will be providing.
   */
  void My_HG_Size_Edge_Weights(int* num_edges, int* ierr);
  
  /** My_HG_Edge_Weights() is a query function used for hypergraph partitioning
      only. Zoltan calls this query to get hyperedge weights from this
      QueryObject.
   */
  void My_HG_Edge_Weights(int num_gid_entries, int num_lid_entries, int num_edges, int edge_weight_dim,
        ZOLTAN_ID_PTR edge_GID, ZOLTAN_ID_PTR edge_LID, float* edge_weights, int* ierr);

  /** My_Number_Geom() is a query function used for geometric partitioning
      only. Zoltan calls this query to get the dimension of the geometric
      coordinates.
   */
  int My_Number_Geom(int *ierr);

  /** My_Geom_Multi() is a query function used for geometric partitioning
      only. Zoltan calls this query to get a list of coordinates from the QueryObject.
   */
  void My_Geom_Multi(int num_gid_entries, int num_lid_entries,
        int num_obj, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int num_dim,
        double *geom_vec, int *ierr);

 public:

  /** Constructor
   */
  QueryObject( Teuchos::RCP<const ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> graph,
	       Teuchos::RCP<const ::Isorropia::Tpetra::CostDescriber<LocalOrdinal, GlobalOrdinal, Node>> costs,
	       int inputType);

#if 0 // EEP
  /** Constructor
   */
  QueryObject( Teuchos::RCP<const ::Tpetra::RowMatrix<LocalOrdinal, GlobalOrdinal, Node>> matrix,
	       Teuchos::RCP<const ::Isorropia::Tpetra::CostDescriber<LocalOrdinal, GlobalOrdinal, Node>> costs,
	       int inputType);

  /** Constructor
   */
  QueryObject( Teuchos::RCP<const ::Tpetra::MultiVector<LocalOrdinal, GlobalOrdinal, Node>> coords,
               Teuchos::RCP<const ::Tpetra::MultiVector<LocalOrdinal, GlobalOrdinal, Node>> weights);


  /** Constructor
   */
  QueryObject(Teuchos::RCP<const ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> graph,
	      Teuchos::RCP<const ::Isorropia::Tpetra::CostDescriber<LocalOrdinal, GlobalOrdinal, Node>> costs, 
              Teuchos::RCP<const ::Tpetra::MultiVector<LocalOrdinal, GlobalOrdinal, Node>> coords,
              Teuchos::RCP<const ::Tpetra::MultiVector<LocalOrdinal, GlobalOrdinal, Node>> weights,
              int inputType);


  /** Constructor
   */
  QueryObject(Teuchos::RCP<const ::Tpetra::RowMatrix<LocalOrdinal, GlobalOrdinal, Node>> matrix,
	      Teuchos::RCP<const ::Isorropia::Tpetra::CostDescriber<LocalOrdinal, GlobalOrdinal, Node>> costs,
              Teuchos::RCP<const ::Tpetra::MultiVector<LocalOrdinal, GlobalOrdinal, Node>> coords,
              Teuchos::RCP<const ::Tpetra::MultiVector<LocalOrdinal, GlobalOrdinal, Node>> weights,
              int inputType);


  /** Constructor
   */
  QueryObject( Teuchos::RCP<const ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> inputMap,
	       int inputType);
#endif



  /** Destructor
   */
  virtual ~QueryObject();

  /** input_type_ == hgraph_input_.
      This indicates that the matrix or graph represents a hypergraph.  Columns
      represent hyperedges, and row (vertex) partitioning is to be performed.
    */

  static const int hgraph_input_ = 1;

  /** input_type_ == hgraph2d_finegrain_input_.
      This indicates that the matrix or graph represents a hypergraph.  Columns
      represent hyperedges, and non-zeroes are to be partitioned.
    */
  static const int hgraph2d_finegrain_input_ = 2;

  /** input_type_ == graph_input_.
      This indicates that the square symmetric matrix or graph represents a graph
      in the sense that row/column IDs are vertices and non-zeroes represent
      edges.  The vertices are to be partitioned.
    */
  static const int graph_input_ = 3;

  /** input_type_ == geometric_input_.
      This indicates that the Tpetra::MultiVector represents geometric
      coordinates.  The MultiVector should have 1, 2 or 3 vectors,
      representing 1, 2 or 3 dimensional coordinates.  The coordinates
      are to be partitioned.
    */
  static const int geometric_input_ = 4;

  /** input_type_ == hgraph_graph_input_
      This indicates that the Tpetra::MultiVector represents a hypergraph
      and graph (see above).  This is necessary for hierarchical partitioning
      with both hypergraph and graph methods.
  */
  static const int hgraph_graph_input_ = 5;

  /** input_type_ == hgraph_geometric_input_
      This indicates that the Tpetra::MultiVector represents a hypergraph
      and graph (see above).  This is necessary for hierarchical partitioning
      with both hypergraph and geometric methods.
  */
  static const int hgraph_geometric_input_ = 6;

  /** input_type_ == graph_geometric_input_
      This indicates that the Tpetra::MultiVector represents graph and
      has geometric coordinates (see above).  
      This is necessary for hierarchical partitioning                                                                                                            
      with both graph and geometric methods.
  */
  static const int graph_geometric_input_ = 7;

  /** input_type_ == hgraph_graph_geometric_input_
      This indicates that the Tpetra::MultiVector represents a hypergraph
      and graph and has geometric coordinates(see above).  This is necessary for 
      hierarchical partitioning using hypergraph, graph, and geometric methods.
  */
  static const int hgraph_graph_geometric_input_ = 8;

  /** input_type_ == simple_input_
      This indicates that a simple method (block, cyclic, or random) will be used.
  */
  static const int simple_input_ = 9;


  /** input_type_ == unspecified_input_.
      This value is the "unset" state for the input_type_ instance variable.
    */
  static const int unspecified_input_ = 10;



  /** The input_type_ indicates how the object to be partitioned is to
      be interpreted - as a graph or a hypergraph for row partitioning, 
      as a hypergraph for fine-grain partitioning, or as a list of coordinates for geometric
      partitioning.
    */
  int input_type_;

  /** Return the map associated with the object to be partitioned.
   */
  const ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> &RowMap(void){ return *rowMap_;};

  /** Return true if any of the processes in the application have defined
      vertex weights.
    */
  bool haveVertexWeights();

  /** Return true if any of the processes in the application have defined
      graph edge weights.
    */
  bool haveGraphEdgeWeights();

  /** Return true if any of the processes in the application have defined
      hypergraph edge weights.
    */
  bool haveHypergraphEdgeWeights();

  // General query functions

  /** The interface to a particular QueryObject's My_Number_Objects query function.
   */
  static int Number_Objects(void *data, int *ierr);
  
  /** The interface to a particular QueryObject's My_Object_List query function.
   */
  static void Object_List  ( void * data, int num_gid_entries, int num_lid_entries,
		     ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
		     int weight_dim, float * object_weights, int * ierr );

  // Query functions for graph partitioning only
  
  /** The interface to a particular QueryObject's My_Number_Edges_Multi query function.
   */
  static void Number_Edges_Multi  ( void * data, int num_gid_entries, int num_lid_entries,
	       int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
	       int *num_edges, int * ierr );

  /** The interface to a particular QueryObject's My_Edges_Multi query function.
   */
  static void Edge_List_Multi( void * data, int num_gid_entries, int num_lid_entries,
	       int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
	       int *num_edges, ZOLTAN_ID_PTR neighbor_global_ids, int * neighbor_procs,
	       int weight_dim, float * edge_weights, int * ierr );

  // Query functions for hypergraph partitioning only
  
  /** The interface to a particular QueryObject's My_HG_Size_CS query function.
   */
  static void HG_Size_CS ( void * data, int* num_lists, int* num_pins, int* format,
			  int * ierr );
  /** The interface to a particular QueryObject's My_HG_CS query function.
   */
  static void HG_CS ( void * data, int num_gid_entries, int num_row_or_col, int num_pins,
	   int format, ZOLTAN_ID_PTR vtxedge_GID, int* vtxedge_ptr, ZOLTAN_ID_PTR pin_GID,
				       int * ierr );
  /** The interface to a particular QueryObject's My_HG_Size_Edge_Weights query function.
   */
  static void HG_Size_Edge_Weights(void * data, int* num_edges, int* ierr);
  
  /** The interface to a particular QueryObject's My_HG_Edge_Weights query function.
   */
  static void HG_Edge_Weights(void * data,
        int num_gid_entries, int num_lid_entries, int num_edges, int edge_weight_dim,
        ZOLTAN_ID_PTR edge_GID, ZOLTAN_ID_PTR edge_LID, float* edge_weights, int* ierr);

  /** The interface to a particular QueryObject's My_Number_Geom query function.
   */
  static int Number_Geom(void *data, int *ierr);

  /** The interface to a particular QueryObject's My_Geom_Multi query function.
   */
  static void Geom_Multi(void *data, int num_gid_entries, int num_lid_entries,
        int num_obj, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int num_dim,
        double *geom_vec, int *ierr);


  //M.M.W. need to add hierarchical query functions

};

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::QueryObject( Teuchos::RCP< const ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > graph,
	   Teuchos::RCP< const Isorropia::Tpetra::CostDescriber<LocalOrdinal, GlobalOrdinal, Node> > costs,
                                     int inputType) 
  : haveGraph_(true) ,
    graph_(graph),
    //matrix_(0),
    //coords_(0),
    rowMap_(/*&*/(graph->getRowMap().getRawPtr())), // EEP
    colMap_(/*&*/(graph->getColMap().getRawPtr())), // EEP
    costs_(costs),
    //weights_(0),
    input_type_(inputType) 
{
  myProc_ = graph->getComm()->getRank();
  base_ = rowMap_->getIndexBase();

  // If graph
  if (input_type_ == graph_input_)
  {

    // graph queries need to know processes owning my column entries
    fill_procmap();

    if (graph->getLocalNumRows()/*MyDiagonals()*/ > 0){ // EEP
      // In graph partitioning, we need to omit graph diagonal entries
      // (self edges) from the query functions.

      int nRows = rowMap_->getLocalNumElements(); // EEP
      int *rowGIDs = nullptr; // rowMap_->MyGlobalElements(); EEP___

      for (int i=0; i < nRows; i++){

	int numEntries;
	int *idx;

	//graph->ExtractMyRowView(i, numEntries, idx); // EEP___

	for (int j=0; j<numEntries; j++){
	  if (rowGIDs[i] == colMap_->getGlobalElement(idx[j])){ // EEP Epetra_Map::GID()
	    graph_self_edges_.insert(rowGIDs[i]);
	  }
	}
      }
    }
  }
}

#if 0 // EEP
QueryObject::QueryObject( Teuchos::RCP<const Epetra_RowMatrix> matrix,
	     Teuchos::RCP<const Isorropia::Epetra::CostDescriber> costs,
                                 int inputType) 
  : haveGraph_(false),
    graph_(0),
    matrix_(matrix),
    coords_(0),
    rowMap_((const Epetra_BlockMap*)&(matrix->RowMatrixRowMap())),
    colMap_((const Epetra_BlockMap*)&(matrix->RowMatrixColMap())),
    costs_(costs),
    weights_(0),
    input_type_(inputType) 
{
  myProc_ = matrix->Comm().MyPID();
  base_ = rowMap_->IndexBase();

  // If graph or hierarchical with graph and hypergraph 
  if (input_type_ == graph_input_ || input_type_ == hgraph_graph_input_)
  {
    // graph queries need to know processes owning my column entries
    fill_procmap();

    if (matrix->NumMyDiagonals() > 0){
      // In graph partitioning, we need to omit graph diagonal entries
      // (self edges) from the query functions.

      int *rowGIDs;
      int nRows;

      Epetra_Vector diagonal(matrix->RowMatrixRowMap());
      nRows = rowMap_->NumMyElements();
      rowGIDs = rowMap_->MyGlobalElements();
      matrix->ExtractDiagonalCopy(diagonal);

      for (int i=0; i < nRows; i++){
	if (diagonal[i] != 0){
	  graph_self_edges_.insert(rowGIDs[i]);
	}
      }
    }
  }
  else if(input_type_ != hgraph2d_finegrain_input_) // otherwise must be hypergraph
  {
    input_type_ = hgraph_input_;
  }
}

// For geometric partitioning
QueryObject::QueryObject( Teuchos::RCP<const Epetra_MultiVector> coords,
                          Teuchos::RCP<const Epetra_MultiVector> weights)
  : haveGraph_(false),
    graph_(0),
    matrix_(0),
    coords_(coords),
    rowMap_(&(coords->Map())),
    costs_(0),
    weights_(weights)
{
  myProc_ = rowMap_->Comm().MyPID();
  base_ = rowMap_->IndexBase();
  input_type_ = geometric_input_;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
QueryObject::QueryObject( Teuchos::RCP<const Epetra_BlockMap> input_map,
                                     int inputType) 
  : haveGraph_(false) ,
    graph_(0),
    matrix_(0),
    coords_(0),
    rowMap_(input_map.getRawPtr()), // MMW: this is dangerous but so are all of the above rowMap_ assignments
    costs_(0),
    weights_(0),
    input_type_(inputType) 
{
  myProc_ = rowMap_->Comm().MyPID();
  base_ = rowMap_->IndexBase();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// For hierarchical with graph or hgraph and geometric
////////////////////////////////////////////////////////////////////////////////
QueryObject::QueryObject(Teuchos::RCP<const Epetra_CrsGraph> graph,
	                 Teuchos::RCP<const Isorropia::Epetra::CostDescriber> costs, 
                         Teuchos::RCP<const Epetra_MultiVector> coords,
                         Teuchos::RCP<const Epetra_MultiVector> weights, int inputType)
: haveGraph_(true) ,
  graph_(graph),
  matrix_(0),
  coords_(coords),
  rowMap_(&(graph->RowMap())),
  colMap_(&(graph->ColMap())),
  costs_(costs),
  weights_(weights),
  input_type_(inputType) 
{
  myProc_ = graph->Comm().MyPID();
  base_ = rowMap_->IndexBase();

  // if graph
  if (input_type_ == graph_geometric_input_ ||
      input_type_ == hgraph_graph_geometric_input_)
  {

    // graph queries need to know processes owning my column entries
    fill_procmap();

    if (graph->NumMyDiagonals() > 0){
      // In graph partitioning, we need to omit graph diagonal entries
      // (self edges) from the query functions.

      int nRows = rowMap_->NumMyElements();
      int *rowGIDs = rowMap_->MyGlobalElements();

      for (int i=0; i < nRows; i++){

	int numEntries;
	int *idx;

	graph->ExtractMyRowView(i, numEntries, idx);

	for (int j=0; j<numEntries; j++){
	  if (rowGIDs[i] == colMap_->GID(idx[j])){
	    graph_self_edges_.insert(rowGIDs[i]);
	  }
	}
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// For hierarchical with graph or hgraph and geometric
////////////////////////////////////////////////////////////////////////////////
QueryObject::QueryObject(Teuchos::RCP<const Epetra_RowMatrix> matrix,
	                 Teuchos::RCP<const Isorropia::Epetra::CostDescriber> costs,
                         Teuchos::RCP<const Epetra_MultiVector> coords,
                         Teuchos::RCP<const Epetra_MultiVector> weights,
                         int inputType) 
  : haveGraph_(false),graph_(0), matrix_(matrix),coords_(coords),
    rowMap_((const Epetra_BlockMap*)&(matrix->RowMatrixRowMap())),
    colMap_((const Epetra_BlockMap*)&(matrix->RowMatrixColMap())),
    costs_(costs), weights_(weights), input_type_(inputType) 
{
  myProc_ = matrix->Comm().MyPID();
  base_ = rowMap_->IndexBase();

  // If graph or hierarchical with graph and hypergraph 
  if (input_type_ == graph_geometric_input_ || input_type_ == hgraph_graph_geometric_input_)
  {
    // graph queries need to know processes owning my column entries
    fill_procmap();

    if (matrix->NumMyDiagonals() > 0){
      // In graph partitioning, we need to omit graph diagonal entries
      // (self edges) from the query functions.

      int *rowGIDs;
      int nRows;

      Epetra_Vector diagonal(matrix->RowMatrixRowMap());
      nRows = rowMap_->NumMyElements();
      rowGIDs = rowMap_->MyGlobalElements();
      matrix->ExtractDiagonalCopy(diagonal);

      for (int i=0; i < nRows; i++){
	if (diagonal[i] != 0){
	  graph_self_edges_.insert(rowGIDs[i]);
	}
      }
    }
  }
}
#endif // EEP
////////////////////////////////////////////////////////////////////////////////


template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::~QueryObject()
{
}

// Create a map required by graph queries

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::fill_procmap()
{
  int num_colmap_elems = colMap_->getLocalNumElements();

  int *colIDs = new int [num_colmap_elems];
  int *procIDs = new int [num_colmap_elems];
  int *tmp = new int [num_colmap_elems];

#if 0 // EEP___
  colMap_->MyGlobalElements(colIDs);

  rowMap_->RemoteIDList(num_colmap_elems, colIDs, procIDs, tmp);
#endif
  
  delete [] tmp;

  int i;
  for(i=0; i<num_colmap_elems; ++i) {

    // map from global column ID to process owning the row
    // with that global ID (this is a square matrix)

    procmap_[colIDs[i]] = procIDs[i];
  }

  delete [] colIDs;
  delete [] procIDs;
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
bool
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::haveVertexWeights()
{
  return costs_->haveGlobalVertexWeights();
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
bool
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::haveGraphEdgeWeights()
{
  return costs_->haveGlobalGraphEdgeWeights();
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
bool
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::haveHypergraphEdgeWeights()
{
  return costs_->haveGlobalHypergraphEdgeWeights();
}

// Static query functions.  These will call the query function of the
// appropriate QueryObject object.

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::Number_Objects(void *data, int *ierr)
{
  int numObj = 0;

  QueryObject *zq = (QueryObject *)data;

  if (zq){
    numObj = zq->My_Number_Objects(ierr);
  }
  else{
    *ierr = ZOLTAN_FATAL;
  }

  return numObj;
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::Object_List  ( void * data,
		   int num_gid_entries, int num_lid_entries,
		   ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
		   int weight_dim, float * object_weights, int * ierr )
{
  QueryObject *zq = (QueryObject *)data;

  if (zq)
  {
    zq->My_Object_List(num_gid_entries, num_lid_entries,
		    global_ids, local_ids, weight_dim, object_weights, ierr);
  }
  else
  {
    *ierr = ZOLTAN_FATAL;
  }

  return;
}

#if 0 // EEP
void QueryObject::Number_Edges_Multi  ( void * data,
	     int num_gid_entries, int num_lid_entries,
	     int num_obj,
	     ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
	     int *num_edges, int * ierr )
{
  QueryObject *zq = (QueryObject *)data;

  if (zq){
    zq->My_Number_Edges_Multi(num_gid_entries, num_lid_entries, num_obj,
	     global_ids, local_ids, num_edges, ierr );
  }
  else{
    *ierr = ZOLTAN_FATAL;
  }

  return;
}
void QueryObject::Edge_List_Multi( void * data,
	       int num_gid_entries, int num_lid_entries, int num_obj,
	       ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *num_edges,
	       ZOLTAN_ID_PTR neighbor_global_ids, int * neighbor_procs,
	       int weight_dim, float * edge_weights, int * ierr )
{
  QueryObject *zq = (QueryObject *)data;

  if (zq){
    zq->My_Edge_List_Multi(num_gid_entries, num_lid_entries, num_obj,
	       global_ids, local_ids, num_edges,
	       neighbor_global_ids,  neighbor_procs,
	       weight_dim, edge_weights, ierr );
  }
  else{
    *ierr = ZOLTAN_FATAL;
  }

  return;
}
#endif // EEP

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::HG_Size_CS ( void * data,
	 int* num_lists, int* num_pins, int* format, int * ierr )
{
  QueryObject *zq = (QueryObject *)data;

  if (zq)
  {
    zq->My_HG_Size_CS(num_lists, num_pins, format, ierr );
  }
  else
  {
    *ierr = ZOLTAN_FATAL;
  }

}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::HG_CS ( void * data,
	    int num_gid_entries, int num_row_or_col, int num_pins, int format,
	    ZOLTAN_ID_PTR vtxedge_GID, int* vtxedge_ptr, ZOLTAN_ID_PTR pin_GID,
				     int * ierr )
{
  QueryObject *zq = (QueryObject *)data;

  if (zq){
    zq->My_HG_CS(num_gid_entries, num_row_or_col, num_pins, format,
	    vtxedge_GID, vtxedge_ptr, pin_GID, ierr );
  }
  else{
    *ierr = ZOLTAN_FATAL;
  }
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::HG_Size_Edge_Weights(void * data,
			    int* num_edges, int* ierr)
{
  QueryObject *zq = (QueryObject *)data;

  if (zq){
    zq->My_HG_Size_Edge_Weights(num_edges, ierr);
  }
  else{
    *ierr = ZOLTAN_FATAL;
  }
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::HG_Edge_Weights(void * data,
      int num_gid_entries, int num_lid_entries, int num_edges, int edge_weight_dim,
      ZOLTAN_ID_PTR edge_GID, ZOLTAN_ID_PTR edge_LID, float* edge_weights, int* ierr)
{
  QueryObject *zq = (QueryObject *)data;

  if (zq){
    zq->My_HG_Edge_Weights(num_gid_entries, num_lid_entries, num_edges, 
         edge_weight_dim, edge_GID, edge_LID, edge_weights, ierr);
  }
  else{
    *ierr = ZOLTAN_FATAL;
  }
}

#if 0 // EEP
int QueryObject::Number_Geom(void *data, int *ierr)
{
  int dim=0;
  QueryObject *zq = (QueryObject *)data;

  if (zq){
    dim = zq->My_Number_Geom(ierr);
  }
  else{
    *ierr = ZOLTAN_FATAL;
  }
  return dim;
}

void QueryObject::Geom_Multi(void *data, int num_gid_entries, int num_lid_entries,
        int num_obj, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int num_dim,
        double *geom_vec, int *ierr)
{
  QueryObject *zq = (QueryObject *)data;

  if (zq){
    zq->My_Geom_Multi(num_gid_entries, num_lid_entries,
        num_obj, gids, lids, num_dim, geom_vec, ierr);
  }
  else{
    *ierr = ZOLTAN_FATAL;
  }
}
#endif // EEP

// Member general query functions.

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::My_Number_Objects(int * ierr )
{
  *ierr = ZOLTAN_OK;
#if 0 // EEP
  if (input_type_ == geometric_input_) 
  {
    return coords_->MyLength(); // probably can use rowMap_ instead
  }
  else if (input_type_ == hgraph2d_finegrain_input_)
  {
    if (!haveGraph_)
    {
      return matrix_->NumMyNonzeros();
    }
    else
    {
      return graph_->NumMyNonzeros();
    }
  }
  else
#endif // EEP
  {
    return rowMap_->getLocalNumElements();  // graph or hypergraph or simple partitioning methods
  }
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::My_Object_List(int num_gid_entries, int num_lid_entries,
		  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
		  int weight_dim, float * object_weights, int * ierr )
{
  *ierr = ZOLTAN_OK;
  int ngids = 0;
  int *ibuf;

  //M.M.W. need to add hierarchical support

  if ((input_type_ == geometric_input_) || (input_type_ == hgraph_input_) || (input_type_ == graph_input_) || (input_type_==simple_input_)) 
  {
    ngids = rowMap_->getLocalNumElements();

    if (sizeof(ZOLTAN_ID_TYPE) != sizeof(int)){
      ibuf = nullptr; // rowMap_->MyGlobalElements(); // EEP___

      for (int i=0; i < ngids; i++){
        global_ids[i] = (ZOLTAN_ID_TYPE)ibuf[i];
      }
    }
    else{
      //rowMap_->MyGlobalElements( ((int *) global_ids) ); // EEP___
    }
  }
#if 0 // EEP
  else if(input_type_ == hgraph2d_finegrain_input_)
  {
    // if statement for now, maybe set function pointer for different function instead

    int numRows = rowMap_->NumMyElements();

    object_weights = 0;

    if (!haveGraph_) //matrix
    {
      ngids = 0;
      for (int rowNum=0; rowNum<numRows; rowNum++)
      {
        int rowSize;
        matrix_->NumMyRowEntries(rowNum,rowSize);

        //MMW need to make this memory allocation more efficient
        int *tmprowCols = new int[rowSize];
        double *tmprowVals = new double[rowSize];
        int numEntries;

        matrix_->ExtractMyRowCopy (rowNum, rowSize, numEntries, tmprowVals, tmprowCols);

        for(int colIndx=0; colIndx<rowSize; colIndx++)
        {
          global_ids[ngids] = (ZOLTAN_ID_TYPE)rowMap_->GID(rowNum);
          ngids++;
          global_ids[ngids]= (ZOLTAN_ID_TYPE)colMap_->GID(tmprowCols[colIndx]);
          ngids++;
        }

	delete [] tmprowCols;
	delete [] tmprowVals;
      }
    }
    else // graph    
    {
      ngids = 0;
      for (int rowNum=0; rowNum<numRows; rowNum++)
      {
        int rowSize;
        rowSize = graph_->NumMyIndices(rowNum);

        int *tmprowCols = new int[rowSize];
        int numEntries;

        graph_->ExtractMyRowCopy (rowNum, rowSize, numEntries, tmprowCols);

        for(int colIndx=0;colIndx<rowSize;colIndx++)
        {
          global_ids[ngids] = (ZOLTAN_ID_TYPE)rowMap_->GID(rowNum);
          ngids++;
          global_ids[ngids]= (ZOLTAN_ID_TYPE)colMap_->GID(tmprowCols[colIndx]);
          ngids++;
        }

        delete [] tmprowCols;
      }
    }

    for (int indx=0; indx<ngids/2; indx++)
    {
      local_ids[indx] = indx;
    }

    return;
  } // fine-grain hypergraph case
#endif // EEP

  if (ngids < 1)
  {
    return;
  }

  for (int i=0; i<ngids; i++){
    local_ids[i] = (ZOLTAN_ID_TYPE)i;
  }
#if 0 // EEP weight_dim is 0 here, during the epetra example run
  if (weight_dim >= 1) // Note we only supply 1-D weights
  {          
    float *to_wgts = object_weights;
    if (input_type_ == geometric_input_){
      double *wgts=NULL;
      int ld=0;
      if (weights_.get()){
        weights_->ExtractView(&wgts, &ld);
        if (ld < ngids){
          *ierr = ZOLTAN_FATAL;
          std::cerr << "Proc:" << myProc_ << " Error: ";
          std::cerr << "My_Object_List: not enough object weights" << std::endl;
          return;
        }
        for (int i=0; i < ngids; i++){
          *to_wgts= static_cast<float>(wgts[i]);
          to_wgts += weight_dim;
        }
      }
    }
    else if (costs_->haveVertexWeights())
    {
      std::map<int, float> weightMap;
      std::map<int, float>::iterator curr;
      std::map<int, float>::iterator end = weightMap.end();
  
      int mapSize = costs_->getVertexWeights(weightMap);
      if (mapSize != ngids)
      {
        *ierr = ZOLTAN_FATAL;
        std::cout << "Proc:" << myProc_ << " Error: ";
        std::cout << "QueryObject::My_Object_List, number of vertex weights" << std::endl;
        return;
      }
      for (int i=0; i < ngids; i++){
        curr = weightMap.find(global_ids[i]);
        if (curr == end){
          *ierr = ZOLTAN_FATAL;
          std::cout << "Proc:" << myProc_ << " Error: ";
          std::cout << "QueryObject::My_Object_List, missing vertex weight" << std::endl;
          return;
        }
        *to_wgts = curr->second;
        to_wgts += weight_dim;
      }
    }
  }
#endif // EEP
  return;
}

#if 0 // EEP
// member graph query functions

void QueryObject::My_Number_Edges_Multi(int num_gid_entries, int num_lid_entries,
	     int num_obj,
	     ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
	     int *num_edges, int * ierr )
{
  *ierr = ZOLTAN_OK;

  int num_indices, rc;

  if (!haveGraph_){
    for (int i=0; i<num_obj; i++){
      rc = matrix_->NumMyRowEntries((int)local_ids[i], num_indices);
      if (rc != 0){
	*ierr = ZOLTAN_FATAL;
	std::cout << "Proc:" << myProc_ << " Error: ";
	std::cout << "QueryObject::My_Number_Edges_Multi, local id not in matrix" << std::endl;
	return;
      }
      num_edges[i] = num_indices;
      if (graph_self_edges_.size() > 0){
	std::set<int>::iterator it = graph_self_edges_.find((int)global_ids[i]);
	if (it != graph_self_edges_.end()){
	  num_edges[i]--;
	}
      }
    }
  }
  else{
    for (int i=0; i<num_obj; i++){
      num_edges[i] = graph_->NumMyIndices((int)local_ids[i]);
      if (num_edges[i] < 0){
	*ierr = ZOLTAN_FATAL;
	std::cout << "Proc:" << myProc_ << " Error: ";
	std::cout << "QueryObject::My_Number_Edges_Multi, local id not in graph" << std::endl;
	return;
      }
      if (graph_self_edges_.size() > 0){
	std::set<int>::iterator it = graph_self_edges_.find((int)global_ids[i]);
	if (it != graph_self_edges_.end()){
	  num_edges[i]--;
	}
      }
    }
  }
  return;
}

void QueryObject::My_Edge_List_Multi(int num_gid_entries, int num_lid_entries, int num_obj,
	       ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *num_edges,
	       ZOLTAN_ID_PTR neighbor_global_ids, int * neighbor_procs,
	       int weight_dim, float * edge_weights, int * ierr )
{
  int *gids=NULL;
  int num_indices, rc, tmpSize;
  double *tmp = NULL;

  *ierr = ZOLTAN_OK;

  if (num_obj < 1) {
    return;
  }

  int *nborProc = neighbor_procs;
  float *wgt = edge_weights;
  ZOLTAN_ID_PTR nborID = neighbor_global_ids;
  bool use_weights = false;
  std::map<int, float> wgtMap;
  std::map<int, float>::iterator wgtIter;
  std::map<int, int>::iterator procIter;

  if ((weight_dim >= 1) && costs_->haveGraphEdgeWeights()){
    use_weights = true;
  }

  if (!haveGraph_)
    tmpSize = matrix_->MaxNumEntries();
  else
    tmpSize = graph_->MaxNumIndices();

  if (tmpSize > 0){
    if (!haveGraph_)
      tmp = new double [tmpSize];

    gids = new int [tmpSize];

    if ((!haveGraph_ && !tmp) || !gids){
      *ierr = ZOLTAN_MEMERR;
      return;
    }
  }

  for (int i=0; (i<num_obj) && (*ierr == ZOLTAN_OK); i++){

    int self_edge = 0;
    if (graph_self_edges_.size() > 0){
      std::set<int>::iterator it = graph_self_edges_.find((int)global_ids[i]);
      if (it != graph_self_edges_.end()){
	self_edge = 1;
      }
    }

    if (!haveGraph_){
      rc = matrix_->ExtractMyRowCopy((int)local_ids[i], tmpSize, num_indices, tmp, gids);
    }
    else{
      rc = graph_->ExtractMyRowCopy((int)local_ids[i], tmpSize,  num_indices, gids);
    }

    if ((rc < 0) || (num_indices != (num_edges[i] + self_edge))){
      *ierr = ZOLTAN_FATAL;
      std::cout << "Proc:" << myProc_ << " Error: ";
      std::cout << "QueryObject::My_Edge_List_Multi, extracting row" << std::endl;
      break;
    }

    for (int j=0; j<num_indices; j++){
      gids[j] = colMap_->GID(gids[j]); // convert to global IDs
    }

    if (use_weights){
      rc = costs_->getGraphEdgeWeights((int)global_ids[i], wgtMap);
      if (rc != num_edges[i]){
	*ierr = ZOLTAN_FATAL;
	std::cout << "Proc:" << myProc_ << " Error: ";
	std::cout << "QueryObject::My_Edge_List_Multi, getting weights" << std::endl;
	break;
      }
    }
    for (int j=0; j < num_indices; j++){

      if (self_edge && (gids[j] == (int)global_ids[i])) continue;  // skip self edge

      *nborID++ = (ZOLTAN_ID_TYPE)gids[j];                   // Global ID of neighbor

      procIter = procmap_.find(gids[j]);
      if (procIter == procmap_.end()){
	*ierr = ZOLTAN_FATAL;
	std::cout << "Proc:" << myProc_ << " Error: ";
	std::cout << "QueryObject::My_Edge_List_Multi, process ID map" << std::endl;
	break;
      }

      *nborProc++ = procIter->second;        // Process owning neighbor

      if (use_weights){
        wgtIter = wgtMap.find(gids[j]); 
        if (wgtIter == wgtMap.end()){
          *ierr = ZOLTAN_FATAL;
          std::cout << "Proc:" << myProc_ << " Error: ";
          std::cout << "QueryObject::My_Edge_List_Multi, weight map" << std::endl;
          break;
        }
        *wgt++ = (float)wgtIter->second;    // edge weight
      }
    }
  }
  if (gids) delete [] gids;
  if (tmp) delete [] tmp;

  return;
}
#endif // EEP

// member hypergraph query functions

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::My_HG_Size_CS(int* num_lists, int* num_pins, int* format, int * ierr )
{
  *ierr = ZOLTAN_OK;

  *format = ZOLTAN_COMPRESSED_VERTEX;   // We will return row (vertex) lists

#if 0 // EEP
  if(input_type_ == hgraph2d_finegrain_input_) // 2D fine-grain hypergraph
  {
    if (haveGraph_) // graph
    {
       *num_lists = graph_->NumMyNonzeros(); 
    }
    else // matrix
    {
      *num_lists = matrix_->NumMyNonzeros(); 
    }
    *num_pins = 2 * (*num_lists);
  }
  else // 1D hypergraph
  {
#endif // EEP
    if (haveGraph_)
    {
      *num_lists = graph_->getLocalNumRows();       // Number of rows
      *num_pins = graph_->getLocalNumEntries();    // Total nonzeros in these rows // EEP
    }
    else
    {
#if 0 // EEP
      *num_lists = matrix_->NumMyRows();       // Number of rows
      *num_pins = matrix_->NumMyNonzeros();    // Total nonzeros in these rows
#else
      throw std::runtime_error("Invalid situation in QueryObject<>::My_HG_Size_CS");
#endif // EEP
    }
    //}



  return;
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::My_HG_CS (int num_gid_entries, int num_row_or_col, int num_pins, int format,
	    ZOLTAN_ID_PTR vtxedge_GID, int* vtxedge_ptr, ZOLTAN_ID_PTR pin_GID,
				     int * ierr )
{
  int num_indices, rc, npins;
  double *tmp=NULL;


  *ierr = ZOLTAN_OK;
#if 0 // EEP
  /////////////////////////////////////////////////////
  // if fine-grain hypergraph, call this other function
  /////////////////////////////////////////////////////
  if(input_type_ == hgraph2d_finegrain_input_)
  {
    My_FGHG_CS(num_gid_entries, num_row_or_col, num_pins, format, 
               vtxedge_GID, vtxedge_ptr, pin_GID, ierr);
    return;
  }
  /////////////////////////////////////////////////////
#endif // EEP

  int maxrow = 0;

  if (haveGraph_)
  {
    npins = graph_->getLocalNumEntries(); // NumMyNonzeros(); // EEP
    maxrow = graph_->getLocalMaxNumRowEntries(); // MaxNumIndices(); // EEP
  }
#if 0 // EEP
  else
  {
    npins = matrix_->NumMyNonzeros();
    maxrow = matrix_->MaxNumEntries();
    if (maxrow > 0){
      tmp = new double [maxrow];
      if (!tmp){
	*ierr = ZOLTAN_MEMERR;
	return;
      }
    }
  }
#endif // EEP

  if ((format != ZOLTAN_COMPRESSED_VERTEX) ||
      (num_row_or_col != rowMap_->getLocalNumElements()) ||
      (num_pins != npins)){
    *ierr = ZOLTAN_FATAL;
    std::cout << "Proc:" << myProc_ << " Error: ";
    std::cout << "QueryObject::My_HG_CS, bad arguments" << std::endl;
    if (tmp) delete [] tmp;
    return;
  }

  int pin_start_pos = 0;
  unsigned int *gids = NULL;

  if (sizeof(ZOLTAN_ID_TYPE) != sizeof(unsigned int)){
    gids = new unsigned int [maxrow];
    if (maxrow && !gids){
      *ierr = ZOLTAN_MEMERR;
      return;
    }
  }
  else{
    gids = (unsigned int *)pin_GID;
  }

  for (int i=0; i<num_row_or_col; i++){ 
    vtxedge_GID[i] = (ZOLTAN_ID_TYPE)rowMap_->getGlobalElement(i);
    vtxedge_ptr[i] = pin_start_pos;

    if (sizeof(ZOLTAN_ID_TYPE) != sizeof(int)){
      if (haveGraph_){
        rc = 0; // graph_->ExtractMyRowCopy(i, maxrow, num_indices, (int *)gids); // EEP___
      }
#if 0 // EEP
      else{
        rc = matrix_->ExtractMyRowCopy(i, maxrow, num_indices, tmp,  (int *)gids);
      }
#endif // EEP
      if (rc == 0){
        for (int j=pin_start_pos, count=0; count < num_indices; j++, count++){
          pin_GID[j] = (ZOLTAN_ID_TYPE)colMap_->getGlobalElement(gids[count]);
          if (pin_GID[i] < base_){
     	    *ierr = ZOLTAN_FATAL;
    	    std::cout << "Proc:" << myProc_ << " Error: ";
    	    std::cout << "QueryObject::My_HG_CS, local ID not in column map" << std::endl;
    	    return;
          }
        }
      }
      else{
        *ierr = ZOLTAN_FATAL;
        std::cout << "Proc:" << myProc_ << " Error: ";
        std::cout << "QueryObject::My_HG_CS, extracting row" << std::endl;
        break;
      }
    }
    else{
      if (haveGraph_){
        rc = 0; // graph_->ExtractMyRowCopy(i, npins, num_indices,(int *)gids + pin_start_pos); // EEP___
      }
#if 0 // EEP
      else{
        rc = matrix_->ExtractMyRowCopy(i, npins, num_indices, tmp, (int *)gids + pin_start_pos);
      }
#endif // EEP

      if (rc == 0){
        for (int i=pin_start_pos; i<pin_start_pos+num_indices; i++){
          gids[i] = (unsigned int)colMap_->getGlobalElement(gids[i]); // convert to global IDs
          if (gids[i] < base_){
     	    *ierr = ZOLTAN_FATAL;
    	    std::cout << "Proc:" << myProc_ << " Error: ";
    	    std::cout << "QueryObject::My_HG_CS, local ID not in column map" << std::endl;
    	    return;
          }
        }
      }
      else{
        *ierr = ZOLTAN_FATAL;
        std::cout << "Proc:" << myProc_ << " Error: ";
        std::cout << "QueryObject::My_HG_CS, extracting row" << std::endl;
        break;
      }
    }

    pin_start_pos += num_indices;
    npins -= num_indices;
  }

  if (gids && (gids != (unsigned int *)pin_GID)){
    delete [] gids;
  }

  if (tmp)
  {
    delete [] tmp;
  }
}

#if 0 // EEP
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void QueryObject::My_FGHG_CS (int num_gid_entries, int num_nonzeros, int num_pins, int format,
	    ZOLTAN_ID_PTR vtx_GID, int* vtx_ptr, ZOLTAN_ID_PTR pin_GID,
				     int * ierr )
{
  int nverts, npins;

   *ierr = ZOLTAN_OK;

   //////////////////////////////////////////
   // Sanity checks of function input
   //////////////////////////////////////////
   if (haveGraph_)  // graph
   {
      nverts = graph_->NumMyNonzeros();
   }
   else   //matrix
   {
      nverts = matrix_->NumMyNonzeros();
   }
   npins = 2*nverts;

   if ((format != ZOLTAN_COMPRESSED_VERTEX) ||
       (num_nonzeros != nverts) || (num_pins != npins))
   {
      *ierr = ZOLTAN_FATAL;
      std::cout << "Proc:" << myProc_ << " Error: ";
      std::cout << "QueryObject::My_HG_CS, bad arguments" << std::endl;
      return;
   }
   //////////////////////////////////////////

   //////////////////////////////////////////
   // Allocate temporary memory needed for row extraction
   //////////////////////////////////////////
   int *tmprowCols=0;
   double *tmprowVals=0;

   int maxrow=0;

   if (haveGraph_)  // graph
   {
      maxrow = graph_->MaxNumIndices();

      if (maxrow > 0)
      {
         tmprowCols = new int [maxrow];
         if (!tmprowCols)
         {
 	   *ierr = ZOLTAN_MEMERR;
 	   return;
         }
      }
   }
   else   //matrix
   {
      maxrow = matrix_->MaxNumEntries();
      if (maxrow > 0)
      {
         tmprowCols = new int [maxrow];
         if (!tmprowCols)
         {
 	   *ierr = ZOLTAN_MEMERR;
 	   return;
         }
         tmprowVals = new double [maxrow];
         if (!tmprowVals)
         {
 	   *ierr = ZOLTAN_MEMERR;
 	   return;
         }
      }
   }
   //////////////////////////////////////////

   //////////////////////////////////////////
   int numRows = rowMap_->NumMyElements();
   int num_indices = 0;

   vtx_ptr[0] = 0;
   int nzindx=0;
   int rc=0;
   for (int rowi=0; rowi<numRows; rowi++)
   {

      if (haveGraph_)
      {
         rc = graph_->ExtractMyRowCopy(rowi, maxrow, num_indices, tmprowCols);
      }
      else
      {
         rc = matrix_->ExtractMyRowCopy(rowi, maxrow, num_indices, tmprowVals, tmprowCols);
      }

      if (rc != 0)
      {
         *ierr = ZOLTAN_FATAL;
         std::cout << "Proc:" << myProc_ << " Error: ";
         std::cout << "QueryObject::My_FGHG_CS, extracting row" << std::endl;
         return;
      }

      for(int indx=0; indx<num_indices; indx++)
      {
        vtx_GID[2*nzindx] = (ZOLTAN_ID_TYPE)rowMap_->GID(rowi);
        vtx_GID[2*nzindx+1] = (ZOLTAN_ID_TYPE)colMap_->GID(tmprowCols[indx]);

        // First pin for vertex (row pin)
        pin_GID[4*nzindx] = 0;
        pin_GID[4*nzindx+1] = (ZOLTAN_ID_TYPE)rowMap_->GID(rowi);

        // Second pin for vertex (column pin)
        pin_GID[4*nzindx+2] = 1;
        pin_GID[4*nzindx+3] = (ZOLTAN_ID_TYPE)colMap_->GID(tmprowCols[indx]);

        nzindx++;
        vtx_ptr[nzindx]= 2*nzindx;
      }

   }
   //////////////////////////////////////////


    if (tmprowCols) 
    {
      delete [] tmprowCols;
    }
    if (tmprowVals) 
    {
      delete [] tmprowVals;
    }
}
////////////////////////////////////////////////////////////////////////////////
#endif // EEP

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::My_HG_Size_Edge_Weights( int* num_edges, int* ierr)
{
  *num_edges = costs_->getNumHypergraphEdgeWeights();

  *ierr = ZOLTAN_OK;
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
QueryObject<LocalOrdinal, GlobalOrdinal, Node>::My_HG_Edge_Weights(
      int num_gid_entries, int num_lid_entries, int num_edges, int edge_weight_dim,
      ZOLTAN_ID_PTR edge_GID, ZOLTAN_ID_PTR edge_LID, float* edge_weights, int* ierr)
{
  *ierr = ZOLTAN_OK;

  if ((num_edges == 0) || (edge_weight_dim == 0)) return;

  if (num_edges != costs_->getNumHypergraphEdgeWeights()){
    *ierr = ZOLTAN_FATAL;
    std::cerr << "Proc:" << myProc_ << " Error: ";
    std::cerr << "QueryObject::My_HG_Edge_Weights, bad arguments" << std::endl;
    return;
  }

  if (num_lid_entries > 0){
    for (int i=0; i<num_edges; i++){
      edge_LID[i] = (ZOLTAN_ID_TYPE)i;
    }
  }

  if (sizeof(ZOLTAN_ID_TYPE) != sizeof(int)){
    int *ibuf = new int [num_edges];
    if (num_edges && !ibuf){
      *ierr = ZOLTAN_MEMERR;
      return;
    }
    costs_->getHypergraphEdgeWeights(num_edges, ibuf, edge_weights);
    for (int i=0; i < num_edges; i++){
      edge_GID[i] = (ZOLTAN_ID_TYPE)ibuf[i];
    }
    delete [] ibuf;
    ibuf = NULL;
  }
  else{
    costs_->getHypergraphEdgeWeights(num_edges, (int *)edge_GID, edge_weights);
  }
}

#if 0 // EEP
// member geometric query functions

int QueryObject::My_Number_Geom(int *ierr)
{
  int dim = 0;
  *ierr = ZOLTAN_FATAL;

  if (coords_.get()){
    dim = coords_->NumVectors();
    *ierr = ZOLTAN_OK;
  }
  else{
    std::cerr << "in My_Number_Geom, no MultiVector present" << std::endl;
  }

  return dim;
}

void QueryObject::My_Geom_Multi(int num_gid_entries, int num_lid_entries,
        int num_obj, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int num_dim,
        double *geom_vec, int *ierr)
{
  if ((num_obj < 1) || (num_dim < 1)){
    *ierr = ZOLTAN_OK;
    return;
  }

  if (coords_.get() && 
      (coords_->MyLength() >= num_obj) && 
      (coords_->NumVectors() >= num_dim)){

    double *coords;
    int stride;

    coords_->ExtractView(&coords, &stride);

    double **v = new double * [num_dim];
    for (int j=0; j < num_dim; j++){
      v[j] = coords + (stride * j);
    }

    for (int i=0; i<num_obj; i++){
      for (int j=0; j<num_dim; j++){
        *geom_vec++ = v[j][i];
      }
    }

    delete [] v;
    *ierr = ZOLTAN_OK;
  }
  else{
    std::cerr << "Error in My_Geom_Multi" << std::endl;
    *ierr = ZOLTAN_FATAL;
  }
}
#endif // EEP

} //namespace ZoltanLib
} //namespace Tpetra
} //namespace Isorropia

#endif //ISORROPIA_TPETRA_ZOLTAN_QUERYOBJECT_H

#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

