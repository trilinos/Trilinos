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

#ifndef ISORROPIA_EPETRA_ZOLTAN_QUERYOBJECT_H
#define ISORROPIA_EPETRA_ZOLTAN_QUERYOBJECT_H

#include "Isorropia_ConfigDefs.hpp"

#include <Teuchos_RCP.hpp>

#include <zoltan_cpp.h>

#include <set>
#include <map>

class Epetra_BlockMap;
class Epetra_CrsGraph;
class Epetra_RowMatrix;
class Epetra_MultiVector;

namespace Isorropia {

namespace Epetra {

  class CostDescriber;

/** The ZoltanLib namespace within the Epetra namespace contains the
    classes and functions that use the Zoltan library to partition an
    Epetra object.
*/



namespace ZoltanLib {

/** QueryObject is a class that contains the query functions required
    by the Zoltan library.

    These methods are not part of the Isorropia API (except to Zoltan).
    They are called by Isorropia itself and by Zoltan.

    For a better understanding of Zoltan's query functions, see the
    Zoltan User's Guide at http://www.cs.sandia.gov/zoltan
 */
class QueryObject
{

  /** haveGraph is true if we have CrsGraph, and not a CrsMatrix or
      a MultiVector.
    */
  const bool haveGraph_;

  /** The CrsGraph.  The QueryObject must be constructed with one of
      an Epetra_CrsGraph, an Epetra_RowMatrix or an Epetra_MultiVector.
    */
  Teuchos::RCP<const Epetra_CrsGraph> graph_;

  /** The CrsMatrix. 
      The QueryObject must be constructed with one of
      an Epetra_CrsGraph, an Epetra_RowMatrix or an Epetra_MultiVector.
    */
  Teuchos::RCP<const Epetra_RowMatrix> matrix_;

  /** The MultiVector containing 1, 2 or 3 dimensional coordinates.  If
      supplied, we will perform geometric partitioning.
      The QueryObject must be constructed with one of
      an Epetra_CrsGraph, an Epetra_RowMatrix or an Epetra_MultiVector.
    */
  Teuchos::RCP<const Epetra_MultiVector> coords_;

  /** The graph or matrix row map, or the MultiVector map
    */
  const Epetra_BlockMap *rowMap_;

  /** The graph or matrix column map
    */
  const Epetra_BlockMap *colMap_;

  /** The CostDescriber contains optional vertex and/or edge weights for
      graph and hypergraph partitioning.
    */

  Teuchos::RCP<const Isorropia::Epetra::CostDescriber> costs_;

  /** The MultiVector contains optional object (point) weights for
      geometric partitioning.  Zoltan currently will use only the
      weights in the first vector (1 dimensional point weights).
    */

  Teuchos::RCP<const Epetra_MultiVector> weights_;

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
  QueryObject( Teuchos::RCP<const Epetra_CrsGraph> graph,
	       Teuchos::RCP<const Isorropia::Epetra::CostDescriber> costs,
	       int inputType);


  /** Constructor
   */
  QueryObject( Teuchos::RCP<const Epetra_RowMatrix> matrix,
	       Teuchos::RCP<const Isorropia::Epetra::CostDescriber> costs,
	       int inputType);

  /** Constructor
   */
  QueryObject( Teuchos::RCP<const Epetra_MultiVector> coords,
               Teuchos::RCP<const Epetra_MultiVector> weights);


  /** Constructor
   */
  QueryObject(Teuchos::RCP<const Epetra_CrsGraph> graph,
	      Teuchos::RCP<const Isorropia::Epetra::CostDescriber> costs, 
              Teuchos::RCP<const Epetra_MultiVector> coords,
              Teuchos::RCP<const Epetra_MultiVector> weights,
              int inputType);


  /** Constructor
   */
  QueryObject(Teuchos::RCP<const Epetra_RowMatrix> matrix,
	      Teuchos::RCP<const Isorropia::Epetra::CostDescriber> costs,
              Teuchos::RCP<const Epetra_MultiVector> coords,
              Teuchos::RCP<const Epetra_MultiVector> weights,
              int inputType);


  /** Constructor
   */
  QueryObject( Teuchos::RCP<const Epetra_BlockMap> inputMap,
	       int inputType);




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
      This indicates that the Epetra_MultiVector represents geometric
      coordinates.  The MultiVector should have 1, 2 or 3 vectors,
      representing 1, 2 or 3 dimensional coordinates.  The coordinates
      are to be partitioned.
    */
  static const int geometric_input_ = 4;

  /** input_type_ == hgraph_graph_input_
      This indicates that the Epetra_MultiVector represents a hypergraph
      and graph (see above).  This is necessary for hierarchical partitioning
      with both hypergraph and graph methods.
  */
  static const int hgraph_graph_input_ = 5;

  /** input_type_ == hgraph_geometric_input_
      This indicates that the Epetra_MultiVector represents a hypergraph
      and graph (see above).  This is necessary for hierarchical partitioning
      with both hypergraph and geometric methods.
  */
  static const int hgraph_geometric_input_ = 6;

  /** input_type_ == graph_geometric_input_
      This indicates that the Epetra_MultiVector represents graph and
      has geometric coordinates (see above).  
      This is necessary for hierarchical partitioning                                                                                                            
      with both graph and geometric methods.
  */
  static const int graph_geometric_input_ = 7;

  /** input_type_ == hgraph_graph_geometric_input_
      This indicates that the Epetra_MultiVector represents a hypergraph
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
  const Epetra_BlockMap &RowMap(void){ return *rowMap_;};

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

} //namespace ZoltanLib
} //namespace Epetra
} //namespace Isorropia

#endif //ISORROPIA_EPETRA_ZOLTAN_QUERYOBJECT_H

#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

