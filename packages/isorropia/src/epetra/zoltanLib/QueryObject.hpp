//@HEADER
// ***********************************************************************
// 
//            Isorropia: Partitioning and Load Balancing Package
//              Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// ***********************************************************************
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
 */
class QueryObject
{
  /** The CrsGraph.  If supplied, we will perform graph or hypergraph partitioning.
    */
  Teuchos::RCP<const Epetra_CrsGraph> graph_;

  /** The CrsMatrix.  If supplied, we will perform graph or hypergraph partitioning.
    */
  Teuchos::RCP<const Epetra_RowMatrix> matrix_;

  /** The MultiVector containing 1, 2 or 3 dimensional coordinates.  If
      supplied, we will perform geometric partitioning.
    */
  Teuchos::RCP<const Epetra_MultiVector> coords_;

  /** The graph or matrix row map
    */
  const Epetra_BlockMap *rowMap_;

  /** The graph or matrix column map
    */
  const Epetra_BlockMap *colMap_;

  /** The multivector map
    */
  const Epetra_BlockMap *map_;

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

  /** haveGraph is true if we have CrsGraph, and not a CrsMatrix or
      a MultiVector.
    */
  const bool haveGraph_;

  int myProc_;
  int base_;
  void fill_procmap();

  /** My_Number_Objects() returns the number of rows currently
      assigned to this process.  (The rows are interpreted as
      graph vertices for Graph partitioning, and as hypergraph
      vertices for hypergraph partitioning.)
   */
  int My_Number_Objects(int *ierr);

  /** My_ObjectList() returns to Zoltan the global ID and weight of the
      rows currently assigned to this process.
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

  /** to do
   */
  void My_HG_Size_CS (int* num_lists, int* num_pins, int* format, 
                          int * ierr );
  /** to do
   */
  void My_HG_CS (int num_gid_entries, int num_row_or_col, int num_pins, 
           int format, ZOLTAN_ID_PTR vtxedge_GID, int* vtxedge_ptr, ZOLTAN_ID_PTR pin_GID,
                                       int * ierr );
  /** to do
   */
  void My_HG_Size_Edge_Weights(int* num_edges, int* ierr);
  /** to do
   */
  void My_HG_Edge_Weights(int num_gid_entries, int num_lid_entries, int num_edges, int edge_weight_dim,
        ZOLTAN_ID_PTR edge_GID, ZOLTAN_ID_PTR edge_LID, float* edge_weights, int* ierr);

  int My_Number_Geom(int *ierr);

  void My_Geom_Multi(int num_gid_entries, int num_lid_entries,
        int num_obj, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int num_dim,
        double *geom_vec, int *ierr);

  /** Object which queries describe - does it represent a graph,
      a hypergraph, or a collection of geometric coordinates?
    */

  static const int hypergraph_input_ = 1;
  static const int graph_input_ = 2;
  static const int geometric_input_ = 3;
  static const int unspecified_input_ = 4;

  int input_type_;

 public:

  /** Constructor
   */
  QueryObject( Teuchos::RCP<const Epetra_CrsGraph> graph,
	       Teuchos::RCP<const Isorropia::Epetra::CostDescriber> costs,
               bool isHypergraph);


  /** Constructor
   */
  QueryObject( Teuchos::RCP<const Epetra_RowMatrix> matrix,
	       Teuchos::RCP<const Isorropia::Epetra::CostDescriber> costs,
               bool isHypergraph);

  /** Constructor
   */
  QueryObject( Teuchos::RCP<const Epetra_MultiVector> coords,
	       Teuchos::RCP<const Epetra_MultiVector> weights);

  /** Destructor
   */
  virtual ~QueryObject();

  /** to do
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
  /** to do
   */
  static int Number_Objects(void *data, int *ierr);
  /** to do
   */
  static void Object_List  ( void * data, int num_gid_entries, int num_lid_entries,
                     ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
                     int weight_dim, float * object_weights, int * ierr );

  // Query functions for graph partitioning only
  /** to do
   */
  static void Number_Edges_Multi  ( void * data, int num_gid_entries, int num_lid_entries,
               int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
               int *num_edges, int * ierr );
  /** to do
   */
  static void Edge_List_Multi( void * data, int num_gid_entries, int num_lid_entries, 
               int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, 
               int *num_edges, ZOLTAN_ID_PTR neighbor_global_ids, int * neighbor_procs,
               int weight_dim, float * edge_weights, int * ierr );

  // Query functions for hypergraph partitioning only
  /** to do
   */
  static void HG_Size_CS ( void * data, int* num_lists, int* num_pins, int* format, 
                          int * ierr );
  /** to do
   */
  static void HG_CS ( void * data, int num_gid_entries, int num_row_or_col, int num_pins, 
           int format, ZOLTAN_ID_PTR vtxedge_GID, int* vtxedge_ptr, ZOLTAN_ID_PTR pin_GID,
                                       int * ierr );
  /** to do
   */
  static void HG_Size_Edge_Weights(void * data, int* num_edges, int* ierr);
  /** to do
   */
  static void HG_Edge_Weights(void * data,
        int num_gid_entries, int num_lid_entries, int num_edges, int edge_weight_dim,
        ZOLTAN_ID_PTR edge_GID, ZOLTAN_ID_PTR edge_LID, float* edge_weights, int* ierr);

  static int Number_Geom(void *data, int *ierr);

  static void Geom_Multi(void *data, int num_gid_entries, int num_lid_entries,
        int num_obj, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int num_dim,
        double *geom_vec, int *ierr);

};

} //namespace ZoltanLib
} //namespace Epetra
} //namespace Isorropia

#endif //ISORROPIA_EPETRA_ZOLTAN_QUERYOBJECT_H

