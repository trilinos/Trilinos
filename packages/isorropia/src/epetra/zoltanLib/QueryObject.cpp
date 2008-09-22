// @HEADER
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
// @HEADER

#include <QueryObject.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <Epetra_CrsGraph.h>
#include <Epetra_Comm.h>
#include <Epetra_BlockMap.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_RowMatrix.h>

#include <map>

// if non-zero, print out query function responses
//#define DEBUG_QUERIES 1
// process to print out responses, or -1 for all
//#define DEBUG_PROC -1

namespace Isorropia{

namespace Epetra {

namespace ZoltanLib{

QueryObject::QueryObject( Teuchos::RefCountPtr<const Epetra_CrsGraph> graph,
	   Teuchos::RefCountPtr<const Isorropia::Epetra::CostDescriber> costs,
				     const std::string &inputType)
  : graph_(graph),
    matrix_(0),
    rowMap_(&(graph->RowMap())),
    colMap_(&(graph->ColMap())),
    costs_(costs),
    inputType_(inputType),
    haveGraph_(true)
{
  myProc_ = graph->Comm().MyPID();
  base_ = rowMap_->IndexBase();
  int rc = 0;

  if (inputType == "GRAPH" )
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

	rc = graph->ExtractMyRowView(i, numEntries, idx);

	for (int j=0; j<numEntries; j++){
	  if (rowGIDs[i] == colMap_->GID(idx[j])){
	    graph_self_edges_.insert(rowGIDs[i]);
	  }
	}
      }
    }
  }
}

QueryObject::QueryObject( Teuchos::RefCountPtr<const Epetra_RowMatrix> matrix,
	     Teuchos::RefCountPtr<const Isorropia::Epetra::CostDescriber> costs,
				 const std::string &inputType)
  : graph_(0),
    matrix_(matrix),
    rowMap_((const Epetra_BlockMap*)&(matrix->RowMatrixRowMap())),
    colMap_((const Epetra_BlockMap*)&(matrix->RowMatrixColMap())),
    costs_(costs),
    inputType_(inputType),
    haveGraph_(false)
{
  myProc_ = matrix->Comm().MyPID();
  base_ = rowMap_->IndexBase();

  if (inputType == "GRAPH")
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

QueryObject::~QueryObject()
{
}

// Create a map required by graph queries

void QueryObject::fill_procmap()
{
  int num_colmap_elems = colMap_->NumMyElements();

  int *colIDs = new int [num_colmap_elems];
  int *procIDs = new int [num_colmap_elems];
  int *tmp = new int [num_colmap_elems];

  colMap_->MyGlobalElements(colIDs);

  rowMap_->RemoteIDList(num_colmap_elems, colIDs, procIDs, tmp);

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

bool QueryObject::haveVertexWeights()
{
  return costs_->haveGlobalVertexWeights();
}
bool QueryObject::haveGraphEdgeWeights()
{
  return costs_->haveGlobalGraphEdgeWeights();
}
bool QueryObject::haveHypergraphEdgeWeights()
{
  return costs_->haveGlobalHypergraphEdgeWeights();
}

// Static query functions.  These will call the query function of the
// appropriate QueryObject object.

int QueryObject::Number_Objects(void *data, int *ierr)
{
  int numObj = 0;

  QueryObject *zq = (QueryObject *)data;

  if (zq)
  {
    numObj = zq->My_Number_Objects(ierr);
  }
  else{
    *ierr = ZOLTAN_FATAL;
  }

  return numObj;
}
void QueryObject::Object_List  ( void * data,
		   int num_gid_entries, int num_lid_entries,
		   ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
		   int weight_dim, float * object_weights, int * ierr )
{
  QueryObject *zq = (QueryObject *)data;

  if (zq){
    zq->My_Object_List(num_gid_entries, num_lid_entries,
		    global_ids, local_ids, weight_dim, object_weights, ierr);
  }
  else{
    *ierr = ZOLTAN_FATAL;
  }

  return;
}
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
void QueryObject::HG_Size_CS ( void * data,
	 int* num_lists, int* num_pins, int* format, int * ierr )
{
  QueryObject *zq = (QueryObject *)data;

  if (zq){
    zq->My_HG_Size_CS(num_lists, num_pins, format, ierr );
  }
  else{
    *ierr = ZOLTAN_FATAL;
  }
}
void QueryObject::HG_CS ( void * data,
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
void QueryObject::HG_Size_Edge_Weights(void * data,
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
void QueryObject::HG_Edge_Weights(void * data,
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

// Member general query functions.

int QueryObject::My_Number_Objects(int * ierr )
{
  *ierr = ZOLTAN_OK;

#if DEBUG_QUERIES
  if ((DEBUG_PROC < 0) || (myProc_ == DEBUG_PROC)){
    std::cout << myProc_ << ": in My_Number_Objects, return " << rowMap_->NumMyElements() << std::endl;
  }
#endif

  // if statement for now, maybe set function pointer for different function instead
  if(inputType_ == "HGRAPH2D_FINEGRAIN")
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

  return rowMap_->NumMyElements();
}

void QueryObject::My_Object_List(int num_gid_entries, int num_lid_entries,
		  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
		  int weight_dim, float * object_weights, int * ierr )
{
  *ierr = ZOLTAN_OK;

  int numRows = rowMap_->NumMyElements();

  if (numRows < 1){
#if DEBUG_QUERIES
    if ((DEBUG_PROC < 0) || (myProc_ == DEBUG_PROC)){
      std::cout << myProc_ << ": in My_Object_List, return due to no objects" << std::endl;
    }
#endif
    return;
  }

  // if statement for now, maybe set function pointer for different function instead
  if(inputType_ == "HGRAPH2D_FINEGRAIN")
  {
    object_weights = 0;

    if (!haveGraph_) //matrix
    {
      int globIDindx = 0;
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
	  global_ids[globIDindx] = rowMap_->GID(rowNum);
	  globIDindx++;
	  global_ids[globIDindx]= colMap_->GID(tmprowCols[colIndx]);
	  globIDindx++;
	}

	delete [] tmprowCols;
	delete [] tmprowVals;
      }
    }
    else // graph
    {
      int globIDindx = 0;
      for (int rowNum=0; rowNum<numRows; rowNum++)
      {
	int rowSize;
	rowSize = graph_->NumMyIndices(rowNum);

	int *tmprowCols = new int[rowSize];
	int numEntries;

	graph_->ExtractMyRowCopy (rowNum, rowSize, numEntries, tmprowCols);

	for(int colIndx=0;colIndx<rowSize;colIndx++)
	{
	  global_ids[globIDindx] = rowMap_->GID(rowNum);
	  globIDindx++;
	  global_ids[globIDindx]= colMap_->GID(tmprowCols[colIndx]);
	  globIDindx++;
	}

	delete [] tmprowCols;
      }
    }
    return;
  } // fine-grain hypergraph case

  rowMap_->MyGlobalElements( ((int *) global_ids) );
  for (int i=0; i<numRows; i++)
  {
    local_ids[i] = i;
  }

  if ((weight_dim >= 1) && costs_->haveVertexWeights()){
    std::map<int, float> weightMap;
    std::map<int, float>::iterator curr;
    std::map<int, float>::iterator end = weightMap.end();

    int mapSize = costs_->getVertexWeights(weightMap);
    if (mapSize != numRows){
      *ierr = ZOLTAN_FATAL;
      std::cout << "Proc:" << myProc_ << " Error: ";
      std::cout << "QueryObject::My_Object_List, number of vertex weights" << std::endl;
      return;
    }
    for (int i=0; i < numRows; i++){
      curr = weightMap.find(global_ids[i]);
      if (curr == end){
	*ierr = ZOLTAN_FATAL;
	std::cout << "Proc:" << myProc_ << " Error: ";
	std::cout << "QueryObject::My_Object_List, missing vertex weight" << std::endl;
	return;
      }
      object_weights[i] = curr->second;
    }
  }
#if DEBUG_QUERIES
  if ((DEBUG_PROC < 0) || (myProc_ == DEBUG_PROC)){
    std::ostringstream msg;

    msg << myProc_ << ": in My_Object_List, num_obj " << numRows << ", weight_dim " << weight_dim;
    msg << std::endl;

    for (int i=0; i<numRows; i++){
      if (weight_dim){
	msg << "   obj gid " << global_ids[i] << ", weight " << object_weights[i] << std::endl;
      }
      else{
	msg << "   obj gid " << global_ids[i] << std::endl;
      }
    }
    std::string s = msg.str();
    std::cout << s;
  }
#endif
  return;
}

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
      rc = matrix_->NumMyRowEntries(local_ids[i], num_indices);
      if (rc != 0){
	*ierr = ZOLTAN_FATAL;
	std::cout << "Proc:" << myProc_ << " Error: ";
	std::cout << "QueryObject::My_Number_Edges_Multi, local id not in matrix" << std::endl;
	return;
      }
      num_edges[i] = num_indices;
      if (graph_self_edges_.size() > 0){
	std::set<int>::iterator it = graph_self_edges_.find(global_ids[i]);
	if (it != graph_self_edges_.end()){
	  num_edges[i]--;
	}
      }
    }
  }
  else{
    for (int i=0; i<num_obj; i++){
      num_edges[i] = graph_->NumMyIndices(local_ids[i]);
      if (num_edges[i] < 0){
	*ierr = ZOLTAN_FATAL;
	std::cout << "Proc:" << myProc_ << " Error: ";
	std::cout << "QueryObject::My_Number_Edges_Multi, local id not in graph" << std::endl;
	return;
      }
      if (graph_self_edges_.size() > 0){
	std::set<int>::iterator it = graph_self_edges_.find(global_ids[i]);
	if (it != graph_self_edges_.end()){
	  num_edges[i]--;
	}
      }
    }
  }
#if DEBUG_QUERIES
  if ((DEBUG_PROC < 0) || (myProc_ == DEBUG_PROC)){
    std::ostringstream msg;

    msg << myProc_ << ": in My_Number_Edges_Multi, num_objs " << num_obj << std::endl;
    for (int i=0; i<num_obj; i++){
      msg << "   obj gid " << global_ids[i] << ", num edges " << num_edges[i] << std::endl;
    }
    std::string s = msg.str();
    std::cout << s;
  }
#endif
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
#if DEBUG_QUERIES
    if ((DEBUG_PROC < 0) || (myProc_ == DEBUG_PROC)){
      std::cout << myProc_ << ": in My_Edge_List_Multi, return due to no objects" << std::endl;
    }
#endif
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
      std::set<int>::iterator it = graph_self_edges_.find(global_ids[i]);
      if (it != graph_self_edges_.end()){
	self_edge = 1;
      }
    }

    if (!haveGraph_){
      rc = matrix_->ExtractMyRowCopy(local_ids[i], tmpSize, num_indices, tmp, gids);
    }
    else{
      rc = graph_->ExtractMyRowCopy(local_ids[i], tmpSize,  num_indices, gids);
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
      rc = costs_->getGraphEdgeWeights(global_ids[i], wgtMap);
      if (rc != num_edges[i]){
	*ierr = ZOLTAN_FATAL;
	std::cout << "Proc:" << myProc_ << " Error: ";
	std::cout << "QueryObject::My_Edge_List_Multi, getting weights" << std::endl;
	break;
      }
    }
    for (int j=0; j < num_indices; j++){

      if (self_edge && (gids[j] == (int)global_ids[i])) continue;  // skip self edge

      *nborID++ = gids[j];                   // Global ID of neighbor

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
#if DEBUG_QUERIES
  if ((DEBUG_PROC < 0) || (myProc_ == DEBUG_PROC)){
    std::ostringstream msg;

    int k = 0;
    msg << myProc_ << ": in My_Edge_List_Multi, num_objs " << num_obj << std::endl;
    for (int i=0; i<num_obj; i++){
      msg << "   obj gid " << global_ids[i] << ", num edges " << num_edges[i] << std::endl;
      for (int j=0; j<num_edges[i]; j++,k++){
	msg << "      n'bor gid " << neighbor_global_ids[k] << ", n'bor proc " << neighbor_procs[k];
	for (int w=0; w < weight_dim ; w++) {
	  msg << " edge weight " << edge_weights[k*weight_dim + w];
	}
	msg << std::endl;
      }
    }
    std::string s = msg.str();
    std::cout << s;
  }
#endif
  if (gids) delete [] gids;
  if (tmp) delete [] tmp;

  return;
}

// member hypergraph query functions

void QueryObject::My_HG_Size_CS(int* num_lists, int* num_pins, int* format, int * ierr )
{
  *ierr = ZOLTAN_OK;

  *format = ZOLTAN_COMPRESSED_VERTEX;   // We will return row (vertex) lists

  if(inputType_ == "HGRAPH2D_FINEGRAIN")
  {
    if (!haveGraph_)
    {
      *num_lists = matrix_->NumMyNonzeros();
      *num_pins = 2*matrix_->NumMyNonzeros();
    }
    else
    {
      *num_lists = graph_->NumMyNonzeros();
      *num_pins = 2*graph_->NumMyNonzeros();
    }
    return;
  }

  if (haveGraph_)
  {
    *num_lists = graph_->NumMyRows();       // Number of rows
    *num_pins = graph_->NumMyNonzeros();    // Total nonzeros in these rows
  }
  else
  {
    *num_lists = matrix_->NumMyRows();       // Number of rows
    *num_pins = matrix_->NumMyNonzeros();    // Total nonzeros in these rows
  }

#if DEBUG_QUERIES
  if (myProc_ == DEBUG_PROC){
    std::cout << "in My_HG_Size_CS" << std::endl;
    std::cout << "  num_lists and num_pins " << *num_lists << " " << *num_pins << std::endl;
  }
#endif

  return;
}

void QueryObject::My_HG_CS (int num_gid_entries, int num_row_or_col, int num_pins, int format,
	    ZOLTAN_ID_PTR vtxedge_GID, int* vtxedge_ptr, ZOLTAN_ID_PTR pin_GID,
				     int * ierr )
{
  int num_indices, rc, npins;
  double *tmp=NULL;


  *ierr = ZOLTAN_OK;


  /////////////////////////////////////////////////////////////
  if(inputType_ == "HGRAPH2D_FINEGRAIN")
  {
    int vIdx=0;
    vtxedge_ptr[0] = 0;
    int nzCnt = 1;
    int gNumRows;

    if ((format != ZOLTAN_COMPRESSED_VERTEX) ||  (num_gid_entries != 2) )
    {
      *ierr = ZOLTAN_FATAL;
      std::cout << "Proc:" << myProc_ << " Error: ";
      std::cout << "QueryObject::My_HG_CS, bad arguments" << std::endl;
      return;
    }


    int maxrow;
    int *tmpCols;

    if (haveGraph_)
    {
      gNumRows = graph_->NumGlobalRows();
      maxrow=graph_->MaxNumIndices();
    }
    else
    {
      gNumRows = matrix_->NumGlobalRows();
      maxrow = matrix_->MaxNumEntries();
      tmp = new double [maxrow];
      if (!tmp)
      {
	*ierr = ZOLTAN_MEMERR;
	return;
      }
    }

    tmpCols = new int[maxrow];
    if (!tmpCols)
    {
      *ierr = ZOLTAN_MEMERR;
      return;
    }

    for (int rownum=0; rownum<num_row_or_col; rownum++)
    {

      if (haveGraph_)
      {
	rc = graph_->ExtractMyRowCopy(rownum, num_pins, num_indices, tmpCols);
      }
      else
      {
	rc = matrix_->ExtractMyRowCopy(rownum, num_pins, num_indices, tmp, tmpCols);
      }
      if (rc != 0)
      {
	*ierr = ZOLTAN_FATAL;
	std::cout << "Proc:" << myProc_ << " Error: ";
	std::cout << "QueryObject::My_HG_CS, extracting row" << std::endl;
	break;
      }

      for(int nzIdx=0;nzIdx < num_indices; nzIdx++)
      {
	vtxedge_GID[vIdx] = rowMap_->GID(rownum);
	vIdx++;
	vtxedge_GID[vIdx] = colMap_->GID(tmpCols[nzIdx]);
	vIdx++;

	// Each vertex belongs to exactly 2 hyperedges.
	// Row hyperedge
	vtxedge_GID[vtxedge_ptr[nzCnt]] = rowMap_->GID(rownum);
	// Column hyperedge
	vtxedge_GID[vtxedge_ptr[nzCnt]+1] = gNumRows + colMap_->GID(tmpCols[nzIdx]);

	vtxedge_ptr[nzCnt] = vtxedge_ptr[nzCnt-1] + 2;
	nzCnt++;
      }

    }

    if (tmp) delete [] tmp;
    if (tmpCols) delete [] tmpCols;

    return;

  }
  /////////////////////////////////////////////////////////////



  if (haveGraph_)
  {
    npins = graph_->NumMyNonzeros();
  }
  else
  {
    npins = matrix_->NumMyNonzeros();
    int maxrow = matrix_->MaxNumEntries();
    if (maxrow > 0){
      tmp = new double [maxrow];
      if (!tmp){
	*ierr = ZOLTAN_MEMERR;
	return;
      }
    }
  }

  if ((format != ZOLTAN_COMPRESSED_VERTEX) ||
      (num_row_or_col != rowMap_->NumMyElements()) ||
      (num_pins != npins))
  {
    *ierr = ZOLTAN_FATAL;
    std::cout << "Proc:" << myProc_ << " Error: ";
    std::cout << "QueryObject::My_HG_CS, bad arguments" << std::endl;
    if (tmp) delete [] tmp;
    return;
  }

  int pin_start_pos = 0;
  int *gids = (int *)pin_GID;

  for (int i=0; i<num_row_or_col; i++)
  {
    vtxedge_GID[i] = rowMap_->GID(i);
    vtxedge_ptr[i] = pin_start_pos;

    if (haveGraph_)
    {
      rc = graph_->ExtractMyRowCopy(i, npins, num_indices, gids + pin_start_pos);
    }
    else{
      rc = matrix_->ExtractMyRowCopy(i, npins, num_indices, tmp, gids + pin_start_pos);
    }

    if (rc != 0){
      *ierr = ZOLTAN_FATAL;
      std::cout << "Proc:" << myProc_ << " Error: ";
      std::cout << "QueryObject::My_HG_CS, extracting row" << std::endl;
      break;
    }

    for (int i=pin_start_pos; i<pin_start_pos+num_indices; i++){
      gids[i] = colMap_->GID(gids[i]); // convert to global IDs
      if (gids[i] < base_){
	*ierr = ZOLTAN_FATAL;
	std::cout << "Proc:" << myProc_ << " Error: ";
	std::cout << "QueryObject::My_HG_CS, local ID not in column map" << std::endl;
	break;
      }
    }

    pin_start_pos += num_indices;
    npins -= num_indices;
  }
#if DEBUG_QUERIES
  if (myProc_ == DEBUG_PROC){
    std::cout << "in My_HG_CS" << std::endl;
    std::cout << "  num rows and num pins " << num_row_or_col << " " << num_pins << std::endl;
    int idx = 0;
    int len = 0;
    for (int i=0; i < num_row_or_col; i++){
      if (i == (num_row_or_col - 1)) len = num_pins - idx;
      else                           len = vtxedge_ptr[i+1] - vtxedge_ptr[i];

      std::cout <<   "  vtx/row " << vtxedge_GID[i] << std::endl;
      std::cout <<   "    " ;
      for (int j=0; j < len; j++){
	std::cout << pin_GID[idx+j] << " ";
      }
      std::cout << std::endl;
      idx += len;
    }
  }
#endif

  if (tmp)
  {
    delete [] tmp;
  }
}

void QueryObject::My_HG_Size_Edge_Weights( int* num_edges, int* ierr)
{
  *num_edges = costs_->getNumHypergraphEdgeWeights();

#if DEBUG_QUERIES
 if (myProc_ == DEBUG_PROC){
   std::cout << "in My_HG_Size_Edge_Weights " << *num_edges << std::endl;
 }
#endif
  *ierr = ZOLTAN_OK;
}

void QueryObject::My_HG_Edge_Weights(
      int num_gid_entries, int num_lid_entries, int num_edges, int edge_weight_dim,
      ZOLTAN_ID_PTR edge_GID, ZOLTAN_ID_PTR edge_LID, float* edge_weights, int* ierr)
{
  *ierr = ZOLTAN_OK;

  if ((num_edges == 0) || (edge_weight_dim == 0)) return;

  if (num_edges != costs_->getNumHypergraphEdgeWeights()){
    *ierr = ZOLTAN_FATAL;
    std::cout << "Proc:" << myProc_ << " Error: ";
    std::cout << "QueryObject::My_HG_Edge_Weights, bad arguments" << std::endl;
    return;
  }

  if (num_lid_entries > 0){
    for (int i=0; i<num_edges; i++){
      edge_LID[i] = i;
    }
  }

#if DEBUG_QUERIES
  if (myProc_ == DEBUG_PROC){
    std::cout << "in My_HG_Edge_Weights, num edges " << num_edges << std::endl;
    for (int j=0; j < num_edges; j++)
      std::cout << "Edge " << edge_GID[j] << " (local " << edge_LID[j] << " weight " << edge_weights[j] << std::endl;
}
#endif

  costs_->getHypergraphEdgeWeights(num_edges, (int *)edge_GID, edge_weights);
}

} // end namespace ZoltanLib
} // end namespace Epetra
} // end namespace Isorropia
