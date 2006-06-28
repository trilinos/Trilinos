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
//Questions? Contact Alan Williams (william@sandia.gov)
//                or Erik Boman    (egboman@sandia.gov)
// 
// ***********************************************************************
// @HEADER

#include <Isorropia_ZoltanQuery.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_BlockMap.h>
#include <Epetra_Comm.h>

#include <algorithm>
// #include <fstream>
// #include <sstream>

Isorropia::ZoltanQuery::ZoltanQuery( Teuchos::RefCountPtr<const Epetra_CrsGraph> graph,
                                     Teuchos::RefCountPtr<const Epetra_CrsGraph> tgraph,
                                     bool localEdgesOnly )
  : matrix_(),
    tmatrix_(),
    graph_(graph),
    tgraph_(tgraph),
    map_(&(graph->RowMap())),
  localEdgesOnly_(localEdgesOnly)
{
  int localProc = map_->Comm().MyPID();

  const Epetra_BlockMap& colmap = graph_->ColMap();
  int num_colmap_elems = colmap.NumMyElements();

  std::vector<int> elems(num_colmap_elems*2);
  colmap.MyGlobalElements(&elems[0]);
  map_->RemoteIDList(num_colmap_elems, &elems[0],
		      &elems[num_colmap_elems], 0);

  int i;
  for(i=0; i<num_colmap_elems; ++i) {
    int proc = elems[num_colmap_elems+i];
    if (proc != localProc) {
      procmap_.insert(std::pair<int,int>(elems[i], proc));
    }
  }

  int num_local_rows = map_->NumMyElements();
  ugraph_.resize(num_local_rows);
  int maxNumIndices = graph_->MaxNumIndices();
  std::vector<int> indices(maxNumIndices);
  int numIndices;

  for(i=0; i<num_local_rows; ++i) {
    std::set<int>& rowset = ugraph_[i];
    int row = map_->GID(i);
    graph_->ExtractGlobalRowCopy(row, maxNumIndices,
				numIndices, &indices[0]);
    for(int j=0; j<numIndices; ++j) {
      if (indices[j] != row) {
	rowset.insert(indices[j]);
      }
    }
  }

  if (tgraph.get() != 0) {
    tmap_ = &(tgraph->RowMap());
    const Epetra_BlockMap& tcolmap = tgraph->ColMap();
    int num_tcolmap_elems = tcolmap.NumMyElements();

    elems.resize(num_tcolmap_elems*2);

    tcolmap.MyGlobalElements(&elems[0]);
    map_->RemoteIDList(num_tcolmap_elems, &elems[0],
		       &elems[num_tcolmap_elems], 0);

    std::map<int,int>::iterator iter_end = procmap_.end();

    for(i=0; i<num_tcolmap_elems; ++i) {
      int elem = elems[i];
      int proc = elems[num_tcolmap_elems+i];
      if (proc != localProc) {
	std::map<int,int>::iterator iter = procmap_.find(elem);
	if (iter == iter_end) {
	  procmap_.insert(std::pair<int,int>(elem, proc));
	}
      }
    }

    maxNumIndices = tgraph_->MaxNumIndices();
    indices.resize(maxNumIndices);

    for(i=0; i<num_local_rows; ++i) {
      std::set<int>& rowset = ugraph_[i];
      int row = tmap_->GID(i);
      tgraph_->ExtractGlobalRowCopy(row, maxNumIndices,
				    numIndices, &indices[0]);
      for(int j=0; j<numIndices; ++j) {
	if (indices[j] != row) {
	  rowset.insert(indices[j]);
	}
      }
    }
  }
}

Isorropia::ZoltanQuery::ZoltanQuery( Teuchos::RefCountPtr<const Epetra_RowMatrix> matrix,
                                     Teuchos::RefCountPtr<const Epetra_RowMatrix> tmatrix,
                                     bool localEdgesOnly )
  : matrix_(matrix),
    tmatrix_(tmatrix),
    graph_(),
    tgraph_(),
    map_((const Epetra_BlockMap*)&(matrix->RowMatrixRowMap())),
    localEdgesOnly_(localEdgesOnly)
{
  int localProc = map_->Comm().MyPID();

  const Epetra_BlockMap& colmap = (const Epetra_BlockMap&)matrix_->RowMatrixColMap();
  int num_colmap_elems = colmap.NumMyElements();

  std::vector<int> elems(num_colmap_elems*2);
  colmap.MyGlobalElements(&elems[0]);
  map_->RemoteIDList(num_colmap_elems, &elems[0],
		     &elems[num_colmap_elems], 0);

  int i;
  for(i=0; i<num_colmap_elems; ++i) {
    int proc = elems[num_colmap_elems+i];
    if (proc != localProc) {
      procmap_.insert(std::pair<int,int>(elems[i], proc));
    }
  }

  int num_local_rows = map_->NumMyElements();
  ugraph_.resize(num_local_rows);
  int maxNumIndices = matrix_->MaxNumEntries();
  std::vector<int> indices(maxNumIndices);
  std::vector<double> coefs(maxNumIndices);
  int numIndices;

  for(i=0; i<num_local_rows; ++i) {
    std::set<int>& rowset = ugraph_[i];
    int row = map_->GID(i);
    matrix_->ExtractMyRowCopy(i, maxNumIndices,
			      numIndices, &coefs[0], &indices[0]);
    for(int j=0; j<numIndices; ++j) {
      int globalj = colmap.GID(indices[j]);
      if (globalj != row) {
	rowset.insert(globalj);
      }
    }
  }

  if (tmatrix.get() != 0) {
    tmap_ = (const Epetra_BlockMap*)&(tmatrix->RowMatrixRowMap());
    const Epetra_BlockMap& tcolmap =
      (const Epetra_BlockMap&)tmatrix->RowMatrixColMap();
    int num_tcolmap_elems = tcolmap.NumMyElements();

    elems.resize(num_tcolmap_elems*2);

    tcolmap.MyGlobalElements(&elems[0]);
    map_->RemoteIDList(num_tcolmap_elems, &elems[0],
		       &elems[num_tcolmap_elems], 0);

    std::map<int,int>::iterator iter_end = procmap_.end();

    for(i=0; i<num_tcolmap_elems; ++i) {
      int elem = elems[i];
      int proc = elems[num_tcolmap_elems+i];
      if (proc != localProc) {
	std::map<int,int>::iterator iter = procmap_.find(elem);
	if (iter == iter_end) {
	  procmap_.insert(std::pair<int,int>(elem, proc));
	}
      }
    }

    maxNumIndices = tmatrix_->MaxNumEntries();
    indices.resize(maxNumIndices);
    coefs.resize(maxNumIndices);

    for(i=0; i<num_local_rows; ++i) {
      std::set<int>& rowset = ugraph_[i];
      int row = tmap_->GID(i);
      tmatrix_->ExtractMyRowCopy(i, maxNumIndices,
				 numIndices, &coefs[0], &indices[0]);
      for(int j=0; j<numIndices; ++j) {
	int globalj = tcolmap.GID(indices[j]);
	if (globalj != row) {
	  rowset.insert(globalj);
	}
      }
    }
  }
}

//General Functions
int Isorropia::ZoltanQuery::Number_Objects        ( void * data,
                                                    int * ierr )
{
  *ierr = ZOLTAN_OK;

  return map_->NumMyElements();
}

void Isorropia::ZoltanQuery::Object_List  ( void * data,
                                        int num_gid_entries,
                                        int num_lid_entries,
                                        ZOLTAN_ID_PTR global_ids,
                                        ZOLTAN_ID_PTR local_ids,
                                        int weight_dim,
                                        float * object_weights,
                                        int * ierr )
{
  *ierr = ZOLTAN_OK;

  int rows = map_->NumMyElements();

  map_->MyGlobalElements( ((int *) global_ids) );

  int Index = map_->IndexBase();
  for( int i = 0; i < rows; i++, Index++ ) {
    local_ids[i] = Index;
  }
}

//Graph Based Functions
int Isorropia::ZoltanQuery::Number_Edges  ( void * data,
                                        int num_gid_entities,
                                        int num_lid_entities,
                                        ZOLTAN_ID_PTR global_id,
                                        ZOLTAN_ID_PTR local_id,
                                        int * ierr )
{
  int row = *global_id;
  int LocalRow = map_->LID(row);

  if( LocalRow >= 0 && LocalRow == (int)*local_id )
  {
    *ierr = ZOLTAN_OK;

    int NumIndices = ugraph_[LocalRow].size();

    return NumIndices;
  }
  else
  {
    *ierr = ZOLTAN_FATAL;
    return -1;
  }
}

void Isorropia::ZoltanQuery::Edge_List    ( void * data,
                                        int num_gid_entities,
                                        int num_lid_entities,
                                        ZOLTAN_ID_PTR global_id,
                                        ZOLTAN_ID_PTR local_id,
                                        ZOLTAN_ID_PTR neighbor_global_ids,
                                        int * neighbor_procs,
                                        int weight_dim,
                                        float * edge_weights,
                                        int * ierr )
{
  int row = *global_id;
  int LocalRow = map_->LID(row);

  if (LocalRow < 0 || LocalRow != (int)*local_id) {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  int localProc = map_->Comm().MyPID();

  std::set<int>& rowset = ugraph_[LocalRow];

  std::set<int>::const_iterator
    iter = rowset.begin(),
    iter_end = rowset.end();

  int offset = 0;
  for(; iter != iter_end; ++iter) {
    int index = *iter;
    neighbor_global_ids[offset] = index;

    if (map_->MyGID(index)) {
      neighbor_procs[offset] = localProc;
    }
    else {
      std::map<int,int>::const_iterator
	iter = procmap_.find(index),
	iter_end = procmap_.end();
      if (iter == iter_end) {
	*ierr = ZOLTAN_FATAL;
	return;
      }

      neighbor_procs[offset] = (*iter).second;
    }

    ++offset;
  }

//   std::ostringstream osstr;
//   osstr << "ZQ2.out."<<localProc;
//   static std::ofstream* ofs =
//     new std::ofstream(osstr.str().c_str(), std::ios::out);
//   *ofs << localProc << " row " << row << ": ";
//   for(int i=0; i<offset; ++i) {
//     *ofs << "("<<neighbor_global_ids[i]<<","<<neighbor_procs[i]<<") ";
//   }
//   *ofs << std::endl;

  *ierr = ZOLTAN_OK;
}

void Isorropia::ZoltanQuery::HG_Size_CS ( void * data,
                                          int* num_lists,
                                          int* num_pins,
                                          int* format,
					  int * ierr )
{
  *num_lists = tmap_->NumMyElements();
  if (tgraph_.get() != 0) {
    *num_pins = tgraph_->NumMyNonzeros();
  }
  else {
    *num_pins = tmatrix_->NumMyNonzeros();
  }

  *format = ZOLTAN_COMPRESSED_EDGE;

  *ierr = ZOLTAN_OK;
}

void Isorropia::ZoltanQuery::HG_CS ( void * data,
                                     int num_gid_entries,
                                     int num_row_or_col,
                                     int num_pins,
                                     int format,
                                     ZOLTAN_ID_PTR vtxedge_GID,
                                     int* vtxedge_ptr,
                                     ZOLTAN_ID_PTR pin_GID,
                                     int * ierr )
{
  if (num_gid_entries != 1) {
    std::cout << "ZoltanQuery::HG_CS: num_gid_entries="<<num_gid_entries
        << ", not yet allowed by this implementation!"<< std::endl;
     *ierr = ZOLTAN_FATAL;
     return;
  }

  int numMyRows = 0;
  int numMyNonzeros = 0;
  if (tgraph_.get() != 0) {
    numMyRows = tgraph_->NumMyRows();
    numMyNonzeros = tgraph_->NumMyNonzeros();
  }
  else {
    numMyRows = tmatrix_->NumMyRows();
    numMyNonzeros = tmatrix_->NumMyNonzeros();
  }

  if (num_row_or_col != numMyRows) {
    std::cout << "ZoltanQuery::HG_CS: num_row_or_col (= "<<num_row_or_col
        << ") != NumMyRows (= " << numMyRows << std::endl;
    *ierr = ZOLTAN_FATAL;
    return;
  }

  if (num_pins != numMyNonzeros) {
    std::cout << "ZoltanQuery::HG_CS: num_pins (= "<<num_pins
        << ") != NumMyNonzeros (= " << numMyNonzeros << std::endl;
    *ierr = ZOLTAN_FATAL;
    return;
  }

  int* rows = (int*)vtxedge_GID;
  tmap_->MyGlobalElements(rows);

  int offset = 0;
  if (tgraph_.get() != 0) {
    for(int i=0; i<num_row_or_col; ++i) {
      vtxedge_ptr[i] = offset;

      int rowlen = tgraph_->NumMyIndices(i);

      int checkNumIndices;
      int* indices = (int*)(&(pin_GID[offset]));
      tgraph_->ExtractMyRowCopy(i, rowlen, checkNumIndices, indices);

      offset += rowlen;
    }
  }
  else {
    std::vector<double> coefs(tmatrix_->MaxNumEntries());
    for(int i=0; i<num_row_or_col; ++i) {
      vtxedge_ptr[i] = offset;

      int rowlen;
      tmatrix_->NumMyRowEntries(i, rowlen);

      int checkNumIndices;
      int* indices = (int*)(&(pin_GID[offset]));
      tmatrix_->ExtractMyRowCopy(i, rowlen, checkNumIndices, &coefs[0], indices);

      offset += rowlen;
    }
  }

  *ierr = ZOLTAN_OK;
}
