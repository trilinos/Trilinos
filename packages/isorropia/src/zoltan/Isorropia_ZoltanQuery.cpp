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
#include <Epetra_BlockMap.h>
#include <Epetra_Comm.h>

#include <algorithm>

Isorropia::ZoltanQuery::ZoltanQuery( const Epetra_CrsGraph & graph,
                                     const Epetra_CrsGraph * tgraph,
                                     bool localEdgesOnly )
: graph_(graph),
  tgraph_(tgraph),
  map_(graph.RowMap()),
  localEdgesOnly_(localEdgesOnly)
{
  if (tgraph != 0) {
    tmap_ = &(tgraph->RowMap());
  }

  int numMyRows = map_.NumMyElements();
  int maxRows;
  map_.Comm().MaxAll( &numMyRows, &maxRows, 1 );

  LBProc_.resize( numMyRows );

  int numIndices;
  int maxNumIndices = graph_.MaxNumIndices();
  vector<int> indexList( maxNumIndices );
  for( int i = 0; i < numMyRows; ++i )
  {
    graph_.ExtractGlobalRowCopy( graph_.GRID(i),
                                 maxNumIndices,
                                 numIndices,
                                 &indexList[0] );
    LBProc_[i].resize( numIndices );
    map_.RemoteIDList( numIndices, &indexList[0], &LBProc_[i][0], 0 );
  }

  for( int i = numMyRows; i < maxRows; ++i )
    map_.RemoteIDList( numIndices, &indexList[0], &LBProc_[numMyRows-1][0], 0 );

  if( tgraph_ )
  {
    LBProc_Trans_.resize( numMyRows );

    maxNumIndices = tgraph_->MaxNumIndices();
    indexList.resize(maxNumIndices);
    for( int i = 0; i < numMyRows; ++i )
    {
      tgraph_->ExtractGlobalRowCopy( tgraph_->GRID(i),
                                     maxNumIndices,
                                     numIndices,
                                     &indexList[0] );
      LBProc_Trans_[i].resize( numIndices );
      tmap_->RemoteIDList( numIndices, &indexList[0], &LBProc_Trans_[i][0], 0 );
    }

    for( int i = numMyRows; i < maxRows; ++i )
      tmap_->RemoteIDList( numIndices, &indexList[0], &LBProc_Trans_[numMyRows-1][0], 0 );
  }

  map_.Comm().Barrier();
}

//General Functions
int Isorropia::ZoltanQuery::Number_Objects        ( void * data,
                                                    int * ierr )
{
  *ierr = ZOLTAN_OK;

  return map_.NumMyElements();
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

  int rows = map_.NumMyElements();

  map_.MyGlobalElements( ((int *) global_ids) );

  int Index = map_.IndexBase();
  for( int i = 0; i < rows; i++, Index++ )
    local_ids[i] = Index;
}

//Graph Based Functions
int Isorropia::ZoltanQuery::Number_Edges  ( void * data,
                                        int num_gid_entities,
                                        int num_lid_entities,
                                        ZOLTAN_ID_PTR global_id,
                                        ZOLTAN_ID_PTR local_id,
                                        int * ierr )
{
  int LocalRow = map_.LID( *global_id );

  if( LocalRow != -1 && LocalRow == (int)*local_id )
  {
    *ierr = ZOLTAN_OK;

    int NumIndices = graph_.NumMyIndices( LocalRow );
    int IndiceCountReturn;

    vector<int> nbr_edges( NumIndices );
    assert( graph_.ExtractGlobalRowCopy( ((int) *global_id),
                                         NumIndices,
                                         IndiceCountReturn,
                                         &(nbr_edges[0]) ) == 0 );
    assert( NumIndices == IndiceCountReturn );
    sort( nbr_edges.begin(), nbr_edges.end() );

    bool self = false;
    for(int i = 0; i < NumIndices; ++i )
      if( nbr_edges[i] == ((int) *global_id) )
      { self = true; break; }

    int nonLocalEdges = 0;
    if( localEdgesOnly_ )
      for( int i = 0; i < NumIndices; ++i )
        if( !graph_.MyGRID(nbr_edges[i]) ) ++nonLocalEdges;

    if( tgraph_ )
    {
      int tNumIndices = tgraph_->NumMyIndices( LocalRow );
      vector<int> t_nbr_edges( tNumIndices );

      assert( tgraph_->ExtractGlobalRowCopy( ((int) *global_id),
                                           tNumIndices,
                                           IndiceCountReturn,
                                           &(t_nbr_edges[0]) ) == 0 );
      assert( tNumIndices == IndiceCountReturn );

      for( int i = 0; i < tNumIndices; ++i )
        if( !binary_search( nbr_edges.begin(), nbr_edges.end(), t_nbr_edges[i] ) )
        {
          ++NumIndices;
          if( localEdgesOnly_ && !graph_.MyGRID(t_nbr_edges[i]) ) ++nonLocalEdges;
        }
    }

//cout << "Indices Cnt: " << ((int)*global_id) << " " << ((int)*local_id) << " " << NumIndices-(self?1:0)-nonLocalEdges  << endl;
    return NumIndices - (self?1:0) - nonLocalEdges;

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
  int NumIndices = graph_.NumMyIndices( ((int) *local_id) );

  int IndiceCountReturn;

  if( NumIndices != -1 )
  {
    vector<int> nbr_edges( NumIndices );
    assert( graph_.ExtractGlobalRowCopy( ((int) *global_id), 
                                         NumIndices,
                                         IndiceCountReturn,
                                         &(nbr_edges[0]) ) == 0 );
    assert( NumIndices == IndiceCountReturn );

    int ii = 0;
    for( int i = 0; i < NumIndices; ++i )
      if( nbr_edges[ i ] != ((int) *global_id) )
        if( !localEdgesOnly_ || graph_.MyGRID(nbr_edges[i]) )
        {
          neighbor_global_ids[ ii ] = nbr_edges[ i ];
          neighbor_procs[ ii ] = LBProc_[(int) *local_id][ i ];
          ++ii;
        }

    if( tgraph_ )
    {
      sort( nbr_edges.begin(), nbr_edges.end() );

      int tNumIndices = tgraph_->NumMyIndices( ((int) *local_id) );
      vector<int> t_nbr_edges( tNumIndices );

      assert( tgraph_->ExtractGlobalRowCopy( ((int) *global_id),
                                           tNumIndices,
                                           IndiceCountReturn,
                                           &(t_nbr_edges[0]) ) == 0 );
      assert( tNumIndices == IndiceCountReturn );

      for( int i = 0; i < tNumIndices; ++i )
        if( !binary_search( nbr_edges.begin(), nbr_edges.end(), t_nbr_edges[i] ) )
          if( !localEdgesOnly_ || graph_.MyGRID(t_nbr_edges[i]) )
          {
            neighbor_global_ids[ii] = t_nbr_edges[i];
            neighbor_procs[ii] = LBProc_Trans_[(int) *local_id][i];
            ++ii;
          }
    }


/*
cout << "Edge List: " << ((int) *global_id) << " " << ((int) *local_id) << endl;
cout << "NumIndices: " << NumIndices << " " << "ii: " << ii << endl;
for( int i = 0; i < ii; ++i )
  cout << " " << ((int *) neighbor_global_ids)[i] << " " <<
        neighbor_procs[i] << endl;
cout << endl;
*/

    *ierr = ZOLTAN_OK;
  }
  else
    *ierr = ZOLTAN_FATAL;

}

void Isorropia::ZoltanQuery::HG_Size_CS ( void * data,
                                          int* num_lists,
                                          int* num_pins,
                                          int* format,
					  int * ierr )
{
  *num_lists = tmap_->NumMyElements();
  *num_pins = tgraph_->NumMyNonzeros();

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

  if (num_row_or_col != tgraph_->NumMyRows()) {
    std::cout << "ZoltanQuery::HG_CS: num_row_or_col (= "<<num_row_or_col
        << ") != NumMyRows (= " << graph_.NumMyRows() << std::endl;
    *ierr = ZOLTAN_FATAL;
    return;
  }

  if (num_pins != tgraph_->NumMyNonzeros()) {
    std::cout << "ZoltanQuery::HG_CS: num_pins (= "<<num_pins
        << ") != NumMyNonzeros (= " << graph_.NumMyNonzeros() << std::endl;
    *ierr = ZOLTAN_FATAL;
    return;
  }

  int* rows = (int*)vtxedge_GID;
  tmap_->MyGlobalElements(rows);

  int offset = 0;
  for(int i=0; i<num_row_or_col; ++i) {
    vtxedge_ptr[i] = offset;

    int rowlen = tgraph_->NumMyIndices(i);

    int checkNumIndices;
    int* indices = (int*)(&(pin_GID[offset]));
    tgraph_->ExtractMyRowCopy(i, rowlen, checkNumIndices, indices);

    offset += rowlen;
  }

  *ierr = ZOLTAN_OK;
}

