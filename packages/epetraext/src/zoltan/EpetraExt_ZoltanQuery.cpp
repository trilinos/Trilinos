// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <EpetraExt_ZoltanQuery.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_BlockMap.h>
#include <Epetra_Comm.h>

#include <algorithm>

EPETRAEXT_DEPRECATED
EpetraExt::ZoltanQuery::ZoltanQuery( const Epetra_CrsGraph & graph,
                                     const Epetra_CrsGraph * tgraph,
                                     bool localEdgesOnly )
: graph_(graph),
  tgraph_(tgraph),
  localEdgesOnly_(localEdgesOnly)
{
  int numMyRows = graph_.NumMyRows();
  int maxRows;
  graph_.Comm().MaxAll( &numMyRows, &maxRows, 1 );

  LBProc_.resize( numMyRows );

  int numIndices;
  int maxNumIndices = graph_.MaxNumIndices();
  std::vector<int> indexList( maxNumIndices );
  for( int i = 0; i < numMyRows; ++i )
  {
    graph_.ExtractGlobalRowCopy( graph_.GRID(i),
                                 maxNumIndices,
                                 numIndices,
                                 &indexList[0] );
    LBProc_[i].resize( numIndices );
    graph_.RowMap().RemoteIDList( numIndices, &indexList[0], &LBProc_[i][0], 0 );
  }

  for( int i = numMyRows; i < maxRows; ++i )
    graph_.RowMap().RemoteIDList( numIndices, &indexList[0], &LBProc_[numMyRows-1][0], 0 );

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
      tgraph_->RowMap().RemoteIDList( numIndices, &indexList[0], &LBProc_Trans_[i][0], 0 );
    }

    for( int i = numMyRows; i < maxRows; ++i )
      tgraph_->RowMap().RemoteIDList( numIndices, &indexList[0], &LBProc_Trans_[numMyRows-1][0], 0 );
  }

  graph_.Comm().Barrier();
}

//General Functions
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanQuery::Number_Objects        ( void * data,
                                                    int * ierr )
{
  *ierr = ZOLTAN_OK;

  return graph_.NumMyRows();
}

EPETRAEXT_DEPRECATED
void EpetraExt::ZoltanQuery::Object_List  ( void * data,
                                        int num_gid_entries,
                                        int num_lid_entries,
                                        ZOLTAN_ID_PTR global_ids,
                                        ZOLTAN_ID_PTR local_ids,
                                        int weight_dim,
                                        float * object_weights,
                                        int * ierr )
{
  *ierr = ZOLTAN_OK;

  int rows = graph_.NumMyRows();

  graph_.RowMap().MyGlobalElements( ((int *) global_ids) );

  int Index = graph_.RowMap().IndexBase();
  for( int i = 0; i < rows; i++, Index++ )
    local_ids[i] = Index;
}

//Graph Based Functions
EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanQuery::Number_Edges  ( void * data,
                                        int num_gid_entities,
                                        int num_lid_entities,
                                        ZOLTAN_ID_PTR global_id,
                                        ZOLTAN_ID_PTR local_id,
                                        int * ierr )
{
  int LocalRow = graph_.LRID( *global_id );

  if( LocalRow != -1 && LocalRow == *local_id )
  {
    *ierr = ZOLTAN_OK;

    int NumIndices = graph_.NumMyIndices( LocalRow );
    int IndiceCountReturn;

    std::vector<int> nbr_edges( NumIndices );
    int flag = graph_.ExtractGlobalRowCopy( ((int) *global_id),
                                         NumIndices,
                                         IndiceCountReturn,
                                         &(nbr_edges[0]) );
    assert( flag == 0 );
    assert( NumIndices == IndiceCountReturn );
    std::sort( nbr_edges.begin(), nbr_edges.end() );

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
      std::vector<int> t_nbr_edges( tNumIndices );

      flag = tgraph_->ExtractGlobalRowCopy( ((int) *global_id),
                                           tNumIndices,
                                           IndiceCountReturn,
                                           &(t_nbr_edges[0]) );
      assert( flag == 0 );
      assert( tNumIndices == IndiceCountReturn );

      for( int i = 0; i < tNumIndices; ++i )
        if( !std::binary_search( nbr_edges.begin(), nbr_edges.end(), t_nbr_edges[i] ) )
        {
          ++NumIndices;
          if( localEdgesOnly_ && !graph_.MyGRID(t_nbr_edges[i]) ) ++nonLocalEdges;
        }
    }

//std::cout << "Indices Cnt: " << ((int)*global_id) << " " << ((int)*local_id) << " " << NumIndices-(self?1:0)-nonLocalEdges  << std::endl;
    return NumIndices - (self?1:0) - nonLocalEdges;

  }
  else
  {
    *ierr = ZOLTAN_FATAL;
    return -1;
  }
}

EPETRAEXT_DEPRECATED
void EpetraExt::ZoltanQuery::Edge_List    ( void * data,
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
    std::vector<int> nbr_edges( NumIndices );
    int flag = graph_.ExtractGlobalRowCopy( ((int) *global_id), 
                                         NumIndices,
                                         IndiceCountReturn,
                                         &(nbr_edges[0]) );
    assert( flag == 0 );
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
      std::sort( nbr_edges.begin(), nbr_edges.end() );

      int tNumIndices = tgraph_->NumMyIndices( ((int) *local_id) );
      std::vector<int> t_nbr_edges( tNumIndices );

      flag = tgraph_->ExtractGlobalRowCopy( ((int) *global_id),
                                           tNumIndices,
                                           IndiceCountReturn,
                                           &(t_nbr_edges[0]) );
      assert( flag == 0 );
      assert( tNumIndices == IndiceCountReturn );

      for( int i = 0; i < tNumIndices; ++i )
        if( !std::binary_search( nbr_edges.begin(), nbr_edges.end(), t_nbr_edges[i] ) )
          if( !localEdgesOnly_ || graph_.MyGRID(t_nbr_edges[i]) )
          {
            neighbor_global_ids[ii] = t_nbr_edges[i];
            neighbor_procs[ii] = LBProc_Trans_[(int) *local_id][i];
            ++ii;
          }
    }


/*
std::cout << "Edge List: " << ((int) *global_id) << " " << ((int) *local_id) << std::endl;
std::cout << "NumIndices: " << NumIndices << " " << "ii: " << ii << std::endl;
for( int i = 0; i < ii; ++i )
  std::cout << " " << ((int *) neighbor_global_ids)[i] << " " <<
        neighbor_procs[i] << std::endl;
std::cout << std::endl;
*/

    *ierr = ZOLTAN_OK;
  }
  else
    *ierr = ZOLTAN_FATAL;

}

EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanQuery::Number_HG_Edges ( void * data,
					  int * ierr )
{
  int num = graph_.NumMyRows();
  std::cout << "NRows: " << num << std::endl;

  *ierr = ZOLTAN_OK;
  return num;
}

EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanQuery::Number_HG_Pins ( void * data,
				         int * ierr )
{
  int num = graph_.NumMyEntries();
  std::cout << "NNZ: " << num << std::endl;

  *ierr = ZOLTAN_OK;
  return num;
}

EPETRAEXT_DEPRECATED
int EpetraExt::ZoltanQuery::HG_Edge_List   ( void * data,
                                         int num_gid_entries,
                                         int ewgt_dim,
                                         int nedge,
                                         int maxsize,
                                         int * edge_sizes,
                                         ZOLTAN_ID_PTR edge_verts,
                                         int * edge_procs,
                                         float * edge_weights )

{
  int NumHEs = graph_.NumMyRows();
  int maxEntries = graph_.MaxNumIndices();

  int numIndices;
  std::vector<int> indices( maxEntries );

  std::cout << "nedge: " << nedge << std::endl;
  std::cout << "maxsize: " << maxsize << std::endl;

  int loc = 0;
  for( int i = 0; i < NumHEs; ++i )
  {
    int flag = graph_.ExtractGlobalRowCopy( graph_.GRID(i), maxEntries, numIndices, &indices[0] );
    assert( flag == 0 );

    edge_sizes[i] = numIndices;
    for( int j = 0; j < numIndices; ++j )
    {
      edge_verts[loc] = indices[j];
      edge_procs[loc] = LBProc_[i][j];
      ++loc;
    }
  }

  std::cout << "last_loc: " << loc << std::endl;

  return ZOLTAN_OK;
}

