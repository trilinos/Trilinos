
#include <Epetra_ZoltanQuery.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_BlockMap.h>
#include <Epetra_Comm.h>

Epetra_ZoltanQuery::Epetra_ZoltanQuery( const Epetra_CrsGraph & graph )
: graph_(graph)
{
  int numMyRows = graph_.NumMyRows();
  int maxRows;
  graph_.Comm().MaxAll( &numMyRows, &maxRows, 1 );

  LBProc_.resize( numMyRows );

  int numIndices;
  int maxNumIndices = graph_.MaxNumIndices();
  vector<int> indexList( maxNumIndices );
  for( int i = 0; i < numMyRows; ++i )
  {
    graph_.ExtractGlobalRowCopy( graph_.GCID(i),
                                 maxNumIndices,
                                 numIndices,
                                 &indexList[0] );
    LBProc_[i].resize( numIndices );
    graph_.RowMap().RemoteIDList( numIndices, &indexList[0], &LBProc_[i][0], 0 );
  }

  for( int i = numMyRows; i < maxRows; ++i )
    graph_.RowMap().RemoteIDList( numIndices, &indexList[0], &LBProc_[numMyRows-1][0], 0 );

  graph_.Comm().Barrier();
}

//General Functions
int Epetra_ZoltanQuery::Number_Objects        ( void * data,
                                                int * ierr )
{
  *ierr = LB_OK;

  return graph_.NumMyRows();
}

void Epetra_ZoltanQuery::Object_List  ( void * data,
                                        int num_gid_entries,
                                        int num_lid_entries,
                                        LB_ID_PTR global_ids,
                                        LB_ID_PTR local_ids,
                                        int weight_dim,
                                        float * object_weights,
                                        int * ierr )
{
  *ierr = LB_OK;

  int rows = graph_.NumMyRows();

  graph_.RowMap().MyGlobalElements( ((int *) global_ids) );

  int Index = graph_.RowMap().IndexBase();
  for( int i = 0; i < rows; i++, Index++ )
    local_ids[i] = Index;
}

int Epetra_ZoltanQuery::First_Object  ( void * data,
                                        int num_gid_entries,
                                        int num_lid_entries,
                                        LB_ID_PTR first_global_id,
                                        LB_ID_PTR first_local_id,
                                        int weight_dim,
                                        float * first_weight,
                                        int * ierr )
{ 
  cout << "Error: Epetra_ZoltanQuery::First_Object( void *, int, int, "
        << "LB_ID_PTR, LB_ID_PTR, int, float *, int * ) must be implemented."
        << endl;
  
  *ierr = LB_FATAL;
  
  return 0;
}

int Epetra_ZoltanQuery::Next_Object   ( void * data,
                                        int num_gid_entries,
                                        int num_lid_entries,
                                        LB_ID_PTR global_id,
                                        LB_ID_PTR local_id,
                                        LB_ID_PTR next_global_id,
                                        LB_ID_PTR next_local_id,
                                        int weight_dim,
                                        float * next_weight,
                                        int * ierr )
{ 
  cout << "Error: Epetra_ZoltanQuery::Next_Object( void *, int, int, " 
        << "LB_ID_PTR, LB_ID_PTR, LB_ID_PTR, LB_ID_PTR, int, float *, int * ) "
        << "must be implemented."
        << endl;
  
  *ierr = LB_FATAL;
  
  return 0;
}

int Epetra_ZoltanQuery::Number_Border_Objects ( void * data,
                                                int number_neighbor_procs,
                                                int * ierr )
{ 
  cout << "Error: Epetra_ZoltanQuery::Number_Border_Objects( void *, "
        << "int, int * ) must be implemented."
        << endl;
  
  *ierr = LB_FATAL;
  
  return 0;
}

void Epetra_ZoltanQuery::Border_Object_List   ( void * data,
                                                int num_gid_entities,
                                                int num_lid_entities,
                                                int number_neighbor_procs,
                                                LB_ID_PTR global_ids,
                                                LB_ID_PTR local_ids,
                                                int weight_dim,
                                                float * object_weights,
                                                int * ierr )
{ 
  cout << "Error: Epetra_ZoltanQuery::Border_Object_List( void *, int, "
        << "int, int, LB_ID_PTR, LB_ID_PTR, int, float *, int * ) must be "
        << "implemented." << endl;
  
  *ierr = LB_FATAL;
}

int Epetra_ZoltanQuery::First_Border_Object   ( void * data,
                                                int num_gid_entities,
                                                int num_lid_entities,
                                                int number_neighbor_procs,
                                                LB_ID_PTR first_global_id,
                                                LB_ID_PTR first_local_id,
                                                int weight_dim,
                                                float * first_weight,
                                                int * ierr )
{ 
  cout << "Error: Epetra_ZoltanQuery::First_Border_Object( void *, "
        << "int, int, int, LB_ID_PTR, LB_ID_PTR, int, float *, int * ) must be "
        << "implemented." << endl;
  
  *ierr = LB_FATAL;
  
  return 0;
}

int Epetra_ZoltanQuery::Next_Border_Object    ( void * data,
                                                int num_gid_entities,
                                                int num_lid_entities,
                                                LB_ID_PTR global_id,
                                                LB_ID_PTR local_id,
                                                int number_neighbor_procs,
                                                LB_ID_PTR next_global_id,
                                                LB_ID_PTR next_local_id,
                                                int weight_dim,
                                                float * next_weight,
                                                int * ierr )
{
  cout << "Error: Epetra_ZoltanQuery::Next_Border_Object( void *, "
        << "int, int, LB_ID_PTR, LB_ID_PTR, int, LB_GID *, LB_LID *, int, "
        << "float *, int * ) must be "
        << "implemented." << endl;

  *ierr = LB_FATAL;

  return 0;
}

//Geometry Based Functions
int Epetra_ZoltanQuery::Number_Geometry_Objects       ( void * data,
                                                        int * ierr )
{
  cout << "Error: Epetra_ZoltanQuery::Number_Geometry_Objects( void *, "
        << "int * ) must be implemented." << endl;

  *ierr = LB_FATAL;

  return 0;
}

void Epetra_ZoltanQuery::Geometry_Values      ( void * data,
                                                int num_gid_entities,
                                                int num_lid_entities,
                                                LB_ID_PTR global_id,
                                                LB_ID_PTR local_id,
                                                double * geometry_vector,
                                                int * ierr )
{
  cout << "Error: Epetra_ZoltanQuery::Geometry_Values( void *, int, int, "
        << "LB_ID_PTR, LB_ID_PTR, double *, int * ) must be implemented."
        << endl;

  *ierr = LB_FATAL;
}

//Graph Based Functions
int Epetra_ZoltanQuery::Number_Edges  ( void * data,
                                        int num_gid_entities,
                                        int num_lid_entities,
                                        LB_ID_PTR global_id,
                                        LB_ID_PTR local_id,
                                        int * ierr )
{
  int LocalRow = graph_.LRID( *global_id );

  if( LocalRow != -1 && LocalRow == *local_id )
  {
    *ierr = LB_OK;

    int NumIndices = graph_.NumMyIndices( LocalRow );
    int IndiceCountReturn;

    vector<int> nbr_edges( NumIndices );

    assert( graph_.ExtractGlobalRowCopy( ((int) *global_id),
                                         NumIndices,
                                         IndiceCountReturn,
                                         &(nbr_edges[0]) ) == 0 );

    assert( NumIndices == IndiceCountReturn );

    for(int i = 0; i < NumIndices; i++ )
      if( nbr_edges[i] == ((int) *global_id) )
      {
        NumIndices--;
        break;
      }

    return NumIndices;

  }
  else
  {
    *ierr = LB_FATAL;
    return -1;
  }
}

void Epetra_ZoltanQuery::Edge_List    ( void * data,
                                        int num_gid_entities,
                                        int num_lid_entities,
                                        LB_ID_PTR global_id,
                                        LB_ID_PTR local_id,
                                        LB_ID_PTR neighbor_global_ids,
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
    for( int i = 0; i < NumIndices; i++ )
    {
      if( nbr_edges[ i ] != ((int) *global_id) )
      {
        neighbor_global_ids[ ii ] = nbr_edges[ i ];
        neighbor_procs[ ii ] = LBProc_[(int) *local_id][ i ];

        ii++;
      }
    }

/*
cout << "Edge List: " << ((int) *global_id) << " " << ((int) *local_id) << endl;
cout << "NumIndices: " << NumIndices << " " << "ii: " << ii << endl;
for( int i = 0; i < (NumIndices); i++ )
  cout << " " << ((int *) neighbor_global_ids)[i] << " " <<
        neighbor_procs[i] << endl;
cout << endl;
*/

    *ierr = LB_OK;
  }
  else
    *ierr = LB_FATAL;

}

//Tree Based Functions
int Epetra_ZoltanQuery::Number_Coarse_Objects ( void * data,
                                                int * ierr )
{
  cout << "Error: Epetra_ZoltanQuery::Number_Coarse_Objects( void *, "
        << "int * ) must be implemented." << endl;

  *ierr = LB_FATAL;

  return 0;
}

void Epetra_ZoltanQuery::Coarse_Object_List   ( void * data,
                                                int num_gid_entities,
                                                int num_lid_entities,
                                                LB_ID_PTR global_ids,
                                                LB_ID_PTR local_ids,
                                                int * assigned,
                                                int * number_vertices,
                                                int * vertices,
                                                int * in_order,
                                                int * in_vertex,
                                                int * out_vertex,
                                                int * ierr )
{
  cout << "Error: Epetra_ZoltanQuery::Coarse_Object_List( void *, "
        << "int, int, LB_ID_PTR, LB_ID_PTR, int *, int *, int *, int *, "
        << "int *, int * ) "
        << "must be implemented." << endl;

  *ierr = LB_FATAL;
}

int Epetra_ZoltanQuery::First_Coarse_Object   ( void * data,
                                                int num_gid_entities,
                                                int num_lid_entities,
                                                LB_ID_PTR first_global_id,
                                                LB_ID_PTR first_local_id,
                                                int * assigned,
                                                int * number_vertices,
                                                int * vertices,
                                                int * in_order,
                                                int * in_vertex,
                                                int * out_vertex,
                                                int * ierr )
{
  cout << "Error: Epetra_ZoltanQuery::First_Coarse_Object( void *, "
        << "int, int, LB_ID_PTR, LB_ID_PTR, int *, int *, int *, int *, "
        << "int *, int * ) "
        << "must be implemented." << endl;

  *ierr = LB_FATAL;

  return 0;
}

int Epetra_ZoltanQuery::Next_Coarse_Object    ( void * data,
                                                int num_gid_entities,
                                                int num_lid_entities,
                                                LB_ID_PTR global_id,
                                                LB_ID_PTR local_id,
                                                LB_ID_PTR next_global_id,
                                                LB_ID_PTR next_local_id,
                                                int * assigned,
                                                int * number_vertices,
                                                int * vertices,
                                                int * in_vertex,
                                                int * out_vertex,
                                                int * ierr )
{
  cout << "Error: Epetra_ZoltanQuery::Next_Coarse_Object( void *, "
        << "int, int, LB_ID_PTR, LB_ID_PTR, LB_ID_PTR, LB_ID_PTR, int *, "
        << "int *, int *, int *, int * ) "
        << "must be implemented." << endl;

  *ierr = LB_FATAL;

  return 0;
} 

int Epetra_ZoltanQuery::Number_Children       ( void * data,
                                                int num_gid_entities,
                                                int num_lid_entities,
                                                LB_ID_PTR global_id,
                                                LB_ID_PTR local_id,
                                                int * ierr )
{
  cout << "Error: Epetra_ZoltanQuery::Number_Children( void *, int, int, "
        << "LB_ID_PTR, LB_ID_PTR, int * ) "
        << "must be implemented." << endl;

  *ierr = LB_FATAL;

  return 0;
}

void Epetra_ZoltanQuery::Child_List   ( void * data,
                                        int num_gid_entities, 
                                        int num_lid_entities,
                                        LB_ID_PTR parent_global_id,
                                        LB_ID_PTR parent_local_id,
                                        LB_ID_PTR child_global_ids,
                                        LB_ID_PTR child_local_ids,
                                        int * assigned,
                                        int * number_vertices,
                                        int * vertices,
                                        LB_REF_TYPE * reference_type,
                                        int * in_vertex,
                                        int * out_vertex, 
                                        int * ierr  ) 
{
  cout << "Error: Epetra_ZoltanQuery::Child_List( void *, int, int, "
        << "LB_ID_PTR, LB_ID_PTR, LB_ID_PTR, LB_ID_PTR, int *, int *, int *, "
        << "LB_REF_TYPE *, int *, int *, int * ) must be implemented."
        << endl;

  *ierr = LB_FATAL;
} 

void Epetra_ZoltanQuery::Child_Weight ( void * data,
                                        int num_gid_entities,
                                        int num_lid_entities,
                                        LB_ID_PTR global_id,
                                        LB_ID_PTR local_id,
                                        int weight_dim,
                                        float * object_weight,
                                        int * ierr )
{ 
  cout << "Error: Epetra_ZoltanQuery::Child_Weight( void *, int, int, "
        << "LB_ID_PTR, LB_ID_PTR, int, float *, int * ) must be implemented."
        << endl;
  
  *ierr = LB_FATAL;
}



