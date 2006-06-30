//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_QueryObject.C$
//
// Purpose        : Base Class for dynamic versions of query
//                  functions to be used by Zoltan.  The application
//                  requires a class derived from this base to
//                  be instantiated and registered with the
//                  Zoltan_LoadBalance object.
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra
//
// Creation Date  : 08/04/2000
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision$
//
// Revision Date  : $Date$
//
// Current Owner  : $Author$
//-------------------------------------------------------------------------

#include <IZoltan_QueryObject.h>

  //General Functions
int Zoltan::QueryObject::Number_Objects        (	void * data,
						int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Number_Objects( void *, "
	<< "int * ) must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

void Zoltan::QueryObject::Object_List  (	void * data,
					int num_gid_entries,
					int num_lid_entries,
					ZOLTAN_ID_PTR global_ids, 
					ZOLTAN_ID_PTR local_ids,
					int weight_dim,
					float * object_weights,
					int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Object_List( void *, int, int, "
	<< "ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int, float *, int * ) must be implemented." 
	<< std::endl;

  *ierr = ZOLTAN_FATAL;
}

int Zoltan::QueryObject::First_Object  (	void * data,
					int num_gid_entries,
					int num_lid_entries,
					ZOLTAN_ID_PTR first_global_id, 
					ZOLTAN_ID_PTR first_local_id,
					int weight_dim,
					float * first_weight,
					int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::First_Object( void *, int, int, "
	<< "ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int, float *, int * ) must be implemented." 
	<< std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

int Zoltan::QueryObject::Next_Object   ( void * data,
					int num_gid_entries,
					int num_lid_entries,
					ZOLTAN_ID_PTR global_id, 
					ZOLTAN_ID_PTR local_id,
					ZOLTAN_ID_PTR next_global_id, 
					ZOLTAN_ID_PTR next_local_id,
					int weight_dim,
					float * next_weight,
					int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Next_Object( void *, int, int, "
	<< "ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int, float *, int * ) "
        << "must be implemented." 
	<< std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

int Zoltan::QueryObject::Number_Border_Objects (	void * data,
						int number_neighbor_procs,
						int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Number_Border_Objects( void *, "
	<< "int, int * ) must be implemented." 
	<< std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

void Zoltan::QueryObject::Border_Object_List   (	void * data,
						int num_gid_entities,
						int num_lid_entities,
						int number_neighbor_procs,
						ZOLTAN_ID_PTR global_ids,
						ZOLTAN_ID_PTR local_ids,
						int weight_dim,
						float * object_weights,
						int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Border_Object_List( void *, int, "
        << "int, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int, float *, int * ) must be "
	<< "implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;
}

int Zoltan::QueryObject::First_Border_Object   (	void * data,
						int num_gid_entities,
						int num_lid_entities,
						int number_neighbor_procs,
						ZOLTAN_ID_PTR first_global_id,
						ZOLTAN_ID_PTR first_local_id,
						int weight_dim,
						float * first_weight,
						int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::First_Border_Object( void *, "
        << "int, int, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int, float *, int * ) must be "
	<< "implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

int Zoltan::QueryObject::Next_Border_Object    (	void * data,
						int num_gid_entities,
						int num_lid_entities,
						ZOLTAN_ID_PTR global_id,
						ZOLTAN_ID_PTR local_id,
						int number_neighbor_procs,
						ZOLTAN_ID_PTR next_global_id,
						ZOLTAN_ID_PTR next_local_id,
						int weight_dim,
						float * next_weight,
						int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Next_Border_Object( void *, "
	<< "int, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int, ZOLTAN_GID *, ZOLTAN_LID *, int, "
	<< "float *, int * ) must be "
	<< "implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

  //Geometry Based Functions
int Zoltan::QueryObject::Number_Geometry_Objects       (	void * data,
							int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Number_Geometry_Objects( void *, "
	<< "int * ) must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

void Zoltan::QueryObject::Geometry_Values      (	void * data,
						int num_gid_entities,
						int num_lid_entities,
						ZOLTAN_ID_PTR global_id, 
						ZOLTAN_ID_PTR local_id,
						double * geometry_vector,
						int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Geometry_Values( void *, int, int, "
	<< "ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, double *, int * ) must be implemented." 
	<< std::endl;

  *ierr = ZOLTAN_FATAL;
}

  //Graph Based Functions
int Zoltan::QueryObject::Number_Edges  (	void * data,
					int num_gid_entities,
					int num_lid_entities,
					ZOLTAN_ID_PTR global_id, 
					ZOLTAN_ID_PTR local_id,
					int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Number_Edges( void *, int, int, "
	<< "ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int * ) must be implemented." 
	<< std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

void Zoltan::QueryObject::Edge_List    (	void * data,
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
  std::cout << "Error: Zoltan::QueryObject::Edge_List( void *, int, int, "
	<< "ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int, float *, int * ) "
	<< "must be implemented." 
	<< std::endl;

  *ierr = ZOLTAN_FATAL;
}

  //Hypergraph Based Functions
void Zoltan::QueryObject::HG_Size_CS(void* data,
                          int* num_lists,
                          int* num_pins,
                          int* format,
                          int* ierr)
{
  std::cout << "Error: Zoltan::QueryObject::HG_Size_CS( void *, int*, int*, "
	<< " int*, int* ) "
	<< "must be implemented." 
	<< std::endl;

  *ierr = ZOLTAN_FATAL;
}

void Zoltan::QueryObject::HG_CS(void* data,
                     int num_gid_entries,
                     int num_row_or_col,
                     int num_pins,
                     int format,
                     ZOLTAN_ID_PTR vtxedge_GID,
                     int* vtxedge_ptr,
                     ZOLTAN_ID_PTR pin_GID,
                     int* ierr)
{
  std::cout << "Error: Zoltan::QueryObject::HG_CS( void *, int, int, "
        << " int, int, ZOLTAN_ID_PTR, int*, ZOLTAN_ID_PTR, int* ) "
        << "must be implemented."
        << std::endl;

  *ierr = ZOLTAN_FATAL;
}

void Zoltan::QueryObject::HG_Size_Edge_Weights(void * data,
				    int* num_edges,
				    int* ierr)
{
  std::cout << "Error: Zoltan::QueryObject::HG_Size_Edge_Weights( void *, int*, "
        << " int* ) must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;
}

void Zoltan::QueryObject::HG_Edge_Weights(void * data,
			      int num_gid_entries,
			      int num_lid_entries,
			      int num_edges,
			      int edge_weight_dim,
			      ZOLTAN_ID_PTR edge_GID,
			       ZOLTAN_ID_PTR edge_LID,
			      float* edge_weights,
			      int* ierr)
{
  std::cout << "Error: Zoltan::QueryObject::HG_Edge_Weights( void *, int*, "
        << " int* ) must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;
}

  //Tree Based Functions
int Zoltan::QueryObject::Number_Coarse_Objects (	void * data,
						int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Number_Coarse_Objects( void *, "
	<< "int * ) must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

void Zoltan::QueryObject::Coarse_Object_List   (	void * data,
						int num_gid_entities,
						int num_lid_entities,
						ZOLTAN_ID_PTR global_ids,
						ZOLTAN_ID_PTR local_ids,
						int * assigned,
						int * number_vertices,
						ZOLTAN_ID_PTR vertices,
						int * in_order,
						ZOLTAN_ID_PTR in_vertex,
						ZOLTAN_ID_PTR out_vertex,
						int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Coarse_Object_List( void *, "
	<< "int, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, int *, int *, int *, "
	<< "int *, int * ) "
	<< "must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;
}

int Zoltan::QueryObject::First_Coarse_Object   ( void * data,
						int num_gid_entities,
						int num_lid_entities,
						ZOLTAN_ID_PTR first_global_id,
						ZOLTAN_ID_PTR first_local_id,
						int * assigned,
						int * number_vertices,
						ZOLTAN_ID_PTR vertices,
						int * in_order,
						ZOLTAN_ID_PTR in_vertex,
						ZOLTAN_ID_PTR out_vertex,
						int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::First_Coarse_Object( void *, "
	<< "int, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, int *, int *, int *, "
	<< "int *, int * ) "
	<< "must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

int Zoltan::QueryObject::Next_Coarse_Object    ( void * data,
						int num_gid_entities,
						int num_lid_entities,
						ZOLTAN_ID_PTR global_id,
						ZOLTAN_ID_PTR local_id,
						ZOLTAN_ID_PTR next_global_id,
						ZOLTAN_ID_PTR next_local_id,
						int * assigned,
						int * number_vertices, 
						ZOLTAN_ID_PTR vertices,
						ZOLTAN_ID_PTR in_vertex,
						ZOLTAN_ID_PTR out_vertex,
						int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Next_Coarse_Object( void *, "
	<< "int, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, "
	<< "int *, int *, int *, int * ) "
	<< "must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

int Zoltan::QueryObject::Number_Children       (	void * data,
						int num_gid_entities,
						int num_lid_entities,
						ZOLTAN_ID_PTR global_id,
						ZOLTAN_ID_PTR local_id,
						int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Number_Children( void *, int, int, "
	<< "ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int * ) "
	<< "must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

void Zoltan::QueryObject::Child_List   (	void * data,
					int num_gid_entities,
					int num_lid_entities,
					ZOLTAN_ID_PTR parent_global_id,
					ZOLTAN_ID_PTR parent_local_id,
					ZOLTAN_ID_PTR child_global_ids,
					ZOLTAN_ID_PTR child_local_ids,
					int * assigned,
					int * number_vertices,
					ZOLTAN_ID_PTR vertices,
					ZOLTAN_REF_TYPE * reference_type,
					ZOLTAN_ID_PTR in_vertex,
					ZOLTAN_ID_PTR out_vertex,
					int * ierr  )
{
  std::cout << "Error: Zoltan::QueryObject::Child_List( void *, int, int, "
	<< "ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, int *, int *, "
	<< "ZOLTAN_REF_TYPE *, int *, int *, int * ) must be implemented."
	<< std::endl;

  *ierr = ZOLTAN_FATAL;
}

void Zoltan::QueryObject::Child_Weight (	void * data,
					int num_gid_entities,
					int num_lid_entities,
					ZOLTAN_ID_PTR global_id,
					ZOLTAN_ID_PTR local_id,
					int weight_dim,
					float * object_weight,
					int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Child_Weight( void *, int, int, "
	<< "ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int, float *, int * ) must be implemented."
	<< std::endl;

  *ierr = ZOLTAN_FATAL;
}
