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

#include <iostream>

#include <Zoltan_QueryObject.h>

  //General Functions
int Zoltan_QueryObject::Number_Objects        (	void * data,
						int * ierr )
{
  std::cout << "Error: Zoltan_QueryObject::Number_Objects( void *, "
	<< "int * ) must be implemented." << std::endl;

  *ierr = LB_FATAL;

  return 0;
}

void Zoltan_QueryObject::Object_List  (	void * data,
					int num_gid_entries,
					int num_lid_entries,
					LB_ID_PTR global_ids, 
					LB_ID_PTR local_ids,
					int weight_dim,
					float * object_weights,
					int * ierr )
{
  std::cout << "Error: Zoltan_QueryObject::Object_List( void *, int, int, "
	<< "LB_ID_PTR, LB_ID_PTR, int, float *, int * ) must be implemented." 
	<< std::endl;

  *ierr = LB_FATAL;
}

int Zoltan_QueryObject::First_Object  (	void * data,
					int num_gid_entries,
					int num_lid_entries,
					LB_ID_PTR first_global_id, 
					LB_ID_PTR first_local_id,
					int weight_dim,
					float * first_weight,
					int * ierr )
{
  std::cout << "Error: Zoltan_QueryObject::First_Object( void *, int, int, "
	<< "LB_ID_PTR, LB_ID_PTR, int, float *, int * ) must be implemented." 
	<< std::endl;

  *ierr = LB_FATAL;

  return 0;
}

int Zoltan_QueryObject::Next_Object   ( void * data,
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
  std::cout << "Error: Zoltan_QueryObject::Next_Object( void *, int, int, "
	<< "LB_ID_PTR, LB_ID_PTR, LB_ID_PTR, LB_ID_PTR, int, float *, int * ) "
        << "must be implemented." 
	<< std::endl;

  *ierr = LB_FATAL;

  return 0;
}

int Zoltan_QueryObject::Number_Border_Objects (	void * data,
						int number_neighbor_procs,
						int * ierr )
{
  std::cout << "Error: Zoltan_QueryObject::Number_Border_Objects( void *, "
	<< "int, int * ) must be implemented." 
	<< std::endl;

  *ierr = LB_FATAL;

  return 0;
}

void Zoltan_QueryObject::Border_Object_List   (	void * data,
						int num_gid_entities,
						int num_lid_entities,
						int number_neighbor_procs,
						LB_ID_PTR global_ids,
						LB_ID_PTR local_ids,
						int weight_dim,
						float * object_weights,
						int * ierr )
{
  std::cout << "Error: Zoltan_QueryObject::Border_Object_List( void *, int, "
        << "int, int, LB_ID_PTR, LB_ID_PTR, int, float *, int * ) must be "
	<< "implemented." << std::endl;

  *ierr = LB_FATAL;
}

int Zoltan_QueryObject::First_Border_Object   (	void * data,
						int num_gid_entities,
						int num_lid_entities,
						int number_neighbor_procs,
						LB_ID_PTR first_global_id,
						LB_ID_PTR first_local_id,
						int weight_dim,
						float * first_weight,
						int * ierr )
{
  std::cout << "Error: Zoltan_QueryObject::First_Border_Object( void *, "
        << "int, int, int, LB_ID_PTR, LB_ID_PTR, int, float *, int * ) must be "
	<< "implemented." << std::endl;

  *ierr = LB_FATAL;

  return 0;
}

int Zoltan_QueryObject::Next_Border_Object    (	void * data,
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
  std::cout << "Error: Zoltan_QueryObject::Next_Border_Object( void *, "
	<< "int, int, LB_ID_PTR, LB_ID_PTR, int, LB_GID *, LB_LID *, int, "
	<< "float *, int * ) must be "
	<< "implemented." << std::endl;

  *ierr = LB_FATAL;

  return 0;
}

  //Geometry Based Functions
int Zoltan_QueryObject::Number_Geometry_Objects       (	void * data,
							int * ierr )
{
  std::cout << "Error: Zoltan_QueryObject::Number_Geometry_Objects( void *, "
	<< "int * ) must be implemented." << std::endl;

  *ierr = LB_FATAL;

  return 0;
}

void Zoltan_QueryObject::Geometry_Values      (	void * data,
						int num_gid_entities,
						int num_lid_entities,
						LB_ID_PTR global_id, 
						LB_ID_PTR local_id,
						double * geometry_vector,
						int * ierr )
{
  std::cout << "Error: Zoltan_QueryObject::Geometry_Values( void *, int, int, "
	<< "LB_ID_PTR, LB_ID_PTR, double *, int * ) must be implemented." 
	<< std::endl;

  *ierr = LB_FATAL;
}

  //Graph Based Functions
int Zoltan_QueryObject::Number_Edges  (	void * data,
					int num_gid_entities,
					int num_lid_entities,
					LB_ID_PTR global_id, 
					LB_ID_PTR local_id,
					int * ierr )
{
  std::cout << "Error: Zoltan_QueryObject::Number_Edges( void *, int, int, "
	<< "LB_ID_PTR, LB_ID_PTR, int * ) must be implemented." 
	<< std::endl;

  *ierr = LB_FATAL;

  return 0;
}

void Zoltan_QueryObject::Edge_List    (	void * data,
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
  std::cout << "Error: Zoltan_QueryObject::Edge_List( void *, int, int, "
	<< "LB_ID_PTR, LB_ID_PTR, LB_ID_PTR, int, float *, int * ) "
	<< "must be implemented." 
	<< std::endl;

  *ierr = LB_FATAL;
}

  //Tree Based Functions
int Zoltan_QueryObject::Number_Coarse_Objects (	void * data,
						int * ierr )
{
  std::cout << "Error: Zoltan_QueryObject::Number_Coarse_Objects( void *, "
	<< "int * ) must be implemented." << std::endl;

  *ierr = LB_FATAL;

  return 0;
}

void Zoltan_QueryObject::Coarse_Object_List   (	void * data,
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
  std::cout << "Error: Zoltan_QueryObject::Coarse_Object_List( void *, "
	<< "int, int, LB_ID_PTR, LB_ID_PTR, int *, int *, int *, int *, "
	<< "int *, int * ) "
	<< "must be implemented." << std::endl;

  *ierr = LB_FATAL;
}

int Zoltan_QueryObject::First_Coarse_Object   ( void * data,
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
  std::cout << "Error: Zoltan_QueryObject::First_Coarse_Object( void *, "
	<< "int, int, LB_ID_PTR, LB_ID_PTR, int *, int *, int *, int *, "
	<< "int *, int * ) "
	<< "must be implemented." << std::endl;

  *ierr = LB_FATAL;

  return 0;
}

int Zoltan_QueryObject::Next_Coarse_Object    ( void * data,
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
  std::cout << "Error: Zoltan_QueryObject::Next_Coarse_Object( void *, "
	<< "int, int, LB_ID_PTR, LB_ID_PTR, LB_ID_PTR, LB_ID_PTR, int *, "
	<< "int *, int *, int *, int * ) "
	<< "must be implemented." << std::endl;

  *ierr = LB_FATAL;

  return 0;
}

int Zoltan_QueryObject::Number_Children       (	void * data,
						int num_gid_entities,
						int num_lid_entities,
						LB_ID_PTR global_id,
						LB_ID_PTR local_id,
						int * ierr )
{
  std::cout << "Error: Zoltan_QueryObject::Number_Children( void *, int, int, "
	<< "LB_ID_PTR, LB_ID_PTR, int * ) "
	<< "must be implemented." << std::endl;

  *ierr = LB_FATAL;

  return 0;
}

void Zoltan_QueryObject::Child_List   (	void * data,
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
  std::cout << "Error: Zoltan_QueryObject::Child_List( void *, int, int, "
	<< "LB_ID_PTR, LB_ID_PTR, LB_ID_PTR, LB_ID_PTR, int *, int *, int *, "
	<< "LB_REF_TYPE *, int *, int *, int * ) must be implemented."
	<< std::endl;

  *ierr = LB_FATAL;
}

void Zoltan_QueryObject::Child_Weight (	void * data,
					int num_gid_entities,
					int num_lid_entities,
					LB_ID_PTR global_id,
					LB_ID_PTR local_id,
					int weight_dim,
					float * object_weight,
					int * ierr )
{
  std::cout << "Error: Zoltan_QueryObject::Child_Weight( void *, int, int, "
	<< "LB_ID_PTR, LB_ID_PTR, int, float *, int * ) must be implemented."
	<< std::endl;

  *ierr = LB_FATAL;
}
