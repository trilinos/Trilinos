//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_QueryFunctions.C$
//
// Purpose        : Static methods which are directly registered with
//                  Zoltan.  They us the static container to access
//                  the dynamic object methods.
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, Parallel Computational Sciences
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

#include <Zoltan_QueryFunctions.h>
#include <Zoltan_QueryContainer.h>
#include <Zoltan_QueryObject.h>

  //General Functions
int Zoltan_QueryFunctions::Number_Objects     (	void * data,
						int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  return obj_ptr->Number_Objects( data, ierr );
}

void Zoltan_QueryFunctions::Object_List       (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						LB_ID_PTR global_ids,
						LB_ID_PTR local_ids, 
						int weight_dim,
						float * object_weights,
						int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  obj_ptr->Object_List( data, num_gid_entries, num_lid_entries, global_ids,
	local_ids, weight_dim, object_weights, ierr );
}

int Zoltan_QueryFunctions::First_Object       (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						LB_ID_PTR first_global_id,
						LB_ID_PTR first_local_id, 
						int weight_dim,
						float * first_weight,
						int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  return obj_ptr->First_Object( data, num_gid_entries, num_lid_entries,
	first_global_id, first_local_id, weight_dim, first_weight, ierr );
}

int Zoltan_QueryFunctions::Next_Object        (	void * data,
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
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  return obj_ptr->Next_Object( data, num_gid_entries, num_lid_entries,
	global_id, local_id, next_global_id, next_local_id,
	weight_dim, next_weight, ierr );
}

int Zoltan_QueryFunctions::Number_Border_Objects    ( void * data,
						      int number_neighbor_procs,
						      int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  return obj_ptr->Number_Border_Objects( data, number_neighbor_procs, ierr );
}

void Zoltan_QueryFunctions::Border_Object_List     ( void * data,
						     int num_gid_entries,
						     int num_lid_entries,
						     int number_neighbor_procs, 
						     LB_ID_PTR global_ids,
						     LB_ID_PTR local_ids,
						     int weight_dim, 
						     float * object_weights,
						     int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  obj_ptr->Border_Object_List( data, num_gid_entries, num_lid_entries,
	number_neighbor_procs, global_ids, local_ids, weight_dim,
	object_weights, ierr );
}


int Zoltan_QueryFunctions::First_Border_Object     ( void * data,
						     int num_gid_entries,
						     int num_lid_entries,
						     int number_neighbor_procs, 
						     LB_ID_PTR first_global_id,
						     LB_ID_PTR first_local_id,
						     int weight_dim,
						     float * first_weight,
						     int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  return obj_ptr->First_Border_Object( data, num_gid_entries, num_lid_entries,
	number_neighbor_procs, first_global_id, first_local_id, weight_dim,
	first_weight, ierr );
}

int Zoltan_QueryFunctions::Next_Border_Object      ( void * data,
						     int num_gid_entries,
						     int num_lid_entries,
						     LB_ID_PTR global_id,
						     LB_ID_PTR local_id,
						     int number_neighbor_procs,
						     LB_ID_PTR next_global_id,
						     LB_ID_PTR next_local_id,
						     int weight_dim,
						     float * next_weight,
						     int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  return obj_ptr->Next_Border_Object( data, num_gid_entries, num_lid_entries,
	global_id, local_id, number_neighbor_procs, next_global_id,
	next_local_id, weight_dim, next_weight, ierr );
}

  //Geometry Based Functions
int Zoltan_QueryFunctions::Number_Geometry_Objects    (	void * data,
							int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  return obj_ptr->Number_Geometry_Objects( data, ierr );
}

void Zoltan_QueryFunctions::Geometry_Values   (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						LB_ID_PTR global_id,
						LB_ID_PTR local_id,
						double * geometry_vector,
						int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  obj_ptr->Geometry_Values( data, num_gid_entries, num_lid_entries, global_id,
	local_id, geometry_vector, ierr );
}

  //Graph Based Functions
int Zoltan_QueryFunctions::Number_Edges       (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						LB_ID_PTR global_id,
						LB_ID_PTR local_id,
						int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  return obj_ptr->Number_Edges( data, num_gid_entries, num_lid_entries,
	global_id, local_id, ierr );
}

void Zoltan_QueryFunctions::Edge_List (	void * data,
					int num_gid_entries,
					int num_lid_entries,
					LB_ID_PTR global_id,
					LB_ID_PTR local_id,
					LB_ID_PTR neighbor_global_ids,
					int * neighbor_procs,
					int weight_dim,
					float * edge_weights,
					int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  obj_ptr->Edge_List( data, num_gid_entries, num_lid_entries, global_id,
	local_id, neighbor_global_ids, neighbor_procs, weight_dim,
	edge_weights, ierr );
}

  //Tree Based Functions
int Zoltan_QueryFunctions::Number_Coarse_Objects      (	void * data,
							int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  return obj_ptr->Number_Coarse_Objects( data, ierr );
}

void Zoltan_QueryFunctions::Coarse_Object_List        (	void * data,
							int num_gid_entries,
							int num_lid_entries,
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
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  obj_ptr->Coarse_Object_List( data, num_gid_entries, num_lid_entries,
	global_ids, local_ids, assigned, number_vertices, vertices,
	in_order, in_vertex, out_vertex, ierr );
}

int Zoltan_QueryFunctions::First_Coarse_Object      ( void * data,
						      int num_gid_entries,
						      int num_lid_entries,
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
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  return obj_ptr->First_Coarse_Object( data, num_gid_entries, num_lid_entries,
	first_global_id, first_local_id, assigned, number_vertices, vertices,
	in_order, in_vertex, out_vertex, ierr );
}

int Zoltan_QueryFunctions::Next_Coarse_Object       ( void * data,
						      int num_gid_entries,
						      int num_lid_entries,
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
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  return obj_ptr->Next_Coarse_Object( data, num_gid_entries, num_lid_entries,
	global_id, local_id, next_global_id, next_local_id,
	assigned, number_vertices, vertices, in_vertex, out_vertex,
	ierr );
}

int Zoltan_QueryFunctions::Number_Children    (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						LB_ID_PTR global_id,
						LB_ID_PTR local_id,
						int * ierr)
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  return obj_ptr->Number_Children( data, num_gid_entries, num_lid_entries,
	global_id, local_id, ierr );
}

void Zoltan_QueryFunctions::Child_List        (	void * data,
						int num_gid_entries,
						int num_lid_entries,
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
						int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  obj_ptr->Child_List( data, num_gid_entries, num_lid_entries,
	parent_global_id, parent_local_id,
	child_global_ids, child_local_ids, assigned, number_vertices,
	vertices, reference_type, in_vertex, out_vertex, ierr );
}

void Zoltan_QueryFunctions::Child_Weight      (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						LB_ID_PTR global_id,
						LB_ID_PTR local_id,
						int weight_dim,
						float * object_weight,
						int * ierr )
{
  Zoltan_QueryObject * obj_ptr = Zoltan_QueryContainer::getQueryObject(
	Zoltan_QueryContainer::getQueryID() );

  obj_ptr->Child_Weight( data, num_gid_entries, num_lid_entries,
	global_id, local_id, weight_dim, object_weight, ierr );
}

