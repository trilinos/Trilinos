/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

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
EPETRAEXT_DEPRECATED
int Zoltan::QueryFunctions::Number_Objects     (	void * data,
						int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  return obj_ptr->Number_Objects( data, ierr );
}

EPETRAEXT_DEPRECATED
void Zoltan::QueryFunctions::Object_List       (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						ZOLTAN_ID_PTR global_ids,
						ZOLTAN_ID_PTR local_ids, 
						int weight_dim,
						float * object_weights,
						int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  obj_ptr->Object_List( data, num_gid_entries, num_lid_entries, global_ids,
	local_ids, weight_dim, object_weights, ierr );
}

EPETRAEXT_DEPRECATED
int Zoltan::QueryFunctions::First_Object       (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						ZOLTAN_ID_PTR first_global_id,
						ZOLTAN_ID_PTR first_local_id, 
						int weight_dim,
						float * first_weight,
						int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  return obj_ptr->First_Object( data, num_gid_entries, num_lid_entries,
	first_global_id, first_local_id, weight_dim, first_weight, ierr );
}

EPETRAEXT_DEPRECATED
int Zoltan::QueryFunctions::Next_Object        (	void * data,
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
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  return obj_ptr->Next_Object( data, num_gid_entries, num_lid_entries,
	global_id, local_id, next_global_id, next_local_id,
	weight_dim, next_weight, ierr );
}

EPETRAEXT_DEPRECATED
int Zoltan::QueryFunctions::Number_Border_Objects    ( void * data,
						      int number_neighbor_procs,
						      int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  return obj_ptr->Number_Border_Objects( data, number_neighbor_procs, ierr );
}

EPETRAEXT_DEPRECATED
void Zoltan::QueryFunctions::Border_Object_List     ( void * data,
						     int num_gid_entries,
						     int num_lid_entries,
						     int number_neighbor_procs, 
						     ZOLTAN_ID_PTR global_ids,
						     ZOLTAN_ID_PTR local_ids,
						     int weight_dim, 
						     float * object_weights,
						     int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  obj_ptr->Border_Object_List( data, num_gid_entries, num_lid_entries,
	number_neighbor_procs, global_ids, local_ids, weight_dim,
	object_weights, ierr );
}


EPETRAEXT_DEPRECATED
int Zoltan::QueryFunctions::First_Border_Object     ( void * data,
						     int num_gid_entries,
						     int num_lid_entries,
						     int number_neighbor_procs, 
						     ZOLTAN_ID_PTR first_global_id,
						     ZOLTAN_ID_PTR first_local_id,
						     int weight_dim,
						     float * first_weight,
						     int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  return obj_ptr->First_Border_Object( data, num_gid_entries, num_lid_entries,
	number_neighbor_procs, first_global_id, first_local_id, weight_dim,
	first_weight, ierr );
}

EPETRAEXT_DEPRECATED
int Zoltan::QueryFunctions::Next_Border_Object      ( void * data,
						     int num_gid_entries,
						     int num_lid_entries,
						     ZOLTAN_ID_PTR global_id,
						     ZOLTAN_ID_PTR local_id,
						     int number_neighbor_procs,
						     ZOLTAN_ID_PTR next_global_id,
						     ZOLTAN_ID_PTR next_local_id,
						     int weight_dim,
						     float * next_weight,
						     int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  return obj_ptr->Next_Border_Object( data, num_gid_entries, num_lid_entries,
	global_id, local_id, number_neighbor_procs, next_global_id,
	next_local_id, weight_dim, next_weight, ierr );
}

  //Geometry Based Functions
EPETRAEXT_DEPRECATED
int Zoltan::QueryFunctions::Number_Geometry_Objects    (	void * data,
							int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  return obj_ptr->Number_Geometry_Objects( data, ierr );
}

EPETRAEXT_DEPRECATED
void Zoltan::QueryFunctions::Geometry_Values   (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						ZOLTAN_ID_PTR global_id,
						ZOLTAN_ID_PTR local_id,
						double * geometry_vector,
						int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  obj_ptr->Geometry_Values( data, num_gid_entries, num_lid_entries, global_id,
	local_id, geometry_vector, ierr );
}

  //Graph Based Functions
EPETRAEXT_DEPRECATED
int Zoltan::QueryFunctions::Number_Edges       (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						ZOLTAN_ID_PTR global_id,
						ZOLTAN_ID_PTR local_id,
						int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  return obj_ptr->Number_Edges( data, num_gid_entries, num_lid_entries,
	global_id, local_id, ierr );
}

EPETRAEXT_DEPRECATED
void Zoltan::QueryFunctions::Edge_List (	void * data,
					int num_gid_entries,
					int num_lid_entries,
					ZOLTAN_ID_PTR global_id,
					ZOLTAN_ID_PTR local_id,
					ZOLTAN_ID_PTR neighbor_global_ids,
					int * neighbor_procs,
					int weight_dim,
					float * edge_weights,
					int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  obj_ptr->Edge_List( data, num_gid_entries, num_lid_entries, global_id,
	local_id, neighbor_global_ids, neighbor_procs, weight_dim,
	edge_weights, ierr );
}

  //Tree Based Functions
EPETRAEXT_DEPRECATED
int Zoltan::QueryFunctions::Number_Coarse_Objects      (	void * data,
							int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  return obj_ptr->Number_Coarse_Objects( data, ierr );
}

EPETRAEXT_DEPRECATED
void Zoltan::QueryFunctions::Coarse_Object_List        (	void * data,
							int num_gid_entries,
							int num_lid_entries,
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
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  obj_ptr->Coarse_Object_List( data, num_gid_entries, num_lid_entries,
	global_ids, local_ids, assigned, number_vertices, vertices,
	in_order, in_vertex, out_vertex, ierr );
}

EPETRAEXT_DEPRECATED
int Zoltan::QueryFunctions::First_Coarse_Object      ( void * data,
						      int num_gid_entries,
						      int num_lid_entries,
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
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  return obj_ptr->First_Coarse_Object( data, num_gid_entries, num_lid_entries,
	first_global_id, first_local_id, assigned, number_vertices, vertices,
	in_order, in_vertex, out_vertex, ierr );
}

EPETRAEXT_DEPRECATED
int Zoltan::QueryFunctions::Next_Coarse_Object       ( void * data,
						      int num_gid_entries,
						      int num_lid_entries,
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
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  return obj_ptr->Next_Coarse_Object( data, num_gid_entries, num_lid_entries,
	global_id, local_id, next_global_id, next_local_id,
	assigned, number_vertices, vertices, in_vertex, out_vertex,
	ierr );
}

EPETRAEXT_DEPRECATED
int Zoltan::QueryFunctions::Number_Children    (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						ZOLTAN_ID_PTR global_id,
						ZOLTAN_ID_PTR local_id,
						int * ierr)
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  return obj_ptr->Number_Children( data, num_gid_entries, num_lid_entries,
	global_id, local_id, ierr );
}

EPETRAEXT_DEPRECATED
void Zoltan::QueryFunctions::Child_List        (	void * data,
						int num_gid_entries,
						int num_lid_entries,
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
						int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  obj_ptr->Child_List( data, num_gid_entries, num_lid_entries,
	parent_global_id, parent_local_id,
	child_global_ids, child_local_ids, assigned, number_vertices,
	vertices, reference_type, in_vertex, out_vertex, ierr );
}

EPETRAEXT_DEPRECATED
void Zoltan::QueryFunctions::Child_Weight      (	void * data,
						int num_gid_entries,
						int num_lid_entries,
						ZOLTAN_ID_PTR global_id,
						ZOLTAN_ID_PTR local_id,
						int weight_dim,
						float * object_weight,
						int * ierr )
{
  Zoltan::QueryObject * obj_ptr = Zoltan::QueryContainer::getQueryObject(
	Zoltan::QueryContainer::getQueryID() );

  obj_ptr->Child_Weight( data, num_gid_entries, num_lid_entries,
	global_id, local_id, weight_dim, object_weight, ierr );
}

