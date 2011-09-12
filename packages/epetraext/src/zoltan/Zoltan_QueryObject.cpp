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

#include <Zoltan_QueryObject.h>

  //General Functions
EPETRAEXT_DEPRECATED
int Zoltan::QueryObject::Number_Objects        (	void * data,
						int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Number_Objects( void *, "
	<< "int * ) must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

EPETRAEXT_DEPRECATED
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

EPETRAEXT_DEPRECATED
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

EPETRAEXT_DEPRECATED
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

EPETRAEXT_DEPRECATED
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

EPETRAEXT_DEPRECATED
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

EPETRAEXT_DEPRECATED
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

EPETRAEXT_DEPRECATED
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
EPETRAEXT_DEPRECATED
int Zoltan::QueryObject::Number_Geometry_Objects       (	void * data,
							int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Number_Geometry_Objects( void *, "
	<< "int * ) must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

EPETRAEXT_DEPRECATED
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
EPETRAEXT_DEPRECATED
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

EPETRAEXT_DEPRECATED
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

  //Tree Based Functions
EPETRAEXT_DEPRECATED
int Zoltan::QueryObject::Number_Coarse_Objects (	void * data,
						int * ierr )
{
  std::cout << "Error: Zoltan::QueryObject::Number_Coarse_Objects( void *, "
	<< "int * ) must be implemented." << std::endl;

  *ierr = ZOLTAN_FATAL;

  return 0;
}

EPETRAEXT_DEPRECATED
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

EPETRAEXT_DEPRECATED
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

EPETRAEXT_DEPRECATED
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

EPETRAEXT_DEPRECATED
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

EPETRAEXT_DEPRECATED
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

EPETRAEXT_DEPRECATED
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
