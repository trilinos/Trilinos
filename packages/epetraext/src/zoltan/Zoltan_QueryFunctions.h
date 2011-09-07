//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_QueryFunctions.h$
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

#ifndef ZOLTAN_QUERYFUNCTIONS_H_
#define ZOLTAN_QUERYFUNCTIONS_H_

#include "EpetraExt_ConfigDefs.h"

#include <zoltan.h>

namespace Zoltan {

class EPETRAEXT_DEPRECATED QueryFunctions
{

public:

  //General Functions
  static int Number_Objects   (	void * data,
				int * ierr );

  static void Object_List     (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_ids,
				ZOLTAN_ID_PTR local_ids, 
				int weight_dim,
				float * object_weights,
				int * ierr );

  static int First_Object     ( void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR first_global_id,
				ZOLTAN_ID_PTR first_local_id, 
				int weight_dim,
				float * first_weight,
				int * ierr );

  static int Next_Object      ( void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id,
				ZOLTAN_ID_PTR local_id, 
				ZOLTAN_ID_PTR next_global_id,
				ZOLTAN_ID_PTR next_local_id,
				int weight_dim,
				float * next_weight,
				int * ierr );

  static int Number_Border_Objects    (	void * data,
					int number_neighbor_procs,
					int * ierr );

  static void Border_Object_List      ( void * data,
					int num_gid_entries,
					int num_lid_entries,
					int number_neighbor_procs, 
					ZOLTAN_ID_PTR global_ids,
					ZOLTAN_ID_PTR local_ids,
					int weight_dim, 
					float * object_weights,
					int * ierr );

  static int First_Border_Object      ( void * data,
					int num_gid_entries,
					int num_lid_entries,
					int number_neighbor_procs, 
					ZOLTAN_ID_PTR first_global_id,
					ZOLTAN_ID_PTR first_local_id,
					int weight_dim,
					float * first_weight,
					int * ierr );

  static int Next_Border_Object       ( void * data,
					int num_gid_entries,
					int num_lid_entries,
					ZOLTAN_ID_PTR global_id,
					ZOLTAN_ID_PTR local_id,
					int number_neighbor_procs,
					ZOLTAN_ID_PTR next_global_id,
					ZOLTAN_ID_PTR next_local_id,
					int weight_dim,
					float * next_weight,
					int * ierr );

  //Geometry Based Functions
  static int Number_Geometry_Objects  ( void * data,
					int * ierr );

  static void Geometry_Values ( void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id, 
				ZOLTAN_ID_PTR local_id,
				double * geometry_vector,
				int * ierr );

  //Graph Based Functions
  static int Number_Edges     ( void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id, 
				ZOLTAN_ID_PTR local_id,
				int * ierr );

  static void Edge_List       ( void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id, 
				ZOLTAN_ID_PTR local_id,
				ZOLTAN_ID_PTR neighbor_global_ids,
				int * neighbor_procs,
				int weight_dim,
				float * edge_weights,
				int * ierr );

  //Tree Based Functions
  static int Number_Coarse_Objects    (	void * data,
					int * ierr );

  static void Coarse_Object_List      ( void * data,
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
					int * ierr );

  static int First_Coarse_Object      (	void * data,
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
					int * ierr );

  static int Next_Coarse_Object       ( void * data,
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
					int * ierr );

  static int Number_Children  (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id, 
				ZOLTAN_ID_PTR local_id,
				int * ierr);

  static void Child_List      ( void * data,
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
				int * ierr );

  static void Child_Weight    (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id, 
				ZOLTAN_ID_PTR local_id,
				int weight_dim,
				float * object_weight,
				int * ierr );

};

} //namespace Zoltan

#endif
