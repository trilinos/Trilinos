//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_QueryObject.h$
//
// Purpose        : Base Class for dynamic versions of query
//		    functions to be used by Zoltan.  The application
//		    requires a class derived from this base to
//		    be instantiated and registered with the
//		    Zoltan_LoadBalance object.
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

#ifndef ZOLTAN_QUERYOBJECT_H_
#define ZOLTAN_QUERYOBJECT_H_

// #include <lbi_const.h>
#include <zoltan.h>

class Zoltan_QueryObject
{

public:

  //General Functions
  virtual int Number_Objects  (	void * data,
				int * ierr );

  virtual void Object_List    (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_ids, 
				ZOLTAN_ID_PTR local_ids,
				int weight_dim,
				float * object_weights,
				int * ierr );

  virtual int First_Object    (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR first_global_id, 
				ZOLTAN_ID_PTR first_local_id,
				int weight_dim,
				float * first_weight,
				int * ierr );

  virtual int Next_Object     (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id, 
				ZOLTAN_ID_PTR local_id,
				ZOLTAN_ID_PTR next_global_id, 
				ZOLTAN_ID_PTR next_local_id,
				int weight_dim,
				float * next_weight,
				int * ierr );

  virtual int Number_Border_Objects   (	void * data,
					int number_neighbor_procs,
					int * ierr );

  virtual void Border_Object_List     (	void * data,
					int num_gid_entries,
					int num_lid_entries,
					int number_neighbor_procs,
					ZOLTAN_ID_PTR global_ids,
					ZOLTAN_ID_PTR local_ids,
					int weight_dim,
					float * object_weights,
					int * ierr );

  virtual int First_Border_Object     (	void * data,
					int num_gid_entries,
					int num_lid_entries,
					int number_neighbor_procs,
					ZOLTAN_ID_PTR first_global_id,
					ZOLTAN_ID_PTR first_local_id,
					int weight_dim,
					float * first_weight,
					int * ierr );

  virtual int Next_Border_Object      ( void * data,
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
  virtual int Number_Geometry_Objects (	void * data,
					int * ierr );

  virtual void Geometry_Values        ( void * data,
					int num_gid_entries,
					int num_lid_entries,
					ZOLTAN_ID_PTR global_id, 
					ZOLTAN_ID_PTR local_id,
					double * geometry_vector,
					int * ierr );

  //Graph Based Functions
  virtual int Number_Edges    (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id,
				ZOLTAN_ID_PTR local_id,
				int * ierr );

  virtual void Edge_List      (	void * data,
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
  virtual int Number_Coarse_Objects   (	void * data,
					int * ierr );

  virtual void Coarse_Object_List     (	void * data,
					int num_gid_entries,
					int num_lid_entries,
					ZOLTAN_ID_PTR global_ids,
					ZOLTAN_ID_PTR local_ids,
					int * assigned,
					int * number_vertices,
					int * vertices,
					int * in_order,
					int * in_vertex,
					int * out_vertex,
					int * ierr );

  virtual int First_Coarse_Object     ( void * data,
					int num_gid_entries,
					int num_lid_entries,
					ZOLTAN_ID_PTR first_global_id,
					ZOLTAN_ID_PTR first_local_id,
					int * assigned,
					int * number_vertices,
					int * vertices,
					int * in_order,
					int * in_vertex,
					int * out_vertex,
					int * ierr );

  virtual int Next_Coarse_Object      ( void * data,
					int num_gid_entries,
					int num_lid_entries,
					ZOLTAN_ID_PTR global_id,
					ZOLTAN_ID_PTR local_id,
					ZOLTAN_ID_PTR next_global_id,
					ZOLTAN_ID_PTR next_local_id,
					int * assigned,
					int * number_vertices,
					int * vertices,
					int * in_vertex,
					int * out_vertex,
					int * ierr );

  virtual int Number_Children (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id,
				ZOLTAN_ID_PTR local_id,
				int * ierr );

  virtual void Child_List     (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR parent_global_id,
				ZOLTAN_ID_PTR parent_local_id,
				ZOLTAN_ID_PTR child_global_id,
				ZOLTAN_ID_PTR child_local_id,
				int * assigned,
				int * number_vertices,
				int * vertices,
				ZOLTAN_REF_TYPE * reference_type,
				int * in_vertex,
				int * out_vertex,
				int * ierr  );

  virtual void Child_Weight   ( void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id,
				ZOLTAN_ID_PTR local_id,
				int weight_dim,
				float * object_weight,
				int * ierr );

};

#endif
