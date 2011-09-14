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

#ifndef ZOLTAN_QUERYOBJECT_H_
#define ZOLTAN_QUERYOBJECT_H_

#include "EpetraExt_ConfigDefs.h"

#include <zoltan.h>

namespace Zoltan {

//! Zoltan::QueryObject: A base class from which the user can derive an application specific support class for Zoltan's query callback functions.

/*! As with Zoltan, the user only need implement those methods used by Zoltan during
    their application execution.  If Zoltan calls an unimplemented method, a fatal error
    will be generated.
*/

class EPETRAEXT_DEPRECATED QueryObject
{

public:

  //@{ \name Virtual General Query Methods.

  //! Supports ZOLTAN_NUM_OBJ_FN_TYPE
  virtual int Number_Objects  (	void * data,
				int * ierr );

  //! Supports ZOLTAN_OBJ_LIST_FN_TYPE
  virtual void Object_List    (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_ids, 
				ZOLTAN_ID_PTR local_ids,
				int weight_dim,
				float * object_weights,
				int * ierr );

  //! Supports ZOLTAN_FIRST_OBJ_FN_TYPE
  virtual int First_Object    (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR first_global_id, 
				ZOLTAN_ID_PTR first_local_id,
				int weight_dim,
				float * first_weight,
				int * ierr );

  //! Supports ZOLTAN_NEXT_OBJ_FN_TYPE
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

  //! Supports ZOLTAN_BORDER_OBJ_FN_TYPE
  virtual int Number_Border_Objects   (	void * data,
					int number_neighbor_procs,
					int * ierr );

  //! Supports ZOLTAN_BORDER_OBJ_LIST_FN_TYPE
  virtual void Border_Object_List     (	void * data,
					int num_gid_entries,
					int num_lid_entries,
					int number_neighbor_procs,
					ZOLTAN_ID_PTR global_ids,
					ZOLTAN_ID_PTR local_ids,
					int weight_dim,
					float * object_weights,
					int * ierr );

  //! Supports ZOLTAN_FIRST_BORDER_OBJ_FN_TYPE
  virtual int First_Border_Object     (	void * data,
					int num_gid_entries,
					int num_lid_entries,
					int number_neighbor_procs,
					ZOLTAN_ID_PTR first_global_id,
					ZOLTAN_ID_PTR first_local_id,
					int weight_dim,
					float * first_weight,
					int * ierr );

  //! Supports ZOLTAN_NEXT_BORDER_OBJ_FN_TYPE
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

  //@}

  //@{ \name Geometry Based Functions
 
  //! Supports ZOLTAN_NUM_GEOM_OBJ_FN_TYPE
  virtual int Number_Geometry_Objects (	void * data,
					int * ierr );

  //! Supports ZOLTAN_GEOM_VALUES_FN_TYPE
  virtual void Geometry_Values        ( void * data,
					int num_gid_entries,
					int num_lid_entries,
					ZOLTAN_ID_PTR global_id, 
					ZOLTAN_ID_PTR local_id,
					double * geometry_vector,
					int * ierr );

  //@}

  //@{ \name Graph Based Functions
 
  //! Supports ZOLTAN_NUM_EDGES_FN_TYPE
  virtual int Number_Edges    (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id,
				ZOLTAN_ID_PTR local_id,
				int * ierr );

  //! Supports ZOLTAN_EDGE_LIST_FN_TYPE
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

  //@}

  //@{ Tree Based Functions
 
  //! Supports ZOLTAN_NUM_COARSE_OBJ_FN_TYPE
  virtual int Number_Coarse_Objects   (	void * data,
					int * ierr );

  //! Supports ZOLTAN_COARSE_OBJ_LIST_FN_TYPE
  virtual void Coarse_Object_List     (	void * data,
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

  //! Supports ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE
  virtual int First_Coarse_Object     ( void * data,
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

  //! Supports ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE
  virtual int Next_Coarse_Object      ( void * data,
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

  //! Supports ZOLTAN_NUM_CHILD_FN_TYPE
  virtual int Number_Children (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id,
				ZOLTAN_ID_PTR local_id,
				int * ierr );

  //! Supports ZOLTAN_CHILD_LIST_FN_TYPE
  virtual void Child_List     (	void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR parent_global_id,
				ZOLTAN_ID_PTR parent_local_id,
				ZOLTAN_ID_PTR child_global_id,
				ZOLTAN_ID_PTR child_local_id,
				int * assigned,
				int * number_vertices,
				ZOLTAN_ID_PTR vertices,
				ZOLTAN_REF_TYPE * reference_type,
				ZOLTAN_ID_PTR in_vertex,
				ZOLTAN_ID_PTR out_vertex,
				int * ierr  );

  //! Supports ZOLTAN_CHILD_WEIGHT_FN_TYPE
  virtual void Child_Weight   ( void * data,
				int num_gid_entries,
				int num_lid_entries,
				ZOLTAN_ID_PTR global_id,
				ZOLTAN_ID_PTR local_id,
				int weight_dim,
				float * object_weight,
				int * ierr );

  //@}

};

} //namespace Zoltan

#endif
