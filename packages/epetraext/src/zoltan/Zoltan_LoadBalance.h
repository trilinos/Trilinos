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

#ifndef ZOLTAN_LOADBALANCE_H_
#define ZOLTAN_LOADBALANCE_H_

#include <EpetraExt_ConfigDefs.h>

#include <zoltan.h>

#include <string>

namespace Zoltan {

class QueryObject;
class MigrationObject;

//! Zoltan::LoadBalance: A class for interfacing the Load Balancing functions of the Zoltan library in a C++/Object Oriented environment.

/*! The Zoltan::LoadBalance class is a wrapper for the C functions at the top level
 interface of Zoltan.  The primary differences include the removal the Zoltan_Struct
 parameter since this object is now stored in the class and the use of user derived
 classes from Zoltan::QueryObject and Zoltan::MigrationObject for support of the
 "callback" functionality used by Zoltan.
*/

class EPETRAEXT_DEPRECATED LoadBalance
{

public:

 //@{ \name Constructors/Destructors.

 //! Constructor
 /*! This constructor replaces the Zoltan_Initialize call. Params are the same.
 */
 LoadBalance   ( int argc = 0,
                 char ** argv = 0,
                 float * ver = 0 );

 //! Destructor
 ~LoadBalance();

 //@}

 //@{ \name General Load Balance Functionality
 
 //! Replaces Zoltan_Create
 int Create( MPI_Comm communicator );

 //! Register UserDerived/ApplicationSpecific Query Object to support callbacks
 int Set_QueryObject( QueryObject * query_obj_ptr );

 //! Register UserDerived/ApplicationSpecific Migration Object to support callbacks
 int Set_MigrationObject( MigrationObject * migration_obj_ptr );

 //! Replaces Zoltan_Set_Param
 int Set_Param( const std::string & param, const std::string & value );

 //! Replaces Zoltan_LB_Balance
 int Balance  ( int * changes,
                int * num_gid_entries,
                int * num_lid_entries,
		int * num_import,
		ZOLTAN_ID_PTR * import_global_ids,
		ZOLTAN_ID_PTR * import_local_ids,
		int ** import_procs,
		int * num_export,
		ZOLTAN_ID_PTR * export_global_ids,
		ZOLTAN_ID_PTR * export_local_ids,
		int ** export_procs );

#ifdef ZOLTAN_ORDER
 int Order    ( int * num_gid_entries,
                int * num_lid_entries,
                int num_objs,
		ZOLTAN_ID_PTR global_ids,
		ZOLTAN_ID_PTR local_ids,
                int * rank,
                int * iperm );
#endif

 //! Replaces Zoltan_LB_Eval
 void Evaluate( int print_stats,
		int * num_objects,
		float * object_weights,
    		int * num_cut,
    		float * cut_weights,
		int * num_boundary_objects,
		int * num_adj_procs );

 //! Replaces Zoltan_LB_Free_Data
 int Free_Data( ZOLTAN_ID_PTR * import_global_ids,
		ZOLTAN_ID_PTR * import_local_ids,
		int ** import_procs,
		ZOLTAN_ID_PTR * export_global_ids,
		ZOLTAN_ID_PTR * export_local_ids,
		int ** export_procs );

 //@}

 //@{ \name Support for direct access to Zoltan callback functionality

#ifdef ZOLTAN_OLD_CALLBACK

 //! Old style callback support
 int Set_CallBack_Fn  ( ZOLTAN_FN_TYPE fn_type,
			void * fn_ptr,
			void * data = 0 );

#else /* ZOLTAN_OLD_CALLBACK */

 //! Individual callback support

 //!
 int Set_Num_Edges_Fn        ( ZOLTAN_NUM_EDGES_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Edge_List_Fn        ( ZOLTAN_EDGE_LIST_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Num_Geom_Fn         ( ZOLTAN_NUM_GEOM_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Geom_Fn             ( ZOLTAN_GEOM_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Num_Obj_Fn          ( ZOLTAN_NUM_OBJ_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Obj_List_Fn         ( ZOLTAN_OBJ_LIST_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_First_Obj_Fn        ( ZOLTAN_FIRST_OBJ_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Next_Obj_Fn         ( ZOLTAN_NEXT_OBJ_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Num_Border_Obj_Fn   ( ZOLTAN_NUM_BORDER_OBJ_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Border_Obj_List_Fn  ( ZOLTAN_BORDER_OBJ_LIST_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_First_Border_Obj_Fn ( ZOLTAN_FIRST_BORDER_OBJ_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Next_Border_Obj_Fn  ( ZOLTAN_NEXT_BORDER_OBJ_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Num_Coarse_Obj_Fn   ( ZOLTAN_NUM_COARSE_OBJ_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Coarse_Obj_List_Fn  ( ZOLTAN_COARSE_OBJ_LIST_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_First_Coarse_Obj_Fn ( ZOLTAN_FIRST_COARSE_OBJ_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Next_Coarse_Obj_Fn  ( ZOLTAN_NEXT_COARSE_OBJ_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Num_Child_Fn        ( ZOLTAN_NUM_CHILD_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Child_List_Fn       ( ZOLTAN_CHILD_LIST_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Child_Weight_Fn     ( ZOLTAN_CHILD_WEIGHT_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Pre_Migrate_Fn      ( ZOLTAN_PRE_MIGRATE_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Mid_Migrate_Fn      ( ZOLTAN_MID_MIGRATE_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Post_Migrate_Fn     ( ZOLTAN_POST_MIGRATE_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Obj_Size_Fn         ( ZOLTAN_OBJ_SIZE_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Pack_Obj_Fn         ( ZOLTAN_PACK_OBJ_FN * fn_ptr,
                              void * data = 0 );
 ///
 int Set_Unpack_Obj_Fn       ( ZOLTAN_UNPACK_OBJ_FN * fn_ptr,
                              void * data = 0 );

#endif /* ZOLTAN_OLD_CALLBACK */

 //@}

 //@{ \name Decomposition Augmentation

 //! Replaces Zoltan_LB_Point_Assign
 int Point_Assign     ( double * coords,
			int * proc );

 //! Replaces Zoltan_LB_Box_Assign
 int Box_Assign       ( double xmin,
			double ymin,
			double zmin,
			double xmax,
			double ymax,
			double zmax,
			int * procs,
			int * numprocs );

 //@}
 
 //@{ \name Migration Funtionality

 //! Replaces Zoltan_Compute_Destinations
 int Compute_Destinations     ( int num_import,
				ZOLTAN_ID_PTR import_global_ids,
				ZOLTAN_ID_PTR import_local_ids,
				int * import_procs,
				int * num_export,
				ZOLTAN_ID_PTR * export_global_ids,
				ZOLTAN_ID_PTR * export_local_ids,
				int ** export_procs );

 //! Replaces Zoltan_Help_Migrate
 int Help_Migrate     ( int num_import,
			ZOLTAN_ID_PTR import_global_ids,
			ZOLTAN_ID_PTR import_local_ids,
			int * import_procs,
			int num_export,
			ZOLTAN_ID_PTR export_global_ids,
			ZOLTAN_ID_PTR export_local_ids,
			int * export_procs );

 //@}

 //@{ \name Extra
 
 //! Direct access to underlying Zoltan_Struct Object
 Zoltan_Struct * Return_Zoltan_Struct()
 { return LoadBalance_Ptr_; }

 //@}

private:

 static int ObjectCount;

 int ObjectID;

 Zoltan_Struct * LoadBalance_Ptr_; 

 QueryObject * QueryObject_Ptr_;
 MigrationObject * MigrationObject_Ptr_;

};

} //namespace Zoltan

#endif
