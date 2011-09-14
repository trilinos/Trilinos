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
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_LoadBalance.C$
//
// Purpose        : C++ wrapper object for Zoltan Load Balance Routines.
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
// Revision Number: $$
//
// Revision Date  : $$
//
// Current Owner  : $$
//-------------------------------------------------------------------------

#include <Zoltan_LoadBalance.h>
#include <Zoltan_QueryContainer.h>
#include <Zoltan_QueryFunctions.h>
#include <Zoltan_MigrationContainer.h>
#include <Zoltan_MigrationFunctions.h>

int Zoltan::LoadBalance::ObjectCount = 0;

EPETRAEXT_DEPRECATED
Zoltan::LoadBalance::LoadBalance( int argc,
                                  char ** argv,
                                  float * ver )
	: ObjectID( 0 ),
	  LoadBalance_Ptr_( 0 ),
	  QueryObject_Ptr_( 0 ),
	  MigrationObject_Ptr_( 0 )
{
  ObjectID = ObjectCount++;

  int tmpReturn = Zoltan_Initialize( argc, argv, ver );
}

EPETRAEXT_DEPRECATED
Zoltan::LoadBalance::~LoadBalance()
{
  Zoltan_Destroy( &LoadBalance_Ptr_ );
}

 //Load Balance Calls
EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Create( MPI_Comm communicator )
{
  LoadBalance_Ptr_ = Zoltan_Create( communicator );

  if( !LoadBalance_Ptr_ )
    return ZOLTAN_FATAL;
  else
    return ZOLTAN_OK;
}

#ifdef ZOLTAN_OLD_CALLBACK

//Old callback support
//--------------------------------

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_CallBack_Fn       ( ZOLTAN_FN_TYPE fn_type,
	 					 void * fn_ptr,
						 void * data )
{
  return Zoltan_Set_Fn( LoadBalance_Ptr_, fn_type, fn_ptr, data );
}

#else /* ZOLTAN_OLD_CALLBACK */

//Individual callback registration
//--------------------------------

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Num_Edges_Fn      ( ZOLTAN_NUM_EDGES_FN * fn_ptr,
						 void * data )
{
  return Zoltan_Set_Num_Edges_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Edge_List_Fn      (	ZOLTAN_EDGE_LIST_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Edge_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Num_Geom_Fn      (	ZOLTAN_NUM_GEOM_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Num_Geom_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Geom_Fn      (	ZOLTAN_GEOM_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Geom_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Num_Obj_Fn      (	ZOLTAN_NUM_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Num_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Obj_List_Fn      (	ZOLTAN_OBJ_LIST_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Obj_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_First_Obj_Fn      (	ZOLTAN_FIRST_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_First_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Next_Obj_Fn      (	ZOLTAN_NEXT_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Next_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Num_Border_Obj_Fn (	ZOLTAN_NUM_BORDER_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Num_Border_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Border_Obj_List_Fn ( ZOLTAN_BORDER_OBJ_LIST_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Border_Obj_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_First_Border_Obj_Fn ( ZOLTAN_FIRST_BORDER_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_First_Border_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Next_Border_Obj_Fn ( ZOLTAN_NEXT_BORDER_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Next_Border_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Num_Coarse_Obj_Fn (	ZOLTAN_NUM_COARSE_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Num_Coarse_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Coarse_Obj_List_Fn ( ZOLTAN_COARSE_OBJ_LIST_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Coarse_Obj_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_First_Coarse_Obj_Fn ( ZOLTAN_FIRST_COARSE_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_First_Coarse_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Next_Coarse_Obj_Fn ( ZOLTAN_NEXT_COARSE_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Next_Coarse_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Num_Child_Fn      (	ZOLTAN_NUM_CHILD_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Num_Child_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Child_List_Fn ( ZOLTAN_CHILD_LIST_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Child_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Child_Weight_Fn ( ZOLTAN_CHILD_WEIGHT_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Child_Weight_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Pre_Migrate_Fn    (	ZOLTAN_PRE_MIGRATE_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Pre_Migrate_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Mid_Migrate_Fn    (	ZOLTAN_MID_MIGRATE_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Mid_Migrate_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Post_Migrate_Fn   (	ZOLTAN_POST_MIGRATE_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Post_Migrate_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Obj_Size_Fn       (	ZOLTAN_OBJ_SIZE_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Obj_Size_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Pack_Obj_Fn       (	ZOLTAN_PACK_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Pack_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Unpack_Obj_Fn     (	ZOLTAN_UNPACK_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Unpack_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

#endif /* ZOLTAN_OLD_CALLBACK */

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_QueryObject( Zoltan::QueryObject * query_obj_ptr )
{
  Zoltan::QueryContainer::registerQueryObject( ObjectID, query_obj_ptr );

#ifdef ZOLTAN_OLD_CALLBACK

  //LB_Set_Fn all Cstyle Functions
  Set_CallBack_Fn( ZOLTAN_NUM_EDGES_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Number_Edges),
	0 );
  Set_CallBack_Fn( ZOLTAN_EDGE_LIST_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Edge_List),
	0 );
  Set_CallBack_Fn( ZOLTAN_NUM_GEOM_FN_TYPE,
    reinterpret_cast<void *> (Zoltan::QueryFunctions::Number_Geometry_Objects),
	0 );
  Set_CallBack_Fn( ZOLTAN_GEOM_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Geometry_Values),
	0 );
  Set_CallBack_Fn( ZOLTAN_NUM_OBJ_FN_TYPE, 
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Number_Objects),
	0 );
  Set_CallBack_Fn( ZOLTAN_OBJ_LIST_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Object_List),
	0 );
  Set_CallBack_Fn( ZOLTAN_FIRST_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::First_Object),
	0 );
  Set_CallBack_Fn( ZOLTAN_NEXT_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Next_Object),
	0 );
  Set_CallBack_Fn( ZOLTAN_NUM_BORDER_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Number_Border_Objects),
	0 );
  Set_CallBack_Fn( ZOLTAN_BORDER_OBJ_LIST_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Border_Object_List),
	0 );
  Set_CallBack_Fn( ZOLTAN_FIRST_BORDER_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::First_Border_Object),
	0 );
  Set_CallBack_Fn( ZOLTAN_NEXT_BORDER_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Next_Border_Object),
	0 );
  Set_CallBack_Fn( ZOLTAN_NUM_COARSE_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Number_Coarse_Objects),
	0 );
  Set_CallBack_Fn( ZOLTAN_COARSE_OBJ_LIST_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Coarse_Object_List),
	0 );
  Set_CallBack_Fn( ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::First_Coarse_Object),
	0 );
  Set_CallBack_Fn( ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Next_Coarse_Object),
	0 );
  Set_CallBack_Fn( ZOLTAN_NUM_CHILD_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Number_Children),
	0 );
  Set_CallBack_Fn( ZOLTAN_CHILD_LIST_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Child_List),
	0 );
  Set_CallBack_Fn( ZOLTAN_CHILD_WEIGHT_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::QueryFunctions::Child_Weight),
	0 );

#else /* ZOLTAN_OLD_CALLBACK */

  Set_Num_Edges_Fn        ( Zoltan::QueryFunctions::Number_Edges,
                            0 );
  Set_Edge_List_Fn        ( Zoltan::QueryFunctions::Edge_List,
                            0 );
  Set_Num_Geom_Fn         ( Zoltan::QueryFunctions::Number_Geometry_Objects,
                            0 );
  Set_Geom_Fn             ( Zoltan::QueryFunctions::Geometry_Values,
                            0 );
  Set_Num_Obj_Fn          ( Zoltan::QueryFunctions::Number_Objects,
                            0 );
  Set_Obj_List_Fn         ( Zoltan::QueryFunctions::Object_List,
                            0 );
  Set_First_Obj_Fn        ( Zoltan::QueryFunctions::First_Object,
                            0 );
  Set_Next_Obj_Fn         ( Zoltan::QueryFunctions::Next_Object,
                            0 );
  Set_Num_Border_Obj_Fn   ( Zoltan::QueryFunctions::Number_Border_Objects,
                            0 );
  Set_Border_Obj_List_Fn  ( Zoltan::QueryFunctions::Border_Object_List,
                            0 );
  Set_First_Border_Obj_Fn ( Zoltan::QueryFunctions::First_Border_Object,
                            0 );
  Set_Next_Border_Obj_Fn  ( Zoltan::QueryFunctions::Next_Border_Object,
                            0 );
  Set_Num_Coarse_Obj_Fn   ( Zoltan::QueryFunctions::Number_Coarse_Objects,
                            0 );
  Set_Coarse_Obj_List_Fn  ( Zoltan::QueryFunctions::Coarse_Object_List,
                            0 );
  Set_First_Coarse_Obj_Fn ( Zoltan::QueryFunctions::First_Coarse_Object,
                            0 );
  Set_Next_Coarse_Obj_Fn  ( Zoltan::QueryFunctions::Next_Coarse_Object,
                            0 );
  Set_Num_Child_Fn        ( Zoltan::QueryFunctions::Number_Children,
                            0 );
  Set_Child_List_Fn       ( Zoltan::QueryFunctions::Child_List,
                            0 );
  Set_Child_Weight_Fn     ( Zoltan::QueryFunctions::Child_Weight,
                            0 );

#endif /* ZOLTAN_OLD_CALLBACK */

  return ZOLTAN_OK;
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_MigrationObject( 
			Zoltan::MigrationObject * migration_obj_ptr )
{
  Zoltan::MigrationContainer::registerMigrationObject( ObjectID,
	migration_obj_ptr );

#ifdef ZOLTAN_OLD_CALLBACK

  //LB_Set_Fn all Cstyle Functions
  Set_CallBack_Fn( ZOLTAN_OBJ_SIZE_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::MigrationFunctions::Object_Size),
        0 );
  Set_CallBack_Fn( ZOLTAN_PRE_MIGRATE_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::MigrationFunctions::Pre_Migrate),
        0 );
  Set_CallBack_Fn( ZOLTAN_MID_MIGRATE_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::MigrationFunctions::Mid_Migrate),
        0 );
  Set_CallBack_Fn( ZOLTAN_POST_MIGRATE_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::MigrationFunctions::Post_Migrate),
        0 );
  Set_CallBack_Fn( ZOLTAN_PACK_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::MigrationFunctions::Pack_Object),
        0 );
  Set_CallBack_Fn( ZOLTAN_UNPACK_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan::MigrationFunctions::Unpack_Object),
        0 );

#else /* ZOLTAN_OLD_CALLBACK */

  Set_Obj_Size_Fn        ( Zoltan::MigrationFunctions::Object_Size, 0 );
  Set_Pre_Migrate_Fn     ( Zoltan::MigrationFunctions::Pre_Migrate, 0 );
  Set_Mid_Migrate_Fn     ( Zoltan::MigrationFunctions::Mid_Migrate, 0 );
  Set_Post_Migrate_Fn    ( Zoltan::MigrationFunctions::Post_Migrate, 0 );
  Set_Pack_Obj_Fn        ( Zoltan::MigrationFunctions::Pack_Object, 0 );
  Set_Unpack_Obj_Fn      ( Zoltan::MigrationFunctions::Unpack_Object, 0 );

#endif /* ZOLTAN_OLD_CALLBACK */

  return ZOLTAN_OK;
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Set_Param( const std::string & param, const std::string & value )
{
  return Zoltan_Set_Param( LoadBalance_Ptr_,
                           const_cast<char*>(param.c_str()),
                           const_cast<char*>(value.c_str()) );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Balance(        int * changes,
					int * num_gid_entries,
					int * num_lid_entries,
					int * num_import,
					ZOLTAN_ID_PTR * import_global_ids,
					ZOLTAN_ID_PTR * import_local_ids,
					int ** import_procs,
					int * num_export,
					ZOLTAN_ID_PTR * export_global_ids,
					ZOLTAN_ID_PTR * export_local_ids,
					int ** export_procs )
{
  Zoltan::QueryContainer::setQueryID( ObjectID );
  Zoltan::MigrationContainer::setMigrationID( ObjectID );

  return Zoltan_LB_Balance( LoadBalance_Ptr_, changes, num_gid_entries,
	num_lid_entries, num_import, import_global_ids,
	import_local_ids, import_procs, num_export, export_global_ids,
	export_local_ids, export_procs );
}

#ifdef ZOLTAN_ORDER
EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Order(		int * num_gid_entries,
					int * num_lid_entries,
                                        int num_objs,
					ZOLTAN_ID_PTR global_ids,
					ZOLTAN_ID_PTR local_ids,
                                        int * rank,
                                        int * iperm )
{
  Zoltan::QueryContainer::setQueryID( ObjectID );

  return Zoltan_Order( LoadBalance_Ptr_, num_gid_entries,
	num_lid_entries, num_objs, global_ids, local_ids,
        rank, iperm, NULL );
}
#endif

EPETRAEXT_DEPRECATED
void Zoltan::LoadBalance::Evaluate     ( int print_stats,
					int * num_objects,
					float * object_weights,
                                        int * num_cuts,
					float * cut_weights,
					int * num_boundary_objects,
					int * num_adj_procs )
{
  Zoltan::QueryContainer::setQueryID( ObjectID );

//  HKT 6/24/2009 Commented out call to Zoltan because interface under construction.
//  Zoltan_LB_Eval( LoadBalance_Ptr_, print_stats, num_objects, object_weights,
//	num_cuts, cut_weights, num_boundary_objects, num_adj_procs );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Free_Data     ( ZOLTAN_ID_PTR * import_global_ids,
					ZOLTAN_ID_PTR * import_local_ids,
					int ** import_procs,
					ZOLTAN_ID_PTR * export_global_ids,
					ZOLTAN_ID_PTR * export_local_ids,
					int ** export_procs )
{
  return Zoltan_LB_Free_Data( import_global_ids, import_local_ids, import_procs,
	export_global_ids, export_local_ids, export_procs );
}

 //Decomposition Augmentation
EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Point_Assign  ( double * coords,
					int * proc )
{
  return Zoltan_LB_Point_Assign( LoadBalance_Ptr_, coords, proc );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Box_Assign    ( double xmin,
					double ymin,
					double zmin,
					double xmax,
					double ymax,
					double zmax,
					int * procs,
					int * numprocs )
{
  return Zoltan_LB_Box_Assign( LoadBalance_Ptr_, xmin, ymin, zmin, xmax, ymax, zmax,
	procs, numprocs );
}

 //Migration Functions
EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Compute_Destinations  ( int num_import,
						ZOLTAN_ID_PTR import_global_ids,
						ZOLTAN_ID_PTR import_local_ids,
						int * import_procs,
						int * num_export,
						ZOLTAN_ID_PTR * export_global_ids,
						ZOLTAN_ID_PTR * export_local_ids,
						int ** export_procs )
{
  return Zoltan_Compute_Destinations( LoadBalance_Ptr_, num_import,
	import_global_ids, import_local_ids, import_procs, num_export,
	export_global_ids, export_local_ids, export_procs );
}

EPETRAEXT_DEPRECATED
int Zoltan::LoadBalance::Help_Migrate  ( int num_import,
					ZOLTAN_ID_PTR import_global_ids,
					ZOLTAN_ID_PTR import_local_ids,
					int * import_procs,
					int num_export,
					ZOLTAN_ID_PTR export_global_ids,
					ZOLTAN_ID_PTR export_local_ids,
					int * export_procs )
{
  Zoltan::MigrationContainer::setMigrationID( ObjectID );

  return Zoltan_Help_Migrate( LoadBalance_Ptr_, num_import,
	import_global_ids, import_local_ids, import_procs,
	num_export, export_global_ids, export_local_ids, export_procs );
}

