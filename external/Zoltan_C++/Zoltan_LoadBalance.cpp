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

int Zoltan_LoadBalance::ObjectCount = 0;

Zoltan_LoadBalance::Zoltan_LoadBalance( int argc,
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

Zoltan_LoadBalance::~Zoltan_LoadBalance()
{
  Zoltan_Destroy( &LoadBalance_Ptr_ );
}

 //Load Balance Calls
int Zoltan_LoadBalance::Create( MPI_Comm communicator )
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

int Zoltan_LoadBalance::Set_CallBack_Fn       ( ZOLTAN_FN_TYPE fn_type,
						void * fn_ptr,
						void * data )
{
  return Zoltan_Set_Fn( LoadBalance_Ptr_, fn_type, fn_ptr, data );
}

#else /* ZOLTAN_OLD_CALLBACK */

//Individual callback registration
//--------------------------------

int Zoltan_LoadBalance::Set_Num_Edges_Fn      (	ZOLTAN_NUM_EDGES_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Num_Edges_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Edge_List_Fn      (	ZOLTAN_EDGE_LIST_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Edge_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Num_Geom_Fn      (	ZOLTAN_NUM_GEOM_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Num_Geom_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Geom_Fn      (	ZOLTAN_GEOM_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Geom_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Num_Obj_Fn      (	ZOLTAN_NUM_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Num_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Obj_List_Fn      (	ZOLTAN_OBJ_LIST_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Obj_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_First_Obj_Fn      (	ZOLTAN_FIRST_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_First_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Next_Obj_Fn      (	ZOLTAN_NEXT_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Next_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Num_Border_Obj_Fn (	ZOLTAN_NUM_BORDER_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Num_Border_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Border_Obj_List_Fn ( ZOLTAN_BORDER_OBJ_LIST_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Border_Obj_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_First_Border_Obj_Fn ( ZOLTAN_FIRST_BORDER_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_First_Border_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Next_Border_Obj_Fn ( ZOLTAN_NEXT_BORDER_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Next_Border_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Num_Coarse_Obj_Fn (	ZOLTAN_NUM_COARSE_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Num_Coarse_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Coarse_Obj_List_Fn ( ZOLTAN_COARSE_OBJ_LIST_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Coarse_Obj_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_First_Coarse_Obj_Fn ( ZOLTAN_FIRST_COARSE_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_First_Coarse_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Next_Coarse_Obj_Fn ( ZOLTAN_NEXT_COARSE_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Next_Coarse_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Num_Child_Fn      (	ZOLTAN_NUM_CHILD_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Num_Child_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Child_List_Fn ( ZOLTAN_CHILD_LIST_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Child_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Child_Weight_Fn ( ZOLTAN_CHILD_WEIGHT_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Child_Weight_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Pre_Migrate_Fn    (	ZOLTAN_PRE_MIGRATE_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Pre_Migrate_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Mid_Migrate_Fn    (	ZOLTAN_MID_MIGRATE_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Mid_Migrate_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Post_Migrate_Fn   (	ZOLTAN_POST_MIGRATE_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Post_Migrate_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Obj_Size_Fn       (	ZOLTAN_OBJ_SIZE_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Obj_Size_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Pack_Obj_Fn       (	ZOLTAN_PACK_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Pack_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

int Zoltan_LoadBalance::Set_Unpack_Obj_Fn     (	ZOLTAN_UNPACK_OBJ_FN * fn_ptr,
						void * data )
{
  return Zoltan_Set_Unpack_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

#endif /* ZOLTAN_OLD_CALLBACK */

int Zoltan_LoadBalance::Set_QueryObject( Zoltan_QueryObject * query_obj_ptr )
{
  Zoltan_QueryContainer::registerQueryObject( ObjectID, query_obj_ptr );

#ifdef ZOLTAN_OLD_CALLBACK

  //LB_Set_Fn all Cstyle Functions
  Set_CallBack_Fn( ZOLTAN_NUM_EDGES_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Number_Edges),
	0 );
  Set_CallBack_Fn( ZOLTAN_EDGE_LIST_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Edge_List),
	0 );
  Set_CallBack_Fn( ZOLTAN_NUM_GEOM_FN_TYPE,
    reinterpret_cast<void *> (Zoltan_QueryFunctions::Number_Geometry_Objects),
	0 );
  Set_CallBack_Fn( ZOLTAN_GEOM_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Geometry_Values),
	0 );
  Set_CallBack_Fn( ZOLTAN_NUM_OBJ_FN_TYPE, 
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Number_Objects),
	0 );
  Set_CallBack_Fn( ZOLTAN_OBJ_LIST_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Object_List),
	0 );
  Set_CallBack_Fn( ZOLTAN_FIRST_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::First_Object),
	0 );
  Set_CallBack_Fn( ZOLTAN_NEXT_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Next_Object),
	0 );
  Set_CallBack_Fn( ZOLTAN_NUM_BORDER_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Number_Border_Objects),
	0 );
  Set_CallBack_Fn( ZOLTAN_BORDER_OBJ_LIST_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Border_Object_List),
	0 );
  Set_CallBack_Fn( ZOLTAN_FIRST_BORDER_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::First_Border_Object),
	0 );
  Set_CallBack_Fn( ZOLTAN_NEXT_BORDER_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Next_Border_Object),
	0 );
  Set_CallBack_Fn( ZOLTAN_NUM_COARSE_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Number_Coarse_Objects),
	0 );
  Set_CallBack_Fn( ZOLTAN_COARSE_OBJ_LIST_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Coarse_Object_List),
	0 );
  Set_CallBack_Fn( ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::First_Coarse_Object),
	0 );
  Set_CallBack_Fn( ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Next_Coarse_Object),
	0 );
  Set_CallBack_Fn( ZOLTAN_NUM_CHILD_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Number_Children),
	0 );
  Set_CallBack_Fn( ZOLTAN_CHILD_LIST_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Child_List),
	0 );
  Set_CallBack_Fn( ZOLTAN_CHILD_WEIGHT_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_QueryFunctions::Child_Weight),
	0 );

#else /* ZOLTAN_OLD_CALLBACK */

  Set_Num_Edges_Fn        ( Zoltan_QueryFunctions::Number_Edges,
                            0 );
  Set_Edge_List_Fn        ( Zoltan_QueryFunctions::Edge_List,
                            0 );
  Set_Num_Geom_Fn         ( Zoltan_QueryFunctions::Number_Geometry_Objects,
                            0 );
  Set_Geom_Fn             ( Zoltan_QueryFunctions::Geometry_Values,
                            0 );
  Set_Num_Obj_Fn          ( Zoltan_QueryFunctions::Number_Objects,
                            0 );
  Set_Obj_List_Fn         ( Zoltan_QueryFunctions::Object_List,
                            0 );
  Set_First_Obj_Fn        ( Zoltan_QueryFunctions::First_Object,
                            0 );
  Set_Next_Obj_Fn         ( Zoltan_QueryFunctions::Next_Object,
                            0 );
  Set_Num_Border_Obj_Fn   ( Zoltan_QueryFunctions::Number_Border_Objects,
                            0 );
  Set_Border_Obj_List_Fn  ( Zoltan_QueryFunctions::Border_Object_List,
                            0 );
  Set_First_Border_Obj_Fn ( Zoltan_QueryFunctions::First_Border_Object,
                            0 );
  Set_Next_Border_Obj_Fn  ( Zoltan_QueryFunctions::Next_Border_Object,
                            0 );
  Set_Num_Coarse_Obj_Fn   ( Zoltan_QueryFunctions::Number_Coarse_Objects,
                            0 );
  Set_Coarse_Obj_List_Fn  ( Zoltan_QueryFunctions::Coarse_Object_List,
                            0 );
  Set_First_Coarse_Obj_Fn ( Zoltan_QueryFunctions::First_Coarse_Object,
                            0 );
  Set_Next_Coarse_Obj_Fn  ( Zoltan_QueryFunctions::Next_Coarse_Object,
                            0 );
  Set_Num_Child_Fn        ( Zoltan_QueryFunctions::Number_Children,
                            0 );
  Set_Child_List_Fn       ( Zoltan_QueryFunctions::Child_List,
                            0 );
  Set_Child_Weight_Fn     ( Zoltan_QueryFunctions::Child_Weight,
                            0 );

#endif /* ZOLTAN_OLD_CALLBACK */

  return ZOLTAN_OK;
}

int Zoltan_LoadBalance::Set_MigrationObject( 
			Zoltan_MigrationObject * migration_obj_ptr )
{
  Zoltan_MigrationContainer::registerMigrationObject( ObjectID,
	migration_obj_ptr );

#ifdef ZOLTAN_OLD_CALLBACK

  //LB_Set_Fn all Cstyle Functions
  Set_CallBack_Fn( ZOLTAN_OBJ_SIZE_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_MigrationFunctions::Object_Size),
        0 );
  Set_CallBack_Fn( ZOLTAN_PRE_MIGRATE_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_MigrationFunctions::Pre_Migrate),
        0 );
  Set_CallBack_Fn( ZOLTAN_MID_MIGRATE_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_MigrationFunctions::Mid_Migrate),
        0 );
  Set_CallBack_Fn( ZOLTAN_POST_MIGRATE_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_MigrationFunctions::Post_Migrate),
        0 );
  Set_CallBack_Fn( ZOLTAN_PACK_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_MigrationFunctions::Pack_Object),
        0 );
  Set_CallBack_Fn( ZOLTAN_UNPACK_OBJ_FN_TYPE,
	reinterpret_cast<void *> (Zoltan_MigrationFunctions::Unpack_Object),
        0 );

#else /* ZOLTAN_OLD_CALLBACK */

  Set_Obj_Size_Fn        ( Zoltan_MigrationFunctions::Object_Size, 0 );
  Set_Pre_Migrate_Fn     ( Zoltan_MigrationFunctions::Pre_Migrate, 0 );
  Set_Mid_Migrate_Fn     ( Zoltan_MigrationFunctions::Mid_Migrate, 0 );
  Set_Post_Migrate_Fn    ( Zoltan_MigrationFunctions::Post_Migrate, 0 );
  Set_Pack_Obj_Fn        ( Zoltan_MigrationFunctions::Pack_Object, 0 );
  Set_Unpack_Obj_Fn      ( Zoltan_MigrationFunctions::Unpack_Object, 0 );

#endif /* ZOLTAN_OLD_CALLBACK */

  return ZOLTAN_OK;
}

int Zoltan_LoadBalance::Set_Param( const std::string & param, const std::string & value )
{
  return Zoltan_Set_Param( LoadBalance_Ptr_,
                           const_cast<char*>(param.c_str()),
                           const_cast<char*>(value.c_str()) );
}

int Zoltan_LoadBalance::Balance(        int * changes,
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
  Zoltan_QueryContainer::setQueryID( ObjectID );
  Zoltan_MigrationContainer::setMigrationID( ObjectID );

  return Zoltan_LB_Balance( LoadBalance_Ptr_, changes, num_gid_entries,
	num_lid_entries, num_import, import_global_ids,
	import_local_ids, import_procs, num_export, export_global_ids,
	export_local_ids, export_procs );
}

#ifdef ZOLTAN_ORDER
int Zoltan_LoadBalance::Order(		int * num_gid_entries,
					int * num_lid_entries,
                                        int num_objs,
					ZOLTAN_ID_PTR global_ids,
					ZOLTAN_ID_PTR local_ids,
                                        int * rank,
                                        int * iperm )
{
  Zoltan_QueryContainer::setQueryID( ObjectID );

  return Zoltan_Order( LoadBalance_Ptr_, num_gid_entries,
	num_lid_entries, num_objs, global_ids, local_ids,
        rank, iperm, NULL );
}
#endif

void Zoltan_LoadBalance::Evaluate     ( int print_stats,
					int * num_objects,
					float * object_weights,
                                        int * num_cuts,
					float * cut_weights,
					int * num_boundary_objects,
					int * num_adj_procs )
{
  Zoltan_LB_Eval( LoadBalance_Ptr_, print_stats, num_objects, object_weights,
	num_cuts, cut_weights, num_boundary_objects, num_adj_procs );
}

int Zoltan_LoadBalance::Free_Data     ( ZOLTAN_ID_PTR * import_global_ids,
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
int Zoltan_LoadBalance::Point_Assign  ( double * coords,
					int * proc )
{
  return Zoltan_LB_Point_Assign( LoadBalance_Ptr_, coords, proc );
}

int Zoltan_LoadBalance::Box_Assign    ( double xmin,
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
int Zoltan_LoadBalance::Compute_Destinations  ( int num_import,
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

int Zoltan_LoadBalance::Help_Migrate  ( int num_import,
					ZOLTAN_ID_PTR import_global_ids,
					ZOLTAN_ID_PTR import_local_ids,
					int * import_procs,
					int num_export,
					ZOLTAN_ID_PTR export_global_ids,
					ZOLTAN_ID_PTR export_local_ids,
					int * export_procs )
{
  Zoltan_MigrationContainer::setMigrationID( ObjectID );

  return Zoltan_Help_Migrate( LoadBalance_Ptr_, num_import,
	import_global_ids, import_local_ids, import_procs,
	num_export, export_global_ids, export_local_ids, export_procs );
}

