//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_LoadBalance.h$
//
// Purpose        : C++ wrapper object for Zoltan Load Balance routines.
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

#ifndef ZOLTAN_LOADBALANCE_H_
#define ZOLTAN_LOADBALANCE_H_

#include <zoltan.h>

#include <string>

class Zoltan_QueryObject;
class Zoltan_MigrationObject;

class Zoltan_LoadBalance
{

public:

 Zoltan_LoadBalance   ( int argc = 0,
			char ** argv = 0,
			float * ver = 0);

 ~Zoltan_LoadBalance();

 //Load Balance Calls
 int Create( MPI_Comm communicator );

#ifdef ZOLTAN_OLD_CALLBACK

 int Set_CallBack_Fn  ( ZOLTAN_FN_TYPE fn_type,
			void * fn_ptr,
			void * data = 0 );

#else /* ZOLTAN_OLD_CALLBACK */

//Individual Callback Function Registration
//-----------------------------------------

int Set_Num_Edges_Fn        ( ZOLTAN_NUM_EDGES_FN * fn_ptr,
                              void * data = 0 );
int Set_Edge_List_Fn        ( ZOLTAN_EDGE_LIST_FN * fn_ptr,
                              void * data = 0 );
int Set_Num_Geom_Fn         ( ZOLTAN_NUM_GEOM_FN * fn_ptr,
                              void * data = 0 );
int Set_Geom_Fn             ( ZOLTAN_GEOM_FN * fn_ptr,
                              void * data = 0 );
int Set_Num_Obj_Fn          ( ZOLTAN_NUM_OBJ_FN * fn_ptr,
                              void * data = 0 );
int Set_Obj_List_Fn         ( ZOLTAN_OBJ_LIST_FN * fn_ptr,
                              void * data = 0 );
int Set_First_Obj_Fn        ( ZOLTAN_FIRST_OBJ_FN * fn_ptr,
                              void * data = 0 );
int Set_Next_Obj_Fn         ( ZOLTAN_NEXT_OBJ_FN * fn_ptr,
                              void * data = 0 );
int Set_Num_Border_Obj_Fn   ( ZOLTAN_NUM_BORDER_OBJ_FN * fn_ptr,
                              void * data = 0 );
int Set_Border_Obj_List_Fn  ( ZOLTAN_BORDER_OBJ_LIST_FN * fn_ptr,
                              void * data = 0 );
int Set_First_Border_Obj_Fn ( ZOLTAN_FIRST_BORDER_OBJ_FN * fn_ptr,
                              void * data = 0 );
int Set_Next_Border_Obj_Fn  ( ZOLTAN_NEXT_BORDER_OBJ_FN * fn_ptr,
                              void * data = 0 );
int Set_Num_Coarse_Obj_Fn   ( ZOLTAN_NUM_COARSE_OBJ_FN * fn_ptr,
                              void * data = 0 );
int Set_Coarse_Obj_List_Fn  ( ZOLTAN_COARSE_OBJ_LIST_FN * fn_ptr,
                              void * data = 0 );
int Set_First_Coarse_Obj_Fn ( ZOLTAN_FIRST_COARSE_OBJ_FN * fn_ptr,
                              void * data = 0 );
int Set_Next_Coarse_Obj_Fn  ( ZOLTAN_NEXT_COARSE_OBJ_FN * fn_ptr,
                              void * data = 0 );
int Set_Num_Child_Fn        ( ZOLTAN_NUM_CHILD_FN * fn_ptr,
                              void * data = 0 );
int Set_Child_List_Fn       ( ZOLTAN_CHILD_LIST_FN * fn_ptr,
                              void * data = 0 );
int Set_Child_Weight_Fn     ( ZOLTAN_CHILD_WEIGHT_FN * fn_ptr,
                              void * data = 0 );
int Set_Pre_Migrate_Fn      ( ZOLTAN_PRE_MIGRATE_FN * fn_ptr,
                              void * data = 0 );
int Set_Mid_Migrate_Fn      ( ZOLTAN_MID_MIGRATE_FN * fn_ptr,
                              void * data = 0 );
int Set_Post_Migrate_Fn     ( ZOLTAN_POST_MIGRATE_FN * fn_ptr,
                              void * data = 0 );
int Set_Obj_Size_Fn         ( ZOLTAN_OBJ_SIZE_FN * fn_ptr,
                              void * data = 0 );
int Set_Pack_Obj_Fn         ( ZOLTAN_PACK_OBJ_FN * fn_ptr,
                              void * data = 0 );
int Set_Unpack_Obj_Fn       ( ZOLTAN_UNPACK_OBJ_FN * fn_ptr,
                              void * data = 0 );

#endif /* ZOLTAN_OLD_CALLBACK */

 int Set_QueryObject( Zoltan_QueryObject * query_obj_ptr );

 int Set_MigrationObject( Zoltan_MigrationObject * migration_obj_ptr );

 int Set_Param( const std::string & param, const std::string & value );

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

 void Evaluate( int print_stats,
		int * num_objects,
		float * object_weights,
    		int * num_cut,
    		float * cut_weights,
		int * num_boundary_objects,
		int * num_adj_procs );

 int Free_Data( ZOLTAN_ID_PTR * import_global_ids,
		ZOLTAN_ID_PTR * import_local_ids,
		int ** import_procs,
		ZOLTAN_ID_PTR * export_global_ids,
		ZOLTAN_ID_PTR * export_local_ids,
		int ** export_procs );

 //Decomposition Augmentation
 int Point_Assign     ( double * coords,
			int * proc );

 int Box_Assign       ( double xmin,
			double ymin,
			double zmin,
			double xmax,
			double ymax,
			double zmax,
			int * procs,
			int * numprocs );

 //Migration Functions
 int Compute_Destinations     ( int num_import,
				ZOLTAN_ID_PTR import_global_ids,
				ZOLTAN_ID_PTR import_local_ids,
				int * import_procs,
				int * num_export,
				ZOLTAN_ID_PTR * export_global_ids,
				ZOLTAN_ID_PTR * export_local_ids,
				int ** export_procs );

 int Help_Migrate     ( int num_import,
			ZOLTAN_ID_PTR import_global_ids,
			ZOLTAN_ID_PTR import_local_ids,
			int * import_procs,
			int num_export,
			ZOLTAN_ID_PTR export_global_ids,
			ZOLTAN_ID_PTR export_local_ids,
			int * export_procs );

 Zoltan_Struct * Return_Zoltan_Struct() { return LoadBalance_Ptr_; }

private:

 static int ObjectCount;

 int ObjectID;

 Zoltan_Struct * LoadBalance_Ptr_; 

 Zoltan_QueryObject * QueryObject_Ptr_;
 Zoltan_MigrationObject * MigrationObject_Ptr_;

};

#endif
