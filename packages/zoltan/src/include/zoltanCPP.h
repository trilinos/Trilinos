// @HEADER
// ************************************************************************
// 
//            Zoltan_CPP: An Object-Oriented Interface To Zoltan
//                    Copyright (2001) Sandia Corporation
// 
// Questions? Contact Robert J. Hoekstra (rjhoeks@sandia.gov)
// 
// ************************************************************************
// @HEADER

#ifndef ZOLTANCPP_H_
#define ZOLTANCPP_H_

#include <zoltan.h>

#include <string>

namespace Zoltan {

//! Zoltan::LoadBalance: A class for interfacing the Load Balancing functions of the Zoltan library in a C++/Object Oriented environment.

/*! The Zoltan::LoadBalance class is a wrapper for the C functions at the top level
 interface of Zoltan.  The primary differences include the removal the Zoltan_Struct
 parameter since this object is now stored in the class
*/

class LoadBalance
{

public:

 //@{ \name Constructors/Destructors.

 //! Constructor
 /*! This constructor replaces the Zoltan_Initialize call. Params are the same.
 */
 LoadBalance   ( int argc = 0,
                 char ** argv = 0,
                 float * ver = 0 )
 : LoadBalance_Ptr_( 0 )
{
  int tmpReturn = Zoltan_Initialize( argc, argv, ver );
}

 //! Destructor
 ~LoadBalance()
{
  Zoltan_Destroy( &LoadBalance_Ptr_ );
}

 //@}

 //@{ \name General Load Balance Functionality
 
 //! Replaces Zoltan_Create
 int Create( MPI_Comm communicator )
{
  LoadBalance_Ptr_ = Zoltan_Create( communicator );
                                                                                                       
  if( !LoadBalance_Ptr_ )
    return ZOLTAN_FATAL;
  else
    return ZOLTAN_OK;
}

 //! Replaces Zoltan_Set_Param
 int Set_Param( const std::string & param, const std::string & value )
{
  return Zoltan_Set_Param( LoadBalance_Ptr_,
                           const_cast<char*>(param.c_str()),
                           const_cast<char*>(value.c_str()) );
}

 //! Replaces Zoltan_Set_Param_Vec
 int Set_Param_Vec( const std::string & param, const std::string & value, int index )
{
  return Zoltan_Set_Param_Vec( LoadBalance_Ptr_,
                               const_cast<char*>(param.c_str()),
                               const_cast<char*>(value.c_str()),
                               index );
}

 //! Replaces Zoltan_LB_Partition
 int Partition ( int * changes,
                 int * num_gid_entries,
                 int * num_lid_entries,
                 int * num_import,
                 ZOLTAN_ID_PTR * import_global_ids,
                 ZOLTAN_ID_PTR * import_local_ids,
                 int ** import_procs,
                 int ** import_to_part,
                 int * num_export,
                 ZOLTAN_ID_PTR * export_global_ids,
                 ZOLTAN_ID_PTR * export_local_ids,
                 int ** export_procs,
                 int ** export_to_part )
{
  return Zoltan_LB_Partiton( LoadBalance_Ptr_, changes,
                             num_gid_entries, num_lid_entries,
                             num_import, import_global_ids, import_local_ids, import_procs, import_to_part,
                             num_export, export_global_ids, export_local_ids, export_procs, export_to_part );
}

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
		int ** export_procs )
{
  return Zoltan_LB_Balance( LoadBalance_Ptr_, changes,
                            num_gid_entries, num_lid_entries,
                            num_import, import_global_ids, import_local_ids, import_procs,
                            num_export, export_global_ids, export_local_ids, export_procs );
}

 int Order    ( int * num_gid_entries,
                int * num_lid_entries,
                int num_objs,
		ZOLTAN_ID_PTR global_ids,
		ZOLTAN_ID_PTR local_ids,
                int * rank,
                int * iperm )
{
  return Zoltan_Order( LoadBalance_Ptr_,
                       num_gid_entries, num_lid_entries,
                       num_objs, global_ids, local_ids,
                       rank, iperm, NULL );
}

 //! Replaces Zoltan_LB_Eval
 void Evaluate( int print_stats,
		int * num_objects,
		float * object_weights,
    		int * num_cut,
    		float * cut_weights,
		int * num_boundary_objects,
		int * num_adj_procs )
{
  Zoltan_LB_Eval( LoadBalance_Ptr_, print_stats,
                  num_objects, object_weights,
                  num_cuts, cut_weights,
                  num_boundary_objects, num_adj_procs );
}


 //! Replaces Zoltan_LB_Free_Part
 int Free_Data( ZOLTAN_ID_PTR * global_ids,
		ZOLTAN_ID_PTR * local_ids,
                int ** procs,
                int ** to_part )
{
  return Zoltan_LB_Free_Part( global_ids, local_ids, procs, to_part );
}

 //! Replaces Zoltan_LB_Free_Data
 int Free_Data( ZOLTAN_ID_PTR * import_global_ids,
		ZOLTAN_ID_PTR * import_local_ids,
		int ** import_procs,
		ZOLTAN_ID_PTR * export_global_ids,
		ZOLTAN_ID_PTR * export_local_ids,
		int ** export_procs )
{
  return Zoltan_LB_Free_Data( import_global_ids, import_local_ids, import_procs,
                              export_global_ids, export_local_ids, export_procs );
}

 //@}

 //@{ \name Support for direct access to Zoltan callback functionality

 //! Old style callback support
 int Set_CallBack_Fn  ( ZOLTAN_FN_TYPE fn_type,
			void (*fn_ptr)(),
			void * data = 0 )
{
  return Zoltan_Set_Fn( LoadBalance_Ptr_, fn_type, fn_ptr, data );
}

 //! Individual callback support

 //!
 int Set_Partition_Multi_Fn  ( ZOLTAN_PARTITION_MULTI_FN * fn_ptr,
                               void * data = 0 )
{
  return Zoltan_Set_Partition_Multi_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 //!
 int Set_Partition_Fn        ( ZOLTAN_PARTITION_FN * fn_ptr,
                               void * data = 0 )
{
  return Zoltan_Set_Partition_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 //!
 int Set_Num_Edges_Multi_Fn ( ZOLTAN_NUM_EDGES_MULTI_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Num_Edges_Multi_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 //!
 int Set_Num_Edges_Fn       ( ZOLTAN_NUM_EDGES_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Num_Edges_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Edge_List_Multi_Fn ( ZOLTAN_EDGE_LIST_MULTI_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Edge_List_Multi_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Edge_List_Fn       ( ZOLTAN_EDGE_LIST_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Edge_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Num_Geom_Fn        ( ZOLTAN_NUM_GEOM_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Num_Geom_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Geom_Multi_Fn      ( ZOLTAN_GEOM_MULTI_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Geom_Multi_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Geom_Fn            ( ZOLTAN_GEOM_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Geom_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Num_Obj_Fn         ( ZOLTAN_NUM_OBJ_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Num_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Obj_List_Fn        ( ZOLTAN_OBJ_LIST_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Obj_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_First_Obj_Fn       ( ZOLTAN_FIRST_OBJ_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_First_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Next_Obj_Fn        ( ZOLTAN_NEXT_OBJ_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Next_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Num_Border_Obj_Fn  ( ZOLTAN_NUM_BORDER_OBJ_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Num_Border_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Border_Obj_List_Fn ( ZOLTAN_BORDER_OBJ_LIST_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Border_Obj_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_First_Border_Obj_Fn( ZOLTAN_FIRST_BORDER_OBJ_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_First_Border_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Next_Border_Obj_Fn ( ZOLTAN_NEXT_BORDER_OBJ_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Next_Border_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Num_Coarse_Obj_Fn  ( ZOLTAN_NUM_COARSE_OBJ_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Num_Coarse_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Coarse_Obj_List_Fn ( ZOLTAN_COARSE_OBJ_LIST_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Coarse_Obj_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_First_Coarse_Obj_Fn( ZOLTAN_FIRST_COARSE_OBJ_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_First_Coarse_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Next_Coarse_Obj_Fn ( ZOLTAN_NEXT_COARSE_OBJ_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Next_Coarse_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Num_Child_Fn       ( ZOLTAN_NUM_CHILD_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Num_Child_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Child_List_Fn      ( ZOLTAN_CHILD_LIST_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Child_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Child_Weight_Fn    ( ZOLTAN_CHILD_WEIGHT_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Child_Weight_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Num_HG_Edges_Fn    ( ZOLTAN_NUM_HG_EDGES_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Num_HG_Edges_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Num_HG_Pins_Fn     ( ZOLTAN_NUM_HG_PINS_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Num_HG_Pins_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_HG_Edges_List_Fn   ( ZOLTAN_HG_EDGES_LIST_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_HG_Edges_List_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Pre_Migrate_PP_Fn  ( ZOLTAN_PRE_MIGRATE_PP_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Pre_Migrate_PP_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Mid_Migrate_PP_Fn  ( ZOLTAN_MID_MIGRATE_PP_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Mid_Migrate_PP_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Post_Migrate_PP_Fn ( ZOLTAN_POST_MIGRATE_PP_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Post_Migrate_PP_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Pre_Migrate_Fn     ( ZOLTAN_PRE_MIGRATE_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Pre_Migrate_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Mid_Migrate_Fn     ( ZOLTAN_MID_MIGRATE_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Mid_Migrate_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Post_Migrate_Fn    ( ZOLTAN_POST_MIGRATE_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Post_Migrate_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Obj_Size_Multi_Fn  ( ZOLTAN_OBJ_SIZE_MULTI_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Obj_Size_Multi_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Obj_Size_Fn        ( ZOLTAN_OBJ_SIZE_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Obj_Size_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Pack_Obj_Multi_Fn  ( ZOLTAN_PACK_OBJ_MULTI_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Pack_Obj_Multi_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Pack_Obj_Fn        ( ZOLTAN_PACK_OBJ_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Pack_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Unpack_Obj_Multi_Fn( ZOLTAN_UNPACK_OBJ_MULTI_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Unpack_Obj_Multi_Fn( LoadBalance_Ptr_, fn_ptr, data );
}

 ///
 int Set_Unpack_Obj_Fn      ( ZOLTAN_UNPACK_OBJ_FN * fn_ptr,
                              void * data = 0 )
{
  return Zoltan_Set_Unpack_Obj_Fn( LoadBalance_Ptr_, fn_ptr, data );
}


 //@}

 //@{ \name Decomposition Augmentation

 //! Replaces Zoltan_LB_Point_PP_Assign
 int Point_PP_Assign     ( double * coords,
                           int * proc,
                           int * part )
{
  return Zoltan_LB_Point_PP_Assign( LoadBalance_Ptr_, coords, proc, part );
}

 //! Replaces Zoltan_LB_Point_Assign
 int Point_Assign     ( double * coords,
			int * proc )
{
  return Zoltan_LB_Point_Assign( LoadBalance_Ptr_, coords, proc );
}


 //! Replaces Zoltan_LB_Box_PP_Assign
 int Box_Assign       ( double xmin,
			double ymin,
			double zmin,
			double xmax,
			double ymax,
			double zmax,
			int * procs,
			int * numprocs,
                        int * parts,
                        int * numparts )
{
  return Zoltan_LB_Box_PP_Assign( LoadBalance_Ptr_,
                                  xmin, ymin, zmin,
                                  xmax, ymax, zmax,
                                  procs, numprocs,
                                  parts, numparts );
}

 //! Replaces Zoltan_LB_Box_Assign
 int Box_Assign       ( double xmin,
			double ymin,
			double zmin,
			double xmax,
			double ymax,
			double zmax,
			int * procs,
			int * numprocs )
{
  return Zoltan_LB_Box_Assign( LoadBalance_Ptr_,
                               xmin, ymin, zmin,
                               xmax, ymax, zmax,
                               procs, numprocs );
}

 //@}
 
 //@{ \name Migration Functionality

 //! Replaces Zoltan_Invert_Lists
 int Invert_Lists             ( int num_import,
				ZOLTAN_ID_PTR import_global_ids,
				ZOLTAN_ID_PTR import_local_ids,
				int * import_procs,
				int * import_to_part,
				int * num_export,
				ZOLTAN_ID_PTR * export_global_ids,
				ZOLTAN_ID_PTR * export_local_ids,
				int ** export_procs,
				int ** export_to_part )
{
  return Zoltan_Invert_Lists( LoadBalance_Ptr_,
                              num_import, import_global_ids, import_local_ids, import_procs, import_to_part,
                              num_export, export_global_ids, export_local_ids, export_procs, export_to_part );
}
 //! Replaces Zoltan_Compute_Destinations
 int Compute_Destinations     ( int num_import,
				ZOLTAN_ID_PTR import_global_ids,
				ZOLTAN_ID_PTR import_local_ids,
				int * import_procs,
				int * num_export,
				ZOLTAN_ID_PTR * export_global_ids,
				ZOLTAN_ID_PTR * export_local_ids,
				int ** export_procs )
{
  return Zoltan_Compute_Destinations( LoadBalance_Ptr_,
                                      num_import, import_global_ids, import_local_ids, import_procs,
                                      num_export, export_global_ids, export_local_ids, export_procs );
}

 //! Replaces Zoltan_Help_Migrate
 int Migrate          ( int num_import,
			ZOLTAN_ID_PTR import_global_ids,
			ZOLTAN_ID_PTR import_local_ids,
			int * import_procs,
			int * import_to_part,
			int num_export,
			ZOLTAN_ID_PTR export_global_ids,
			ZOLTAN_ID_PTR export_local_ids,
			int * export_procs,
			int * export_to_part )
{
  return Zoltan_Migrate( LoadBalance_Ptr_,
                         num_import, import_global_ids, import_local_ids, import_procs, import_to_part,
                         num_export, export_global_ids, export_local_ids, export_procs, export_to_part );
}

 //! Replaces Zoltan_Help_Migrate
 int Help_Migrate     ( int num_import,
			ZOLTAN_ID_PTR import_global_ids,
			ZOLTAN_ID_PTR import_local_ids,
			int * import_procs,
			int num_export,
			ZOLTAN_ID_PTR export_global_ids,
			ZOLTAN_ID_PTR export_local_ids,
			int * export_procs )
{
  return Zoltan_Help_Migrate( LoadBalance_Ptr_,
                              num_import, import_global_ids, import_local_ids, import_procs,
                              num_export, export_global_ids, export_local_ids, export_procs );
}

 //@}

 //@{ \name Partitioning Help Functionality

 //! Replaces Zoltan_LB_Set_Part_Sizes
 int Set_Part_Sizes( int global_num,
                     int len,
                     int * part_ids,
                     int * wgt_idx,
                     float * part_sizes )
{
  return Zoltan_LB_Set_Part_Sizes( LoadBalance_Ptr_,
                                   global_num, len,
                                   part_ids, wgt_idx, part_sizes );
}

 //@}

 //@{ \name File Generation Functionality

 //! Replaces Zoltan_Generate_Files
 int Generate_Files( char * fname,
                     int base_index,
                     int gen_geom,
                     int gen_graph,
                     int gen_hg )
{
  return Zoltan_Generate_Files( LoadBalance_Ptr_,
                                fname, base_index,
                                gen_geom, gen_graph, gen_hg );
}

 //@}

 //@{ \name Extra
 
 //! Direct access to underlying Zoltan_Struct Object
 Zoltan_Struct * Return_Zoltan_Struct()
 { return LoadBalance_Ptr_; }

 //@}

private:

 Zoltan_Struct * LoadBalance_Ptr_; 

};

} //namespace Zoltan

#endif
