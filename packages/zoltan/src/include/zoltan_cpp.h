// *****************************************************************************
// * Zoltan Library for Parallel Applications                                  *
// * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
// * This software is distributed under the GNU Lesser General Public License. *
// * For more info, see the README file in the top-level Zoltan directory.     *
// *****************************************************************************
// *****************************************************************************
// * CVS File Information :
// *    $RCSfile$
// *    $Author$
// *    $Date$
// *    $Revision$
// *****************************************************************************

// ************************************************************************
// C++ class representing a Zoltan_Struct object.
//
// Assumption: Zoltan_Initialize has been called prior to creating
// a Zoltan_Object.
// ************************************************************************

#ifndef ZOLTAN_CPP_H_
#define ZOLTAN_CPP_H_

#include "zoltan.h"
#include "zoltan_comm_cpp.h"
#include "zoltan_dd_cpp.h"
#include <string>

class Zoltan_Object {

public:

  Zoltan_Object (MPI_Comm communicator = MPI_COMM_WORLD) 
  {
  this->ZZ_Ptr = Zoltan_Create(communicator);

  // int fail = (this->ZZ_Ptr == NULL);  should catch this exception
  }

  ~Zoltan_Object()
  {
    // Warning: Zoltan_Destroy calls MPI.  If you are not using the
    // C++ bindings for MPI, along with these C++ bindings for Zoltan,
    // you may be calling MPI_Finalize() before this destructor gets
    // called, unless you allocate and destroy Zoltan_Objects
    // explicitly:
    //
    //   MPI_Init(...);
    //   Zoltan_Initialize(...);
    //   Zoltan_Object *zz = new Zoltan_Object();
    //    ... more code ...
    //   delete zz;
    //   MPI_Finalize();
    //
    // ZoltanObject's created on the stack will be deleted at exit,
    //  after MPI_Finalize().

    Zoltan_Destroy( &ZZ_Ptr );
  }

  int Set_Param( const std::string & param, const std::string & value )
  {
    return Zoltan_Set_Param( ZZ_Ptr,
                             const_cast<char*>(param.c_str()),
                             const_cast<char*>(value.c_str()) );
  }

  //! Replaces Zoltan_Set_Param_Vec
  int Set_Param_Vec( const std::string & param, const std::string & value, 
                     int index )
  {
    return Zoltan_Set_Param_Vec( ZZ_Ptr,
                                 const_cast<char*>(param.c_str()),
                                 const_cast<char*>(value.c_str()),
                                 index );
  }

  int LB_Partition ( int * changes,
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
    return Zoltan_LB_Partition( ZZ_Ptr, changes,
                                num_gid_entries, num_lid_entries,
                                num_import, import_global_ids, import_local_ids,
                                import_procs, import_to_part,
                                num_export, export_global_ids, export_local_ids,
                                export_procs, export_to_part );
  }

  int LB_Balance  ( int * changes,
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
    return Zoltan_LB_Balance( ZZ_Ptr, changes,
                              num_gid_entries, num_lid_entries,
                              num_import, import_global_ids, import_local_ids, 
                              import_procs,
                              num_export, export_global_ids, export_local_ids, 
                              export_procs );
  }

  int LB_Set_Part_Sizes( int global_num,
                         int len,
                         int * part_ids,
                         int * wgt_idx,
                         float * part_sizes )
  {
    return Zoltan_LB_Set_Part_Sizes( ZZ_Ptr,
                                     global_num, len,
                                     part_ids, wgt_idx, part_sizes );
  }

  int Order    ( int * num_gid_entries,
                 int * num_lid_entries,
                 int num_objs,
                 ZOLTAN_ID_PTR global_ids,
                 ZOLTAN_ID_PTR local_ids,
                 int * rank,
                 int * iperm )
  {
    //  Note:  Zoltan_Order_Struct set to NULL.
    return Zoltan_Order( ZZ_Ptr,
                         num_gid_entries, num_lid_entries,
                         num_objs, global_ids, local_ids,
                         rank, iperm, NULL );
  }

  int LB_Eval( int print_stats,
                int * num_objects,
                float * object_weights,
                int * num_cuts,
                float * cut_weights,
                int * num_boundary_objects,
                int * num_adj_procs )
  {
    return Zoltan_LB_Eval( ZZ_Ptr, print_stats,
                    num_objects, object_weights,
                    num_cuts, cut_weights,
                    num_boundary_objects, num_adj_procs );
  }

  int RCB_Box( int part,
                int * ndim, 
                double * xmin, 
                double * ymin, 
                double * zmin, 
                double * xmax, 
                double * ymax, 
                double * zmax )
  {
    return Zoltan_RCB_Box( ZZ_Ptr, part, ndim, xmin, ymin, zmin, xmax, ymax, zmax);
  }

  int LB_Free_Part( ZOLTAN_ID_PTR * global_ids,
                    ZOLTAN_ID_PTR * local_ids,
                    int ** procs,
                    int ** to_part )
  {
    return Zoltan_LB_Free_Part( global_ids, local_ids, procs, to_part );
  }

  int LB_Free_Data( ZOLTAN_ID_PTR * import_global_ids,
                    ZOLTAN_ID_PTR * import_local_ids,
                    int ** import_procs,
                    ZOLTAN_ID_PTR * export_global_ids,
                    ZOLTAN_ID_PTR * export_local_ids,
                    int ** export_procs )
  {
    return Zoltan_LB_Free_Data( import_global_ids, import_local_ids, 
                                import_procs,
                                export_global_ids, export_local_ids, 
                                export_procs );
  }

  int Set_Fn  ( ZOLTAN_FN_TYPE fn_type,
                void (*fn_ptr)(),
                void * data = 0 )
  {
    return Zoltan_Set_Fn( ZZ_Ptr, fn_type, fn_ptr, data );
  }

  // Individual callback support

  ///--------------------------
  int Set_Partition_Multi_Fn  ( ZOLTAN_PARTITION_MULTI_FN * fn_ptr,
                                void * data = 0 )
  {
    return Zoltan_Set_Partition_Multi_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Partition_Fn        ( ZOLTAN_PARTITION_FN * fn_ptr,
                                void * data = 0 )
  {
    return Zoltan_Set_Partition_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Num_Edges_Multi_Fn ( ZOLTAN_NUM_EDGES_MULTI_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Num_Edges_Multi_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Num_Edges_Fn       ( ZOLTAN_NUM_EDGES_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Num_Edges_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Edge_List_Multi_Fn ( ZOLTAN_EDGE_LIST_MULTI_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Edge_List_Multi_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Edge_List_Fn       ( ZOLTAN_EDGE_LIST_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Edge_List_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Num_Geom_Fn        ( ZOLTAN_NUM_GEOM_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Num_Geom_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Geom_Multi_Fn      ( ZOLTAN_GEOM_MULTI_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Geom_Multi_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Geom_Fn            ( ZOLTAN_GEOM_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Geom_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Num_Obj_Fn         ( ZOLTAN_NUM_OBJ_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Num_Obj_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Obj_List_Fn        ( ZOLTAN_OBJ_LIST_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Obj_List_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_First_Obj_Fn       ( ZOLTAN_FIRST_OBJ_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_First_Obj_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Next_Obj_Fn        ( ZOLTAN_NEXT_OBJ_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Next_Obj_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Num_Border_Obj_Fn  ( ZOLTAN_NUM_BORDER_OBJ_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Num_Border_Obj_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Border_Obj_List_Fn ( ZOLTAN_BORDER_OBJ_LIST_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Border_Obj_List_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_First_Border_Obj_Fn( ZOLTAN_FIRST_BORDER_OBJ_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_First_Border_Obj_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Next_Border_Obj_Fn ( ZOLTAN_NEXT_BORDER_OBJ_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Next_Border_Obj_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Num_Coarse_Obj_Fn  ( ZOLTAN_NUM_COARSE_OBJ_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Num_Coarse_Obj_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Coarse_Obj_List_Fn ( ZOLTAN_COARSE_OBJ_LIST_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Coarse_Obj_List_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_First_Coarse_Obj_Fn( ZOLTAN_FIRST_COARSE_OBJ_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_First_Coarse_Obj_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Next_Coarse_Obj_Fn ( ZOLTAN_NEXT_COARSE_OBJ_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Next_Coarse_Obj_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Num_Child_Fn       ( ZOLTAN_NUM_CHILD_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Num_Child_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Child_List_Fn      ( ZOLTAN_CHILD_LIST_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Child_List_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Child_Weight_Fn    ( ZOLTAN_CHILD_WEIGHT_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Child_Weight_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Num_HG_Edges_Fn    ( ZOLTAN_NUM_HG_EDGES_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Num_HG_Edges_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Num_HG_Pins_Fn     ( ZOLTAN_NUM_HG_PINS_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Num_HG_Pins_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_HG_Edge_List_Fn    ( ZOLTAN_HG_EDGE_LIST_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_HG_Edge_List_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Pre_Migrate_PP_Fn  ( ZOLTAN_PRE_MIGRATE_PP_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Pre_Migrate_PP_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Mid_Migrate_PP_Fn  ( ZOLTAN_MID_MIGRATE_PP_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Mid_Migrate_PP_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Post_Migrate_PP_Fn ( ZOLTAN_POST_MIGRATE_PP_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Post_Migrate_PP_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Pre_Migrate_Fn     ( ZOLTAN_PRE_MIGRATE_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Pre_Migrate_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Mid_Migrate_Fn     ( ZOLTAN_MID_MIGRATE_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Mid_Migrate_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Post_Migrate_Fn    ( ZOLTAN_POST_MIGRATE_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Post_Migrate_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Obj_Size_Multi_Fn  ( ZOLTAN_OBJ_SIZE_MULTI_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Obj_Size_Multi_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Obj_Size_Fn        ( ZOLTAN_OBJ_SIZE_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Obj_Size_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Pack_Obj_Multi_Fn  ( ZOLTAN_PACK_OBJ_MULTI_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Pack_Obj_Multi_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Pack_Obj_Fn        ( ZOLTAN_PACK_OBJ_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Pack_Obj_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Unpack_Obj_Multi_Fn( ZOLTAN_UNPACK_OBJ_MULTI_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Unpack_Obj_Multi_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Unpack_Obj_Fn      ( ZOLTAN_UNPACK_OBJ_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Unpack_Obj_Fn( ZZ_Ptr, fn_ptr, data );
  }

  int LB_Point_PP_Assign ( double * coords,
                           int * proc,
                           int * part )
  {
    return Zoltan_LB_Point_PP_Assign( ZZ_Ptr, coords, proc, part );
  }

  int LB_Point_Assign ( double * coords,
                        int * proc )
  {
    return Zoltan_LB_Point_Assign( ZZ_Ptr, coords, proc );
  }

  int LB_Box_PP_Assign ( double xmin,
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
    return Zoltan_LB_Box_PP_Assign( ZZ_Ptr,
                                    xmin, ymin, zmin,
                                    xmax, ymax, zmax,
                                    procs, numprocs,
                                    parts, numparts );
  }

  int LB_Box_Assign ( double xmin,
                      double ymin,
                      double zmin,
                      double xmax,
                      double ymax,
                      double zmax,
                      int * procs,
                      int * numprocs )
  {
    return Zoltan_LB_Box_Assign( ZZ_Ptr,
                                 xmin, ymin, zmin,
                                 xmax, ymax, zmax,
                                 procs, numprocs );
  }

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
    return Zoltan_Invert_Lists( ZZ_Ptr,
                                num_import, import_global_ids, import_local_ids,
                                import_procs, import_to_part,
                                num_export, export_global_ids, export_local_ids,
                                export_procs, export_to_part );
  }

  int Compute_Destinations     ( int num_import,
                                 ZOLTAN_ID_PTR import_global_ids,
                                 ZOLTAN_ID_PTR import_local_ids,
                                 int * import_procs,
                                 int * num_export,
                                 ZOLTAN_ID_PTR * export_global_ids,
                                 ZOLTAN_ID_PTR * export_local_ids,
                                 int ** export_procs )
  {
    return Zoltan_Compute_Destinations( ZZ_Ptr,
                                        num_import, import_global_ids, 
                                        import_local_ids, import_procs,
                                        num_export, export_global_ids, 
                                        export_local_ids, export_procs );
  }

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
    return Zoltan_Migrate( ZZ_Ptr,
                           num_import, import_global_ids, import_local_ids, 
                           import_procs, import_to_part,
                           num_export, export_global_ids, export_local_ids, 
                           export_procs, export_to_part );
  }

  int Help_Migrate     ( int num_import,
                         ZOLTAN_ID_PTR import_global_ids,
                         ZOLTAN_ID_PTR import_local_ids,
                         int * import_procs,
                         int num_export,
                         ZOLTAN_ID_PTR export_global_ids,
                         ZOLTAN_ID_PTR export_local_ids,
                         int * export_procs )
  {
    return Zoltan_Help_Migrate( ZZ_Ptr,
                                num_import, import_global_ids, import_local_ids,
                                import_procs,
                                num_export, export_global_ids, export_local_ids,
                                export_procs );
  }

  int Generate_Files( char * fname,
                      int base_index,
                      int gen_geom,
                      int gen_graph,
                      int gen_hg )
  {
    return Zoltan_Generate_Files( ZZ_Ptr,
                                  fname, base_index,
                                  gen_geom, gen_graph, gen_hg );
  }

private:

  Zoltan_Struct * ZZ_Ptr; 

};

#endif

