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

// @HEADER
// ************************************************************************
// 
//            Zoltan_CPP: An Object-Oriented Interface To Zoltan
//                    Copyright (2001) Sandia Corporation
// 
// Questions? Contact Robert J. Hoekstra (rjhoeks@sandia.gov)
//                 or Karen D. Devine (kddevin@sandia.gov)
// 
// ************************************************************************
// @HEADER

#ifndef ZOLTANCPP_H_
#define ZOLTANCPP_H_

#include "zoltan.h"

#include <string>

namespace Zoltan {

//! Zoltan::Zoltan_Object: A class for interfacing the functions of the Zoltan 
//! library in a C++/Object Oriented environment.

/*! The Zoltan::Zoltan_Object class is a wrapper for the C functions at the 
    top level interface of Zoltan.  The primary differences include the removal
    the Zoltan_Struct parameter since this object is now stored in the class.
*/

//KDD  Questions remaining:
//KDD     How should Zoltan_Initialize be called?
//KDD     Does it need to be in the namespace?
//KDD     Should we track whether it is called and, if not, call it during
//KDD     constructor?

class Zoltan_Object
{

public:

  /////////////////////////////////////////////////////////////////////////
  //@{ \name Constructors/Destructors.

  //! Constructor
  /*! This constructor replaces the Zoltan_Create call. Params are the same.
   */

  Zoltan_Object(MPI_Comm communicator = MPI_COMM_WORLD) 
  {
    ZZ_Ptr = Zoltan_Create(communicator);
    if (ZZ_Ptr == NULL) {
      //KDD How should we handle this error condition?
    }
    ZOS_Ptr = NULL;
  }

  //! Destructor
  ~Zoltan_Object()
  {
    Zoltan_Destroy( &ZZ_Ptr );
  }

  //@}

  /////////////////////////////////////////////////////////////////////////
  //@{ \name General Load Balance Functionality
 
  //! Replaces Zoltan_Set_Param
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

  /////////////////////////////////////////////////////////////////////////
  //! Replaces Zoltan_LB_Partition
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

  //! Replaces Zoltan_LB_Balance
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

  /////////////////////////////////////////////////////////////////////////
  //! Replaces Zoltan_LB_Set_Part_Sizes
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

  /////////////////////////////////////////////////////////////////////////
  int Order    ( int * num_gid_entries,
                 int * num_lid_entries,
                 int num_objs,
                 ZOLTAN_ID_PTR global_ids,
                 ZOLTAN_ID_PTR local_ids,
                 int * rank,
                 int * iperm )
  {
    return Zoltan_Order( ZZ_Ptr,
                         num_gid_entries, num_lid_entries,
                         num_objs, global_ids, local_ids,
                         rank, iperm, ZOS_Ptr );
  }

  /////////////////////////////////////////////////////////////////////////
  //! Replaces Zoltan_LB_Eval
  void LB_Eval( int print_stats,
                int * num_objects,
                float * object_weights,
                int * num_cuts,
                float * cut_weights,
                int * num_boundary_objects,
                int * num_adj_procs )
  {
    Zoltan_LB_Eval( ZZ_Ptr, print_stats,
                    num_objects, object_weights,
                    num_cuts, cut_weights,
                    num_boundary_objects, num_adj_procs );
  }


  /////////////////////////////////////////////////////////////////////////
  //! Replaces Zoltan_RCB_Box
  void RCB_Box( int part,
                int * ndim, 
                double * xmin, 
                double * ymin, 
                double * zmin, 
                double * xmax, 
                double * ymax, 
                double * zmax )
  {
    Zoltan_RCB_Box( ZZ_Ptr, part, ndim, xmin, ymin, zmin, xmax, ymax, zmax);
  }

  /////////////////////////////////////////////////////////////////////////
  //! Replaces Zoltan_LB_Free_Part
  int LB_Free_Part( ZOLTAN_ID_PTR * global_ids,
                    ZOLTAN_ID_PTR * local_ids,
                    int ** procs,
                    int ** to_part )
  {
    return Zoltan_LB_Free_Part( global_ids, local_ids, procs, to_part );
  }

  /////////////////////////////////////////////////////////////////////////
  //! Replaces Zoltan_LB_Free_Data
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

  //@}

  /////////////////////////////////////////////////////////////////////////
  //@{ \name Support for direct access to Zoltan callback functionality

  //! Generic callback support
  int Set_Fn  ( ZOLTAN_FN_TYPE fn_type,
                void (*fn_ptr)(),
                void * data = 0 )
  {
    return Zoltan_Set_Fn( ZZ_Ptr, fn_type, fn_ptr, data );
  }

  //! Individual callback support

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


  //@}

  /////////////////////////////////////////////////////////////////////////
  //@{ \name Decomposition Augmentation

  //! Replaces Zoltan_LB_Point_PP_Assign
  int LB_Point_PP_Assign ( double * coords,
                           int * proc,
                           int * part )
  {
    return Zoltan_LB_Point_PP_Assign( ZZ_Ptr, coords, proc, part );
  }

  //! Replaces Zoltan_LB_Point_Assign
  int LB_Point_Assign ( double * coords,
                        int * proc )
  {
    return Zoltan_LB_Point_Assign( ZZ_Ptr, coords, proc );
  }


  /////////////////////////////////////////////////////////////////////////
  //! Replaces Zoltan_LB_Box_PP_Assign
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

  //! Replaces Zoltan_LB_Box_Assign
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

  //@}
 
  /////////////////////////////////////////////////////////////////////////
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
    return Zoltan_Invert_Lists( ZZ_Ptr,
                                num_import, import_global_ids, import_local_ids,
                                import_procs, import_to_part,
                                num_export, export_global_ids, export_local_ids,
                                export_procs, export_to_part );
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
    return Zoltan_Compute_Destinations( ZZ_Ptr,
                                        num_import, import_global_ids, 
                                        import_local_ids, import_procs,
                                        num_export, export_global_ids, 
                                        export_local_ids, export_procs );
  }

  /////////////////////////////////////////////////////////////////////////
  //! Replaces Zoltan_Migrate
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
    return Zoltan_Help_Migrate( ZZ_Ptr,
                                num_import, import_global_ids, import_local_ids,
                                import_procs,
                                num_export, export_global_ids, export_local_ids,
                                export_procs );
  }

  //@}

  /////////////////////////////////////////////////////////////////////////
  //@{ \name File Generation Functionality

  //! Replaces Zoltan_Generate_Files
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

  //@}

  /////////////////////////////////////////////////////////////////////////
  //@{ \name Extra
 
  //! Direct access to underlying Zoltan_Struct Object
  //  Zoltan_Struct * Return_Zoltan_Struct()
  //  { return ZZ_Ptr; }

  //@}

private:

  Zoltan_Struct * ZZ_Ptr; 
  Zoltan_Order_Struct * ZOS_Ptr;

};

} //namespace Zoltan

#endif
