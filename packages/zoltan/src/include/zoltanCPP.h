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
  }

  //! Destructor
  ~Zoltan_Object()
  {
    // Warning: Zoltan_Destroy calls MPI.  If you are not using the
    // C++ bindings for MPI, along with the C++ bindings for Zoltan,
    // you may be calling MPI_Finalize() before this destructor gets
    // called.  Do this instead:
    //   MPI_Init(...);
    //   Zoltan::Zoltan_Object *zz = NULL;
    //   zz = new Zoltan::Zoltan_Object();
    //    ... more code ...
    //   delete zz;
    //   MPI_Finalize();
    //
    // Zoltan_Object's created on the stack will be deleted at exit,
    //  after MPI_Finalize().

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
    //  Note:  Zoltan_Order_Struct set to NULL.
    return Zoltan_Order( ZZ_Ptr,
                         num_gid_entries, num_lid_entries,
                         num_objs, global_ids, local_ids,
                         rank, iperm, NULL );
  }

  /////////////////////////////////////////////////////////////////////////
  //! Replaces Zoltan_LB_Eval
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


  /////////////////////////////////////////////////////////////////////////
  //! Replaces Zoltan_RCB_Box
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
  // Static methods
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  //@{ \name Distributed Directory Functionality

  static int DD_Create(Zoltan_DD_Directory **dd, MPI_Comm comm, int num_gid,
   int num_lid, int user_length,  int table_length, int debug_level) 
    {
      return Zoltan_DD_Create (dd, comm, num_gid, num_lid, user_length,  
                               table_length, debug_level) ;
    }
  
  static void DD_Destroy (Zoltan_DD_Directory **dd) 
    {
      return Zoltan_DD_Destroy (dd) ;
    }
 
  static int DD_Update (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
    ZOLTAN_ID_PTR lid, ZOLTAN_ID_PTR user, int *partition, int count) 
    {
      return Zoltan_DD_Update (dd, gid, lid, user, partition, count) ;
    }
  
  static int DD_Find (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
    ZOLTAN_ID_PTR lid, ZOLTAN_ID_PTR data, int *partition, int count,
    int *owner) 
    {
      return Zoltan_DD_Find (dd, gid, lid, data, partition, count, owner);
    }
  
  static int DD_Remove (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
    int count)
    {
      return Zoltan_DD_Remove (dd, gid, count);
    }
  
  static int DD_Set_Hash_Fn (Zoltan_DD_Directory *dd,
    unsigned int (*hash) (ZOLTAN_ID_PTR, int, unsigned int))
    {
      return Zoltan_DD_Set_Hash_Fn (dd, hash);
    } 
  
  static void DD_Stats (Zoltan_DD_Directory *dd)
    {
    return Zoltan_DD_Stats (dd) ;
    }
  
  static int DD_Print (Zoltan_DD_Directory *dd)
    {
      return Zoltan_DD_Print (dd) ;
    }

  //@}

  /////////////////////////////////////////////////////////////////////////
  //@{ \name Unstructured Communication Functionality

  // Comment above regarding DD methods applies to Comm methods as well.
  // (C++ wishlist)
 
  static int Comm_Create(ZOLTAN_COMM_OBJ** cobj, int nvals, int *assign, MPI_Comm comm,
    int tag, int *pnvals_recv)
    {
      return Zoltan_Comm_Create(cobj, nvals, assign, comm, tag, pnvals_recv);
    }

  static int Comm_Destroy(ZOLTAN_COMM_OBJ **plan)
    {
      return Zoltan_Comm_Destroy(plan);
    }
      
  static int Comm_Invert_Map( int *lengths_to, int *procs_to, int nsends, int self_msg,
    int **plengths_from, int **pprocs_from, int *pnrecvs, int my_proc,
    int nprocs, int out_of_mem, int tag, MPI_Comm  comm)
    {
      return Zoltan_Comm_Invert_Map( lengths_to, procs_to, nsends, self_msg,
        plengths_from, pprocs_from, pnrecvs, my_proc, nprocs, out_of_mem, 
        tag, comm);
    }
      
  static int Comm_Exchange_Sizes( int *sizes_to, int *procs_to, int  nsends, int  self_msg,
    int *sizes_from, int *procs_from, int  nrecvs, int *total_recv_size,
    int  my_proc, int  tag, MPI_Comm  comm)
    {
      return Zoltan_Comm_Exchange_Sizes(sizes_to, procs_to, nsends, self_msg,
        sizes_from, procs_from, nrecvs, total_recv_size, my_proc, tag, 
        comm);
    }
  
  static int Comm_Resize( ZOLTAN_COMM_OBJ *plan, int *sizes, int tag, int *sum_recv_sizes)
    {
      return Zoltan_Comm_Resize( plan, sizes, tag, sum_recv_sizes);
    }
  
  static int Comm_Do(ZOLTAN_COMM_OBJ *plan, int tag, char *send_data, int nbytes,
    char *recv_data)
    {
      return Zoltan_Comm_Do(plan, tag, send_data, nbytes, recv_data);
    }

  static int Comm_Do_Post( ZOLTAN_COMM_OBJ * plan, int tag, char *send_data, int nbytes,
    char *recv_data)
    {
      return Zoltan_Comm_Do_Post(plan, tag, send_data, nbytes, recv_data);
    }

  static int Comm_Do_Wait(ZOLTAN_COMM_OBJ *plan, int tag, char *send_data,
    int nbytes, char *recv_data)
    {
      return Zoltan_Comm_Do_Wait(plan, tag, send_data, nbytes, recv_data);
    }
  
  static int Comm_Do_Reverse(ZOLTAN_COMM_OBJ *plan, int tag, char *send_data,
    int nbytes, int *sizes, char *recv_data)
    {
      return Zoltan_Comm_Do_Reverse(plan, tag, send_data, nbytes, sizes, 
        recv_data);
    }

  static int Comm_Do_Reverse_Post(ZOLTAN_COMM_OBJ *plan, int tag, char *send_data,
    int nbytes, int *sizes, char *recv_data)
    {
      return Zoltan_Comm_Do_Reverse_Post(plan, tag, send_data, nbytes, sizes, 
        recv_data);
    }

  static int Comm_Do_Reverse_Wait(ZOLTAN_COMM_OBJ *plan, int tag, char *send_data,
    int nbytes, int *sizes, char *recv_data)
    {
      return Zoltan_Comm_Do_Reverse_Wait(plan, tag, send_data, nbytes, sizes, 
        recv_data);
    }
  
  static int Comm_Info( ZOLTAN_COMM_OBJ *plan, int *nsends, int *send_procs,
    int *send_lengths, int *send_nvals, int *send_max_size, int *send_list,
    int *nrecvs, int *recv_procs, int *recv_lengths, int *recv_nvals,
    int *recv_total_size, int *recv_list, int *self_msg)
    {
      return Zoltan_Comm_Info( plan, nsends, send_procs, send_lengths, 
        send_nvals, send_max_size, send_list, nrecvs, recv_procs, recv_lengths, 
        recv_nvals, recv_total_size, recv_list, self_msg);
    }
  
  static int Comm_Invert_Plan(ZOLTAN_COMM_OBJ** plan)
    {
    return Zoltan_Comm_Invert_Plan(plan);
    }

  //@}

  /////////////////////////////////////////////////////////////////////////
  //@{ \name Extra
 
  //! Direct access to underlying Zoltan_Struct Object
  //  Zoltan_Struct * Return_Zoltan_Struct()
  //  { return ZZ_Ptr; }

  //@}

  /////////////////////////////////////////////////////////////////////////
  //@{ \name Print out status of malloc/free calls.  Flag any memory leaks.

  static void Memory_Stats()
     {
     Zoltan_Memory_Stats();
     }
  //@}

  /////////////////////////////////////////////////////////////////////////
  //@{ \name Plauger alignment algorithm

   static int Align(int a)
     {
     return Zoltan_Align(a);
     }

  //@}

private:

  Zoltan_Struct * ZZ_Ptr; 

};

} //namespace Zoltan

#endif
