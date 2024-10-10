// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// ************************************************************************
// C++ class representing a Zoltan_Struct object.
//
// Assumption: Zoltan_Initialize has been called prior to creating
// a Zoltan object.
// ************************************************************************

#ifndef ZOLTAN_CPP_H_
#define ZOLTAN_CPP_H_

#include "zoltan.h"
#include "zoltan_comm_cpp.h"
#include "zoltan_dd_cpp.h"
#include <cstdlib>

/* F90 names must be less than 31 characters, support old name */
#define Set_HG_Size_Edge_Weights_Fn    Set_HG_Size_Edge_Wts_Fn 
#define Set_HG_Edge_Weights_Fn         Set_HG_Edge_Wts_Fn 

#ifdef TFLOP
  #include <string.h>
#else
  #include <string>
#endif

class Zoltan {

public:

  // Constructor

  Zoltan (const MPI_Comm &communicator = zoltan_get_global_comm()) 
  {
    this->ZZ_Ptr = Zoltan_Create(communicator);

    // int fail = (this->ZZ_Ptr == NULL);  should catch this exception
  }

  // Copy constructor

  Zoltan (const Zoltan &zz)
  {
  this->ZZ_Ptr = Zoltan_Copy(zz.ZZ_Ptr);
  }

  // Destructor

  ~Zoltan()
  {
    // Warning: Zoltan_Destroy calls MPI.   
    // Do not call MPI_Finalize() before this destructor gets called. 
    // Ensure that ZoltanObject's created on the stack are deleted 
    // before MPI_Finalize().
    // Alternatively, you can allocate and destroy Zoltans explicitly:
    //
    //   MPI_Init(...);
    //   Zoltan_Initialize(...);
    //   Zoltan *zz = new Zoltan();
    //    ... more code ...
    //   delete zz;
    //   MPI_Finalize();

    Zoltan_Destroy( &(this->ZZ_Ptr) );
  }

  // Copy operator

  Zoltan & operator= (const Zoltan &zz)
  {
    Zoltan_Copy_To(this->ZZ_Ptr, zz.ZZ_Ptr);

    return *this;
  }

  // Wrappers for Zoltan functions

  static int LB_Free_Part( ZOLTAN_ID_PTR *global_ids,
                    ZOLTAN_ID_PTR *local_ids,
                    int **procs,
                    int **to_part )
  {
    return Zoltan_LB_Free_Part( global_ids, local_ids, procs, to_part );
  }

  int Set_Param( const std::string & param, const std::string & value )
  {
    return Zoltan_Set_Param( ZZ_Ptr, param.c_str(), value.c_str() );
  }

  //! Replaces Zoltan_Set_Param_Vec
  int Set_Param_Vec( const std::string & param, const std::string & value, 
                     const int &index )
  {
    return Zoltan_Set_Param_Vec( ZZ_Ptr, param.c_str(), value.c_str(), index);
  }

  int LB_Partition ( int &changes,
                     int &num_gid_entries,
                     int &num_lid_entries,
                     int &num_import,
                     ZOLTAN_ID_PTR &import_global_ids,
                     ZOLTAN_ID_PTR &import_local_ids,
                     int * &import_procs,
                     int * &import_to_part,
                     int &num_export,
                     ZOLTAN_ID_PTR &export_global_ids,
                     ZOLTAN_ID_PTR &export_local_ids,
                     int * &export_procs,
                     int * &export_to_part )
  {
    return Zoltan_LB_Partition( ZZ_Ptr, &changes,
                                &num_gid_entries, &num_lid_entries,
                                &num_import, &import_global_ids, &import_local_ids,
                                &import_procs, &import_to_part,
                                &num_export, &export_global_ids, &export_local_ids,
                                &export_procs, &export_to_part );
  }

  int LB_Set_Part_Sizes( const int &global_num,
                         const int &len,
                         int * part_ids,
                         int * wgt_idx,
                         float * part_sizes )
  {
    return Zoltan_LB_Set_Part_Sizes( ZZ_Ptr,
                                     global_num, len,
                                     part_ids, wgt_idx, part_sizes );
  }

  /* This method is deprecated */
  int Order    ( int num_gid_entries,
                 int num_objs,
                 ZOLTAN_ID_PTR global_ids,
                 int *rank,
                 int * /* iperm */ )
  {
    return Order(  num_gid_entries, num_objs, global_ids,
		  (ZOLTAN_ID_PTR)rank);
  }

  int Order    ( int num_gid_entries,
                 int num_objs,
                 ZOLTAN_ID_PTR global_ids,
                 ZOLTAN_ID_PTR permuted_global_ids)
  {
    return Zoltan_Order( ZZ_Ptr,
                         num_gid_entries, num_objs, global_ids,
                         permuted_global_ids);
  }

  int Order_Get_Num_Blocks() {
    return Zoltan_Order_Get_Num_Blocks(ZZ_Ptr);
  }

  int Order_Get_Block_Bounds(int block_num, int &first, int &last) {
    return Zoltan_Order_Get_Block_Bounds(ZZ_Ptr, block_num, &first, &last);
  }

  int Order_Get_Block_Size(int block_num) {
    return Zoltan_Order_Get_Block_Size(ZZ_Ptr, block_num);
  }

  int Order_Get_Block_Parent(int block_num) {
    return Zoltan_Order_Get_Block_Parent(ZZ_Ptr, block_num);
  }

  int Order_Get_Num_Leaves() {
    return Zoltan_Order_Get_Num_Leaves(ZZ_Ptr);
  }

  void Order_Get_Block_Leaves(int *leaves) {
    return Zoltan_Order_Get_Block_Leaves(ZZ_Ptr, leaves);
  }

  int Color (int &num_gid_entries,
             const int &num_objs,
             ZOLTAN_ID_PTR global_ids,
             int *color_exp )
  {
    return Zoltan_Color(ZZ_Ptr, num_gid_entries,
			num_objs, global_ids, color_exp);
  }
  /* Simpler method to do coloring */
/*   int Color (int *color_exp); */

  int Color_Test (int &num_gid_entries,
             int &num_lid_entries,
             const int &num_objs,
             ZOLTAN_ID_PTR global_ids,
             ZOLTAN_ID_PTR local_ids,
             int *color_exp )
  {
    return Zoltan_Color_Test(ZZ_Ptr, &num_gid_entries, &num_lid_entries,
      num_objs, global_ids, local_ids, color_exp);
  }

  int LB_Eval_Balance(int print_stats, ZOLTAN_BALANCE_EVAL *eval)
  {
    return Zoltan_LB_Eval_Balance(ZZ_Ptr, print_stats, eval);
  }

  int LB_Eval_Graph(int print_stats, ZOLTAN_GRAPH_EVAL *graph)
  {
    return Zoltan_LB_Eval_Graph(ZZ_Ptr, print_stats, graph);
  }

  int LB_Eval_HG(int print_stats, ZOLTAN_HG_EVAL *hg)
  {
    return Zoltan_LB_Eval_HG(ZZ_Ptr, print_stats, hg);
  }

  int LB_Eval( const int &print_stats, ZOLTAN_BALANCE_EVAL *eval, 
               ZOLTAN_GRAPH_EVAL *graph, ZOLTAN_HG_EVAL *hg)
  {
    return Zoltan_LB_Eval( ZZ_Ptr, print_stats, eval, graph, hg); 
  }

  int RCB_Box( const int &part,
                int &ndim, 
                double &xmin, 
                double &ymin, 
                double &zmin, 
                double &xmax, 
                double &ymax, 
                double &zmax )
  {
    return Zoltan_RCB_Box( ZZ_Ptr,part,&ndim,&xmin,&ymin,&zmin,&xmax,&ymax,&zmax);
  }

  int Set_Fn  ( const ZOLTAN_FN_TYPE &fn_type,
                void (*fn_ptr)(),
                void * data = 0 )
  {
    return Zoltan_Set_Fn( ZZ_Ptr, fn_type, fn_ptr, data );
  }

  // Individual callback support

  ///--------------------------
  int Set_Part_Multi_Fn  ( ZOLTAN_PART_MULTI_FN * fn_ptr,
                                void * data = 0 )
  {
    return Zoltan_Set_Part_Multi_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Part_Fn        ( ZOLTAN_PART_FN * fn_ptr,
                                void * data = 0 )
  {
    return Zoltan_Set_Part_Fn( ZZ_Ptr, fn_ptr, data );
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
  int Set_HG_Size_CS_Fn    ( ZOLTAN_HG_SIZE_CS_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_HG_Size_CS_Fn( ZZ_Ptr, fn_ptr, data );
  }
  ///--------------------------
  int Set_HG_CS_Fn    ( ZOLTAN_HG_CS_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_HG_CS_Fn( ZZ_Ptr, fn_ptr, data );
  }
  ///--------------------------
  int Set_HG_Size_Edge_Wts_Fn    ( ZOLTAN_HG_SIZE_EDGE_WTS_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_HG_Size_Edge_Wts_Fn( ZZ_Ptr, fn_ptr, data );
  }
  ///--------------------------
  int Set_HG_Edge_Wts_Fn    ( ZOLTAN_HG_EDGE_WTS_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_HG_Edge_Wts_Fn( ZZ_Ptr, fn_ptr, data );
  }
  ///--------------------------
  int Set_Hier_Num_Levels_Fn( ZOLTAN_HIER_NUM_LEVELS_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Hier_Num_Levels_Fn( ZZ_Ptr, fn_ptr, data );
  }
  ///--------------------------
  int Set_Hier_Part_Fn( ZOLTAN_HIER_PART_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Hier_Part_Fn( ZZ_Ptr, fn_ptr, data );
  }
  ///--------------------------
  int Set_Hier_Method_Fn( ZOLTAN_HIER_METHOD_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Hier_Method_Fn( ZZ_Ptr, fn_ptr, data );
  }
  ///--------------------------
  int Set_Num_Fixed_Obj_Fn    ( ZOLTAN_NUM_FIXED_OBJ_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Num_Fixed_Obj_Fn( ZZ_Ptr, fn_ptr, data );
  }
  ///--------------------------
  int Set_Fixed_Obj_List_Fn    ( ZOLTAN_FIXED_OBJ_LIST_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Fixed_Obj_List_Fn( ZZ_Ptr, fn_ptr, data );
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

  /// Backward compatibility with v3.0
  ///--------------------------
  int Set_Partition_Multi_Fn  ( ZOLTAN_PART_MULTI_FN * fn_ptr,
                                void * data = 0 )
  {
    return Zoltan_Set_Part_Multi_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Partition_Fn        ( ZOLTAN_PART_FN * fn_ptr,
                                void * data = 0 )
  {
    return Zoltan_Set_Part_Fn( ZZ_Ptr, fn_ptr, data );
  }

  ///--------------------------
  int Set_Hier_Partition_Fn( ZOLTAN_HIER_PART_FN * fn_ptr,
                               void * data = 0 )
  {
    return Zoltan_Set_Hier_Part_Fn( ZZ_Ptr, fn_ptr, data );
  }
  ///--------------------------
  /// Version that returns only part assignment
  int LB_Point_PP_Assign ( double * const coords,
                           int &part )
  {
    return Zoltan_LB_Point_PP_Assign( ZZ_Ptr, coords, NULL, &part );
  }

  ///--------------------------
  /// Version that returns part assignment and mapping to proc as well
  int LB_Point_PP_Assign ( double * const coords,
                           int &proc,
                           int &part )
  {
    return Zoltan_LB_Point_PP_Assign( ZZ_Ptr, coords, &proc, &part );
  }

  ///--------------------------
  int LB_Box_PP_Assign ( const double &xmin,
                         const double &ymin,
                         const double &zmin,
                         const double &xmax,
                         const double &ymax,
                         const double &zmax,
                         int * const procs,
                         int &numprocs,
                         int * const parts,
                         int &numparts )
  {
    return Zoltan_LB_Box_PP_Assign( ZZ_Ptr,
                                    xmin, ymin, zmin,
                                    xmax, ymax, zmax,
                                    procs, &numprocs,
                                    parts, &numparts );
  }

  ///--------------------------
  int Invert_Lists             ( const int &num_known,
                                 ZOLTAN_ID_PTR const known_global_ids,
                                 ZOLTAN_ID_PTR const known_local_ids,
                                 int * const known_procs,
                                 int * const known_to_part,
                                 int &num_found,
                                 ZOLTAN_ID_PTR &found_global_ids,
                                 ZOLTAN_ID_PTR &found_local_ids,
                                 int * &found_procs,
                                 int * &found_to_part )
  {
    return Zoltan_Invert_Lists( ZZ_Ptr,
                                num_known, known_global_ids, known_local_ids,
                                known_procs, known_to_part,
                                &num_found, &found_global_ids, &found_local_ids,
                                &found_procs, &found_to_part );
  }

  ///--------------------------
  int Migrate          ( const int &num_import,
                         ZOLTAN_ID_PTR const import_global_ids,
                         ZOLTAN_ID_PTR const import_local_ids,
                         int * const import_procs,
                         int * const import_to_part,
                         const int &num_export,
                         ZOLTAN_ID_PTR const export_global_ids,
                         ZOLTAN_ID_PTR const export_local_ids,
                         int * const export_procs,
                         int * const export_to_part )
  {
    return Zoltan_Migrate( ZZ_Ptr,
                           num_import, import_global_ids, import_local_ids, 
                           import_procs, import_to_part,
                           num_export, export_global_ids, export_local_ids, 
                           export_procs, export_to_part );
  }

  ///--------------------------
  int Generate_Files( std::string & fname,
                      const int &base_index,
                      const int &gen_geom,
                      const int &gen_graph,
                      const int &gen_hg )
  {
    // c_str() is a "const char *", and Zoltan_Generate_Files may re-write
    // the name passed in, so in order to compile we need a non-const
    // "char *" file name to pass to Zoltan_Generate_Files.

#ifdef TFLOP
    char *fn = strdup((char *)fname.c_str());
#else
    char *fn = strdup(fname.c_str());
#endif

    int rc = Zoltan_Generate_Files( ZZ_Ptr, fn, base_index,
                                  gen_geom, gen_graph, gen_hg );

    std::free(fn);

    return rc;
  }

  //---------------------------------------------------
  // Expert functions; see zoltan.h for details

  int Serialize_Size() {
    return Zoltan_Serialize_Size( ZZ_Ptr );
  }

  int Serialize(size_t bufSize, char *buf) {
    return Zoltan_Serialize( ZZ_Ptr, bufSize, buf);
  }

  int Deserialize(size_t bufSize, char *buf) {
    return Zoltan_Deserialize( ZZ_Ptr, bufSize, buf);
  }

  //---------------------------------------------------
  // Access to the "C" language structure, required for use of
  // Zoltan "C" functions.  (Deprecated "C" functions were not
  // included in the C++ interface.)

  Zoltan_Struct *Get_C_Handle(void)
  {
    return ZZ_Ptr;
  }

private:

  Zoltan_Struct * ZZ_Ptr; 

};

#endif
