/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"
#include "all_allo_const.h"
#include "sppr_header"

/*--------------------------------------------------------------------*/
/* procedure name mangling                                            */

#define LOWERCASE   1
#define UPPERCASE   2
#define UNDERSCORE  3
#define UNDERSCORE2 4

#if FMANGLE==LOWERCASE

#define Zfw_Initialize                 zfw_initialize
#define Zfw_Initialize1                zfw_initialize1
#define Zfw_Create                     zfw_create       
#define Zfw_Destroy                    zfw_destroy       
#define Zfw_Memory_Stats               zfw_memory_stats       
#define Zfw_Set_Fn0f                   zfw_set_fn0f
#define Zfw_Set_Fn1f                   zfw_set_fn1f
#define Zfw_Set_Fn2f                   zfw_set_fn2f
#define Zfw_Set_Fn3f                   zfw_set_fn3f
#define Zfw_Set_Fn4f                   zfw_set_fn4f
#define Zfw_Set_Fn5f                   zfw_set_fn5f
#define Zfw_Set_Fn6f                   zfw_set_fn6f
#define Zfw_Set_Fn7f                   zfw_set_fn7f
#define Zfw_Set_Fn8f                   zfw_set_fn8f
#define Zfw_Set_Fn9f                   zfw_set_fn9f
#define Zfw_Set_FnAf                   zfw_set_fnaf
#define Zfw_Set_FnBf                   zfw_set_fnbf
#define Zfw_Set_Fn0s                   zfw_set_fn0s
#define Zfw_Set_Fn1s                   zfw_set_fn1s
#define Zfw_Set_Fn2s                   zfw_set_fn2s
#define Zfw_Set_Fn3s                   zfw_set_fn3s
#define Zfw_Set_Fn4s                   zfw_set_fn4s
#define Zfw_Set_Fn5s                   zfw_set_fn5s
#define Zfw_Set_Fn6s                   zfw_set_fn6s
#define Zfw_Set_Fn7s                   zfw_set_fn7s
#define Zfw_Set_Fn8s                   zfw_set_fn8s
#define Zfw_Set_Fn9s                   zfw_set_fn9s
#define Zfw_Set_FnAs                   zfw_set_fnas
#define Zfw_Set_FnBs                   zfw_set_fnbs
#define Zfw_Set_Param                  zfw_set_param
#define Zfw_Set_Param_Vec              zfw_set_param_vec
#define Zfw_LB_Partition               zfw_lb_partition
#define Zfw_LB_Eval                    zfw_lb_eval
#define Zfw_LB_Set_Part_Sizes          zfw_lb_set_part_sizes
#define Zfw_LB_Point_Assign            zfw_lb_point_assign
#define Zfw_LB_Point_PP_Assign         zfw_lb_point_pp_assign
#define Zfw_LB_Box_Assign              zfw_lb_box_assign
#define Zfw_LB_Box_PP_Assign           zfw_lb_box_pp_assign
#define Zfw_Invert_Lists               zfw_invert_lists
#define Zfw_Compute_Destinations       zfw_compute_destinations
#define Zfw_Migrate                    zfw_migrate  
#define Zfw_Help_Migrate               zfw_help_migrate  
#define Zfw_Order                      zfw_order  
#define Zfw_Generate_Files             zfw_generate_files  
#define Zfw_RCB_Box                    zfw_rcb_box  
#define Zfw_Register_Fort_Malloc       zfw_register_fort_malloc
#define Zfw_Get_Address_int            zfw_get_address_int
#define Zfw_Get_Wgt_Dim                zfw_get_wgt_dim
#define Zfw_Get_Comm_Dim               zfw_get_comm_dim
/* TEMP child_order */
#define Zfw_Reftree_Get_Child_Order         zfw_reftree_get_child_order

#elif FMANGLE==UPPERCASE

#define Zfw_Initialize                 ZFW_INITIALIZE
#define Zfw_Initialize1                ZFW_INITIALIZE1
#define Zfw_Create                     ZFW_CREATE       
#define Zfw_Destroy                    ZFW_DESTROY       
#define Zfw_Memory_Stats               ZFW_MEMORY_STATS  
#define Zfw_Set_Fn0f                   ZFW_SET_FN0F
#define Zfw_Set_Fn1f                   ZFW_SET_FN1F
#define Zfw_Set_Fn2f                   ZFW_SET_FN2F
#define Zfw_Set_Fn3f                   ZFW_SET_FN3F
#define Zfw_Set_Fn4f                   ZFW_SET_FN4F
#define Zfw_Set_Fn5f                   ZFW_SET_FN5F
#define Zfw_Set_Fn6f                   ZFW_SET_FN6F
#define Zfw_Set_Fn7f                   ZFW_SET_FN7F
#define Zfw_Set_Fn8f                   ZFW_SET_FN8F
#define Zfw_Set_Fn9f                   ZFW_SET_FN9F
#define Zfw_Set_FnAf                   ZFW_SET_FNAF
#define Zfw_Set_FnBf                   ZFW_SET_FNBF
#define Zfw_Set_Fn0s                   ZFW_SET_FN0S
#define Zfw_Set_Fn1s                   ZFW_SET_FN1S
#define Zfw_Set_Fn2s                   ZFW_SET_FN2S
#define Zfw_Set_Fn3s                   ZFW_SET_FN3S
#define Zfw_Set_Fn4s                   ZFW_SET_FN4S
#define Zfw_Set_Fn5s                   ZFW_SET_FN5S
#define Zfw_Set_Fn6s                   ZFW_SET_FN6S
#define Zfw_Set_Fn7s                   ZFW_SET_FN7S
#define Zfw_Set_Fn8s                   ZFW_SET_FN8S
#define Zfw_Set_Fn9s                   ZFW_SET_FN9S
#define Zfw_Set_FnAs                   ZFW_SET_FNAS
#define Zfw_Set_FnBs                   ZFW_SET_FNBS
#define Zfw_Set_Param                  ZFW_SET_PARAM
#define Zfw_Set_Param_Vec              ZFW_SET_PARAM_VEC
#define Zfw_LB_Partition               ZFW_LB_PARTITION
#define Zfw_LB_Eval                    ZFW_LB_EVAL
#define Zfw_LB_Set_Part_Sizes          ZFW_LB_SET_PART_SIZES
#define Zfw_LB_Point_Assign            ZFW_LB_POINT_ASSIGN
#define Zfw_LB_Point_PP_Assign         ZFW_LB_POINT_PP_ASSIGN
#define Zfw_LB_Box_Assign              ZFW_LB_BOX_ASSIGN
#define Zfw_LB_Box_PP_Assign           ZFW_LB_BOX_PP_ASSIGN
#define Zfw_Invert_Lists               ZFW_INVERT_LISTS
#define Zfw_Compute_Destinations       ZFW_COMPUTE_DESTINATIONS  
#define Zfw_Migrate                    ZFW_MIGRATE  
#define Zfw_Help_Migrate               ZFW_HELP_MIGRATE  
#define Zfw_Order                      ZFW_ORDER  
#define Zfw_Generate_Files             ZFW_GENERATE_FILES  
#define Zfw_RCB_Box                    ZFW_RCB_BOX  
#define Zfw_Register_Fort_Malloc       ZFW_REGISTER_FORT_MALLOC
#define Zfw_Get_Address_int            ZFW_GET_ADDRESS_INT
#define Zfw_Get_Comm_Dim               ZFW_GET_COMM_DIM
/* TEMP child_order */
#define Zfw_Reftree_Get_Child_Order         ZFW_REFTREE_GET_CHILD_ORDER

#elif FMANGLE==UNDERSCORE

#define Zfw_Initialize                 zfw_initialize_
#define Zfw_Initialize1                zfw_initialize1_
#define Zfw_Create                     zfw_create_
#define Zfw_Destroy                    zfw_destroy_
#define Zfw_Memory_Stats               zfw_memory_stats_
#define Zfw_Set_Fn0f                   zfw_set_fn0f_
#define Zfw_Set_Fn1f                   zfw_set_fn1f_
#define Zfw_Set_Fn2f                   zfw_set_fn2f_
#define Zfw_Set_Fn3f                   zfw_set_fn3f_
#define Zfw_Set_Fn4f                   zfw_set_fn4f_
#define Zfw_Set_Fn5f                   zfw_set_fn5f_
#define Zfw_Set_Fn6f                   zfw_set_fn6f_
#define Zfw_Set_Fn7f                   zfw_set_fn7f_
#define Zfw_Set_Fn8f                   zfw_set_fn8f_
#define Zfw_Set_Fn9f                   zfw_set_fn9f_
#define Zfw_Set_FnAf                   zfw_set_fnaf_
#define Zfw_Set_FnBf                   zfw_set_fnbf_
#define Zfw_Set_Fn0s                   zfw_set_fn0s_
#define Zfw_Set_Fn1s                   zfw_set_fn1s_
#define Zfw_Set_Fn2s                   zfw_set_fn2s_
#define Zfw_Set_Fn3s                   zfw_set_fn3s_
#define Zfw_Set_Fn4s                   zfw_set_fn4s_
#define Zfw_Set_Fn5s                   zfw_set_fn5s_
#define Zfw_Set_Fn6s                   zfw_set_fn6s_
#define Zfw_Set_Fn7s                   zfw_set_fn7s_
#define Zfw_Set_Fn8s                   zfw_set_fn8s_
#define Zfw_Set_Fn9s                   zfw_set_fn9s_
#define Zfw_Set_FnAs                   zfw_set_fnas_
#define Zfw_Set_FnBs                   zfw_set_fnbs_
#define Zfw_Set_Param                  zfw_set_param_
#define Zfw_Set_Param_Vec              zfw_set_param_vec_
#define Zfw_LB_Partition               zfw_lb_partition_
#define Zfw_LB_Eval                    zfw_lb_eval_
#define Zfw_LB_Set_Part_Sizes          zfw_lb_set_part_sizes_
#define Zfw_LB_Point_Assign            zfw_lb_point_assign_
#define Zfw_LB_Point_PP_Assign         zfw_lb_point_pp_assign_
#define Zfw_LB_Box_Assign              zfw_lb_box_assign_
#define Zfw_LB_Box_PP_Assign           zfw_lb_box_pp_assign_
#define Zfw_Invert_Lists               zfw_invert_lists_
#define Zfw_Compute_Destinations       zfw_compute_destinations_
#define Zfw_Migrate                    zfw_migrate_
#define Zfw_Help_Migrate               zfw_help_migrate_  
#define Zfw_Order                      zfw_order_  
#define Zfw_Generate_Files             zfw_generate_files_ 
#define Zfw_RCB_Box                    zfw_rcb_box_
#define Zfw_Register_Fort_Malloc       zfw_register_fort_malloc_
#define Zfw_Get_Address_int            zfw_get_address_int_
#define Zfw_Get_Wgt_Dim                zfw_get_wgt_dim_
#define Zfw_Get_Comm_Dim               zfw_get_comm_dim_
/* TEMP child_order */
#define Zfw_Reftree_Get_Child_Order         zfw_reftree_get_child_order_

#elif FMANGLE==UNDERSCORE2

#define Zfw_Initialize                 zfw_initialize__
#define Zfw_Initialize1                zfw_initialize1__
#define Zfw_Create                     zfw_create__
#define Zfw_Destroy                    zfw_destroy__
#define Zfw_Memory_Stats               zfw_memory_stats__
#define Zfw_Set_Fn0f                   zfw_set_fn0f__
#define Zfw_Set_Fn1f                   zfw_set_fn1f__
#define Zfw_Set_Fn2f                   zfw_set_fn2f__
#define Zfw_Set_Fn3f                   zfw_set_fn3f__
#define Zfw_Set_Fn4f                   zfw_set_fn4f__
#define Zfw_Set_Fn5f                   zfw_set_fn5f__
#define Zfw_Set_Fn6f                   zfw_set_fn6f__
#define Zfw_Set_Fn7f                   zfw_set_fn7f__
#define Zfw_Set_Fn8f                   zfw_set_fn8f__
#define Zfw_Set_Fn9f                   zfw_set_fn9f__
#define Zfw_Set_FnAf                   zfw_set_fnaf__
#define Zfw_Set_FnBf                   zfw_set_fnbf__
#define Zfw_Set_Fn0s                   zfw_set_fn0s__
#define Zfw_Set_Fn1s                   zfw_set_fn1s__
#define Zfw_Set_Fn2s                   zfw_set_fn2s__
#define Zfw_Set_Fn3s                   zfw_set_fn3s__
#define Zfw_Set_Fn4s                   zfw_set_fn4s__
#define Zfw_Set_Fn5s                   zfw_set_fn5s__
#define Zfw_Set_Fn6s                   zfw_set_fn6s__
#define Zfw_Set_Fn7s                   zfw_set_fn7s__
#define Zfw_Set_Fn8s                   zfw_set_fn8s__
#define Zfw_Set_Fn9s                   zfw_set_fn9s__
#define Zfw_Set_FnAs                   zfw_set_fnas__
#define Zfw_Set_FnBs                   zfw_set_fnbs__
#define Zfw_Set_Param                  zfw_set_param__
#define Zfw_Set_Param_Vec              zfw_set_param_vec__
#define Zfw_LB_Partition               zfw_lb_partition__
#define Zfw_LB_Eval                    zfw_lb_eval__
#define Zfw_LB_Set_Part_Sizes          zfw_lb_set_part_sizes__
#define Zfw_LB_Point_Assign            zfw_lb_point_assign__
#define Zfw_LB_Point_PP_Assign         zfw_lb_point_pp_assign__
#define Zfw_LB_Box_Assign              zfw_lb_box_assign__
#define Zfw_LB_Box_PP_Assign           zfw_lb_box_pp_assign__
#define Zfw_Invert_Lists               zfw_invert_lists__
#define Zfw_Compute_Destinations       zfw_compute_destinations__
#define Zfw_Migrate                    zfw_migrate__
#define Zfw_Help_Migrate               zfw_help_migrate__
#define Zfw_Order                      zfw_order__
#define Zfw_Generate_Files             zfw_generate_files__
#define Zfw_RCB_Box                    zfw_rcb_box__
#define Zfw_Register_Fort_Malloc       zfw_register_fort_malloc__
#define Zfw_Get_Address_int            zfw_get_address_int__
#define Zfw_Get_Wgt_Dim                zfw_get_wgt_dim__
#define Zfw_Get_Comm_Dim               zfw_get_comm_dim__
/* TEMP child_order */
#define Zfw_Reftree_Get_Child_Order         zfw_reftree_get_child_order__

#endif /* FMANGLE */

/*--------------------------------------------------------------------*/
/* Variables                                                          */

static struct Zoltan_Struct *Zoltan_Current;
void Zoltan_Reftree_Get_Child_Order(struct Zoltan_Struct *, int *, int *);


/*--------------------------------------------------------------------*/
/* Utilities                                                          */

/* some MPI implementations may require conversion between a Fortran
   communicator and a C communicator.  This routine is used to perform the
   conversion.  It may need different forms for different MPI libraries. */

MPI_Comm Zoltan_comm_f2c(int *f_comm)
{
#ifndef NO_MPI2
/* MPI 2 provides a standard way of doing this */
   return MPI_Comm_f2c((MPI_Fint)(*f_comm));
#else
/* will probably need some special cases here */
/* when in doubt, just return the input */
   return (MPI_Comm)(*f_comm);
#endif
}

/*****************************************************************************/
/* These routines get the address of an array allocated by fortran and
   return it */
#ifdef PTR_64BIT
void Zfw_Get_Address_int(int *addr,
                         long *ret_addr)
{
   if (sizeof(long) != sizeof(int *)) {
     ZOLTAN_PRINT_ERROR(-1, "Zfw_Get_Address_int", 
       "sizeof(long) != sizeof(int *); F90 allocation will not work properly.");
   }
   *ret_addr = (long)addr;
}
#else
void Zfw_Get_Address_int(int *addr,
                         int *ret_addr)
{
   if (sizeof(int) != sizeof(int *)) {
     ZOLTAN_PRINT_ERROR(-1, "Zfw_Get_Address_int", 
       "sizeof(int) != sizeof(int *); F90 allocation will not work properly.");
   }
   *ret_addr = (int)addr;
}
#endif  /* PTR_64BIT */

/*****************************************************************************/
int Zfw_Get_Wgt_Dim(int *addr_lb, int *nbytes)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   return lb->Obj_Weight_Dim;
}

/*****************************************************************************/
int Zfw_Get_Comm_Dim(int *addr_lb, int *nbytes)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   return lb->Edge_Weight_Dim;
}

/*--------------------------------------------------------------------*/
/* Reverse wrappers for callbacks                                     */
/*--------------------------------------------------------------------*/

void Zoltan_Partition_Multi_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, int num_obj,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *parts,
  int *ierr)
{
   Zoltan_Current->Get_Partition_Multi_Fort(data,
                       &num_gid_entries, &num_lid_entries, &num_obj,
                       global_id, local_id, parts, ierr);
}

/*****************************************************************************/

int Zoltan_Partition_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  int *ierr)
{
   return Zoltan_Current->Get_Partition_Fort(data,
                                            &num_gid_entries, &num_lid_entries,
                                            global_id, local_id, ierr);
}

/*****************************************************************************/
int Zoltan_Num_Edges_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  int *ierr)
{
   return Zoltan_Current->Get_Num_Edges_Fort(data,
                                            &num_gid_entries, &num_lid_entries,
                                            global_id, local_id, ierr);
}

/*****************************************************************************/
void Zoltan_Num_Edges_Multi_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, int num_obj,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *num_edges,
  int *ierr)
{
   Zoltan_Current->Get_Num_Edges_Multi_Fort(data,
                                            &num_gid_entries, &num_lid_entries,
                                            &num_obj, global_id, local_id, 
                                            num_edges, ierr);
}

/*****************************************************************************/
void Zoltan_Edge_List_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs,
  int wdim, float *nbor_ewgts, int *ierr)
{
   Zoltan_Current->Get_Edge_List_Fort(data, &num_gid_entries, &num_lid_entries,
                                     global_id, local_id,
                                     nbor_global_id, nbor_procs, &wdim,
                                     nbor_ewgts, ierr);
}

/*****************************************************************************/
void Zoltan_Edge_List_Multi_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, int num_obj,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *num_edges,
  ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs,
  int wdim, float *nbor_ewgts, int *ierr)
{
   Zoltan_Current->Get_Edge_List_Multi_Fort(data,
                                            &num_gid_entries, &num_lid_entries,
                                            &num_obj, global_id, local_id,
                                            num_edges,
                                            nbor_global_id, nbor_procs, &wdim,
                                            nbor_ewgts, ierr);
}

/*****************************************************************************/
int Zoltan_Num_Geom_Fort_Wrapper(void *data, int *ierr)
{
   return Zoltan_Current->Get_Num_Geom_Fort(data,ierr);
}

/*****************************************************************************/
void Zoltan_Geom_Multi_Fort_Wrapper(
  void *data, int num_gid_entries, int num_lid_entries, int num_obj,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int num_dim,
  double *geom_vec, int *ierr)
{
   Zoltan_Current->Get_Geom_Multi_Fort(data, &num_gid_entries, &num_lid_entries,
                                       &num_obj, global_id, local_id,
                                       &num_dim, geom_vec, ierr);
}
/*****************************************************************************/
void Zoltan_Geom_Fort_Wrapper(
  void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  double *geom_vec, int *ierr)
{
   Zoltan_Current->Get_Geom_Fort(data, &num_gid_entries, &num_lid_entries,
                                global_id, local_id, geom_vec, ierr);
}

/*****************************************************************************/
int Zoltan_Num_Obj_Fort_Wrapper(void *data, int *ierr)
{
   return Zoltan_Current->Get_Num_Obj_Fort(data, ierr);
}

/*****************************************************************************/
void Zoltan_Obj_List_Fort_Wrapper(void *data,
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
  int wdim, float *objwgts, int *ierr)
{
   Zoltan_Current->Get_Obj_List_Fort(data, &num_gid_entries, &num_lid_entries,
                                    global_ids, local_ids, &wdim,
                                    objwgts, ierr);
}

/*****************************************************************************/
int Zoltan_First_Obj_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, 
  ZOLTAN_ID_PTR first_global_id,
  ZOLTAN_ID_PTR first_local_id,
  int wdim, float *first_obj_wgt, int *ierr)
{
   return Zoltan_Current->Get_First_Obj_Fort(data, 
                                            &num_gid_entries, &num_lid_entries,
                                            first_global_id,
                                            first_local_id, &wdim,
                                            first_obj_wgt, ierr);
}

/*****************************************************************************/
int Zoltan_Next_Obj_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, 
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  ZOLTAN_ID_PTR next_global_id, ZOLTAN_ID_PTR next_local_id,
  int wdim, float *next_obj_wgt, int *ierr)
{
   return Zoltan_Current->Get_Next_Obj_Fort(data, 
                                           &num_gid_entries, &num_lid_entries, 
                                           global_id, local_id,
                                           next_global_id, next_local_id,
                                           &wdim, next_obj_wgt, ierr);
}

/*****************************************************************************/
int Zoltan_Num_Border_Obj_Fort_Wrapper(void *data, int nbor_proc, int *ierr)
{
   return Zoltan_Current->Get_Num_Border_Obj_Fort(data, &nbor_proc, ierr);
}

/*****************************************************************************/
void Zoltan_Border_Obj_List_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, 
  int nbor_proc,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
  int wdim, float *objwgts, int *ierr)
{
   Zoltan_Current->Get_Border_Obj_List_Fort(data, 
                                           &num_gid_entries, &num_lid_entries, 
                                           &nbor_proc, global_ids,
                                           local_ids, &wdim, objwgts, ierr);
}

/*****************************************************************************/
int Zoltan_First_Border_Obj_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int nbor_proc,
  ZOLTAN_ID_PTR first_global_id,
  ZOLTAN_ID_PTR first_local_id,
  int wdim, float *first_obj_wgt,
  int *ierr)
{
   return Zoltan_Current->Get_First_Border_Obj_Fort(data, 
                                                   &num_gid_entries, 
                                                   &num_lid_entries, 
                                                   &nbor_proc,
                                                   first_global_id,
                                                   first_local_id, &wdim,
                                                   first_obj_wgt, ierr);
}

/*****************************************************************************/
int Zoltan_Next_Border_Obj_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id,
  ZOLTAN_ID_PTR local_id, int nbor_proc,
  ZOLTAN_ID_PTR next_global_id,
  ZOLTAN_ID_PTR next_local_id,
  int wdim, float *next_obj_wgt,
  int *ierr)
{
   return Zoltan_Current->Get_Next_Border_Obj_Fort(data, 
                                                  &num_gid_entries,
                                                  &num_lid_entries,
                                                  global_id, local_id,
                                                  &nbor_proc, next_global_id,
                                                  next_local_id, &wdim,
                                                  next_obj_wgt, ierr);
}

/*****************************************************************************/
int Zoltan_Obj_Size_Fort_Wrapper(void *data, int num_gid_entries,
  int num_lid_entries, ZOLTAN_ID_PTR global_id, 
  ZOLTAN_ID_PTR local_id, int *ierr)
{
   return Zoltan_Current->Get_Obj_Size_Fort(data,
             &num_gid_entries, &num_lid_entries,
             global_id, local_id, ierr);
}

/*****************************************************************************/
void Zoltan_Obj_Size_Multi_Fort_Wrapper(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  int num_ids,
  ZOLTAN_ID_PTR global_ids,
  ZOLTAN_ID_PTR local_ids,
  int *num_bytes,
  int *ierr)
{
   Zoltan_Current->Get_Obj_Size_Multi_Fort(data,
             &num_gid_entries, &num_lid_entries, &num_ids,
             global_ids, local_ids, num_bytes, ierr);
}

/*****************************************************************************/
void Zoltan_Pre_Migrate_PP_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids, int *import_procs, int *import_to_proc,
  int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_proc,
  int *ierr)
{
   Zoltan_Current->Migrate.Pre_Migrate_PP_Fort(data, 
                   &num_gid_entries, &num_lid_entries,
                   &num_import, import_global_ids, import_local_ids, 
                   import_procs, import_to_proc,
                   &num_export, export_global_ids, export_local_ids,
                   export_procs, export_to_proc, ierr);
}

/*****************************************************************************/
void Zoltan_Mid_Migrate_PP_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids, int *import_procs, int *import_to_proc,
  int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_proc,
  int *ierr)
{
   Zoltan_Current->Migrate.Mid_Migrate_PP_Fort(data,
                   &num_gid_entries, &num_lid_entries,
                   &num_import, import_global_ids, import_local_ids,
                   import_procs, import_to_proc,
                   &num_export, export_global_ids, export_local_ids,
                   export_procs, export_to_proc, ierr);
}

/*****************************************************************************/
void Zoltan_Post_Migrate_PP_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids, int *import_procs, int *import_to_proc,
  int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_proc,
  int *ierr)
{
   Zoltan_Current->Migrate.Post_Migrate_PP_Fort(data,
                   &num_gid_entries, &num_lid_entries,
                   &num_import, import_global_ids, import_local_ids,
                   import_procs, import_to_proc,
                   &num_export, export_global_ids, export_local_ids,
                   export_procs, export_to_proc, ierr);

}

/*****************************************************************************/
void Zoltan_Pre_Migrate_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids, int *import_procs,
  int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs,
  int *ierr)
{
   Zoltan_Current->Migrate.Pre_Migrate_Fort(data, 
                                           &num_gid_entries,
                                           &num_lid_entries,
                                           &num_import,
                                           import_global_ids,
                                           import_local_ids, import_procs,
                                           &num_export, export_global_ids,
                                           export_local_ids, export_procs,
                                           ierr);
}

/*****************************************************************************/
void Zoltan_Mid_Migrate_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids, int *import_procs,
  int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs,
  int *ierr)
{
   Zoltan_Current->Migrate.Mid_Migrate_Fort(data,
                                           &num_gid_entries,
                                           &num_lid_entries,
                                           &num_import,
                                           import_global_ids,
                                           import_local_ids, import_procs,
                                           &num_export, export_global_ids,
                                           export_local_ids, export_procs,
                                           ierr);
}

/*****************************************************************************/
void Zoltan_Post_Migrate_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids, int *import_procs,
  int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs,
  int *ierr)
{
   Zoltan_Current->Migrate.Post_Migrate_Fort(data,
                                            &num_gid_entries, &num_lid_entries, 
                                            &num_import,
                                            import_global_ids,
                                            import_local_ids, import_procs,
                                            &num_export, export_global_ids,
                                            export_local_ids, export_procs,
                                            ierr);
}

/*****************************************************************************/
void Zoltan_Pack_Obj_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  int dest_proc, int size, char *buf, int *ierr)
{
   Zoltan_Current->Pack_Obj_Fort(data, 
                                        &num_gid_entries, &num_lid_entries, 
                                        global_id, local_id,
                                        &dest_proc, &size, buf, ierr);
}

/*****************************************************************************/
void Zoltan_Pack_Obj_Multi_Fort_Wrapper(
  void *data,
  int num_gid_entries,
  int num_lid_entries,
  int num_ids,
  ZOLTAN_ID_PTR global_ids,
  ZOLTAN_ID_PTR local_ids,
  int *dest_proc,
  int *size,
  int *index,
  char *buffer,
  int *ierr)
{
  int factor = sizeof(int) / sizeof(char);
  int i;
 
  /* Convert index array from indices into char * to indices into int *. */
  /* Add 1 for F90 one-based indexing. */
  for (i = 0; i < num_ids; i++) {
    /* Sanity check */
    if (index[i] % factor != 0) {
      ZOLTAN_PRINT_ERROR(-1, "Zoltan_Pack_Obj_Multi_Fort_Wrapper", 
                         "Alignment problem in index array.");
      
      *ierr = ZOLTAN_FATAL;
      return;
    }
    index[i] = index[i]/factor + 1;
  }
     
  Zoltan_Current->Pack_Obj_Multi_Fort(data, &num_gid_entries, &num_lid_entries,
                                      &num_ids, global_ids, local_ids,
                                      dest_proc, size, index, buffer, ierr);

  /* Restore index array to original condition. */
  for (i = 0; i < num_ids; i++) 
    index[i] = (index[i] - 1) * factor;
}

/*****************************************************************************/
void Zoltan_Unpack_Obj_Fort_Wrapper(void *data, int num_gid_entries,
                                ZOLTAN_ID_PTR global_id, int size,
                                char *buf, int *ierr)
{
   Zoltan_Current->Unpack_Obj_Fort(data, &num_gid_entries, 
                                          global_id, &size, buf, ierr);
}

/*****************************************************************************/
void Zoltan_Unpack_Obj_Multi_Fort_Wrapper(
  void *data,
  int num_gid_entries,
  int num_ids,
  ZOLTAN_ID_PTR global_ids,
  int *size,
  int *index,
  char *buffer,
  int *ierr)
{
  int factor = sizeof(int) / sizeof(char);
  int i;
 
  /* Convert index array from indices into char * to indices into int *. */
  /* Add 1 for F90 one-based indexing. */
  for (i = 0; i < num_ids; i++) {
    /* Sanity check */
    if (index[i] % factor != 0) {
      ZOLTAN_PRINT_ERROR(-1, "Zoltan_Pack_Obj_Multi_Fort_Wrapper", 
                         "Alignment problem in index array.");
      
      *ierr = ZOLTAN_FATAL;
      return;
    }
    index[i] = index[i]/factor + 1;
  }
     
  Zoltan_Current->Unpack_Obj_Multi_Fort(data, &num_gid_entries, &num_ids,
                                        global_ids, size, index, buffer, ierr);

  /* Restore index array to original condition. */
  for (i = 0; i < num_ids; i++) 
    index[i] = (index[i] - 1) * factor;
}


/*****************************************************************************/
int Zoltan_Num_Coarse_Obj_Fort_Wrapper(void *data, int *ierr)
{
   return Zoltan_Current->Get_Num_Coarse_Obj_Fort(data, ierr);
}

/*****************************************************************************/
void Zoltan_Coarse_Obj_List_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_ids,
  ZOLTAN_ID_PTR local_ids, int *assigned, int *num_vert,
  ZOLTAN_ID_PTR vertices, int *in_order, ZOLTAN_ID_PTR in_vertex,
  ZOLTAN_ID_PTR out_vertex, int *ierr)
{
   Zoltan_Current->Get_Coarse_Obj_List_Fort(data, 
                                           &num_gid_entries, &num_lid_entries,
                                           global_ids, local_ids,
                                           assigned, num_vert, vertices,
                                           in_order, in_vertex, out_vertex,
                                           ierr);
}

/*****************************************************************************/
int Zoltan_First_Coarse_Obj_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, 
  ZOLTAN_ID_PTR global_id,
  ZOLTAN_ID_PTR local_id, int *assigned,
  int *num_vert, ZOLTAN_ID_PTR vertices,
  int *in_order, ZOLTAN_ID_PTR in_vertex,
  ZOLTAN_ID_PTR out_vertex, int *ierr)
{
   return Zoltan_Current->Get_First_Coarse_Obj_Fort(data, 
                                                   &num_gid_entries, 
                                                   &num_lid_entries,
                                                   global_id, local_id,
                                                   assigned, num_vert, vertices,
                                                   in_order, in_vertex,
                                                   out_vertex, ierr);
}

/*****************************************************************************/
int Zoltan_Next_Coarse_Obj_Fort_Wrapper(void *data, int num_gid_entries, 
  int num_lid_entries, ZOLTAN_ID_PTR global_id,
  ZOLTAN_ID_PTR local_id, 
  ZOLTAN_ID_PTR next_global_id, 
  ZOLTAN_ID_PTR next_local_id,
  int *assigned,
  int *num_vert, ZOLTAN_ID_PTR vertices,
  ZOLTAN_ID_PTR in_vertex, ZOLTAN_ID_PTR out_vertex, int *ierr)
{
   return Zoltan_Current->Get_Next_Coarse_Obj_Fort(data, &num_gid_entries,
                                                  &num_lid_entries,
                                                  global_id, local_id,
                                                  next_global_id, next_local_id,
                                                  assigned, num_vert, vertices,
                                                  in_vertex, out_vertex, ierr);
}

/*****************************************************************************/
int Zoltan_Num_Child_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, 
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  int *ierr)
{
   return Zoltan_Current->Get_Num_Child_Fort(data, 
                                            &num_gid_entries, &num_lid_entries,
                                            global_id, local_id, ierr);
}

/*****************************************************************************/
void Zoltan_Child_List_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries, 
  ZOLTAN_ID_PTR parent_gid,
  ZOLTAN_ID_PTR parent_lid, ZOLTAN_ID_PTR child_gids,
  ZOLTAN_ID_PTR child_lids, int *assigned,
  int *num_vert, ZOLTAN_ID_PTR vertices,
  ZOLTAN_REF_TYPE *ref_type, ZOLTAN_ID_PTR in_vertex,
  ZOLTAN_ID_PTR out_vertex, int *ierr)
{
   Zoltan_Current->Get_Child_List_Fort(data, &num_gid_entries, &num_lid_entries,
                                      parent_gid, parent_lid,
                                      child_gids, child_lids, assigned,
                                      num_vert, vertices,
                                      ref_type, in_vertex, out_vertex, ierr);
}

/*****************************************************************************/
void Zoltan_Child_Weight_Fort_Wrapper(void *data, 
  int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  int wgt_dim, float *obj_wgt, int *ierr)
{
   Zoltan_Current->Get_Child_Weight_Fort(data, 
                                        &num_gid_entries, &num_lid_entries,
                                        global_id, local_id, &wgt_dim,
                                        obj_wgt, ierr);
}

/*****************************************************************************/
/*--------------------------------------------------------------------*/
/* C wrapper functions                                                */
/*--------------------------------------------------------------------*/

/*****************************************************************************/
int Zfw_Initialize(float *ver)
{
   int myArgc;
   char **myArgv;
   int result;
   myArgc = 1;
   myArgv = (char **) ZOLTAN_MALLOC((myArgc+1)*sizeof(char *));
   myArgv[0] = "unknown";
   myArgv[1] = NULL;
   result = Zoltan_Initialize(myArgc,myArgv,ver);
   ZOLTAN_FREE(&myArgv);
   return result;
}

/*****************************************************************************/
int Zfw_Initialize1(int *argc, int *argv, int *starts, float *ver)
{
   int i, j, result;
   char **myArgv;
   myArgv = (char **) ZOLTAN_MALLOC(((*argc)+1)*sizeof(char *));
   for (i=0; i<(*argc); i++) {
      myArgv[i] = (char *) ZOLTAN_MALLOC((starts[i+1]-starts[i]+1)*sizeof(char));
      for (j=0; j<starts[i+1]-starts[i]; j++) {
         myArgv[i][j] = (char) argv[starts[i]+j-1];
      }
      myArgv[i][starts[i+1]-starts[i]] = '\0';
   }
   myArgv[*argc] = NULL;
   result = Zoltan_Initialize(*argc,myArgv,ver);
   for (i=0; i<(*argc); i++) 
     ZOLTAN_FREE(&(myArgv[i]));
   ZOLTAN_FREE(&myArgv);
   return result;
}

/*****************************************************************************/
void Zfw_Create(int *f_communicator, int *addr_lb, int *nbytes)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   MPI_Comm c_communicator;
   c_communicator = Zoltan_comm_f2c(f_communicator);
   lb = Zoltan_Create(c_communicator);
   lb->Fortran = 1;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {addr_lb[i] = (int)*p; p++;}
}

/*****************************************************************************/
void Zfw_Destroy(int *addr_lb, int *nbytes)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Destroy(&lb);
}

/*****************************************************************************/
void Zfw_Memory_Stats()
{
   Zoltan_Memory_Stats();
}

/*****************************************************************************/
int Zfw_Set_Fn(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
               void *data)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;
   switch(*type) {
   case ZOLTAN_PARTITION_MULTI_FN_TYPE:
      lb->Get_Partition_Multi_Fort = (ZOLTAN_PARTITION_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Partition_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_PARTITION_FN_TYPE:
      lb->Get_Partition_Fort = (ZOLTAN_PARTITION_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Partition_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_EDGES_MULTI_FN_TYPE:
      lb->Get_Num_Edges_Multi_Fort = (ZOLTAN_NUM_EDGES_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Edges_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_EDGE_LIST_MULTI_FN_TYPE:
      lb->Get_Edge_List_Multi_Fort = (ZOLTAN_EDGE_LIST_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Edge_List_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_EDGES_FN_TYPE:
      lb->Get_Num_Edges_Fort = (ZOLTAN_NUM_EDGES_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Edges_Fort_Wrapper, data);
      break;
   case ZOLTAN_EDGE_LIST_FN_TYPE:
      lb->Get_Edge_List_Fort = (ZOLTAN_EDGE_LIST_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Edge_List_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_GEOM_FN_TYPE:
      lb->Get_Num_Geom_Fort = (ZOLTAN_NUM_GEOM_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Geom_Fort_Wrapper, data);
      break;
   case ZOLTAN_GEOM_MULTI_FN_TYPE:
      lb->Get_Geom_Multi_Fort = (ZOLTAN_GEOM_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Geom_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_GEOM_FN_TYPE:
      lb->Get_Geom_Fort = (ZOLTAN_GEOM_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Geom_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_OBJ_FN_TYPE:
      lb->Get_Num_Obj_Fort = (ZOLTAN_NUM_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_OBJ_LIST_FN_TYPE:
      lb->Get_Obj_List_Fort = (ZOLTAN_OBJ_LIST_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Obj_List_Fort_Wrapper, data);
      break;
   case ZOLTAN_FIRST_OBJ_FN_TYPE:
      lb->Get_First_Obj_Fort = (ZOLTAN_FIRST_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_First_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_NEXT_OBJ_FN_TYPE:
      lb->Get_Next_Obj_Fort = (ZOLTAN_NEXT_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Next_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_BORDER_OBJ_FN_TYPE:
      lb->Get_Num_Border_Obj_Fort = (ZOLTAN_NUM_BORDER_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Border_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_BORDER_OBJ_LIST_FN_TYPE:
      lb->Get_Border_Obj_List_Fort = (ZOLTAN_BORDER_OBJ_LIST_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Border_Obj_List_Fort_Wrapper, data);
      break;
   case ZOLTAN_FIRST_BORDER_OBJ_FN_TYPE:
      lb->Get_First_Border_Obj_Fort = (ZOLTAN_FIRST_BORDER_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_First_Border_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_NEXT_BORDER_OBJ_FN_TYPE:
      lb->Get_Next_Border_Obj_Fort = (ZOLTAN_NEXT_BORDER_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Next_Border_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_PRE_MIGRATE_PP_FN_TYPE:
      lb->Migrate.Pre_Migrate_PP_Fort = (ZOLTAN_PRE_MIGRATE_PP_FORT_FN *)fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Pre_Migrate_PP_Fort_Wrapper, data);
      break;
   case ZOLTAN_MID_MIGRATE_PP_FN_TYPE:
      lb->Migrate.Mid_Migrate_PP_Fort = (ZOLTAN_MID_MIGRATE_PP_FORT_FN *)fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Mid_Migrate_PP_Fort_Wrapper, data);
      break;
   case ZOLTAN_POST_MIGRATE_PP_FN_TYPE:
      lb->Migrate.Post_Migrate_PP_Fort =(ZOLTAN_POST_MIGRATE_PP_FORT_FN*)fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Post_Migrate_PP_Fort_Wrapper, data);
      break;
   case ZOLTAN_PRE_MIGRATE_FN_TYPE:
      lb->Migrate.Pre_Migrate_Fort = (ZOLTAN_PRE_MIGRATE_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Pre_Migrate_Fort_Wrapper, data);
      break;
   case ZOLTAN_MID_MIGRATE_FN_TYPE:
      lb->Migrate.Mid_Migrate_Fort = (ZOLTAN_MID_MIGRATE_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Mid_Migrate_Fort_Wrapper, data);
      break;
   case ZOLTAN_POST_MIGRATE_FN_TYPE:
      lb->Migrate.Post_Migrate_Fort = (ZOLTAN_POST_MIGRATE_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Post_Migrate_Fort_Wrapper, data);
      break;
   case ZOLTAN_OBJ_SIZE_FN_TYPE:
      lb->Get_Obj_Size_Fort = (ZOLTAN_OBJ_SIZE_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Obj_Size_Fort_Wrapper, data);
      break;
   case ZOLTAN_PACK_OBJ_FN_TYPE:
      lb->Pack_Obj_Fort = (ZOLTAN_PACK_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Pack_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_UNPACK_OBJ_FN_TYPE:
      lb->Unpack_Obj_Fort = (ZOLTAN_UNPACK_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Unpack_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE:
      lb->Get_Obj_Size_Multi_Fort = (ZOLTAN_OBJ_SIZE_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Obj_Size_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_PACK_OBJ_MULTI_FN_TYPE:
      lb->Pack_Obj_Multi_Fort = (ZOLTAN_PACK_OBJ_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Pack_Obj_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE:
      lb->Unpack_Obj_Multi_Fort = (ZOLTAN_UNPACK_OBJ_MULTI_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Unpack_Obj_Multi_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_COARSE_OBJ_FN_TYPE:
      lb->Get_Num_Coarse_Obj_Fort = (ZOLTAN_NUM_COARSE_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Coarse_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_COARSE_OBJ_LIST_FN_TYPE:
      lb->Get_Coarse_Obj_List_Fort = (ZOLTAN_COARSE_OBJ_LIST_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Coarse_Obj_List_Fort_Wrapper, data);
      break;
   case ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE:
      lb->Get_First_Coarse_Obj_Fort = (ZOLTAN_FIRST_COARSE_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_First_Coarse_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE:
      lb->Get_Next_Coarse_Obj_Fort = (ZOLTAN_NEXT_COARSE_OBJ_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Next_Coarse_Obj_Fort_Wrapper, data);
      break;
   case ZOLTAN_NUM_CHILD_FN_TYPE:
      lb->Get_Num_Child_Fort = (ZOLTAN_NUM_CHILD_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Num_Child_Fort_Wrapper, data);
      break;
   case ZOLTAN_CHILD_LIST_FN_TYPE:
      lb->Get_Child_List_Fort = (ZOLTAN_CHILD_LIST_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Child_List_Fort_Wrapper, data);
      break;
   case ZOLTAN_CHILD_WEIGHT_FN_TYPE:
      lb->Get_Child_Weight_Fort = (ZOLTAN_CHILD_WEIGHT_FORT_FN *) fn;
      return Zoltan_Set_Fn(lb, *type, 
               (void (*)())Zoltan_Child_Weight_Fort_Wrapper, data);
      break;

   default:
      return Zoltan_Set_Fn(lb, *type, (void (*)())NULL, data);
      break;
   }
}

/*****************************************************************************/
int Zfw_Set_Fn0f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)())
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)NULL);
}

/*****************************************************************************/
int Zfw_Set_Fn1f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                  int *data)
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn2f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 float *data)
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn3f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 double *data)
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn4f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_1) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn5f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_2) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn6f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_3) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn7f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_4) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn8f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_1) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn9f(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_2) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_FnAf(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_3) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_FnBf(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_4) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn0s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)())
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)NULL);
}

/*****************************************************************************/
int Zfw_Set_Fn1s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                  int *data)
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn2s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                  float *data)
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn3s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                  double *data)
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn4s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_1) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn5s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_2) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn6s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_3) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn7s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(Zoltan_User_Data_4) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn8s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_1) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Fn9s(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_2) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_FnAs(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_3) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_FnBs(int *addr_lb, int *nbytes, ZOLTAN_FN_TYPE *type, void (*fn)(),
                 void *data)
/* data is type(LB_User_Data_4) */
{
   return Zfw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

/*****************************************************************************/
int Zfw_Set_Param(int *addr_lb, int *nbytes, int *int_param_name,
                   int *param_name_len, int *int_new_value, int *new_value_len)
{
   struct Zoltan_Struct *lb;
   char *param_name, *new_value;
   unsigned char *p;
   int i, result;
   param_name = (char *)ZOLTAN_MALLOC(*param_name_len+1);
   new_value = (char *)ZOLTAN_MALLOC(*new_value_len+1);
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;
   for (i=0; i<(*param_name_len); i++) param_name[i] = (char)int_param_name[i];
   param_name[*param_name_len] = '\0';
   for (i=0; i<(*new_value_len); i++) new_value[i] = (char)int_new_value[i];
   new_value[*new_value_len] = '\0';
   result = Zoltan_Set_Param(lb, param_name, new_value);
   ZOLTAN_FREE(&param_name);
   ZOLTAN_FREE(&new_value);
   return result;
}

/*****************************************************************************/
int Zfw_Set_Param_Vec(int *addr_lb, int *nbytes, int *int_param_name,
                   int *param_name_len, int *int_new_value, int *new_value_len,
                   int index)
{
   struct Zoltan_Struct *lb;
   char *param_name, *new_value;
   unsigned char *p;
   int i, result;
   param_name = (char *)ZOLTAN_MALLOC(*param_name_len+1);
   new_value = (char *)ZOLTAN_MALLOC(*new_value_len+1);
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;
   for (i=0; i<(*param_name_len); i++) param_name[i] = (char)int_param_name[i];
   param_name[*param_name_len] = '\0';
   for (i=0; i<(*new_value_len); i++) new_value[i] = (char)int_new_value[i];
   new_value[*new_value_len] = '\0';
   result = Zoltan_Set_Param_Vec(lb, param_name, new_value, index);
   ZOLTAN_FREE(&param_name);
   ZOLTAN_FREE(&new_value);
   return result;
}

/*****************************************************************************/
int Zfw_LB_Partition(int *addr_lb, int *nbytes, int *changes, 
  int *num_gid_entries, int *num_lid_entries,
  int *num_import,
  ZOLTAN_ID_PTR *import_global_ids, ZOLTAN_ID_PTR *import_local_ids,
  int **import_procs, int **import_to_part, int *num_export,
  ZOLTAN_ID_PTR *export_global_ids, ZOLTAN_ID_PTR *export_local_ids,
  int **export_procs, int **export_to_part
#ifdef PGI
/* PGI uses hidden arguments when it passes pointers */
   ,int *imp_gid_hide, int *imp_lid_hide, int *imp_proc_hide,
    int *imp_to_part_hide,
    int *exp_gid_hide, int *exp_lid_hide, int *exp_proc_hide,
    int *exp_to_part_hide
#endif
#ifdef FUJITSU
/* Fujitsu and Lahey use a hidden argument for every argument */
/* TEMP need to verify this with Fujitsu or Lahey */
   ,int *addr_lb_hide, int *nbytes_hide, int *changes_hide,
    int *num_gid_entries_hide, int *num_lid_entries_hide,
    int *num_import_hide, int *imp_gid_hide, int *imp_lid_hide,
    int *imp_proc_hide, int *imp_to_part_hide,
    int *num_export_hide, int *exp_gid_hide,
    int *exp_lid_hide, int *exp_proc_hide, int *exp_to_part_hide
#endif
)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
#if defined (PGI) || defined (FUJITSU)
#define F90LB_TEMP 3
#else
#define F90LB_TEMP 2
#endif
   ZOLTAN_ID_PTR temp_imp_gid[F90LB_TEMP], temp_exp_gid[F90LB_TEMP];
   ZOLTAN_ID_PTR temp_imp_lid[F90LB_TEMP], temp_exp_lid[F90LB_TEMP];
   int *temp_imp_proc[F90LB_TEMP], *temp_exp_proc[F90LB_TEMP];
   int *temp_imp_to_part[F90LB_TEMP], *temp_exp_to_part[F90LB_TEMP];
#undef F90LB_TEMP

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */

   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;

/* put the address of the Fortran pointer into temp_*[1] to be passed to
   Fortran for allocation.  The address of the allocated space will be
   in temp_*[0] so it can be used by C without messing up the Fortran pointer*/

   temp_imp_gid[1] = (ZOLTAN_ID_PTR)import_global_ids;
   temp_imp_lid[1] = (ZOLTAN_ID_PTR)import_local_ids;
   temp_imp_proc[1] = (int *)import_procs;
   temp_imp_to_part[1] = (int *)import_to_part;
   temp_exp_gid[1] = (ZOLTAN_ID_PTR)export_global_ids;
   temp_exp_lid[1] = (ZOLTAN_ID_PTR)export_local_ids;
   temp_exp_proc[1] = (int *)export_procs;
   temp_exp_to_part[1] = (int *)export_to_part;

/* for PGI and FUJITSU, put the hidden argument in temp_*[2] */

#if defined (PGI) || defined (FUJITSU)
   temp_imp_gid[2] = (ZOLTAN_ID_PTR)imp_gid_hide;
   temp_imp_lid[2] = (ZOLTAN_ID_PTR)imp_lid_hide;
   temp_imp_proc[2] = (int *)imp_proc_hide;
   temp_imp_to_part[2] = (int *)imp_to_part_hide;
   temp_exp_gid[2] = (ZOLTAN_ID_PTR)exp_gid_hide;
   temp_exp_lid[2] = (ZOLTAN_ID_PTR)exp_lid_hide;
   temp_exp_proc[2] = (int *)exp_proc_hide;
   temp_exp_to_part[2] = (int *)exp_to_part_hide;
#endif

/* call Zoltan_LB_Partition */

   return Zoltan_LB_Partition(lb, changes, num_gid_entries, num_lid_entries, 
                     num_import, temp_imp_gid, temp_imp_lid,
                     temp_imp_proc, temp_imp_to_part,
                     num_export, temp_exp_gid, temp_exp_lid,
                     temp_exp_proc, temp_exp_to_part);
}

/*****************************************************************************/
int Zfw_LB_Eval(int *addr_lb, int *nbytes, int *print_stats,
                int *nobj, float *obj_wgt, int *ncuts, float *cut_wgt,
                int *nboundary, int *nadj,
                int *is_nobj, int *is_obj_wgt, int *is_ncuts, int *is_cut_wgt,
                int *is_nboundary, int *is_nadj)
{
   struct Zoltan_Struct *lb;
   int *loc_nobj, *loc_ncuts, *loc_nboundary, *loc_nadj;
   float *loc_obj_wgt, *loc_cut_wgt;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;
   if (*is_nobj) {loc_nobj = nobj;} else {loc_nobj = NULL;}
   if (*is_ncuts) {loc_ncuts = ncuts;} else {loc_ncuts = NULL;}
   if (*is_obj_wgt) {loc_obj_wgt = obj_wgt;} else {loc_obj_wgt = NULL;}
   if (*is_cut_wgt) {loc_cut_wgt = cut_wgt;} else {loc_cut_wgt = NULL;}
   if (*is_nboundary) {loc_nboundary = nboundary;} else {loc_nboundary = NULL;}
   if (*is_nadj) {loc_nadj = nadj;} else {loc_nadj = NULL;}

   return  Zoltan_LB_Eval(lb, *print_stats, loc_nobj, loc_obj_wgt, loc_ncuts, loc_cut_wgt,
           loc_nboundary, loc_nadj);
}

/*****************************************************************************/
int Zfw_LB_Set_Part_Sizes(int *addr_lb, int *nbytes, int *global_part, int *len,
                          int *partids, int *wgtidx, float *partsizes)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;

   return Zoltan_LB_Set_Part_Sizes(lb, *global_part, *len, partids, wgtidx,
                                   partsizes);
}

/*****************************************************************************/
int Zfw_LB_Point_Assign(int *addr_lb, int *nbytes, double *coords, int *proc)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;

   return Zoltan_LB_Point_Assign(lb, coords, proc);
}

/*****************************************************************************/
int Zfw_LB_Point_PP_Assign(int *addr_lb, int *nbytes, double *coords, int *proc,
                           int *part)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;

   return Zoltan_LB_Point_PP_Assign(lb, coords, proc, part);
}

/*****************************************************************************/
int Zfw_LB_Box_Assign(int *addr_lb, int *nbytes, double *xmin, double *ymin,
                     double *zmin, double *xmax, double *ymax, double *zmax,
                     int *procs, int *numprocs)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;

   return Zoltan_LB_Box_Assign(lb, *xmin, *ymin, *zmin, *xmax, *ymax, *zmax, 
                               procs, numprocs);
}

/*****************************************************************************/
int Zfw_LB_Box_PP_Assign(int *addr_lb, int *nbytes, double *xmin, double *ymin,
                     double *zmin, double *xmax, double *ymax, double *zmax,
                     int *procs, int *numprocs, int *parts, int *numparts)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;

   return Zoltan_LB_Box_PP_Assign(lb, *xmin, *ymin, *zmin, *xmax, *ymax, *zmax,
                                  procs, numprocs, parts, numparts);
}

/*****************************************************************************/
int Zfw_Invert_Lists(int *addr_lb, int *nbytes, 
  int *num_gid_entries, int *num_lid_entries, int *num_input,
  ZOLTAN_ID_PTR input_global_ids, ZOLTAN_ID_PTR input_local_ids,
  int *input_procs, int *input_to_part, int *num_output,
  ZOLTAN_ID_PTR *output_global_ids, ZOLTAN_ID_PTR *output_local_ids,
  int **output_procs, int **output_to_part
#ifdef PGI
  ,int *output_gid_hide, int *output_lid_hide, int *output_proc_hide, 
   int *output_to_part_hide
#endif
#ifdef FUJITSU
 ,int *addr_lb_hide, int *nbytes_hide,
  int *num_gid_entries_hide, int *num_lid_entries_hide,
  int *num_input_hide,
  int *input_global_ids_hide, int *input_local_ids_hide,
  int *input_procs_hide, int *input_to_part_hide,
  int *num_output_hide,
  int *output_gid_hide, int *output_lid_hide, 
  int *output_proc_hide, int *output_to_part_hide
#endif
)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
#if defined (PGI) || defined(FUJITSU)
#define F90LB_TEMP 3
#else
#define F90LB_TEMP 2
#endif
   ZOLTAN_ID_PTR temp_output_gid[F90LB_TEMP];
   ZOLTAN_ID_PTR temp_output_lid[F90LB_TEMP];
   int *temp_output_proc[F90LB_TEMP];
   int *temp_output_to_part[F90LB_TEMP];
#undef F90LB_TEMP

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */

   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;

/* put the address of the Fortran pointer into temp_*[1] to be passed to
   Fortran for allocation.  The address of the allocated space will be
   in temp_*[0] so it can be used by C without messing up the Fortran pointer*/

   temp_output_gid[1] = (ZOLTAN_ID_PTR)output_global_ids;
   temp_output_lid[1] = (ZOLTAN_ID_PTR)output_local_ids;
   temp_output_proc[1] = (int *)output_procs;
   temp_output_to_part[1] = (int *)output_to_part;

/* for PGI and FUJITSU, put the hidden argument in temp_*[2] */

#if defined (PGI) || defined(FUJITSU)
   temp_output_gid[2] = (ZOLTAN_ID_PTR)output_gid_hide;
   temp_output_lid[2] = (ZOLTAN_ID_PTR)output_lid_hide;
   temp_output_proc[2] = (int *)output_proc_hide;
   temp_output_to_part[2] = (int *)output_to_part_hide;
#endif

/* call Zoltan_Invert_Lists */

   return Zoltan_Invert_Lists(lb, 
                     *num_input, input_global_ids,
                     input_local_ids, input_procs, input_to_part,
                     num_output, temp_output_gid, temp_output_lid,
                     temp_output_proc, temp_output_to_part);
}

/*****************************************************************************/
int Zfw_Compute_Destinations(int *addr_lb, int *nbytes, 
  int *num_gid_entries, int *num_lid_entries, int *num_input,
  ZOLTAN_ID_PTR input_global_ids, ZOLTAN_ID_PTR input_local_ids,
  int *input_procs, int *num_output,
  ZOLTAN_ID_PTR *output_global_ids, ZOLTAN_ID_PTR *output_local_ids,
  int **output_procs
#ifdef PGI
  ,int *output_gid_hide, int *output_lid_hide, int *output_proc_hide
#endif
#ifdef FUJITSU
 ,int *addr_lb_hide, int *nbytes_hide,
  int *num_gid_entries_hide, int *num_lid_entries_hide,
  int *num_input_hide,
  int *input_global_ids_hide, int *input_local_ids_hide,
  int *input_procs_hide, int *num_output_hide,
  int *output_gid_hide, int *output_lid_hide, int *output_proc_hide
#endif
)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
#if defined (PGI) || defined(FUJITSU)
#define F90LB_TEMP 3
#else
#define F90LB_TEMP 2
#endif
   ZOLTAN_ID_PTR temp_output_gid[F90LB_TEMP];
   ZOLTAN_ID_PTR temp_output_lid[F90LB_TEMP];
   int *temp_output_proc[F90LB_TEMP];
#undef F90LB_TEMP

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */

   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;

/* put the address of the Fortran pointer into temp_*[1] to be passed to
   Fortran for allocation.  The address of the allocated space will be
   in temp_*[0] so it can be used by C without messing up the Fortran pointer*/

   temp_output_gid[1] = (ZOLTAN_ID_PTR)output_global_ids;
   temp_output_lid[1] = (ZOLTAN_ID_PTR)output_local_ids;
   temp_output_proc[1] = (int *)output_procs;

/* for PGI and FUJITSU, put the hidden argument in temp_*[2] */

#if defined (PGI) || defined(FUJITSU)
   temp_output_gid[2] = (ZOLTAN_ID_PTR)output_gid_hide;
   temp_output_lid[2] = (ZOLTAN_ID_PTR)output_lid_hide;
   temp_output_proc[2] = (int *)output_proc_hide;
#endif

/* call Zoltan_Compute_Destinations */

   return Zoltan_Compute_Destinations(lb, 
                     *num_input, input_global_ids,
                     input_local_ids, input_procs, 
                     num_output, temp_output_gid, temp_output_lid,
                     temp_output_proc);
}


/*****************************************************************************/
int Zfw_Migrate(int *addr_lb, int *nbytes, 
 int *num_import,
 ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
 int *import_procs, int *import_to_part, int *num_export,
 ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids,
 int *export_procs, int *export_to_part)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;
   return Zoltan_Migrate(lb,
                         *num_import,import_global_ids,import_local_ids,
                         import_procs,import_to_part,
                         *num_export,export_global_ids,
                         export_local_ids,export_procs,export_to_part);
}

/*****************************************************************************/
int Zfw_Help_Migrate(int *addr_lb, int *nbytes, 
 int *num_import,
 ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
 int *import_procs, int *num_export,
 ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids,
 int *export_procs)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;
   return Zoltan_Help_Migrate(lb,
                          *num_import,import_global_ids,import_local_ids,
                          import_procs,*num_export,export_global_ids,
                          export_local_ids,export_procs);
}

/*****************************************************************************/
int Zfw_Order(
 int *addr_lb, int *nbytes,
 int *num_gid_entries, int *num_lid_entries,
 int *num_obj,
 ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
 int *rank, int *iperm)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   int ierr;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;
   ierr = Zoltan_Order(lb,num_gid_entries,num_lid_entries,*num_obj,
                       gids, lids, rank, iperm, NULL);
   return ierr;
}

/*****************************************************************************/
int Zfw_Generate_Files(int *addr_lb, int *nbytes, int *int_filename,
                   int *filename_len, int *base_index, int *gen_geom,
                   int *gen_graph, int *gen_hg)
{
   struct Zoltan_Struct *lb;
   char *filename;
   unsigned char *p;
   int i, result;
   filename = (char *)ZOLTAN_MALLOC(*filename_len+1);
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;
   for (i=0; i<(*filename_len); i++) filename[i] = (char)int_filename[i];
   filename[*filename_len] = '\0';
   result = Zoltan_Generate_Files(lb, filename, *base_index, *gen_geom,
                                  *gen_graph, *gen_hg);
   ZOLTAN_FREE(&filename);
   return result;
}

/*****************************************************************************/
int Zfw_RCB_Box(int *addr_lb, int *nbytes, int *part, int *ndim,
                double *xmin, double *ymin, double *zmin, 
                double *xmax, double *ymax, double *zmax)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Current = lb;

   return Zoltan_RCB_Box(lb, *part, ndim, xmin, ymin, zmin, xmax, ymax, zmax);
}
/*****************************************************************************/
void Zfw_Register_Fort_Malloc(ZOLTAN_FORT_MALLOC_INT_FN *fort_malloc_int,
                                    ZOLTAN_FORT_FREE_INT_FN *fort_free_int)
{
   Zoltan_Register_Fort_Malloc(fort_malloc_int,fort_free_int);
}

/*****************************************************************************/
/* TEMP child_order */
void Zfw_Reftree_Get_Child_Order(
  int *addr_lb, 
  int *nbytes, 
  int *order, 
  int *ierr)
{
   struct Zoltan_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   Zoltan_Reftree_Get_Child_Order(lb,order,ierr);
}
/* end TEMP child_order */

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
