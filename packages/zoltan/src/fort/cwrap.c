/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lb_const.h"
#include "all_allo_const.h"
#include "sppr_header"

/*--------------------------------------------------------------------*/
/* procedure name mangling                                            */

#define LOWERCASE   1
#define UPPERCASE   2
#define UNDERSCORE  3
#define UNDERSCORE2 4

#if FMANGLE==LOWERCASE

#define LB_fw_Initialize                 lb_fw_initialize
#define LB_fw_Initialize1                lb_fw_initialize1
#define LB_fw_Create                     lb_fw_create       
#define LB_fw_Destroy                    lb_fw_destroy       
#define LB_fw_Memory_Stats               lb_fw_Memory_Stats       
#define LB_fw_Set_Fn0f                   lb_fw_set_fn0f
#define LB_fw_Set_Fn1f                   lb_fw_set_fn1f
#define LB_fw_Set_Fn2f                   lb_fw_set_fn2f
#define LB_fw_Set_Fn3f                   lb_fw_set_fn3f
#define LB_fw_Set_Fn4f                   lb_fw_set_fn4f
#define LB_fw_Set_Fn5f                   lb_fw_set_fn5f
#define LB_fw_Set_Fn6f                   lb_fw_set_fn6f
#define LB_fw_Set_Fn7f                   lb_fw_set_fn7f
#define LB_fw_Set_Fn0s                   lb_fw_set_fn0s
#define LB_fw_Set_Fn1s                   lb_fw_set_fn1s
#define LB_fw_Set_Fn2s                   lb_fw_set_fn2s
#define LB_fw_Set_Fn3s                   lb_fw_set_fn3s
#define LB_fw_Set_Fn4s                   lb_fw_set_fn4s
#define LB_fw_Set_Fn5s                   lb_fw_set_fn5s
#define LB_fw_Set_Fn6s                   lb_fw_set_fn6s
#define LB_fw_Set_Fn7s                   lb_fw_set_fn7s
#define LB_fw_Set_Method                 lb_fw_set_method
#define LB_fw_Set_Param                  lb_fw_set_param
#define LB_fw_Balance                    lb_fw_balance
#define LB_fw_Eval                       lb_fw_eval
#define LB_fw_Point_Assign               lb_fw_point_assign
#define LB_fw_Box_Assign                 lb_fw_box_assign
#define LB_fw_Compute_Destinations       lb_fw_compute_destinations
#define LB_fw_Help_Migrate               lb_fw_help_migrate  
#define LB_fw_Register_Fort_Malloc       lb_fw_register_fort_malloc
#define LB_fw_Get_Address_int            lb_fw_get_address_int
#define LB_fw_Get_Wgt_Dim                lb_fw_get_wgt_dim
#define LB_fw_Get_Comm_Dim               lb_fw_get_comm_dim
/* TEMP child_order */
#define LB_fw_Get_Child_Order            lb_fw_get_child_order

#elif FMANGLE==UPPERCASE

#define LB_fw_Initialize                 LB_FW_INITIALIZE
#define LB_fw_Initialize1                LB_FW_INITIALIZE1
#define LB_fw_Create                     LB_FW_CREATE       
#define LB_fw_Destroy                    LB_FW_DESTROY       
#define LB_fw_Memory_Stats               LB_FW_MEMORY_STATS  
#define LB_fw_Set_Fn0f                   LB_FW_SET_FN0F
#define LB_fw_Set_Fn1f                   LB_FW_SET_FN1F
#define LB_fw_Set_Fn2f                   LB_FW_SET_FN2F
#define LB_fw_Set_Fn3f                   LB_FW_SET_FN3F
#define LB_fw_Set_Fn4f                   LB_FW_SET_FN4F
#define LB_fw_Set_Fn5f                   LB_FW_SET_FN5F
#define LB_fw_Set_Fn6f                   LB_FW_SET_FN6F
#define LB_fw_Set_Fn7f                   LB_FW_SET_FN7F
#define LB_fw_Set_Fn0s                   LB_FW_SET_FN0S
#define LB_fw_Set_Fn1s                   LB_FW_SET_FN1S
#define LB_fw_Set_Fn2s                   LB_FW_SET_FN2S
#define LB_fw_Set_Fn3s                   LB_FW_SET_FN3S
#define LB_fw_Set_Fn4s                   LB_FW_SET_FN4S
#define LB_fw_Set_Fn5s                   LB_FW_SET_FN5S
#define LB_fw_Set_Fn6s                   LB_FW_SET_FN6S
#define LB_fw_Set_Fn7s                   LB_FW_SET_FN7S
#define LB_fw_Set_Method                 LB_FW_SET_METHOD
#define LB_fw_Set_Param                  LB_FW_SET_PARAM
#define LB_fw_Balance                    LB_FW_BALANCE
#define LB_fw_Eval                       LB_FW_EVAL
#define LB_fw_Point_Assign               LB_FW_POINT_ASSIGN
#define LB_fw_Box_Assign                 LB_FW_BOX_ASSIGN
#define LB_fw_Compute_Destinations       LB_FW_COMPUTE_DESTINATIONS  
#define LB_fw_Help_Migrate               LB_FW_HELP_MIGRATE  
#define LB_fw_Register_Fort_Malloc       LB_FW_REGISTER_FORT_MALLOC
#define LB_fw_Get_Address_int            LB_FW_GET_ADDRESS_INT
#define LB_fw_Get_Comm_Dim               LB_FW_GET_COMM_DIM
/* TEMP child_order */
#define LB_fw_Get_Child_Order            LB_FW_GET_CHILD_ORDER

#elif FMANGLE==UNDERSCORE

#define LB_fw_Initialize                 lb_fw_initialize_
#define LB_fw_Initialize1                lb_fw_initialize1_
#define LB_fw_Create                     lb_fw_create_
#define LB_fw_Destroy                    lb_fw_destroy_
#define LB_fw_Memory_Stats               lb_fw_memory_stats_
#define LB_fw_Set_Fn0f                   lb_fw_set_fn0f_
#define LB_fw_Set_Fn1f                   lb_fw_set_fn1f_
#define LB_fw_Set_Fn2f                   lb_fw_set_fn2f_
#define LB_fw_Set_Fn3f                   lb_fw_set_fn3f_
#define LB_fw_Set_Fn4f                   lb_fw_set_fn4f_
#define LB_fw_Set_Fn5f                   lb_fw_set_fn5f_
#define LB_fw_Set_Fn6f                   lb_fw_set_fn6f_
#define LB_fw_Set_Fn7f                   lb_fw_set_fn7f_
#define LB_fw_Set_Fn0s                   lb_fw_set_fn0s_
#define LB_fw_Set_Fn1s                   lb_fw_set_fn1s_
#define LB_fw_Set_Fn2s                   lb_fw_set_fn2s_
#define LB_fw_Set_Fn3s                   lb_fw_set_fn3s_
#define LB_fw_Set_Fn4s                   lb_fw_set_fn4s_
#define LB_fw_Set_Fn5s                   lb_fw_set_fn5s_
#define LB_fw_Set_Fn6s                   lb_fw_set_fn6s_
#define LB_fw_Set_Fn7s                   lb_fw_set_fn7s_
#define LB_fw_Set_Method                 lb_fw_set_method_
#define LB_fw_Set_Param                  lb_fw_set_param_
#define LB_fw_Balance                    lb_fw_balance_
#define LB_fw_Eval                       lb_fw_eval_
#define LB_fw_Point_Assign               lb_fw_point_assign_
#define LB_fw_Box_Assign                 lb_fw_box_assign_
#define LB_fw_Compute_Destinations       lb_fw_compute_destinations_
#define LB_fw_Help_Migrate               lb_fw_help_migrate_
#define LB_fw_Register_Fort_Malloc       lb_fw_register_fort_malloc_
#define LB_fw_Get_Address_int            lb_fw_get_address_int_
#define LB_fw_Get_Wgt_Dim                lb_fw_get_wgt_dim_
#define LB_fw_Get_Comm_Dim               lb_fw_get_comm_dim_
/* TEMP child_order */
#define LB_fw_Get_Child_Order            lb_fw_get_child_order_

#elif FMANGLE==UNDERSCORE2

#define LB_fw_Initialize                 lb_fw_initialize__
#define LB_fw_Initialize1                lb_fw_initialize1__
#define LB_fw_Create                     lb_fw_create__
#define LB_fw_Destroy                    lb_fw_destroy__
#define LB_fw_Memory_Stats               lb_fw_memory_stats__
#define LB_fw_Set_Fn0f                   lb_fw_set_fn0f__
#define LB_fw_Set_Fn1f                   lb_fw_set_fn1f__
#define LB_fw_Set_Fn2f                   lb_fw_set_fn2f__
#define LB_fw_Set_Fn3f                   lb_fw_set_fn3f__
#define LB_fw_Set_Fn4f                   lb_fw_set_fn4f__
#define LB_fw_Set_Fn5f                   lb_fw_set_fn5f__
#define LB_fw_Set_Fn6f                   lb_fw_set_fn6f__
#define LB_fw_Set_Fn7f                   lb_fw_set_fn7f__
#define LB_fw_Set_Fn0s                   lb_fw_set_fn0s__
#define LB_fw_Set_Fn1s                   lb_fw_set_fn1s__
#define LB_fw_Set_Fn2s                   lb_fw_set_fn2s__
#define LB_fw_Set_Fn3s                   lb_fw_set_fn3s__
#define LB_fw_Set_Fn4s                   lb_fw_set_fn4s__
#define LB_fw_Set_Fn5s                   lb_fw_set_fn5s__
#define LB_fw_Set_Fn6s                   lb_fw_set_fn6s__
#define LB_fw_Set_Fn7s                   lb_fw_set_fn7s__
#define LB_fw_Set_Method                 lb_fw_set_method__
#define LB_fw_Set_Param                  lb_fw_set_param__
#define LB_fw_Balance                    lb_fw_balance__
#define LB_fw_Eval                       lb_fw_eval__
#define LB_fw_Point_Assign               lb_fw_point_assign__
#define LB_fw_Box_Assign                 lb_fw_box_assign__
#define LB_fw_Compute_Destinations       lb_fw_compute_destinations__
#define LB_fw_Help_Migrate               lb_fw_help_migrate__
#define LB_fw_Register_Fort_Malloc       lb_fw_register_fort_malloc__
#define LB_fw_Get_Address_int            lb_fw_get_address_int__
#define LB_fw_Get_Wgt_Dim                lb_fw_get_wgt_dim__
#define LB_fw_Get_Comm_Dim               lb_fw_get_comm_dim__
/* TEMP child_order */
#define LB_fw_Get_Child_Order            lb_fw_get_child_order__

#endif /* FMANGLE */

/*--------------------------------------------------------------------*/
/* Variables                                                          */

static struct LB_Struct *LB_Current_lb;

/*--------------------------------------------------------------------*/
/* Utilities                                                          */

/* some MPI implementations may require conversion between a Fortran
   communicator and a C communicator.  This routine is used to perform the
   conversion.  It may need different forms for different MPI libraries. */

MPI_Comm LB_comm_f2c(int *f_comm)
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

/* These routines get the address of an array allocated by fortran and
   return it */
#ifdef PTR_64BIT
void LB_fw_Get_Address_int(int *addr,
                           long *ret_addr)
{
   if (sizeof(long) != sizeof(int *)) {
     LB_PRINT_ERROR(-1, "LB_fw_Get_Address_int", 
       "sizeof(long) != sizeof(int *); F90 allocation will not work properly.");
   }
   *ret_addr = (long)addr;
}
#else
void LB_fw_Get_Address_int(int *addr,
                           int *ret_addr)
{
   if (sizeof(int) != sizeof(int *)) {
     LB_PRINT_ERROR(-1, "LB_fw_Get_Address_int", 
       "sizeof(int) != sizeof(int *); F90 allocation will not work properly.");
   }
   *ret_addr = (int)addr;
}
#endif  /* PTR_64BIT */

int LB_fw_Get_Wgt_Dim(int *addr_lb, int *nbytes)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   return lb->Obj_Weight_Dim;
}

int LB_fw_Get_Comm_Dim(int *addr_lb, int *nbytes)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   return lb->Comm_Weight_Dim;
}

/*--------------------------------------------------------------------*/
/* Reverse wrappers for callbacks                                     */

int LB_Num_Edges_Fort_Wrapper(void *data, 
                              int num_gid_entries, int num_lid_entries,
                              LB_ID_PTR global_id, LB_ID_PTR local_id,
                              int *ierr)
{
   return LB_Current_lb->Get_Num_Edges_Fort(data,
                                            &num_gid_entries, &num_lid_entries,
                                            global_id, local_id, ierr);
}

void LB_Edge_List_Fort_Wrapper(void *data, 
                               int num_gid_entries, int num_lid_entries,
                               LB_ID_PTR global_id, LB_ID_PTR local_id,
                               LB_ID_PTR nbor_global_id, int *nbor_procs,
                               int wdim, float *nbor_ewgts, int *ierr)
{
   LB_Current_lb->Get_Edge_List_Fort(data, &num_gid_entries, &num_lid_entries,
                                     global_id, local_id,
                                     nbor_global_id, nbor_procs, &wdim,
                                     nbor_ewgts, ierr);
}

int LB_Num_Geom_Fort_Wrapper(void *data, int *ierr)
{
   return LB_Current_lb->Get_Num_Geom_Fort(data,ierr);
}

void LB_Geom_Fort_Wrapper(void *data, int num_gid_entries, int num_lid_entries,
                          LB_ID_PTR global_id, LB_ID_PTR local_id,
                          double *geom_vec, int *ierr)
{
   LB_Current_lb->Get_Geom_Fort(data, &num_gid_entries, &num_lid_entries,
                                global_id, local_id, geom_vec, ierr);
}

int LB_Num_Obj_Fort_Wrapper(void *data, int *ierr)
{
   return LB_Current_lb->Get_Num_Obj_Fort(data, ierr);
}

void LB_Obj_List_Fort_Wrapper(void *data,
                              int num_gid_entries, int num_lid_entries,
                              LB_ID_PTR global_ids, LB_ID_PTR local_ids,
                              int wdim, float *objwgts, int *ierr)
{
   LB_Current_lb->Get_Obj_List_Fort(data, &num_gid_entries, &num_lid_entries,
                                    global_ids, local_ids, &wdim,
                                    objwgts, ierr);
}

int LB_First_Obj_Fort_Wrapper(void *data, 
                              int num_gid_entries, int num_lid_entries, 
                              LB_ID_PTR first_global_id,
                              LB_ID_PTR first_local_id,
                              int wdim, float *first_obj_wgt, int *ierr)
{
   return LB_Current_lb->Get_First_Obj_Fort(data, 
                                            &num_gid_entries, &num_lid_entries,
                                            first_global_id,
                                            first_local_id, &wdim,
                                            first_obj_wgt, ierr);
}

int LB_Next_Obj_Fort_Wrapper(void *data, 
                             int num_gid_entries, int num_lid_entries, 
                             LB_ID_PTR global_id, LB_ID_PTR local_id,
                             LB_ID_PTR next_global_id, LB_ID_PTR next_local_id,
                             int wdim, float *next_obj_wgt, int *ierr)
{
   return LB_Current_lb->Get_Next_Obj_Fort(data, 
                                           &num_gid_entries, &num_lid_entries, 
                                           global_id, local_id,
                                           next_global_id, next_local_id,
                                           &wdim, next_obj_wgt, ierr);
}

int LB_Num_Border_Obj_Fort_Wrapper(void *data, int nbor_proc, int *ierr)
{
   return LB_Current_lb->Get_Num_Border_Obj_Fort(data, &nbor_proc, ierr);
}

void LB_Border_Obj_List_Fort_Wrapper(void *data, 
                                     int num_gid_entries, int num_lid_entries, 
                                     int nbor_proc,
                                     LB_ID_PTR global_ids, LB_ID_PTR local_ids,
                                     int wdim, float *objwgts, int *ierr)
{
   LB_Current_lb->Get_Border_Obj_List_Fort(data, 
                                           &num_gid_entries, &num_lid_entries, 
                                           &nbor_proc, global_ids,
                                           local_ids, &wdim, objwgts, ierr);
}

int LB_First_Border_Obj_Fort_Wrapper(void *data, 
                                     int num_gid_entries, int num_lid_entries,
                                     int nbor_proc,
                                     LB_ID_PTR first_global_id,
                                     LB_ID_PTR first_local_id,
                                     int wdim, float *first_obj_wgt,
                                     int *ierr)
{
   return LB_Current_lb->Get_First_Border_Obj_Fort(data, 
                                                   &num_gid_entries, 
                                                   &num_lid_entries, 
                                                   &nbor_proc,
                                                   first_global_id,
                                                   first_local_id, &wdim,
                                                   first_obj_wgt, ierr);
}

int LB_Next_Border_Obj_Fort_Wrapper(void *data, 
                                    int num_gid_entries, int num_lid_entries,
                                    LB_ID_PTR global_id,
                                    LB_ID_PTR local_id, int nbor_proc,
                                    LB_ID_PTR next_global_id,
                                    LB_ID_PTR next_local_id,
                                    int wdim, float *next_obj_wgt,
                                    int *ierr)
{
   return LB_Current_lb->Get_Next_Border_Obj_Fort(data, 
                                                  &num_gid_entries,
                                                  &num_lid_entries,
                                                  global_id, local_id,
                                                  &nbor_proc, next_global_id,
                                                  next_local_id, &wdim,
                                                  next_obj_wgt, ierr);
}

int LB_Obj_Size_Fort_Wrapper(void *data, int num_gid_entries,
            int num_lid_entries, LB_ID_PTR global_id, 
            LB_ID_PTR local_id, int *ierr)
{
   return LB_Current_lb->Migrate.Get_Obj_Size_Fort(data,
             &num_gid_entries, &num_lid_entries,
             global_id, local_id, ierr);
}

void LB_Pre_Migrate_Fort_Wrapper(void *data, 
                                 int num_gid_entries, int num_lid_entries,
                                 int num_import,
                                 LB_ID_PTR import_global_ids,
                                 LB_ID_PTR import_local_ids, int *import_procs,
                                 int num_export, LB_ID_PTR export_global_ids,
                                 LB_ID_PTR export_local_ids, int *export_procs,
                                 int *ierr)
{
   LB_Current_lb->Migrate.Pre_Migrate_Fort(data, 
                                           &num_gid_entries,
                                           &num_lid_entries,
                                           &num_import,
                                           import_global_ids,
                                           import_local_ids, import_procs,
                                           &num_export, export_global_ids,
                                           export_local_ids, export_procs,
                                           ierr);
}

void LB_Mid_Migrate_Fort_Wrapper(void *data, 
                                 int num_gid_entries, int num_lid_entries,
                                 int num_import,
                                 LB_ID_PTR import_global_ids,
                                 LB_ID_PTR import_local_ids, int *import_procs,
                                 int num_export, LB_ID_PTR export_global_ids,
                                 LB_ID_PTR export_local_ids, int *export_procs,
                                 int *ierr)
{
   LB_Current_lb->Migrate.Mid_Migrate_Fort(data,
                                           &num_gid_entries,
                                           &num_lid_entries,
                                           &num_import,
                                           import_global_ids,
                                           import_local_ids, import_procs,
                                           &num_export, export_global_ids,
                                           export_local_ids, export_procs,
                                           ierr);
}

void LB_Post_Migrate_Fort_Wrapper(void *data, 
                                  int num_gid_entries, int num_lid_entries,
                                  int num_import,
                                  LB_ID_PTR import_global_ids,
                                  LB_ID_PTR import_local_ids, int *import_procs,
                                  int num_export, LB_ID_PTR export_global_ids,
                                  LB_ID_PTR export_local_ids, int *export_procs,
                                  int *ierr)
{
   LB_Current_lb->Migrate.Post_Migrate_Fort(data,
                                            &num_gid_entries, &num_lid_entries, 
                                            &num_import,
                                            import_global_ids,
                                            import_local_ids, import_procs,
                                            &num_export, export_global_ids,
                                            export_local_ids, export_procs,
                                            ierr);
}

void LB_Pack_Obj_Fort_Wrapper(void *data, 
                            int num_gid_entries, int num_lid_entries,
                            LB_ID_PTR global_id, LB_ID_PTR local_id,
                            int dest_proc, int size, char *buf, int *ierr)
{
   LB_Current_lb->Migrate.Pack_Obj_Fort(data, 
                                        &num_gid_entries, &num_lid_entries, 
                                        global_id, local_id,
                                        &dest_proc, &size, buf, ierr);
}

void LB_Unpack_Obj_Fort_Wrapper(void *data, int num_gid_entries,
                                LB_ID_PTR global_id, int size,
                                char *buf, int *ierr)
{
   LB_Current_lb->Migrate.Unpack_Obj_Fort(data, &num_gid_entries, 
                                          global_id, &size, buf, ierr);
}

int LB_Num_Coarse_Obj_Fort_Wrapper(void *data, int *ierr)
{
   return LB_Current_lb->Get_Num_Coarse_Obj_Fort(data, ierr);
}

void LB_Coarse_Obj_List_Fort_Wrapper(void *data, 
                           int num_gid_entries, int num_lid_entries,
                           LB_ID_PTR global_ids,
                           LB_ID_PTR local_ids, int *assigned, int *num_vert,
                           int *vertices, int *in_order, int *in_vertex,
                           int *out_vertex, int *ierr)
{
   LB_Current_lb->Get_Coarse_Obj_List_Fort(data, 
                                           &num_gid_entries, &num_lid_entries,
                                           global_ids, local_ids,
                                           assigned, num_vert, vertices,
                                           in_order, in_vertex, out_vertex,
                                           ierr);
}

int LB_First_Coarse_Obj_Fort_Wrapper(void *data, 
                                     int num_gid_entries, int num_lid_entries, 
                                     LB_ID_PTR global_id,
                                     LB_ID_PTR local_id, int *assigned,
                                     int *num_vert, int *vertices,
                                     int *in_order, int *in_vertex,
                                     int *out_vertex, int *ierr)
{
   return LB_Current_lb->Get_First_Coarse_Obj_Fort(data, 
                                                   &num_gid_entries, 
                                                   &num_lid_entries,
                                                   global_id, local_id,
                                                   assigned, num_vert, vertices,
                                                   in_order, in_vertex,
                                                   out_vertex, ierr);
}

int LB_Next_Coarse_Obj_Fort_Wrapper(void *data, int num_gid_entries, 
                                    int num_lid_entries, LB_ID_PTR global_id,
                                    LB_ID_PTR local_id, 
                                    LB_ID_PTR next_global_id, 
                                    LB_ID_PTR next_local_id,
                                    int *assigned,
                                    int *num_vert, int *vertices,
                                    int *in_vertex, int *out_vertex, int *ierr)
{
   return LB_Current_lb->Get_Next_Coarse_Obj_Fort(data, &num_gid_entries,
                                                  &num_lid_entries,
                                                  global_id, local_id,
                                                  next_global_id, next_local_id,
                                                  assigned, num_vert, vertices,
                                                  in_vertex, out_vertex, ierr);
}

int LB_Num_Child_Fort_Wrapper(void *data, 
                              int num_gid_entries, int num_lid_entries, 
                              LB_ID_PTR global_id, LB_ID_PTR local_id,
                              int *ierr)
{
   return LB_Current_lb->Get_Num_Child_Fort(data, 
                                            &num_gid_entries, &num_lid_entries,
                                            global_id, local_id, ierr);
}

void LB_Child_List_Fort_Wrapper(void *data, 
                                int num_gid_entries, int num_lid_entries, 
                                LB_ID_PTR parent_gid,
                                LB_ID_PTR parent_lid, LB_ID_PTR child_gids,
                                LB_ID_PTR child_lids, int *assigned,
                                int *num_vert, int *vertices,
                                LB_REF_TYPE *ref_type, int *in_vertex,
                                int *out_vertex, int *ierr)
{
   LB_Current_lb->Get_Child_List_Fort(data, &num_gid_entries, &num_lid_entries,
                                      parent_gid, parent_lid,
                                      child_gids, child_lids, assigned,
                                      num_vert, vertices,
                                      ref_type, in_vertex, out_vertex, ierr);
}

void LB_Child_Weight_Fort_Wrapper(void *data, 
                                  int num_gid_entries, int num_lid_entries,
                                  LB_ID_PTR global_id, LB_ID_PTR local_id,
                                  int wgt_dim, float *obj_wgt, int *ierr)
{
   LB_Current_lb->Get_Child_Weight_Fort(data, 
                                        &num_gid_entries, &num_lid_entries,
                                        global_id, local_id, &wgt_dim,
                                        obj_wgt, ierr);
}
/*--------------------------------------------------------------------*/
/* C wrapper functions                                                */

int LB_fw_Initialize(float *ver)
{
   int myArgc;
   char **myArgv;
   int result;
   myArgc = 1;
   myArgv = (char **) LB_MALLOC((myArgc+1)*sizeof(char *));
   myArgv[0] = "unknown";
   myArgv[1] = NULL;
   result = LB_Initialize(myArgc,myArgv,ver);
   LB_FREE(&myArgv);
   return result;
}

int LB_fw_Initialize1(int *argc, int *argv, int *starts, float *ver)
{
   int i, j, result;
   char **myArgv;
   myArgv = (char **) LB_MALLOC(((*argc)+1)*sizeof(char *));
   for (i=0; i<(*argc); i++) {
      myArgv[i] = (char *) LB_MALLOC((starts[i+1]-starts[i]+1)*sizeof(char));
      for (j=0; j<starts[i+1]-starts[i]; j++) {
         myArgv[i][j] = (char) argv[starts[i]+j-1];
      }
      myArgv[i][starts[i+1]-starts[i]] = '\0';
   }
   myArgv[*argc] = NULL;
   result = LB_Initialize(*argc,myArgv,ver);
   for (i=0; i<(*argc); i++) 
     LB_FREE(&(myArgv[i]));
   LB_FREE(&myArgv);
   return result;
}

void LB_fw_Create(int *f_communicator, int *addr_lb, int *nbytes)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   MPI_Comm c_communicator;
   c_communicator = LB_comm_f2c(f_communicator);
   lb = LB_Create(c_communicator);
   lb->Fortran = 1;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {addr_lb[i] = (int)*p; p++;}
}

void LB_fw_Destroy(int *addr_lb, int *nbytes)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Destroy(&lb);
}

void LB_fw_Memory_Stats()
{
   LB_Memory_Stats();
}

int LB_fw_Set_Fn(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                 void *data)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;
   switch(*type) {
   case LB_NUM_EDGES_FN_TYPE:
      lb->Get_Num_Edges_Fort = (LB_NUM_EDGES_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Num_Edges_Fort_Wrapper, data);
      break;
   case LB_EDGE_LIST_FN_TYPE:
      lb->Get_Edge_List_Fort = (LB_EDGE_LIST_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Edge_List_Fort_Wrapper, data);
      break;
   case LB_NUM_GEOM_FN_TYPE:
      lb->Get_Num_Geom_Fort = (LB_NUM_GEOM_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Num_Geom_Fort_Wrapper, data);
      break;
   case LB_GEOM_FN_TYPE:
      lb->Get_Geom_Fort = (LB_GEOM_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Geom_Fort_Wrapper, data);
      break;
   case LB_NUM_OBJ_FN_TYPE:
      lb->Get_Num_Obj_Fort = (LB_NUM_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Num_Obj_Fort_Wrapper, data);
      break;
   case LB_OBJ_LIST_FN_TYPE:
      lb->Get_Obj_List_Fort = (LB_OBJ_LIST_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Obj_List_Fort_Wrapper, data);
      break;
   case LB_FIRST_OBJ_FN_TYPE:
      lb->Get_First_Obj_Fort = (LB_FIRST_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_First_Obj_Fort_Wrapper, data);
      break;
   case LB_NEXT_OBJ_FN_TYPE:
      lb->Get_Next_Obj_Fort = (LB_NEXT_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Next_Obj_Fort_Wrapper, data);
      break;
   case LB_NUM_BORDER_OBJ_FN_TYPE:
      lb->Get_Num_Border_Obj_Fort = (LB_NUM_BORDER_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Num_Border_Obj_Fort_Wrapper, data);
      break;
   case LB_BORDER_OBJ_LIST_FN_TYPE:
      lb->Get_Border_Obj_List_Fort = (LB_BORDER_OBJ_LIST_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Border_Obj_List_Fort_Wrapper, data);
      break;
   case LB_FIRST_BORDER_OBJ_FN_TYPE:
      lb->Get_First_Border_Obj_Fort = (LB_FIRST_BORDER_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_First_Border_Obj_Fort_Wrapper, data);
      break;
   case LB_NEXT_BORDER_OBJ_FN_TYPE:
      lb->Get_Next_Border_Obj_Fort = (LB_NEXT_BORDER_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Next_Border_Obj_Fort_Wrapper, data);
      break;
   case LB_PRE_MIGRATE_FN_TYPE:
      lb->Migrate.Pre_Migrate_Fort = (LB_PRE_MIGRATE_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Pre_Migrate_Fort_Wrapper, data);
      break;
   case LB_MID_MIGRATE_FN_TYPE:
      lb->Migrate.Mid_Migrate_Fort = (LB_MID_MIGRATE_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Mid_Migrate_Fort_Wrapper, data);
      break;
   case LB_POST_MIGRATE_FN_TYPE:
      lb->Migrate.Post_Migrate_Fort = (LB_POST_MIGRATE_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Post_Migrate_Fort_Wrapper, data);
      break;
   case LB_OBJ_SIZE_FN_TYPE:
      lb->Migrate.Get_Obj_Size_Fort = (LB_OBJ_SIZE_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Obj_Size_Fort_Wrapper, data);
      break;
   case LB_PACK_OBJ_FN_TYPE:
      lb->Migrate.Pack_Obj_Fort = (LB_PACK_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Pack_Obj_Fort_Wrapper, data);
      break;
   case LB_UNPACK_OBJ_FN_TYPE:
      lb->Migrate.Unpack_Obj_Fort = (LB_UNPACK_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Unpack_Obj_Fort_Wrapper, data);
      break;
   case LB_NUM_COARSE_OBJ_FN_TYPE:
      lb->Get_Num_Coarse_Obj_Fort = (LB_NUM_COARSE_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Num_Coarse_Obj_Fort_Wrapper, data);
      break;
   case LB_COARSE_OBJ_LIST_FN_TYPE:
      lb->Get_Coarse_Obj_List_Fort = (LB_COARSE_OBJ_LIST_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Coarse_Obj_List_Fort_Wrapper, data);
      break;
   case LB_FIRST_COARSE_OBJ_FN_TYPE:
      lb->Get_First_Coarse_Obj_Fort = (LB_FIRST_COARSE_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_First_Coarse_Obj_Fort_Wrapper, data);
      break;
   case LB_NEXT_COARSE_OBJ_FN_TYPE:
      lb->Get_Next_Coarse_Obj_Fort = (LB_NEXT_COARSE_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Next_Coarse_Obj_Fort_Wrapper, data);
      break;
   case LB_NUM_CHILD_FN_TYPE:
      lb->Get_Num_Child_Fort = (LB_NUM_CHILD_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Num_Child_Fort_Wrapper, data);
      break;
   case LB_CHILD_LIST_FN_TYPE:
      lb->Get_Child_List_Fort = (LB_CHILD_LIST_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Child_List_Fort_Wrapper, data);
      break;
   case LB_CHILD_WEIGHT_FN_TYPE:
      lb->Get_Child_Weight_Fort = (LB_CHILD_WEIGHT_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void (*)())LB_Child_Weight_Fort_Wrapper, data);
      break;

   default:
      return LB_Set_Fn(lb, *type, (void (*)())NULL, data);
      break;
   }
}

int LB_fw_Set_Fn0f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)())
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)NULL);
}

int LB_fw_Set_Fn1f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  int *data)
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn2f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  float *data)
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn3f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  double *data)
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn4f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  void *data)
/* data is type(LB_User_Data_1) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn5f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  void *data)
/* data is type(LB_User_Data_2) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn6f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  void *data)
/* data is type(LB_User_Data_3) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn7f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  void *data)
/* data is type(LB_User_Data_4) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn0s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)())
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)NULL);
}

int LB_fw_Set_Fn1s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  int *data)
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn2s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  float *data)
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn3s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  double *data)
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn4s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  void *data)
/* data is type(LB_User_Data_1) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn5s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  void *data)
/* data is type(LB_User_Data_2) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn6s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  void *data)
/* data is type(LB_User_Data_3) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn7s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void (*fn)(),
                  void *data)
/* data is type(LB_User_Data_4) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Method(int *addr_lb, int *nbytes, int *int_str, int *len)
{
   struct LB_Struct *lb;
   char *str;
   unsigned char *p;
   int i, result;
   str = (char *)LB_MALLOC(*len+1);
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;
   for (i=0; i<(*len); i++) str[i] = (char)int_str[i];
   str[*len] = '\0';
   result = LB_Set_Method(lb, str);
   LB_FREE(&str);
   return result;
}

int LB_fw_Set_Param(int *addr_lb, int *nbytes, int *int_param_name,
                   int *param_name_len, int *int_new_value, int *new_value_len)
{
   struct LB_Struct *lb;
   char *param_name, *new_value;
   unsigned char *p;
   int i, result;
   param_name = (char *)LB_MALLOC(*param_name_len+1);
   new_value = (char *)LB_MALLOC(*new_value_len+1);
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;
   for (i=0; i<(*param_name_len); i++) param_name[i] = (char)int_param_name[i];
   param_name[*param_name_len] = '\0';
   for (i=0; i<(*new_value_len); i++) new_value[i] = (char)int_new_value[i];
   new_value[*new_value_len] = '\0';
   result = LB_Set_Param(lb, param_name, new_value);
   LB_FREE(&param_name);
   LB_FREE(&new_value);
   return result;
}

int LB_fw_Balance(int *addr_lb, int *nbytes, int *changes, 
                  int *num_gid_entries, int *num_lid_entries,
                  int *num_import,
                  LB_ID_PTR *import_global_ids, LB_ID_PTR *import_local_ids,
                  int **import_procs, int *num_export,
                  LB_ID_PTR *export_global_ids, LB_ID_PTR *export_local_ids,
                  int **export_procs
#ifdef PGI
/* PGI uses hidden arguments when it passes pointers */
                   ,int *imp_gid_hide, int *imp_lid_hide, int *imp_proc_hide,
                    int *exp_gid_hide, int *exp_lid_hide, int *exp_proc_hide
#endif
#ifdef FUJITSU
/* Fujitsu and Lahey use a hidden argument for every argument */
/* TEMP need to verify this with Fujitsu or Lahey */
                   ,int *addr_lb_hide, int *nbytes_hide, int *changes_hide,
                    int *num_gid_entries_hide, int *num_lid_entries_hide,
                    int *num_import_hide, int *imp_gid_hide, int *imp_lid_hide,
                    int *imp_proc_hide, int *num_export_hide, int *exp_gid_hide,
                    int *exp_lid_hide, int *exp_proc_hide
#endif
                    )
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
#if defined (PGI) || defined (FUJITSU)
#define F90LB_TEMP 3
#else
#define F90LB_TEMP 2
#endif
   LB_ID_PTR temp_imp_gid[F90LB_TEMP], temp_exp_gid[F90LB_TEMP];
   LB_ID_PTR temp_imp_lid[F90LB_TEMP], temp_exp_lid[F90LB_TEMP];
   int *temp_imp_proc[F90LB_TEMP], *temp_exp_proc[F90LB_TEMP];
#undef F90LB_TEMP

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */

   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;

/* put the address of the Fortran pointer into temp_*[1] to be passed to
   Fortran for allocation.  The address of the allocated space will be
   in temp_*[0] so it can be used by C without messing up the Fortran pointer*/

   temp_imp_gid[1] = (LB_ID_PTR)import_global_ids;
   temp_imp_lid[1] = (LB_ID_PTR)import_local_ids;
   temp_imp_proc[1] = (int *)import_procs;
   temp_exp_gid[1] = (LB_ID_PTR)export_global_ids;
   temp_exp_lid[1] = (LB_ID_PTR)export_local_ids;
   temp_exp_proc[1] = (int *)export_procs;

/* for PGI and FUJITSU, put the hidden argument in temp_*[2] */

#if defined (PGI) || defined (FUJITSU)
   temp_imp_gid[2] = (LB_ID_PTR)imp_gid_hide;
   temp_imp_lid[2] = (LB_ID_PTR)imp_lid_hide;
   temp_imp_proc[2] = (int *)imp_proc_hide;
   temp_exp_gid[2] = (LB_ID_PTR)exp_gid_hide;
   temp_exp_lid[2] = (LB_ID_PTR)exp_lid_hide;
   temp_exp_proc[2] = (int *)exp_proc_hide;
#endif

/* call LB_Balance */

   return LB_Balance(lb, changes, num_gid_entries, num_lid_entries, 
                     num_import, temp_imp_gid, temp_imp_lid,
                     temp_imp_proc, num_export, temp_exp_gid, temp_exp_lid,
                     temp_exp_proc);
}


int LB_fw_Eval(int *addr_lb, int *nbytes, int *print_stats,
                int *nobj, float *obj_wgt, int *ncuts, float *cut_wgt,
                int *nboundary, int *nadj,
                int *is_nobj, int *is_obj_wgt, int *is_ncuts, int *is_cut_wgt,
                int *is_nboundary, int *is_nadj)
{
   struct LB_Struct *lb;
   int *loc_nobj, *loc_ncuts, *loc_nboundary, *loc_nadj;
   float *loc_obj_wgt, *loc_cut_wgt;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;
   if (*is_nobj) {loc_nobj = nobj;} else {loc_nobj = NULL;}
   if (*is_ncuts) {loc_ncuts = ncuts;} else {loc_ncuts = NULL;}
   if (*is_obj_wgt) {loc_obj_wgt = obj_wgt;} else {loc_obj_wgt = NULL;}
   if (*is_cut_wgt) {loc_cut_wgt = cut_wgt;} else {loc_cut_wgt = NULL;}
   if (*is_nboundary) {loc_nboundary = nboundary;} else {loc_nboundary = NULL;}
   if (*is_nadj) {loc_nadj = nadj;} else {loc_nadj = NULL;}

   return  LB_Eval(lb, *print_stats, loc_nobj, loc_obj_wgt, loc_ncuts, loc_cut_wgt,
           loc_nboundary, loc_nadj);
}

int LB_fw_Point_Assign(int *addr_lb, int *nbytes, double *coords, int *proc)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;

   return LB_Point_Assign(lb, coords, proc);
}

int LB_fw_Box_Assign(int *addr_lb, int *nbytes, double *xmin, double *ymin,
                     double *zmin, double *xmax, double *ymax, double *zmax,
                     int *procs, int *numprocs)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;

   return LB_Box_Assign(lb, *xmin, *ymin, *zmin, *xmax, *ymax, *zmax, procs,
                        numprocs);
}

int LB_fw_Compute_Destinations(int *addr_lb, int *nbytes, 
                    int *num_gid_entries, int *num_lid_entries, int *num_import,
                    LB_ID_PTR import_global_ids, LB_ID_PTR import_local_ids,
                    int *import_procs, int *num_export,
                    LB_ID_PTR *export_global_ids, LB_ID_PTR *export_local_ids,
                    int **export_procs
#ifdef PGI
                    ,int *exp_gid_hide, int *exp_lid_hide, int *exp_proc_hide
#endif
#ifdef FUJITSU
                   ,int *addr_lb_hide, int *nbytes_hide,
                    int *num_gid_entries_hide, int *num_lid_entries_hide,
                    int *num_import_hide,
                    int *import_global_ids_hide, int *import_local_ids_hide,
                    int *import_procs_hide, int *num_export_hide,
                    int *exp_gid_hide, int *exp_lid_hide, int *exp_proc_hide
#endif
                    )
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
#if defined (PGI) || defined(FUJITSU)
#define F90LB_TEMP 3
#else
#define F90LB_TEMP 2
#endif
   LB_ID_PTR temp_exp_gid[F90LB_TEMP];
   LB_ID_PTR temp_exp_lid[F90LB_TEMP];
   int *temp_exp_proc[F90LB_TEMP];
#undef F90LB_TEMP

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */

   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;

/* put the address of the Fortran pointer into temp_*[1] to be passed to
   Fortran for allocation.  The address of the allocated space will be
   in temp_*[0] so it can be used by C without messing up the Fortran pointer*/

   temp_exp_gid[1] = (LB_ID_PTR)export_global_ids;
   temp_exp_lid[1] = (LB_ID_PTR)export_local_ids;
   temp_exp_proc[1] = (int *)export_procs;

/* for PGI and FUJITSU, put the hidden argument in temp_*[2] */

#if defined (PGI) || defined(FUJITSU)
   temp_exp_gid[2] = (LB_ID_PTR)exp_gid_hide;
   temp_exp_lid[2] = (LB_ID_PTR)exp_lid_hide;
   temp_exp_proc[2] = (int *)exp_proc_hide;
#endif

/* call LB_Compute_Destinations */

   return LB_Compute_Destinations(lb, 
                     *num_import, import_global_ids,
                     import_local_ids, import_procs, 
                     num_export, temp_exp_gid, temp_exp_lid,
                     temp_exp_proc);
}


int LB_fw_Help_Migrate(int *addr_lb, int *nbytes, 
                    int *num_import,
                    LB_ID_PTR import_global_ids, LB_ID_PTR import_local_ids,
                    int *import_procs, int *num_export,
                    LB_ID_PTR export_global_ids, LB_ID_PTR export_local_ids,
                    int *export_procs)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;
   return LB_Help_Migrate(lb,
                          *num_import,import_global_ids,import_local_ids,
                          import_procs,*num_export,export_global_ids,
                          export_local_ids,export_procs);
}

void LB_fw_Register_Fort_Malloc(LB_FORT_MALLOC_INT_FN *fort_malloc_int,
                                LB_FORT_FREE_INT_FN *fort_free_int)
{
   LB_Register_Fort_Malloc(fort_malloc_int,fort_free_int);
}

/* TEMP child_order */
void LB_fw_Get_Child_Order(int *addr_lb, int *nbytes, int *order, int *ierr)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Get_Child_Order(lb,order,ierr);
}
/* end TEMP child_order */
