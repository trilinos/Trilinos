#include "lb_const.h"
#include "all_allo_const.h"

/*--------------------------------------------------------------------*/
/* procedure name mangling                                            */

#define LOWERCASE   1
#define UPPERCASE   2
#define UNDERSCORE  3
#define UNDERSCORE2 4

#if FMANGLE==LOWERCASE

#define LB_fw_Initialize                 lb_fw_initialize
#define LB_fw_Initialize1                lb_fw_initialize1
#define LB_fw_Create_Object              lb_fw_create_object
#define LB_fw_Destroy_Object             lb_fw_destroy_object
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
#define LB_fw_Balance11                  lb_fw_balance11
#define LB_fw_Balance12                  lb_fw_balance12
#define LB_fw_Balance21                  lb_fw_balance21
#define LB_fw_Balance22                  lb_fw_balance22
#define LB_fw_Compute_Destinations11     lb_fw_compute_destinations11
#define LB_fw_Compute_Destinations12     lb_fw_compute_destinations12
#define LB_fw_Compute_Destinations21     lb_fw_compute_destinations21
#define LB_fw_Compute_Destinations22     lb_fw_compute_destinations22
#define LB_fw_Help_Migrate11             lb_fw_help_migrate11
#define LB_fw_Help_Migrate12             lb_fw_help_migrate12
#define LB_fw_Help_Migrate21             lb_fw_help_migrate21
#define LB_fw_Help_Migrate22             lb_fw_help_migrate22
#define LB_fw_Register_Fort_Malloc       lb_fw_register_fort_malloc
#define LB_fw_Get_Address_int            lb_fw_get_address_int
#define LB_fw_Get_Address_GID            lb_fw_get_address_gid
#define LB_fw_Get_Address_LID            lb_fw_get_address_lid

#elif FMANGLE==UPPERCASE

#define LB_fw_Initialize                 LB_FW_INITIALIZE
#define LB_fw_Initialize1                LB_FW_INITIALIZE1
#define LB_fw_Create_Object              LB_FW_CREATE_OBJECT
#define LB_fw_Destroy_Object             LB_FW_DESTROY_OBJECT
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
#define LB_fw_Balance11                  LB_FW_BALANCE11
#define LB_fw_Balance12                  LB_FW_BALANCE12
#define LB_fw_Balance21                  LB_FW_BALANCE21
#define LB_fw_Balance22                  LB_FW_BALANCE22
#define LB_fw_Compute_Destinations11     LB_FW_COMPUTE_DESTINATIONS11
#define LB_fw_Compute_Destinations12     LB_FW_COMPUTE_DESTINATIONS12
#define LB_fw_Compute_Destinations21     LB_FW_COMPUTE_DESTINATIONS21
#define LB_fw_Compute_Destinations22     LB_FW_COMPUTE_DESTINATIONS22
#define LB_fw_Help_Migrate11             LB_FW_HELP_MIGRATE11
#define LB_fw_Help_Migrate12             LB_FW_HELP_MIGRATE12
#define LB_fw_Help_Migrate21             LB_FW_HELP_MIGRATE21
#define LB_fw_Help_Migrate22             LB_FW_HELP_MIGRATE22
#define LB_fw_Register_Fort_Malloc       LB_FW_REGISTER_FORT_MALLOC
#define LB_fw_Get_Address_int            LB_FW_GET_ADDRESS_INT
#define LB_fw_Get_Address_GID            LB_FW_GET_ADDRESS_GID
#define LB_fw_Get_Address_LID            LB_FW_GET_ADDRESS_LID

#elif FMANGLE==UNDERSCORE

#define LB_fw_Initialize                 lb_fw_initialize_
#define LB_fw_Initialize1                lb_fw_initialize1_
#define LB_fw_Create_Object              lb_fw_create_object_
#define LB_fw_Destroy_Object             lb_fw_destroy_object_
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
#define LB_fw_Balance11                  lb_fw_balance11_
#define LB_fw_Balance12                  lb_fw_balance12_
#define LB_fw_Balance21                  lb_fw_balance21_
#define LB_fw_Balance22                  lb_fw_balance22_
#define LB_fw_Compute_Destinations11     lb_fw_compute_destinations11_
#define LB_fw_Compute_Destinations12     lb_fw_compute_destinations12_
#define LB_fw_Compute_Destinations21     lb_fw_compute_destinations21_
#define LB_fw_Compute_Destinations22     lb_fw_compute_destinations22_
#define LB_fw_Help_Migrate11             lb_fw_help_migrate11_
#define LB_fw_Help_Migrate12             lb_fw_help_migrate12_
#define LB_fw_Help_Migrate21             lb_fw_help_migrate21_
#define LB_fw_Help_Migrate22             lb_fw_help_migrate22_
#define LB_fw_Register_Fort_Malloc       lb_fw_register_fort_malloc_
#define LB_fw_Get_Address_int            lb_fw_get_address_int_
#define LB_fw_Get_Address_GID            lb_fw_get_address_gid_
#define LB_fw_Get_Address_LID            lb_fw_get_address_lid_

#elif FMANGLE==UNDERSCORE2

#define LB_fw_Initialize                 lb_fw_initialize__
#define LB_fw_Initialize1                lb_fw_initialize1__
#define LB_fw_Create_Object              lb_fw_create_object__
#define LB_fw_Destroy_Object             lb_fw_destroy_object__
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
#define LB_fw_Balance11                  lb_fw_balance11__
#define LB_fw_Balance12                  lb_fw_balance12__
#define LB_fw_Balance21                  lb_fw_balance21__
#define LB_fw_Balance22                  lb_fw_balance22__
#define LB_fw_Compute_Destinations11     lb_fw_compute_destinations11__
#define LB_fw_Compute_Destinations12     lb_fw_compute_destinations12__
#define LB_fw_Compute_Destinations21     lb_fw_compute_destinations21__
#define LB_fw_Compute_Destinations22     lb_fw_compute_destinations22__
#define LB_fw_Help_Migrate11             lb_fw_help_migrate11__
#define LB_fw_Help_Migrate12             lb_fw_help_migrate12__
#define LB_fw_Help_Migrate21             lb_fw_help_migrate21__
#define LB_fw_Help_Migrate22             lb_fw_help_migrate22__
#define LB_fw_Register_Fort_Malloc       lb_fw_register_fort_malloc__
#define LB_fw_Get_Address_int            lb_fw_get_address_int__
#define LB_fw_Get_Address_GID            lb_fw_get_address_gid__
#define LB_fw_Get_Address_LID            lb_fw_get_address_lid__

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
/* MPI 2 provides a standard way of doing this
   return MPI_Comm_f2c((MPI_Fint)(*f_comm));
#else
/* will probably need some special cases here */
/* when in doubt, just return the input */
   return (MPI_Comm)(*f_comm);
#endif
}

/* These routines get the address of an array allocated by fortran and
   return it */

void LB_fw_Get_Address_int(int *addr, int *ret_addr)
{
   *ret_addr = (int)addr;
}

void LB_fw_Get_Address_GID(LB_GID *addr, int *ret_addr)
{
   *ret_addr = (int)addr;
}

void LB_fw_Get_Address_LID(LB_LID *addr, int *ret_addr)
{
   *ret_addr = (int)addr;
}

/*--------------------------------------------------------------------*/
/* Reverse wrappers for callbacks                                     */

int LB_Num_Edges_Fort_Wrapper(void *data, LB_GID global_id, LB_LID local_id,
                              int *ierr)
{
   return LB_Current_lb->Get_Num_Edges_Fort(data, &global_id, &local_id, ierr);
}

void LB_Edge_List_Fort_Wrapper(void *data, LB_GID global_id, LB_LID local_id,
                               LB_GID *nbor_global_id, int *nbor_procs,
                               int wdim, int *nbor_ewgts, int *ierr)
{
   LB_Current_lb->Get_Edge_List_Fort(data, &global_id, &local_id,
                                     nbor_global_id, nbor_procs, &wdim,
                                     nbor_ewgts, ierr);
}

int LB_Num_Geom_Fort_Wrapper(void *data, int *ierr)
{
   return LB_Current_lb->Get_Num_Geom_Fort(data,ierr);
}

void LB_Geom_Fort_Wrapper(void *data, LB_GID global_id, LB_LID local_id,
                          double *geom_vec, int *ierr)
{
   LB_Current_lb->Get_Geom_Fort(data, &global_id, &local_id, geom_vec, ierr);
}

int LB_Num_Obj_Fort_Wrapper(void *data, int *ierr)
{
   return LB_Current_lb->Get_Num_Obj_Fort(data, ierr);
}

void LB_Obj_List_Fort_Wrapper(void *data, LB_GID *global_ids, LB_LID *local_ids,
                              int wdim, float *objwgts, int *ierr)
{
   LB_Current_lb->Get_Obj_List_Fort(data, global_ids, local_ids, &wdim,
                                    objwgts, ierr);
}

int LB_First_Obj_Fort_Wrapper(void *data, LB_GID *first_global_id,
                              LB_LID *first_local_id,
                              int wdim, float *first_obj_wgt, int *ierr)
{
   return LB_Current_lb->Get_First_Obj_Fort(data, first_global_id,
                                            first_local_id, &wdim,
                                            first_obj_wgt, ierr);
}

int LB_Next_Obj_Fort_Wrapper(void *data, LB_GID global_id, LB_LID local_id,
                             LB_GID *next_global_id, LB_LID *next_local_id,
                             int wdim, float *next_obj_wgt, int *ierr)
{
   return LB_Current_lb->Get_Next_Obj_Fort(data, &global_id, &local_id,
                                           next_global_id, next_local_id,
                                           &wdim, next_obj_wgt, ierr);
}

int LB_Num_Border_Obj_Fort_Wrapper(void *data, int nbor_proc, int *ierr)
{
   return LB_Current_lb->Get_Num_Border_Obj_Fort(data, &nbor_proc, ierr);
}

void LB_Border_Obj_List_Fort_Wrapper(void *data, int nbor_proc,
                                     LB_GID *global_ids, LB_LID *local_ids,
                                     int wdim, float *objwgts, int *ierr)
{
   LB_Current_lb->Get_Border_Obj_List_Fort(data, &nbor_proc, global_ids,
                                           local_ids, &wdim, objwgts, ierr);
}

int LB_First_Border_Obj_Fort_Wrapper(void *data, int nbor_proc,
                                     LB_GID *first_global_id,
                                     LB_LID *first_local_id,
                                     int wdim, float *first_obj_wgt,
                                     int *ierr)
{
   return LB_Current_lb->Get_First_Border_Obj_Fort(data, &nbor_proc,
                                                   first_global_id,
                                                   first_local_id, &wdim,
                                                   first_obj_wgt, ierr);
}

int LB_Next_Border_Obj_Fort_Wrapper(void *data, LB_GID global_id,
                                  LB_LID local_id, int nbor_proc,
                                  LB_GID *next_global_id,
                                  LB_LID *next_local_id,
                                  int wdim, float *next_obj_wgt,
                                  int *ierr)
{
   return LB_Current_lb->Get_Next_Border_Obj_Fort(data, &global_id, &local_id,
                                                  &nbor_proc, next_global_id,
                                                  next_local_id, &wdim,
                                                  next_obj_wgt, ierr);
}

int LB_Obj_Size_Fort_Wrapper(void *data, int *ierr)
{
   return LB_Current_lb->Migrate.Get_Obj_Size_Fort(data,ierr);
}

void LB_Pre_Migrate_Fort_Wrapper(void *data, int num_import,
                                 LB_GID *import_global_ids,
                                 LB_LID *import_local_ids, int *import_procs,
                                 int num_export, LB_GID *export_global_ids,
                                 LB_LID *export_local_ids, int *export_procs,
                                 int *ierr)
{
   LB_Current_lb->Migrate.Pre_Process_Fort(data, &num_import,
                                            import_global_ids,
                                            import_local_ids, import_procs,
                                            &num_export, export_global_ids,
                                            export_local_ids, export_procs,
                                            ierr);
}

void LB_Pack_Obj_Fort_Wrapper(void *data, LB_GID global_id, LB_LID local_id,
                            int dest_proc, int size, char *buf, int *ierr)
{
   LB_Current_lb->Migrate.Pack_Obj_Fort(data, &global_id, &local_id,
                                         &dest_proc, &size, buf, ierr);
}

void LB_Unpack_Obj_Fort_Wrapper(void *data, LB_GID global_id, int size,
                              char *buf, int *ierr)
{
   LB_Current_lb->Migrate.Unpack_Obj_Fort(data, &global_id, &size, buf, ierr);
}

/*--------------------------------------------------------------------*/
/* C wrapper functions                                                */

int LB_fw_Initialize(float *ver)
{
   int myArgc;
   char **myArgv;
   int result;
   myArgc = 1;
   myArgv = (char **) LB_Malloc((myArgc+1)*sizeof(char),__FILE__,__LINE__);
/* TEMP do I have to malloc myArgv[0]? */
   myArgv[0] = "unknown";
   myArgv[1] = NULL;
   result = LB_Initialize(myArgc,myArgv,ver);
/* TEMP should myArgv be freed or kept? */
/* LB_FREE((void **)&myArgv); */
   return result;
}

int LB_fw_Initialize1(int *argc, int *argv, int *starts, float *ver)
{
   int i, j, result;
   char **myArgv;
   myArgv = (char **) LB_Malloc(((*argc)+1)*sizeof(char *),__FILE__,__LINE__);
   for (i=0; i<(*argc); i++) {
      myArgv[i] = (char *) LB_Malloc((starts[i+1]-starts[i]+1)*sizeof(char),
                                     __FILE__,__LINE__);
      for (j=0; j<starts[i+1]-starts[i]; j++) {
         myArgv[i][j] = (char) argv[starts[i]+j-1];
      }
      myArgv[i][starts[i+1]-starts[i]] = '\0';
   }
   myArgv[*argc] = NULL;
   result = LB_Initialize(*argc,myArgv,ver);
   LB_FREE((void **)&myArgv);
   return result;
}

void LB_fw_Create_Object(int *f_communicator, int *addr_lb, int *nbytes)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   MPI_Comm c_communicator;
   c_communicator = LB_comm_f2c(f_communicator);
   lb = LB_Create_Object(c_communicator);
   lb->Fortran = 1;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {addr_lb[i] = (int)*p; p++;}
}

void LB_fw_Destroy_Object(int *addr_lb, int *nbytes)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Destroy_Object(&lb);
}

int LB_fw_Set_Fn(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
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
      return LB_Set_Fn(lb, *type, (void *)LB_Num_Edges_Fort_Wrapper, data);
      break;
   case LB_EDGE_LIST_FN_TYPE:
      lb->Get_Edge_List_Fort = (LB_EDGE_LIST_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_Edge_List_Fort_Wrapper, data);
      break;
   case LB_NUM_GEOM_FN_TYPE:
      lb->Get_Num_Geom_Fort = (LB_NUM_GEOM_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_Num_Geom_Fort_Wrapper, data);
      break;
   case LB_GEOM_FN_TYPE:
      lb->Get_Geom_Fort = (LB_GEOM_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_Geom_Fort_Wrapper, data);
      break;
   case LB_NUM_OBJ_FN_TYPE:
      lb->Get_Num_Obj_Fort = (LB_NUM_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_Num_Obj_Fort_Wrapper, data);
      break;
   case LB_OBJ_LIST_FN_TYPE:
      lb->Get_Obj_List_Fort = (LB_OBJ_LIST_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_Obj_List_Fort_Wrapper, data);
      break;
   case LB_FIRST_OBJ_FN_TYPE:
      lb->Get_First_Obj_Fort = (LB_FIRST_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_First_Obj_Fort_Wrapper, data);
      break;
   case LB_NEXT_OBJ_FN_TYPE:
      lb->Get_Next_Obj_Fort = (LB_NEXT_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_Next_Obj_Fort_Wrapper, data);
      break;
   case LB_NUM_BORDER_OBJ_FN_TYPE:
      lb->Get_Num_Border_Obj_Fort = (LB_NUM_BORDER_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_Num_Border_Obj_Fort_Wrapper, data);
      break;
   case LB_BORDER_OBJ_LIST_FN_TYPE:
      lb->Get_Border_Obj_List_Fort = (LB_BORDER_OBJ_LIST_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_Border_Obj_List_Fort_Wrapper, data);
      break;
   case LB_FIRST_BORDER_OBJ_FN_TYPE:
      lb->Get_First_Border_Obj_Fort = (LB_FIRST_BORDER_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_First_Border_Obj_Fort_Wrapper, data);
      break;
   case LB_NEXT_BORDER_OBJ_FN_TYPE:
      lb->Get_Next_Border_Obj_Fort = (LB_NEXT_BORDER_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_Next_Border_Obj_Fort_Wrapper, data);
      break;
   case LB_PRE_MIGRATE_FN_TYPE:
      lb->Migrate.Pre_Process_Fort = (LB_PRE_MIGRATE_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_Pre_Migrate_Fort_Wrapper, data);
      break;
   case LB_OBJ_SIZE_FN_TYPE:
      lb->Migrate.Get_Obj_Size_Fort = (LB_OBJ_SIZE_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_Obj_Size_Fort_Wrapper, data);
      break;
   case LB_PACK_OBJ_FN_TYPE:
      lb->Migrate.Pack_Obj_Fort = (LB_PACK_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_Pack_Obj_Fort_Wrapper, data);
      break;
   case LB_UNPACK_OBJ_FN_TYPE:
      lb->Migrate.Unpack_Obj_Fort = (LB_UNPACK_OBJ_FORT_FN *) fn;
      return LB_Set_Fn(lb, *type, (void *)LB_Unpack_Obj_Fort_Wrapper, data);
      break;
   default:
      return LB_Set_Fn(lb, *type, (void *)NULL, data);
      break;
   }
}

int LB_fw_Set_Fn0f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn())
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)NULL);
}

int LB_fw_Set_Fn1f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
                  int *data)
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn2f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
                  float *data)
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn3f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
                  double *data)
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn4f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
                  void *data)
/* data is type(LB_User_Data_1) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn5f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
                  void *data)
/* data is type(LB_User_Data_2) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn6f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
                  void *data)
/* data is type(LB_User_Data_3) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn7f(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
                  void *data)
/* data is type(LB_User_Data_4) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn0s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn())
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)NULL);
}

int LB_fw_Set_Fn1s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
                  int *data)
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn2s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
                  float *data)
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn3s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
                  double *data)
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn4s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
                  void *data)
/* data is type(LB_User_Data_1) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn5s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
                  void *data)
/* data is type(LB_User_Data_2) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn6s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
                  void *data)
/* data is type(LB_User_Data_3) */
{
   return LB_fw_Set_Fn(addr_lb, nbytes, type, fn, (void *)data);
}

int LB_fw_Set_Fn7s(int *addr_lb, int *nbytes, LB_FN_TYPE *type, void *fn(),
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
   str = (char *)LB_Malloc(*len+1,__FILE__,__LINE__);
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;
   for (i=0; i<(*len); i++) str[i] = (char)int_str[i];
   str[*len] = '\0';
   result = LB_Set_Method(lb, str);
   LB_FREE((void **)&str);
   return result;
}

int LB_fw_Set_Param(int *addr_lb, int *nbytes, int *int_param_name,
                   int *param_name_len, int *int_new_value, int *new_value_len)
{
   struct LB_Struct *lb;
   char *param_name, *new_value;
   unsigned char *p;
   int i;
   param_name = (char *)LB_Malloc(*param_name_len+1,__FILE__,__LINE__);
   new_value = (char *)LB_Malloc(*new_value_len+1,__FILE__,__LINE__);
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;
   for (i=0; i<(*param_name_len); i++) param_name[i] = (char)int_param_name[i];
   param_name[*param_name_len] = '\0';
   for (i=0; i<(*new_value_len); i++) new_value[i] = (char)int_new_value[i];
   new_value[*new_value_len] = '\0';
   return LB_Set_Param(lb, param_name, new_value);
}

int LB_fw_Balance11(int *addr_lb, int *nbytes, int *changes, int *num_import,
                    LB_GID **import_global_ids, LB_LID **import_local_ids,
                    int **import_procs, int *num_export,
                    LB_GID **export_global_ids, LB_LID **export_local_ids,
                    int **export_procs
#ifdef PGI
/* PGI uses hidden arguments when it passes pointers */
                   ,int *imp_gid_hide, int *imp_lid_hide, int *imp_proc_hide,
                    int *exp_gid_hide, int *exp_lid_hide, int *exp_proc_hide
#endif
                    )
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
#ifdef PGI
#define F90LB_TEMP 3
#else
#define F90LB_TEMP 2
#endif
   LB_GID *temp_imp_gid[F90LB_TEMP], *temp_exp_gid[F90LB_TEMP];
   LB_LID *temp_imp_lid[F90LB_TEMP], *temp_exp_lid[F90LB_TEMP];
   int *temp_imp_proc[F90LB_TEMP], *temp_exp_proc[F90LB_TEMP];
#undef F90LB_TEMP

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */

   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;

/* put the address of the Fortran pointer into temp_*[1] to be passed to
   Fortran for allocation.  The address of the allocated space will be
   in temp_*[0] so it can be used by C without messing up the Fortran pointer*/

   temp_imp_gid[1] = (LB_GID *)import_global_ids;
   temp_imp_lid[1] = (LB_LID *)import_local_ids;
   temp_imp_proc[1] = (int *)import_procs;
   temp_exp_gid[1] = (LB_GID *)export_global_ids;
   temp_exp_lid[1] = (LB_LID *)export_local_ids;
   temp_exp_proc[1] = (int *)export_procs;

/* for PGI, put the hidden argument in temp_*[2] */

#ifdef PGI
   temp_imp_gid[2] = (LB_GID *)imp_gid_hide;
   temp_imp_lid[2] = (LB_LID *)imp_lid_hide;
   temp_imp_proc[2] = (int *)imp_proc_hide;
   temp_exp_gid[2] = (LB_GID *)exp_gid_hide;
   temp_exp_lid[2] = (LB_LID *)exp_lid_hide;
   temp_exp_proc[2] = (int *)exp_proc_hide;
#endif

/* call LB_Balance */

   return LB_Balance(lb, changes, num_import, temp_imp_gid, temp_imp_lid,
                     temp_imp_proc, num_export, temp_exp_gid, temp_exp_lid,
                     temp_exp_proc);
}

int LB_fw_Balance12(int *addr_lb, int *nbytes, int *changes, int *num_import,
                    LB_GID **import_global_ids, LB_LID **import_local_ids,
                    int **import_procs, int *num_export,
                    LB_GID **export_global_ids, LB_LID **export_local_ids,
                    int **export_procs
#ifdef PGI
                   ,int *imp_gid_hide, int *imp_lid_hide, int *imp_proc_hide,
                    int *exp_gid_hide, int *exp_lid_hide, int *exp_proc_hide
#endif
                    )
{return LB_fw_Balance11(addr_lb, nbytes, changes, num_import, import_global_ids,
                    import_local_ids, import_procs, num_export,
                    export_global_ids, export_local_ids, export_procs
#ifdef PGI
                   ,imp_gid_hide, imp_lid_hide, imp_proc_hide,
                    exp_gid_hide, exp_lid_hide, exp_proc_hide
#endif
                    );
}

int LB_fw_Balance21(int *addr_lb, int *nbytes, int *changes, int *num_import,
                    LB_GID **import_global_ids, LB_LID **import_local_ids,
                    int **import_procs, int *num_export,
                    LB_GID **export_global_ids, LB_LID **export_local_ids,
                    int **export_procs
#ifdef PGI
                   ,int *imp_gid_hide, int *imp_lid_hide, int *imp_proc_hide,
                    int *exp_gid_hide, int *exp_lid_hide, int *exp_proc_hide
#endif
                    )
{return LB_fw_Balance11(addr_lb, nbytes, changes, num_import, import_global_ids,
                    import_local_ids, import_procs, num_export,
                    export_global_ids, export_local_ids, export_procs
#ifdef PGI
                   ,imp_gid_hide, imp_lid_hide, imp_proc_hide,
                    exp_gid_hide, exp_lid_hide, exp_proc_hide
#endif
                    );
}

int LB_fw_Balance22(int *addr_lb, int *nbytes, int *changes, int *num_import,
                    LB_GID **import_global_ids, LB_LID **import_local_ids,
                    int **import_procs, int *num_export,
                    LB_GID **export_global_ids, LB_LID **export_local_ids,
                    int **export_procs
#ifdef PGI
                   ,int *imp_gid_hide, int *imp_lid_hide, int *imp_proc_hide,
                    int *exp_gid_hide, int *exp_lid_hide, int *exp_proc_hide
#endif
                    )
{return LB_fw_Balance11(addr_lb, nbytes, changes, num_import, import_global_ids,
                    import_local_ids, import_procs, num_export,
                    export_global_ids, export_local_ids, export_procs
#ifdef PGI
                   ,imp_gid_hide, imp_lid_hide, imp_proc_hide,
                    exp_gid_hide, exp_lid_hide, exp_proc_hide
#endif
                    );
}

int LB_fw_Compute_Destinations11(int *addr_lb, int *nbytes, int *num_import,
                    LB_GID *import_global_ids, LB_LID *import_local_ids,
                    int *import_procs, int *num_export,
                    LB_GID **export_global_ids, LB_LID **export_local_ids,
                    int **export_procs
#ifdef PGI
                    ,int *exp_gid_hide, int *exp_lid_hide, int *exp_proc_hide
#endif
                    )
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
#ifdef PGI
#define F90LB_TEMP 3
#else
#define F90LB_TEMP 2
#endif
   LB_GID *temp_exp_gid[F90LB_TEMP];
   LB_LID *temp_exp_lid[F90LB_TEMP];
   int *temp_exp_proc[F90LB_TEMP];
#undef F90LB_TEMP

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */

   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;

/* put the address of the Fortran pointer into temp_*[1] to be passed to
   Fortran for allocation.  The address of the allocated space will be
   in temp_*[0] so it can be used by C without messing up the Fortran pointer*/

   temp_exp_gid[1] = (LB_GID *)export_global_ids;
   temp_exp_lid[1] = (LB_LID *)export_local_ids;
   temp_exp_proc[1] = (int *)export_procs;

/* for PGI, put the hidden argument in temp_*[2] */

#ifdef PGI
   temp_exp_gid[2] = (LB_GID *)exp_gid_hide;
   temp_exp_lid[2] = (LB_LID *)exp_lid_hide;
   temp_exp_proc[2] = (int *)exp_proc_hide;
#endif

/* call LB_Compute_Destinations */

   return LB_Compute_Destinations(lb, *num_import, import_global_ids,
                     import_local_ids, import_procs, num_export, temp_exp_gid, temp_exp_lid,
                     temp_exp_proc);
}

int LB_fw_Compute_Destinations12(int *addr_lb, int *nbytes, int *num_import,
                    LB_GID *import_global_ids, LB_LID *import_local_ids,
                    int *import_procs, int *num_export,
                    LB_GID **export_global_ids, LB_LID **export_local_ids,
                    int **export_procs
#ifdef PGI
                    ,int *exp_gid_hide, int *exp_lid_hide, int *exp_proc_hide
#endif
                    )
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
#ifdef PGI
#define F90LB_TEMP 3
#else
#define F90LB_TEMP 2
#endif
   LB_GID *temp_exp_gid[F90LB_TEMP];
   LB_LID *temp_exp_lid[F90LB_TEMP];
   int *temp_exp_proc[F90LB_TEMP];
#undef F90LB_TEMP

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */

   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;

/* put the address of the Fortran pointer into temp_*[1] to be passed to
   Fortran for allocation.  The address of the allocated space will be
   in temp_*[0] so it can be used by C without messing up the Fortran pointer*/

   temp_exp_gid[1] = (LB_GID *)export_global_ids;
   temp_exp_lid[1] = (LB_LID *)export_local_ids;
   temp_exp_proc[1] = (int *)export_procs;

/* for PGI, put the hidden argument in temp_*[2] */

#ifdef PGI
   temp_exp_gid[2] = (LB_GID *)exp_gid_hide;
   temp_exp_lid[2] = (LB_LID *)exp_lid_hide;
   temp_exp_proc[2] = (int *)exp_proc_hide;
#endif

/* call LB_Compute_Destinations */

   return LB_Compute_Destinations(lb, *num_import, import_global_ids,
                     import_local_ids, import_procs, num_export, temp_exp_gid, temp_exp_lid,
                     temp_exp_proc);
}

int LB_fw_Compute_Destinations21(int *addr_lb, int *nbytes, int *num_import,
                    LB_GID *import_global_ids, LB_LID *import_local_ids,
                    int *import_procs, int *num_export,
                    LB_GID **export_global_ids, LB_LID **export_local_ids,
                    int **export_procs
#ifdef PGI
                    ,int *exp_gid_hide, int *exp_lid_hide, int *exp_proc_hide
#endif
                    )
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
#ifdef PGI
#define F90LB_TEMP 3
#else
#define F90LB_TEMP 2
#endif
   LB_GID *temp_exp_gid[F90LB_TEMP];
   LB_LID *temp_exp_lid[F90LB_TEMP];
   int *temp_exp_proc[F90LB_TEMP];
#undef F90LB_TEMP

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */

   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;

/* put the address of the Fortran pointer into temp_*[1] to be passed to
   Fortran for allocation.  The address of the allocated space will be
   in temp_*[0] so it can be used by C without messing up the Fortran pointer*/

   temp_exp_gid[1] = (LB_GID *)export_global_ids;
   temp_exp_lid[1] = (LB_LID *)export_local_ids;
   temp_exp_proc[1] = (int *)export_procs;

/* for PGI, put the hidden argument in temp_*[2] */

#ifdef PGI
   temp_exp_gid[2] = (LB_GID *)exp_gid_hide;
   temp_exp_lid[2] = (LB_LID *)exp_lid_hide;
   temp_exp_proc[2] = (int *)exp_proc_hide;
#endif

/* call LB_Compute_Destinations */

   return LB_Compute_Destinations(lb, *num_import, import_global_ids,
                     import_local_ids, import_procs, num_export, temp_exp_gid, temp_exp_lid,
                     temp_exp_proc);
}

int LB_fw_Compute_Destinations22(int *addr_lb, int *nbytes, int *num_import,
                    LB_GID *import_global_ids, LB_LID *import_local_ids,
                    int *import_procs, int *num_export,
                    LB_GID **export_global_ids, LB_LID **export_local_ids,
                    int **export_procs
#ifdef PGI
                    ,int *exp_gid_hide, int *exp_lid_hide, int *exp_proc_hide
#endif
                    )
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
#ifdef PGI
#define F90LB_TEMP 3
#else
#define F90LB_TEMP 2
#endif
   LB_GID *temp_exp_gid[F90LB_TEMP];
   LB_LID *temp_exp_lid[F90LB_TEMP];
   int *temp_exp_proc[F90LB_TEMP];
#undef F90LB_TEMP

/* reconstruct the lb pointer from the nbyte 1-byte integers in addr_lb */

   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;

/* put the address of the Fortran pointer into temp_*[1] to be passed to
   Fortran for allocation.  The address of the allocated space will be
   in temp_*[0] so it can be used by C without messing up the Fortran pointer*/

   temp_exp_gid[1] = (LB_GID *)export_global_ids;
   temp_exp_lid[1] = (LB_LID *)export_local_ids;
   temp_exp_proc[1] = (int *)export_procs;

/* for PGI, put the hidden argument in temp_*[2] */

#ifdef PGI
   temp_exp_gid[2] = (LB_GID *)exp_gid_hide;
   temp_exp_lid[2] = (LB_LID *)exp_lid_hide;
   temp_exp_proc[2] = (int *)exp_proc_hide;
#endif

/* call LB_Compute_Destinations */

   return LB_Compute_Destinations(lb, *num_import, import_global_ids,
                     import_local_ids, import_procs, num_export, temp_exp_gid, temp_exp_lid,
                     temp_exp_proc);
}

int LB_fw_Help_Migrate11(int *addr_lb, int *nbytes, int *num_import,
                    LB_GID *import_global_ids, LB_LID *import_local_ids,
                    int *import_procs, int *num_export,
                    LB_GID *export_global_ids, LB_LID *export_local_ids,
                    int *export_procs)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;
   return LB_Help_Migrate(lb,*num_import,import_global_ids,import_local_ids,
                          import_procs,*num_export,export_global_ids,
                          export_local_ids,export_procs);
}

int LB_fw_Help_Migrate12(int *addr_lb, int *nbytes, int *num_import,
                    LB_GID *import_global_ids, LB_LID *import_local_ids,
                    int *import_procs, int *num_export,
                    LB_GID *export_global_ids, LB_LID *export_local_ids,
                    int *export_procs)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;
   return LB_Help_Migrate(lb,*num_import,import_global_ids,import_local_ids,
                          import_procs,*num_export,export_global_ids,
                          export_local_ids,export_procs);
}

int LB_fw_Help_Migrate21(int *addr_lb, int *nbytes, int *num_import,
                    LB_GID *import_global_ids, LB_LID *import_local_ids,
                    int *import_procs, int *num_export,
                    LB_GID *export_global_ids, LB_LID *export_local_ids,
                    int *export_procs)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;
   return LB_Help_Migrate(lb,*num_import,import_global_ids,import_local_ids,
                          import_procs,*num_export,export_global_ids,
                          export_local_ids,export_procs);
}

int LB_fw_Help_Migrate22(int *addr_lb, int *nbytes, int *num_import,
                    LB_GID *import_global_ids, LB_LID *import_local_ids,
                    int *import_procs, int *num_export,
                    LB_GID *export_global_ids, LB_LID *export_local_ids,
                    int *export_procs)
{
   struct LB_Struct *lb;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lb;
   for (i=0; i<(*nbytes); i++) {*p = (unsigned char)addr_lb[i]; p++;}
   LB_Current_lb = lb;
   return LB_Help_Migrate(lb,*num_import,import_global_ids,import_local_ids,
                          import_procs,*num_export,export_global_ids,
                          export_local_ids,export_procs);
}

void LB_fw_Register_Fort_Malloc(LB_FORT_MALLOC_INT_FN *fort_malloc_int,
                                LB_FORT_MALLOC_GID_FN *fort_malloc_GID,
                                LB_FORT_MALLOC_LID_FN *fort_malloc_LID,
                                LB_FORT_FREE_INT_FN *fort_free_int,
                                LB_FORT_FREE_GID_FN *fort_free_GID,
                                LB_FORT_FREE_LID_FN *fort_free_LID)
{
   LB_Register_Fort_Malloc(fort_malloc_int,fort_malloc_GID,fort_malloc_LID,
                           fort_free_int,fort_free_GID,fort_free_LID);
}
