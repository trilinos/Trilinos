/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
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
#include "ha_const.h"
#include "params_const.h"

/* 
 * These routines build a description of the machine
 * defined by the processes in the MPI communicator
 * and the machine description file.
 */

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/****** Parameters structure for building machine description. *****/

static PARAM_VARS Mach_params[] = {
        { "USE_MACHINE_DESC", NULL, "INT", 0 },
        { "MACHINE_DESC_FILE", NULL, "STRING", 0 },
        { NULL, NULL, NULL, 0 } };

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

int Zoltan_Set_Machine_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
  int status;
  PARAM_UTYPE result;         /* value returned from Check_Param */
  int index;                  /* index returned from Check_Param */

  status = Zoltan_Check_Param(name, val, Mach_params, &result, &index);

  return(status);
}

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

int Zoltan_Build_Machine_Desc(
   ZZ *zz              /* The Zoltan structure.                */
)
{
  char *yo = "Zoltan_Build_Machine_Desc";
  int ierr = ZOLTAN_OK;
  int use_mach_desc;
  char filename[256];

  Zoltan_Bind_Param(Mach_params, "USE_MACHINE_DESC", (void *) &use_mach_desc);
  Zoltan_Bind_Param(Mach_params, "MACHINE_DESC_FILE", (void *) filename);

  use_mach_desc = 0;
  strcpy(filename, MACHINE_DESC_FILE_DEFAULT);

  Zoltan_Assign_Param_Vals(zz->Params, Mach_params, zz->Debug_Level, zz->Proc,
                       zz->Debug_Proc);

  if (use_mach_desc > 0) {
    /* If zz->Machine_Desc already exists, don't rebuild it
     * unless USE_MACHINE_DESC has been set to 2. 
     */

    if ((zz->Machine_Desc == NULL) || (use_mach_desc==2)){
      /* Read machine description from file. 
       * Use Zoltan_Get_Processor_Name to extract the sub-machine
       * on which this Zoltan structure is running. 
       * Broadcast the machine structure to all procs.
       */
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Sorry, heterogeneous load-balancing "
                                  "is still under development!");
      ierr = ZOLTAN_WARN;
    }
  }
  else {
    zz->Machine_Desc = NULL;
  }

  return ierr;
}

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

int Zoltan_Free_Machine_Desc(MachineType **desc)
{
  if (desc && *desc){
    ZOLTAN_FREE(&((*desc)->type)); 
    ZOLTAN_FREE(&((*desc)->xadj)); 
    ZOLTAN_FREE(&((*desc)->adjncy)); 
    ZOLTAN_FREE(&((*desc)->adj_band)); 
    ZOLTAN_FREE(&((*desc)->adj_lat)); 

    ZOLTAN_FREE(desc);
  }

  return ZOLTAN_OK;
}

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

int Zoltan_Copy_Machine_Desc(MachineType **to, MachineType const *from)
{
  MachineType *md = NULL;

  if (*to != NULL) {
    Zoltan_Free_Machine_Desc(to);
  }

  if (!from){
    return ZOLTAN_OK;
  }

  md = *to = (MachineType *)ZOLTAN_MALLOC(sizeof(MachineType));

  if (!md){
    return ZOLTAN_MEMERR;
  }

  *md = *from;

  if (from->ntypes == 1){
    md->type = (int *)malloc(sizeof(int));
    md->type[0] = from->type[0];
  }
  else if (from->ntypes > 1){
    md->type = (int *)malloc(sizeof(int) * from->nnodes);
    memcpy(md->type, from->type, sizeof(int) * from->nnodes);
  }
  else{
    md->type = NULL;
  }

  /*
  ** The following fields are not being used, but should
  ** be copied when they are implemented.
  FUTURE_COPY(xadj);
  FUTURE_COPY(adjncy);
  FUTURE_COPY(adj_band);
  FUTURE_COPY(adj_lat);
  */

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
