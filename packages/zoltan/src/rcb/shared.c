// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_sort.h"
#include "zz_util_const.h"
#include "rcb.h"
#include "rib.h"
#include "par_median_const.h"
#include "all_allo_const.h"
#include "create_proc_list_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* PROTOTYPES */

static int initialize_dot( ZZ *, int , ZOLTAN_ID_PTR , ZOLTAN_ID_PTR , float **, int **,
  struct Dot_Struct *, int *, int , int , int);

static int reallocate_dot_structure(struct Dot_Struct *dots, int newsize);
void Zoltan_Free_And_Reset_Dot_Structure(struct Dot_Struct *dots);

static int send_receive_weights(double *c, int dim, int outgoing, int total, char *sendbuf,
   int firstStaying, int *reorder, int *dotmark, int set, ZOLTAN_COMM_OBJ *cobj, int message_tag);
static int send_receive_doubles(double *c, int outgoing, int total, char *sendbuf,
   int firstStaying, int *reorder, int *dotmark, int set, ZOLTAN_COMM_OBJ *cobj, int message_tag);
static int send_receive_ints(int *c, int outgoing, int total, char *sendbuf,
   int firstStaying, int *reorder, int *dotmark, int set, ZOLTAN_COMM_OBJ *cobj, int message_tag);
static int send_receive_ids(ZOLTAN_ID_TYPE *c, int num_ids, int outgoing, int total, char *sendbuf,
          int firstStaying, int *reorder, int *dotmark, int set, ZOLTAN_COMM_OBJ *cobj, int message_tag);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RB_Build_Structure(
  ZZ *zz,                       /* Zoltan structure */
  ZOLTAN_ID_PTR *global_ids,    /* pointer to array of global IDs; allocated
                                   in this function.  */
  ZOLTAN_ID_PTR *local_ids,     /* pointer to array of local IDs; allocated
                                   in this function.  */
  struct Dot_Struct *dots,       /* pointer to dot info */
  int *num_obj,                 /* number of objects on this processor. */
  int *max_obj,                 /* number of Dots for which storage is 
                                   allocated on this processor. */
  int *num_geom,                /* # values per object used to describe
                                   the geometry.                       */
  int wgtflag,                  /* number of weights per dot. */
  double overalloc,             /* amount to overallocate by when realloc
                                   of dot array must be done.
                                   1.0 = no extra; 1.5 = 50% extra; etc. */
  int use_ids,                  /* true if global and local IDs are to be
                                   kept for RCB or RIB.  In all cases, the 
                                   IDs are allocated and used in Zoltan 
                                   callback functions.  Then, if use_ids is
                                   false, the IDs are deallocated before
                                   this function exits.                */
  int add_unit_weight    /* if wgtflag==0, assume a unit weight for each dot */
)
{
/*
 *  Function to build the geometry-based data structures for 
 *  RCB and RIB.
 */
char *yo = "Zoltan_RB_Build_Structure";
float *objs_wgt = NULL;               /* Array of object weights returned by 
                                         the application.                    */
int *parts = NULL;
int ierr = ZOLTAN_OK;

  /*
   * Allocate space for objects.  Get object info.
   */
  *global_ids = NULL;
  *local_ids = NULL;
  ierr = Zoltan_Get_Obj_List(zz, num_obj, global_ids, local_ids, wgtflag, 
                             &objs_wgt, &parts);
  if (ierr) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Error returned from user function Zoltan_Get_Obj_List.");
    goto End;
  }

  /* Allow extra space for objects that are imported to the processor. */
  *max_obj = (int)(overalloc * *num_obj) + 1;
  *global_ids = ZOLTAN_REALLOC_GID_ARRAY(zz, *global_ids, (*max_obj));
  *local_ids  = ZOLTAN_REALLOC_LID_ARRAY(zz, *local_ids, (*max_obj));

  /* initialize dot frees the objs_wgt array and the parts array */

  ierr = initialize_dot(zz, *num_obj, *global_ids, *local_ids, &objs_wgt, &parts, dots,
                         num_geom, wgtflag, add_unit_weight, *max_obj);

  if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Error returned from user function initialize_dot.");
    goto End;
  }

End:

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    Zoltan_Multifree(__FILE__, __LINE__, 2, global_ids, 
                                            local_ids);
    Zoltan_Free_And_Reset_Dot_Structure(dots);
  }

  if (!use_ids) {
    /* 
     * Do not need IDs in RCB or RIB.  Don't store, manipulate or
     * communicate any IDs.  
     */
    ZOLTAN_FREE(global_ids);
    ZOLTAN_FREE(local_ids);
  }

  return(ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int initialize_dot(
  ZZ *zz, 
  int num_obj,
  ZOLTAN_ID_PTR gid, 
  ZOLTAN_ID_PTR lid, 
  float **wgt,
  int **parts,
  struct Dot_Struct *dots, 
  int *num_geom,
  int wgtflag, int add_unit_weight,
  int max_obj)
{
/*
 *  Function that initializes the dot data structure for RCB and RIB. 
 *  It uses the global ID, coordinates and weight provided by the application.  
 */
int ierr = ZOLTAN_OK;
int i, j, np, fpart, fail, dimcoord;
double *geom_vec = NULL, *coord;
float *obj_weight;
ZOLTAN_ID_PTR id;
char *yo = "initialize_dot";


  memset(dots, 0, sizeof(struct Dot_Struct));

  /* Coordinates *********************************************************/

  ierr = Zoltan_Get_Coordinates(zz, num_obj, gid, lid, num_geom, &geom_vec);
  if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Error returned from Zoltan_Get_Coordinates.");
    goto End;
  }

  fail = 0;
  dimcoord = *num_geom;

  dots->X = (double *)ZOLTAN_MALLOC(max_obj*sizeof(double));
  if (dots->X == NULL) fail = 1;
  if (!fail && (dimcoord > 1)){
    dots->Y = (double *)ZOLTAN_MALLOC(max_obj*sizeof(double));
    if (dots->Y == NULL) fail = 1;
    if (!fail && (dimcoord > 2)){
      dots->Z = (double *)ZOLTAN_MALLOC(max_obj*sizeof(double));
      if (dots->Z == NULL) fail = 1;
    }
  }

  if (fail){
    ierr = ZOLTAN_MEMERR;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "failure to allocate memory for coordinates");
    goto End;
  }

  for (i = 0, coord=geom_vec; i < num_obj; i++) {
    dots->X[i] = *coord++;
    if (dimcoord > 1){
      dots->Y[i] = *coord++;
      if (dimcoord > 2){
        dots->Z[i] = *coord++;
      }
    }
  }

  ZOLTAN_FREE(&geom_vec);

  /* Weights   *********************************************************/

  if ((wgtflag==0) && add_unit_weight){
    dots->uniformWeight = 1.0;
  }
  else if (wgtflag > 0){
    dots->nWeights = wgtflag;
  }

  if (dots->nWeights > 0){

    dots->Weight  = (double *)ZOLTAN_MALLOC(dots->nWeights*max_obj*sizeof(double));
    if (!dots->Weight){
      ierr = ZOLTAN_MEMERR;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "failure to allocate memory for weights");
      goto End;
    }

    obj_weight = *wgt;
    for (i = 0; i < num_obj * dots->nWeights; i++) {
      dots->Weight[i] = *obj_weight++;
    }
    ZOLTAN_FREE(wgt);
  }

  /* The initial part *********************************************************/

  if ((dots->Input_Part = (int *)ZOLTAN_MALLOC(max_obj*sizeof(int))) == NULL){
    ierr = ZOLTAN_MEMERR;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "failure to allocate memory for Input_Part");
    goto End;
  }

  for (i = 0; i < num_obj; i++) {
    dots->Input_Part[i] = (*parts)[i];
  }

  ZOLTAN_FREE(parts);

  /* Object migration sizes **************************************************/

  if ((zz->Get_Obj_Size_Multi) || (zz->Get_Obj_Size)) 
    i = 1;
  else
    i = 0;

  MPI_Allreduce(&i, &j, 1, MPI_INT, MPI_MAX, zz->Communicator);

  if (j == 1){    /* At least one process has defined Obj_Size queries */

    if ((dots->Size = (int *)ZOLTAN_MALLOC(max_obj*sizeof(int))) == NULL){
      ierr = ZOLTAN_MEMERR;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "failure to allocate memory for Size");
      goto End;
    }
  
    if (zz->Get_Obj_Size_Multi) {
      zz->Get_Obj_Size_Multi(zz->Get_Obj_Size_Multi_Data,
                             zz->Num_GID, zz->Num_LID, num_obj,
                             gid, lid, dots->Size, &ierr);
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                        "ZOLTAN_OBJ_SIZE_MULTI function.");
        goto End;
      }
    }
    else if (zz->Get_Obj_Size) {
      for (i = 0; i < num_obj; i++) {
        id = (zz->Num_LID ? &(lid[i*zz->Num_LID]):NULL);
        dots->Size[i] = zz->Get_Obj_Size(zz->Get_Obj_Size_Data,
                                       zz->Num_GID, zz->Num_LID,
                                         &(gid[i*zz->Num_GID]),
                                         id, &ierr);
        if (ierr < 0) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                          "ZOLTAN_OBJ_SIZE function.");
          goto End;
        }
      }
    }
    else if (num_obj){
      for (i=0; i < num_obj; i++){
        dots->Size[i] = 1;  /* default object migration size */
      }
    }
  }

  /* The rest   *********************************************************/

  dots->Proc  = (int *)ZOLTAN_MALLOC(max_obj*sizeof(int));
  if (!dots->Proc) fail = 1;

  if (!fail){
    dots->Part = (int *)ZOLTAN_MALLOC(max_obj*sizeof(int));
    if (!dots->Part) fail = 1;
  }

  if (fail){
    ierr = ZOLTAN_MEMERR;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "failure to allocate memory");
    goto End;
  }

  for (i = 0; i < num_obj; i++) {

    dots->Proc[i] = zz->Proc;

    Zoltan_LB_Proc_To_Part(zz, zz->Proc, &np, &fpart);

    if (fpart >= 0)
      dots->Part[i] = fpart;
    else
      dots->Part[i] = 0;
  }

End:
  
  return(ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RB_Send_Outgoing(
  ZZ *zz,                           /* Load-balancing structure. */
  ZOLTAN_ID_PTR *gidpt,             /* pointer to Global_IDs array. */
  ZOLTAN_ID_PTR *lidpt,             /* pointer to Local_IDs array.  */
  struct Dot_Struct *dotpt,        /* pointer to Dots array. */
  int **dotmark,                    /* which side of median for each dot */
  int *dottop,                      /* dots >= this index are new */
  int *dotnum,                      /* number of dots */
  int *dotmax,                      /* max # of dots arrays can hold */
  int  set,                         /* which part processor is in = 0/1 */
  int *allocflag,                   /* have to re-allocate space */
  double overalloc,                 /* amount to overallocate by when realloc
                                       of dot array must be done.
                                       1.0 = no extra; 1.5 = 50% extra; etc. */
  int stats,                        /* Print timing & count summary? */
  ZOLTAN_GNO_TYPE counters[],         /* diagnostic counts
                                       0 = # of median iterations
                                       1 = # of dots sent
                                       2 = # of dots received
                                       3 = most dots this proc ever owns
                                       4 = most dot memory this proc ever allocs
                                       5 = # of times a previous cut is re-used
                                       6 = # of reallocs of dot array */
  int use_ids,                      /* true if global and local IDs are to be
                                       kept for RCB or RIB.  The IDs must be
                                       communicated if use_ids is true.  */
  MPI_Comm local_comm,
  int proclower,                    /* smallest processor for Tflops_Special */
  int numprocs,                     /* number of processors for Tflops_Special*/
  int partlower,                    /* smallest part # in set 0 */
  int partmid                       /* smallest part # in set 1 */
)
{
/* Routine to determine new processors for outgoing dots. */

  char *yo = "Zoltan_RB_Send_Outgoing";
  int keep, outgoing;               /* message exchange counters */
  int *proc_list = NULL;            /* list of processors to send dots to */
  int i, ierr = ZOLTAN_OK;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* outgoing = number of dots to ship to partner */
  /* dottop = number of dots that have never migrated */
  /* Also, update part assignments. */

  for (i = 0, keep = 0, outgoing = 0; i < *dotnum; i++) {
    if ((*dotmark)[i] != set)
      outgoing++;
    else if (i < *dottop)
      keep++;
    if ((*dotmark)[i]) 
      dotpt->Part[i] = partmid;
    else
      dotpt->Part[i] = partlower;
  }
  *dottop = keep;

  if (outgoing)
    if ((proc_list = (int *) ZOLTAN_MALLOC(outgoing*sizeof(int))) == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

  ierr = Zoltan_RB_Create_Proc_List(zz, set, *dotnum, outgoing, proc_list,
                                    local_comm, proclower, numprocs);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                       "Error returned from Zoltan_RB_Create_Proc_List.");
    goto End;
  }
#if 0
  ierr = Zoltan_RB_Send_Dots(zz, gidpt, lidpt, dotpt, dotmark, proc_list, 
                         outgoing, dotnum, dotmax, set, allocflag, overalloc, 
                         stats, counters, use_ids, local_comm);
#else
  ierr = Zoltan_RB_Send_Dots_less_memory(zz, gidpt, lidpt, dotpt, dotmark, proc_list, 
                         outgoing, dotnum, dotmax, set, allocflag, overalloc, 
                         stats, counters, use_ids, local_comm);
#endif

End:

  if (outgoing) ZOLTAN_FREE(&proc_list);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RB_Send_To_Part(
  ZZ *zz,                           /* Load-balancing structure. */
  ZOLTAN_ID_PTR *gidpt,             /* pointer to Global_IDs array. */
  ZOLTAN_ID_PTR *lidpt,             /* pointer to Local_IDs array.  */
  struct Dot_Struct *dotpt,        /* pointer to Dots array. */
  int **dotmark,                    /* which side of median for each dot */
  int *dottop,                      /* dots >= this index are new */
  int *dotnum,                      /* number of dots */
  int *dotmax,                      /* max # of dots arrays can hold */
  int *allocflag,                   /* have to re-allocate space */
  double overalloc,                 /* amount to overallocate by when realloc
                                       of dot array must be done.
                                       1.0 = no extra; 1.5 = 50% extra; etc. */
  int stats,                        /* Print timing & count summary? */
  ZOLTAN_GNO_TYPE counters[],         /* diagnostic counts
                                       0 = # of median iterations
                                       1 = # of dots sent
                                       2 = # of dots received
                                       3 = most dots this proc ever owns
                                       4 = most dot memory this proc ever allocs
                                       5 = # of times a previous cut is re-used
                                       6 = # of reallocs of dot array */
  int use_ids                      /* true if global and local IDs are to be
                                       kept for RCB or RIB.  The IDs must be
                                       communicated if use_ids is true.  */
)
{
/* When parallel partitioning is done, send dots that are on the wrong
 * processor for their part to the correct processor.
 * This situation arises when a processor has zero parts assigned to
 * it, yet has participated in the parallel partitioning and has, as a result,
 * stored some dots.  
 * (e.g., three processors, two parts, NUM_LOCAL_PARTS = 1 on 
 * procs 1 and 2.  Procs 0 and 1 are in set 0 during parallel partitioning, so
 * Proc 0 may have some dots after the parallel partitioning.  Those dots 
 * must be sent to proc 1.
 * NOTE:  This routine changes values in dotmark.
 */
char *yo = "Zoltan_RB_Send_To_Part";
int outtop, outgoing;               /* message exchange counters */
int *proc_list = NULL;            /* list of processors to send dots to */
int i, ierr = ZOLTAN_OK;
int proc = zz->Proc;
int tmp;
int num_gid = zz->Num_GID;
int set = 0;

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (zz->LB.PartDist == NULL)
    return ierr;  /* Check is not needed for uniform k == p */
  
  if (*dotnum > 0) {
    if (!(proc_list = (int *) ZOLTAN_MALLOC(*dotnum * sizeof(int)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }
  
  outtop = 0;
  outgoing = 0;
  for (i = 0; i < *dotnum; i++) {
    tmp = Zoltan_LB_Part_To_Proc(zz, dotpt->Part[i],
                                     (use_ids ? &((*gidpt)[i*num_gid]) : NULL));
    if (tmp != proc) {
      (*dotmark)[i] = 1;
      proc_list[outgoing] = tmp;
      outgoing++;
      if (i < *dottop) outtop++;
    }
    else {
      (*dotmark)[i] = 0;
    }
  }
  *dottop -= outtop;

#if 0
  ierr = Zoltan_RB_Send_Dots(zz, gidpt, lidpt, dotpt, dotmark, proc_list,
                         outgoing, dotnum, dotmax, set, allocflag, overalloc, 
                         stats, counters, use_ids, zz->Communicator);
#else
  ierr = Zoltan_RB_Send_Dots_less_memory(zz, gidpt, lidpt, dotpt, dotmark, proc_list,
                         outgoing, dotnum, dotmax, set, allocflag, overalloc, 
                         stats, counters, use_ids, zz->Communicator);
#endif

End:
  ZOLTAN_FREE(&proc_list);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ierr);
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#if 0
int Zoltan_RB_Send_Dots(
  ZZ *zz,                           /* Load-balancing structure. */
  ZOLTAN_ID_PTR *gidpt,             /* pointer to Global_IDs array. */
  ZOLTAN_ID_PTR *lidpt,             /* pointer to Local_IDs array.  */
  struct Dot_Struct **dotpt,        /* pointer to Dots array. */
  int **dotmark,                    /* which side of median for each dot */
  int *proc_list,                   /* list of processors to send dots to */
  int outgoing,                     /* message exchange counters */
  int *dotnum,                      /* number of dots */
  int *dotmax,                      /* max # of dots arrays can hold */
  int  set,                         /* which part processor is in = 0/1 */
  int *allocflag,                   /* have to re-allocate space */
  double overalloc,                 /* amount to overallocate by when realloc
                                       of dot array must be done.
                                       1.0 = no extra; 1.5 = 50% extra; etc. */
  int stats,                        /* Print timing & count summary? */
  ZOLTAN_GNO_TYPE counters[],         /* diagnostic counts
                                       0 = # of median iterations
                                       1 = # of dots sent
                                       2 = # of dots received
                                       3 = most dots this proc ever owns
                                       4 = most dot memory this proc ever allocs
                                       5 = # of times a previous cut is re-used
                                       6 = # of reallocs of dot array */
  int use_ids,                      /* true if global and local IDs are to be
                                       kept for RCB or RIB (for LB.Return_Lists 
                                       or high Debug_Levels).  The IDs must be
                                       communicated if use_ids is true.  */
  MPI_Comm local_comm
)
{
/* Routine to send outgoing dots to their new processors. */

  char *yo = "Zoltan_RB_Send_Dots";
  int dotnew;                       /* # of new dots after send/recv */
  int keep, incoming;               /* message exchange counters */
  ZOLTAN_ID_PTR gidbuf = NULL;      /* communication buffer for global IDs. */
  ZOLTAN_ID_PTR lidbuf = NULL;      /* communication buffer for local IDs.  */
  struct Dot_Struct *dotbuf = NULL; /* communication buffer for dots. */
  ZOLTAN_COMM_OBJ *cobj = NULL;     /* pointer for communication object */
  int message_tag = 32760;          /* message tag */
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  int i, ierr = ZOLTAN_OK;

  ZOLTAN_TRACE_ENTER(zz, yo);
  incoming = 0;
  ierr = Zoltan_Comm_Create(&cobj, outgoing, proc_list, local_comm, message_tag,
                        &incoming);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Create.");
    goto End;
  }

  /* check if need to malloc more space */

  dotnew = *dotnum - outgoing + incoming;

  if (dotnew > *dotmax) {
    *allocflag = 1;

    *dotmax = (int) (overalloc * dotnew);
    if (*dotmax < dotnew) *dotmax = dotnew;

    if (use_ids) {
      *gidpt = ZOLTAN_REALLOC_GID_ARRAY(zz, *gidpt, *dotmax);
      *lidpt = ZOLTAN_REALLOC_LID_ARRAY(zz, *lidpt, *dotmax);
      if (!*gidpt || (num_lid_entries && !*lidpt)) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        ierr = ZOLTAN_MEMERR;
        goto End;
      }
    }

    *dotpt = (struct Dot_Struct *) ZOLTAN_REALLOC(*dotpt, *dotmax * sizeof(struct Dot_Struct));
    *dotmark = (int *) ZOLTAN_REALLOC(*dotmark, *dotmax * sizeof(int));
    if (!*dotpt || !*dotmark) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    if (stats) counters[6]++;
  }

  if (stats) {
    counters[1] += outgoing;
    counters[2] += incoming;
    if (dotnew > counters[3])  counters[3] = dotnew;
    if (*dotmax > counters[4]) counters[4] = *dotmax;
  }
    
  /* malloc comm send buffer */

  if (outgoing > 0) {
    if (use_ids) {
      gidbuf = ZOLTAN_MALLOC_GID_ARRAY(zz, outgoing);
      lidbuf = ZOLTAN_MALLOC_LID_ARRAY(zz, outgoing);
    }
    dotbuf = (struct Dot_Struct *)
              ZOLTAN_MALLOC(outgoing * sizeof(struct Dot_Struct));
    if (!dotbuf || (use_ids && (!gidbuf || (num_lid_entries && !lidbuf)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }

  /* fill buffer with dots that are marked for sending */
  /* pack down the unmarked ones */
    
  keep = outgoing = 0;
  for (i = 0; i < *dotnum; i++) {
    if ((*dotmark)[i] != set) {
      if (use_ids) {
        ZOLTAN_SET_GID(zz, &(gidbuf[outgoing*num_gid_entries]), 
                       &((*gidpt)[i*num_gid_entries]));
        ZOLTAN_SET_LID(zz, &(lidbuf[outgoing*num_lid_entries]), 
                       &((*lidpt)[i*num_lid_entries]));
      }
      memcpy((char *) &dotbuf[outgoing], (char *) &((*dotpt)[i]), 
             sizeof(struct Dot_Struct));
      outgoing++;
    }
    else {
      if (use_ids) {
        ZOLTAN_SET_GID(zz, &((*gidpt)[keep*num_gid_entries]), 
                       &((*gidpt)[i*num_gid_entries]));
        ZOLTAN_SET_LID(zz, &((*lidpt)[keep*num_lid_entries]), 
                       &((*lidpt)[i*num_lid_entries]));
      }
      if (keep != i)
        memcpy((char *) &((*dotpt)[keep]), (char *) &((*dotpt)[i]), 
               sizeof(struct Dot_Struct));
      keep++;
    }
  }

  if (use_ids) {
    /* Communicate Global IDs */
    ierr = Zoltan_Comm_Do(cobj, message_tag, (char *) gidbuf, 
                      sizeof(ZOLTAN_ID_TYPE)*num_gid_entries,
                      (char *) &((*gidpt)[keep*num_gid_entries]));
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
      goto End;
    }

    /* Communicate Local IDs, if any */
    if (num_lid_entries) {
      message_tag--;
      ierr = Zoltan_Comm_Do(cobj, message_tag, (char *) lidbuf, 
                        sizeof(ZOLTAN_ID_TYPE)*num_lid_entries,
                        (char *) &((*lidpt)[keep*num_lid_entries]));
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
        goto End;
      }
    }
  }

  /* Communicate Dots */
  message_tag--;
  ierr = Zoltan_Comm_Do(cobj, message_tag, (char *) dotbuf, 
                    sizeof(struct Dot_Struct), (char *) &((*dotpt)[keep]));
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
    goto End;
  }

  ierr = Zoltan_Comm_Destroy(&cobj);

  *dotnum = dotnew;

End:

  if (use_ids) {
    ZOLTAN_FREE(&gidbuf);
    ZOLTAN_FREE(&lidbuf);
  }
  ZOLTAN_FREE(&dotbuf);
    
  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ierr);
}
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#define COMM_DO_ERROR \
    { \
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do."); \
      goto End; \
    } \

int Zoltan_RB_Send_Dots_less_memory(
  ZZ *zz,                           /* Load-balancing structure. */
  ZOLTAN_ID_PTR *gidpt,             /* pointer to Global_IDs array. */
  ZOLTAN_ID_PTR *lidpt,             /* pointer to Local_IDs array.  */
  struct Dot_Struct *dotpt,         /* pointer to Dots info. */
  int **dotmark,                    /* which side of median for each dot */
  int *proc_list,                   /* list of processors to send dots to; 
                                       GETS REORDERED*/
  int outgoing,                     /* message exchange counters */
  int *dotnum,                      /* number of dots */
  int *dotmax,                      /* max # of dots arrays can hold */
  int  set,                         /* which part processor is in = 0/1 */
  int *allocflag,                   /* have to re-allocate space */
  double overalloc,                 /* amount to overallocate by when realloc
                                       of dot array must be done.
                                       1.0 = no extra; 1.5 = 50% extra; etc. */
  int stats,                        /* Print timing & count summary? */
  ZOLTAN_GNO_TYPE counters[],         /* diagnostic counts
                                       0 = # of median iterations
                                       1 = # of dots sent
                                       2 = # of dots received
                                       3 = most dots this proc ever owns
                                       4 = most dot memory this proc ever allocs
                                       5 = # of times a previous cut is re-used
                                       6 = # of reallocs of dot array */
  int use_ids,                      /* true if global and local IDs are to be
                                       kept for RCB or RIB (for LB.Return_Lists 
                                       or high Debug_Levels).  The IDs must be
                                       communicated if use_ids is true.  */
  MPI_Comm local_comm
)
{
/* Routine to send outgoing dots to their new processors.  
 * On very large problems,
 * we run out of memory when allocating the entire Dot_Struct
 * array of outgoing dots.
 * This is slower but uses the minimum amount of memory.
 */

  char *yo = "Zoltan_RB_Send_Dots_less_memory";
  int dotnew;                       /* # of new dots after send/recv */
  int incoming;               /* message exchange counters */
  int bufsize;
  char *sendbuf = NULL;
  int message_tag = 32760;          /* message tag */
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  int i, j, ierr = ZOLTAN_OK;
  int firstStaying;
  int *reorder=NULL;
  ZOLTAN_COMM_OBJ *cobj = NULL;     /* pointer for communication object */
  ZOLTAN_GNO_TYPE verify[4];

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Re-order the data to be sent so that data per proc is contiguous.  This
   * saves memory in Zoltan_Comm_*. The proc_list input array is reordered here.
   */

  if (outgoing > 0){
    reorder = (int *)ZOLTAN_MALLOC(sizeof(int) * outgoing);
    if (outgoing && !reorder){
      MEMORY_ERROR;
    }
  
    firstStaying = *dotnum;       /* location of first dot that is staying */
  
    for (i=0,j=0; i < *dotnum; i++){
      if ((*dotmark)[i] != set) {  /* going */
        reorder[j++] = i;
      }
      else{
        if (i < firstStaying){
          firstStaying = i;
        }
      }
    }
    Zoltan_Comm_Sort_Ints(proc_list, reorder, outgoing);
  }
  else{
    firstStaying = 0;
  }

  /* Create comm object
   * TODO64 - in Zoltan_Comm verify that incoming fits in an integer
   */

  ierr = Zoltan_Comm_Create(&cobj, outgoing, proc_list, local_comm,
                            message_tag, &incoming);

  if (ierr != ZOLTAN_OK) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Create.");
    goto End;
  }

  dotnew = *dotnum - outgoing + incoming;

  if (sizeof(ZOLTAN_GNO_TYPE) > sizeof(int)){
    verify[0] = (ZOLTAN_GNO_TYPE)*dotnum;
    verify[1] = (ZOLTAN_GNO_TYPE)outgoing;
    verify[2] = (ZOLTAN_GNO_TYPE)incoming;
    verify[3] = verify[0] - verify[1] + verify[2];

    if (verify[3] != (ZOLTAN_GNO_TYPE)dotnew){
      ierr = ZOLTAN_FATAL;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
            "Error: Number of incoming objects does not fit in an integer.\n");
      goto End;
    }
  }

  /* check if need to malloc more space */

  if (dotnew > *dotmax) {
    *allocflag = 1;

    *dotmax = (int) (overalloc * dotnew);
    if (*dotmax < dotnew) *dotmax = dotnew;

    if (use_ids) {
      *gidpt = ZOLTAN_REALLOC_GID_ARRAY(zz, *gidpt, *dotmax);
      *lidpt = ZOLTAN_REALLOC_LID_ARRAY(zz, *lidpt, *dotmax);
      if (!*gidpt || (num_lid_entries && !*lidpt)) MEMORY_ERROR;
    }

    ierr = reallocate_dot_structure(dotpt, *dotmax);

    if (ierr != ZOLTAN_OK){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "insufficient memory in reallocate_dot_structure\n");
      goto End;
    }

    *dotmark = (int *) ZOLTAN_REALLOC(*dotmark, *dotmax * sizeof(int));
    if (!*dotmark) MEMORY_ERROR;

    if (stats) counters[6]++;
  }

  if (stats) {
    counters[1] += outgoing;
    counters[2] += incoming;
    if (dotnew > counters[3])  counters[3] = dotnew;
    if (*dotmax > counters[4]) counters[4] = *dotmax;
  }

  /* Figure out the size of the largest single send */

  bufsize = sizeof(double); /* to send coordinates */

  if (dotpt->nWeights > 0){
    bufsize = sizeof(double) * dotpt->nWeights ; /* to send weights */
  }

  if (use_ids){
    j = sizeof(ZOLTAN_ID_TYPE) * num_gid_entries;   /* to send gids */
    if (j > bufsize)
      bufsize = j;
    j = sizeof(ZOLTAN_ID_TYPE) * num_lid_entries;   /* to send lids */
    if (j > bufsize)
      bufsize = j;
  }

  if (outgoing > 0){
    sendbuf = (char *)ZOLTAN_MALLOC(bufsize * outgoing);
    if (!sendbuf) MEMORY_ERROR;
  }

  /***** Send/receive coordinates *****/

  ierr = send_receive_doubles(dotpt->X, outgoing, *dotnum, sendbuf,
                   firstStaying, reorder, *dotmark, set, cobj, message_tag++);

  if ((ierr == ZOLTAN_OK) && dotpt->Y){
    ierr = send_receive_doubles(dotpt->Y, outgoing, *dotnum, sendbuf,
                   firstStaying, reorder, *dotmark, set, cobj, message_tag++);

    if ((ierr == ZOLTAN_OK) && dotpt->Z){
      ierr = send_receive_doubles(dotpt->Z, outgoing, *dotnum, sendbuf,
                     firstStaying, reorder, *dotmark, set, cobj, message_tag++);
    }
  }

  if (ierr != ZOLTAN_OK) COMM_DO_ERROR;

  if (use_ids){

    /***** Send/receive global IDs *****/
    ierr = send_receive_ids(*gidpt, num_gid_entries, outgoing, *dotnum, sendbuf,
                    firstStaying, reorder, *dotmark, set, cobj, message_tag++);

    if (ierr != ZOLTAN_OK) COMM_DO_ERROR;

    if (num_lid_entries){

      /***** Send/receive local IDs *****/
      ierr = send_receive_ids(*lidpt, num_lid_entries, outgoing, *dotnum,
                              sendbuf, firstStaying, reorder, *dotmark, set,
                              cobj, message_tag++);

      if (ierr != ZOLTAN_OK) COMM_DO_ERROR;
    }
  }

  /***** Send/receive weights *****/

  if (dotpt->nWeights > 0){

    ierr = send_receive_weights(dotpt->Weight, dotpt->nWeights,
                                outgoing, *dotnum, sendbuf,
                                firstStaying, reorder, *dotmark, set,
                                cobj, message_tag++);
  
    if (ierr != ZOLTAN_OK) COMM_DO_ERROR;
  }

  /***** Send/receive process owning the dot *****/

  ierr = send_receive_ints(dotpt->Proc, outgoing, *dotnum, sendbuf,
                           firstStaying, reorder, *dotmark, set,
                           cobj, message_tag++);
  
  if (ierr != ZOLTAN_OK) COMM_DO_ERROR;

  /***** Send/receive dot's original part *****/

  ierr = send_receive_ints(dotpt->Input_Part, outgoing, *dotnum, sendbuf,
                           firstStaying, reorder, *dotmark, set,
                           cobj, message_tag++);
  
  if (ierr != ZOLTAN_OK) COMM_DO_ERROR;

  /***** Send/receive dot's new part *****/

  ierr = send_receive_ints(dotpt->Part, outgoing, *dotnum, sendbuf,
                           firstStaying, reorder, *dotmark, set,
                           cobj, message_tag++);
  
  if (ierr != ZOLTAN_OK) COMM_DO_ERROR;

  /***** Send/receive dot's migration size *****/

  if (dotpt->Size){
    ierr = send_receive_ints(dotpt->Size, outgoing, *dotnum, sendbuf,
                             firstStaying, reorder, *dotmark, set,
                             cobj, message_tag++);
  
    if (ierr != ZOLTAN_OK) COMM_DO_ERROR;
  }

  *dotnum = dotnew;

End:

  Zoltan_Comm_Destroy(&cobj);
  ZOLTAN_FREE(&sendbuf);
  ZOLTAN_FREE(&reorder);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Zoltan_RB_Print_All(
                     ZZ *zz, ZOLTAN_ID_PTR global_ids, struct Dot_Struct *dots,
                     int dotnum,
                     int num_import, ZOLTAN_ID_PTR import_global_ids, 
                     int *import_procs)
{
/*
 * Routine to print debugging information.  This routine runs serially
 * over all processors (due to Zoltan_Print_Sync_*) and thus should be used
 * only for debugging.
 */
int kk;
int num_gid_entries = zz->Num_GID;

  Zoltan_Print_Sync_Start(zz->Communicator, TRUE);
  printf("ZOLTAN Proc %d Num_Obj=%d Num_Non_Local=%d\n",
         zz->Proc, dotnum, num_import);
  printf("  Assigned objects:\n");
  for (kk = 0; kk < dotnum; kk++) {
     printf("    Obj:  ");
     ZOLTAN_PRINT_GID(zz, &(global_ids[kk*num_gid_entries]));
     printf("  Orig: %4d\n", dots->Proc[kk]);
  }
  printf("  Non_locals:\n");
  for (kk = 0; kk < num_import; kk++) {
     printf("    Obj:  ");
     ZOLTAN_PRINT_GID(zz, &(import_global_ids[kk*num_gid_entries]));
     printf("     Orig: %4d\n", import_procs[kk]);
  }
  Zoltan_Print_Sync_End(zz->Communicator, TRUE);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RB_Remap(
  ZZ *zz,                           /* Load-balancing structure. */
  ZOLTAN_ID_PTR *gidpt,             /* pointer to Global_IDs array. */
  ZOLTAN_ID_PTR *lidpt,             /* pointer to Local_IDs array.  */
  struct Dot_Struct *dotpt,        /* pointer to Dots array. */
  int *dotnum,                      /* number of dots */
  int *dotmax,                      /* max # of dots arrays can hold */
  int *allocflag,                   /* have to re-allocate space */
  double overalloc,                 /* amount to overallocate by when realloc
                                       of dot array must be done.
                                       1.0 = no extra; 1.5 = 50% extra; etc. */
  int stats,                        /* Print timing & count summary? */
  ZOLTAN_GNO_TYPE counters[],                   /* diagnostic counts
                                       0 = # of median iterations
                                       1 = # of dots sent
                                       2 = # of dots received
                                       3 = most dots this proc ever owns
                                       4 = most dot memory this proc ever allocs
                                       5 = # of times a previous cut is re-used
                                       6 = # of reallocs of dot array */
  int use_ids                       /* true if global and local IDs are to be
                                       kept for RCB or RIB (for LB.Return_Lists 
                                       or high Debug_Levels).  The IDs must be
                                       communicated if use_ids is true.  */
)
{
char *yo = "Zoltan_RB_Remap";
int *old_part = NULL;    /* Array of old part assignments for dots */
int *new_part = NULL;    /* Array of new part assignments for dots;
                            initially determined by partitioning algorithm;
                            may be reset by Zoltan_LB_Remap.  */
int *proc = NULL;        /* Array of processor assignments for dots;
                            initially set to input processors;
                            may be reset by Zoltan_LB_Remap to new processors
                            in Remap. */
int *proc_list = NULL;   /* Temporary array of destinations for dots to be sent
                            due to remap. */
int outgoing;            /* Number of dots to be sent due to remap. */
int new_map;             /* Flag indicating whether a remap has occurred. */
int ierr = ZOLTAN_OK;
int i; 

  old_part = (int *) ZOLTAN_MALLOC(2 * *dotnum * sizeof(int));
  new_part = old_part + *dotnum;
  proc = (int *) ZOLTAN_MALLOC(*dotmax * sizeof(int));
  if (*dotnum && (!old_part || !proc)) {
    /* Can't do remapping, but decomposition should be OK. */
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Memory error.");
    ierr = ZOLTAN_WARN;
    goto End;
  }

  for (i = 0; i < *dotnum; i++) {
    old_part[i] = dotpt->Input_Part[i];
    new_part[i] = dotpt->Part[i];
    proc[i] = dotpt->Proc[i];
  }

  /* Remap parts to reduce data movement. */
  ierr = Zoltan_LB_Remap(zz, &new_map, *dotnum, proc, old_part, new_part, 0);
  if (ierr < 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_LB_Remap");
    goto End;
  }
  
  if (new_map) {
    /* Parts are being remapped; need to move the dots to new procs. */
    for (i = 0; i < *dotnum; i++)
      dotpt->Part[i] = new_part[i];

    if (zz->LB.Return_Lists) {
      /* Call send dots for all dots.  Zoltan_Comm will handle self messages
       * for dots not moving.  Use proc for dotmark and Proc for set; this will
       * cause all dots not remapped to Proc to be sent.
       */
      proc_list = (int *) ZOLTAN_MALLOC(*dotnum * sizeof(int));
      if (*dotnum && !proc_list) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory Error");
        goto End;
      }
      for (outgoing = 0, i = 0; i < *dotnum; i++) 
        if (proc[i] != zz->Proc)
          proc_list[outgoing++] = proc[i];

#if 0
      ierr = Zoltan_RB_Send_Dots(zz, gidpt, lidpt, dotpt, &proc, proc_list,
                                 outgoing, dotnum, dotmax, zz->Proc, allocflag,
                                 overalloc, stats, counters, use_ids,
                                 zz->Communicator);
#else
      ierr = Zoltan_RB_Send_Dots_less_memory(zz, gidpt, lidpt, dotpt, &proc, proc_list,
                                 outgoing, dotnum, dotmax, zz->Proc, allocflag,
                                 overalloc, stats, counters, use_ids,
                                 zz->Communicator);
#endif
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                           "Error returned from Zoltan_RB_Send_Dots");
        goto End;
      }
    }
  }

End:

  ZOLTAN_FREE(&old_part);
  ZOLTAN_FREE(&proc);
  ZOLTAN_FREE(&proc_list);
  return ierr;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RB_Return_Arguments(
  ZZ *zz,                            /* Load-balancing structure */
  ZOLTAN_ID_PTR gidpt,               /* pointer to array of global IDs. */
  ZOLTAN_ID_PTR lidpt,               /* pointer to array of local IDs. */
  struct Dot_Struct *dotpt,          /* pointer to array of Dots. */
  int *num_import,                   /* number of objects to be imported. */
  ZOLTAN_ID_PTR *import_global_ids,  /* global IDs of objects to be imported. */
  ZOLTAN_ID_PTR *import_local_ids,   /* local IDs of objects to be imported. */
  int **import_procs,                /* processors from which objects will be 
                                        imported. */
  int **import_to_part,              /* parts to which objects will be 
                                        imported. */
  int dotnum                         /* number of dots on this processor */
)
{
/*
 * Function to build the return arguments expected by Zoltan.
 * Allocates, fills and returns import_global_ids, import_local_ids, and
 * import_procs.
 */
char *yo = "Zoltan_RB_Return_Arguments";
int i, j;
int ierr = ZOLTAN_OK;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;

  ZOLTAN_TRACE_ENTER(zz, yo);
  /* Compute number of objects to import.  Include those that change only
     part but not processor. */

  *num_import = 0;
  for (i = 0; i < dotnum; i++)  
    if (dotpt->Proc[i] != zz->Proc ||   /* imported from other processors */
        dotpt->Input_Part[i] != dotpt->Part[i])   /* part change only */
      (*num_import)++;

  *import_global_ids = *import_local_ids = NULL;
  *import_procs = *import_to_part = NULL;

  if (*num_import > 0) {
    if (!Zoltan_Special_Malloc(zz, (void **)import_global_ids, *num_import,
                               ZOLTAN_SPECIAL_MALLOC_GID)
     || !Zoltan_Special_Malloc(zz, (void **)import_local_ids, *num_import,
                               ZOLTAN_SPECIAL_MALLOC_LID)
     || !Zoltan_Special_Malloc(zz, (void **)import_procs, *num_import,
                               ZOLTAN_SPECIAL_MALLOC_INT)
     || !Zoltan_Special_Malloc(zz, (void **)import_to_part, *num_import,
                               ZOLTAN_SPECIAL_MALLOC_INT)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    for (i = 0, j = 0; j < dotnum; j++) {
      if (dotpt->Proc[j] != zz->Proc ||    /* imported from other processors */
          dotpt->Input_Part[j] != dotpt->Part[j]) {  /* part change only */
        ZOLTAN_SET_GID(zz, &((*import_global_ids)[i*num_gid_entries]),
                           &(gidpt[j*num_gid_entries]));
        if (num_lid_entries)
          ZOLTAN_SET_LID(zz, &((*import_local_ids)[i*num_lid_entries]),
                             &(lidpt[j*num_lid_entries]));
        (*import_procs)[i] = dotpt->Proc[j];
        (*import_to_part)[i] = dotpt->Part[j];
        i++;
      }
    }
  }


End:
  if (ierr < 0) {
    Zoltan_Special_Free(zz, (void **)import_global_ids,
                        ZOLTAN_SPECIAL_MALLOC_GID);
    Zoltan_Special_Free(zz, (void **)import_local_ids,
                        ZOLTAN_SPECIAL_MALLOC_LID);
    Zoltan_Special_Free(zz, (void **)import_procs,
                        ZOLTAN_SPECIAL_MALLOC_INT);
    Zoltan_Special_Free(zz, (void **)import_to_part,
                        ZOLTAN_SPECIAL_MALLOC_INT);
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RB_Use_IDs(ZZ *zz)
{
/* Function that returns a flag indicating whether or not global and
 * local IDs should be stored, manipulated, and communicated within the
 * RCB and RIB algorithms.  
 * IDs are used for two purposes in RCB and RIB:
 * 1.  To build import lists.  When LB.Return_Lists == NONE, import lists are 
 * not built, so there is no need to carry along the IDs,
 * 2.  To be printed as debugging information.  Only for high debug levels are
 * the IDs printed; for lower levels, the IDs need not be carried along.
 */
  return (zz->LB.Return_Lists || (zz->Debug_Level >= ZOLTAN_DEBUG_ALL));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int Zoltan_RB_Tree_Gatherv(
  ZZ *zz,
  int size,              /* # bytes for rcb_tree or rib_tree */
  int *sendcount,        /* Output:  # bytes to be sent by current proc */
  int *recvcount,        /* Output:  array of # bytes to be received from each
                                     processor. */
  int *displ             /* Output:  array of displacements into recvbuf
                                     for data from each processor. */
) 
{
/* Routine to build arguments to MPI_Allgatherv for consolidating trees
   for KEEP_CUTS=1.  Scans parts array to assign unique processor owners to
   parts; this is important when a part is spread across multiple
   processors.
   Can't do more in common due to different tree types
   for RCB and RIB.  */

int np;       /* number of parts on a processor */
int fp;       /* first part on a processor */
int prev_fp;  /* previous first part.  If fp == prev_fp, part fp
                 is spread across multiple processors.  Include its entry
                 only on the first processor holding fp. */
int ierr = ZOLTAN_OK;
int i, cnt;

  *sendcount = 0;
  prev_fp = -1;
  cnt = 0;
  for (i = 0; i < zz->Num_Proc; i++) {
    recvcount[i] = 0;
    ierr = Zoltan_LB_Proc_To_Part(zz, i, &np, &fp);
    if (ierr >= 0) {
      if (np > 0 && fp != prev_fp) {
        /* processor i has some parts, and the parts start on proc i */
        if (i == zz->Proc) *sendcount = np * size;
        recvcount[i] = np * size;
        displ[i] = cnt * size;
        cnt += np;
        prev_fp = fp;
      }
      else {
        recvcount[i] = 0;
        displ[i] = 0;
      }
    }     
    else if (ierr < 0) 
      break;
  }
  return ierr;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int Zoltan_RB_check_geom_input(
  ZZ *zz,
  struct Dot_Struct *dotpt,
  int dotnum
)
{
/* Routine to check input to geometric methods for consistency. */
  char *yo = "Zoltan_RB_check_geom_input";
  int i, k, count;
  char msg[256];
  int proc = zz->Proc;
  int ierr = ZOLTAN_OK;

  /* Error check the weights. */
  count = 0;
  for (i = 0; i < dotnum * dotpt->nWeights; i++)
    if (dotpt->Weight[i] < 0.0) 
        count++;
  MPI_Allreduce(&count,&k,1,MPI_INT,MPI_SUM,zz->Communicator);
  if (k > 0) {
    if (proc == 0) {
      sprintf(msg, "%d dot weights are < 0.",k);
      ZOLTAN_PRINT_ERROR(proc, yo, msg);
    }
    ierr = ZOLTAN_FATAL;
  }
  return(ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RB_check_geom_output(
  ZZ *zz, 
  struct Dot_Struct *dotpt,
  float *part_sizes,
  int np,               /* number of parts on processor */
  int fp,               /* first part on processor */
  int dotnum,
  int dotorig,
  void *rcbbox_arg)
{
/* Routine to check output of geometric methods for consistency. */

  char *yo = "Zoltan_RB_check_geom_output";
  char msg[1024];
  int dd=0,i,iflag,proc,nprocs,input[2],total[2];
  double *wtsum,tolerance;
  struct rcb_box *rcbbox = (struct rcb_box *) rcbbox_arg;
  int ierr = ZOLTAN_OK;
  int ngp = zz->LB.Num_Global_Parts;
  double *wtpp = NULL;
  double *x, *y, *z;
  int wtdim = (zz->Obj_Weight_Dim > 0 ? zz->Obj_Weight_Dim : 1);

  ZOLTAN_TRACE_ENTER(zz, yo);
  MPI_Comm_rank(zz->Communicator,&proc);
  MPI_Comm_size(zz->Communicator,&nprocs);

  /* check that total # of dots remained the same */

  input[0] = dotorig;  input[1] = dotnum;
  MPI_Allreduce(input,total,2,MPI_INT,MPI_SUM,zz->Communicator);
  if (total[0] != total[1]) {
    if (proc == 0) {
      sprintf(msg, "Points before partitioning = %d, "
                   "Points after partitioning = %d.",
                    total[0],total[1]);
      ZOLTAN_PRINT_WARN(proc, yo, msg);
      ierr = ZOLTAN_WARN;
    }
  }
  
  /* check that result is within Imbalance_Tol of part size target */

  wtpp = (double *) ZOLTAN_MALLOC(2*(1+ngp)*sizeof(double));

  if (dotpt->Weight){
    for (dd = 0; dd < wtdim; dd++) {
      memset(wtpp, 0, 2*(1+ngp)*sizeof(double));
      for (i = 0; i < dotnum; i++) {
        wtpp[ngp] += dotpt->Weight[i*wtdim + dd];
        wtpp[dotpt->Part[i]] += dotpt->Weight[i*wtdim + dd];
      }
      wtsum = wtpp + (1+ngp);
  
      MPI_Allreduce(wtpp, wtsum, ngp+1, MPI_DOUBLE, MPI_SUM, zz->Communicator);
  
      for (i = fp; i < fp + np; i++) {
        tolerance = part_sizes[i*wtdim+dd]*wtsum[ngp]*zz->LB.Imbalance_Tol[dd];
        if (wtsum[i] > tolerance) {
          if (zz->Debug_Level > ZOLTAN_DEBUG_NONE) {
            sprintf(msg, 
                    "Weight of part %d = %f > tolerance %f for weight %d.", 
                    i, wtsum[i], tolerance, dd);
            ZOLTAN_PRINT_WARN(proc, yo, msg);
          }
          ierr = ZOLTAN_WARN;
        }
      }
    }
  }
  else{
      memset(wtpp, 0, 2*(1+ngp)*sizeof(double));
      for (i = 0; i < dotnum; i++) {
        wtpp[ngp] += dotpt->uniformWeight;
        wtpp[dotpt->Part[i]] += dotpt->uniformWeight;
      }
      wtsum = wtpp + (1+ngp);
  
      MPI_Allreduce(wtpp, wtsum, ngp+1, MPI_DOUBLE, MPI_SUM, zz->Communicator);
  
      for (i = fp; i < fp + np; i++) {
        tolerance = part_sizes[i*wtdim]*wtsum[ngp]*zz->LB.Imbalance_Tol[0];
        if (wtsum[i] > tolerance) {
          if (zz->Debug_Level > ZOLTAN_DEBUG_NONE) {
            sprintf(msg, 
                    "Weight of part %d = %f > tolerance %f for weight 0.", 
                    i, wtsum[i], tolerance);
            ZOLTAN_PRINT_WARN(proc, yo, msg);
          }
          ierr = ZOLTAN_WARN;
        }
      }
  }

  ZOLTAN_FREE(&wtpp);
  
  if (zz->LB.Method == RCB) {

    /* check that final set of points is inside RCB box of each proc */
  
    iflag = 0;
    x = dotpt->X;
    y = dotpt->Y;
    z = dotpt->Z;
    for (i = 0; i < dotnum; i++) {
      if (       x[i] < rcbbox->lo[0] || x[i] > rcbbox->hi[0] ||
          (y && (y[i] < rcbbox->lo[1] || y[i] > rcbbox->hi[1])) ||
  	  (z && (z[i] < rcbbox->lo[2] || z[i] > rcbbox->hi[2])) ){
        iflag++;
        dd=i;
      }
    }
    if (iflag > 0) {
      sprintf(msg, "\n%d points are out-of-box on proc %d.\n" 
        "Example (%g, %g, %g) is not in (%g, %g) , (%g, %g), (%g, %g)\n",
        iflag, proc, x[dd], (y ? y[dd] : 0.0), (z ? z[dd] : 0.0),
      rcbbox->lo[0], rcbbox->hi[0] , rcbbox->lo[1], rcbbox->hi[1], rcbbox->lo[2], rcbbox->hi[2]);

      ZOLTAN_PRINT_ERROR(proc, yo, msg);
      ierr = ZOLTAN_FATAL;
    }
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Zoltan_RB_stats(ZZ *zz, double timetotal, struct Dot_Struct *dotpt,
                 int dotnum, float *part_sizes,
                 double *timers, ZOLTAN_GNO_TYPE *counters, int stats,
                 ZOLTAN_GNO_TYPE *reuse_count, void *rcbbox_arg, int reuse)

{
  char *yo = "Zoltan_RB_stats";
  int i,proc,nprocs,print_proc;
  ZOLTAN_GNO_TYPE sum, min, max;
  double ave,rsum,rmin,rmax;
  double weight,wttot,wtmin,wtmax;
  struct rcb_box *rcbbox = (struct rcb_box *) rcbbox_arg;
  int numParts, wdim;
  double move, gmove, bal, max_imbal, ib;
  double *lpartWgt = NULL;
  double *gpartWgt = NULL;
  MPI_Datatype zoltan_gno_mpi_type;

  zoltan_gno_mpi_type = Zoltan_mpi_gno_type();

  MPI_Comm_rank(zz->Communicator,&proc);
  MPI_Comm_size(zz->Communicator,&nprocs);
  print_proc = zz->Debug_Proc;
  
  /* Make certain that timetotal is the same everywhere. */
  MPI_Allreduce(&timetotal,&rmax,1,MPI_DOUBLE,MPI_MAX,zz->Communicator);
  timetotal = rmax;

  if (proc == print_proc) 
    printf("Partitioning total time: %g (secs)\n", timetotal);

  wdim = dotpt->nWeights;

  if (stats) {
    /* EBEB Do we need stats inside RCB? LB_Eval can do this better. */
    if (proc == print_proc) printf("Partitioning Statistics:\n");

    MPI_Barrier(zz->Communicator);

    /* distribution info */
    /* multiple weights not supported. */

    if (wdim){
      for (i = 0, weight = 0.0; i < dotnum; i++) weight += dotpt->Weight[i*wdim];
    }
    else{
      weight = dotnum * dotpt->uniformWeight;
    }
  
    MPI_Allreduce(&weight,&wttot,1,MPI_DOUBLE,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&weight,&wtmin,1,MPI_DOUBLE,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,zz->Communicator);

    if (proc == print_proc) {
      printf(" Total weight of dots = %g\n",wttot);
      printf(" Weight on each proc: ave = %g, max = %g, min = %g\n",
	     wttot/nprocs,wtmax,wtmin);
    }

    MPI_Barrier(zz->Communicator);
    if (stats > 1)  
      printf("    Proc %d has weight = %g\n",proc,weight);

    if (wdim){
      for (i = 0, weight = 0.0; i < dotnum; i++) 
        if (dotpt->Weight[i*wdim] > weight) weight = dotpt->Weight[i*wdim];
    }
    else{
      weight = dotpt->uniformWeight; 
    }
    MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,zz->Communicator);
  
    if (proc == print_proc) 
      printf(" Maximum weight of single dot = %g\n",wtmax);

    MPI_Barrier(zz->Communicator);
    if (stats > 1)  
      printf("    Proc %d max weight = %g\n",proc,weight);
  
    /* counter info */
  
    MPI_Allreduce(&counters[1],&sum,1,zoltan_gno_mpi_type,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&counters[1],&min,1,zoltan_gno_mpi_type,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&counters[1],&max,1,zoltan_gno_mpi_type,MPI_MAX,zz->Communicator);
    ave = ((double) sum)/nprocs;
    if (proc == print_proc) 
      printf(" Send count: ave = %g, min = " ZOLTAN_GNO_SPEC ", max = " ZOLTAN_GNO_SPEC "\n",ave,min,max);
    MPI_Barrier(zz->Communicator);
    if (stats > 1) 
      printf("    Proc %d send count = " ZOLTAN_GNO_SPEC "\n",proc,counters[1]);
    
    MPI_Allreduce(&counters[2],&sum,1,zoltan_gno_mpi_type,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&counters[2],&min,1,zoltan_gno_mpi_type,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&counters[2],&max,1,zoltan_gno_mpi_type,MPI_MAX,zz->Communicator);
    ave = ((double) sum)/nprocs;
    if (proc == print_proc) 
      printf(" Recv count: ave = %g, min = " ZOLTAN_GNO_SPEC ", max = " ZOLTAN_GNO_SPEC "\n",ave,min,max);
    MPI_Barrier(zz->Communicator);
    if (stats > 1) 
      printf("    Proc %d recv count = " ZOLTAN_GNO_SPEC "\n",proc,counters[2]);
  
    if (reuse) {
      MPI_Allreduce(&reuse_count[1],&max,1,zoltan_gno_mpi_type,MPI_MAX,zz->Communicator);
      MPI_Allreduce(&reuse_count[1],&sum,1,zoltan_gno_mpi_type,MPI_SUM,zz->Communicator);
      MPI_Allreduce(&reuse_count[1],&min,1,zoltan_gno_mpi_type,MPI_MIN,zz->Communicator);
      ave = ((double) sum)/nprocs;
      if (proc == print_proc) 
        printf(" Presend count: ave = %g, min = " ZOLTAN_GNO_SPEC ", max = " ZOLTAN_GNO_SPEC "\n",ave,min,max);
      MPI_Barrier(zz->Communicator);
      if (stats > 1) 
        printf("    Proc %d presend count = " ZOLTAN_GNO_SPEC "\n",proc,reuse_count[1]);
    
      MPI_Allreduce(&reuse_count[2],&sum,1,zoltan_gno_mpi_type,MPI_SUM,zz->Communicator);
      MPI_Allreduce(&reuse_count[2],&min,1,zoltan_gno_mpi_type,MPI_MIN,zz->Communicator);
      MPI_Allreduce(&reuse_count[2],&max,1,zoltan_gno_mpi_type,MPI_MAX,zz->Communicator);
      ave = ((double) sum)/nprocs;
      if (proc == print_proc) 
        printf(" Prerecv count: ave = %g, min = " ZOLTAN_GNO_SPEC ", max = " ZOLTAN_GNO_SPEC "\n",ave,min,max);
      MPI_Barrier(zz->Communicator);
      if (stats > 1) 
        printf("    Proc %d prerecv count = " ZOLTAN_GNO_SPEC "\n",proc,reuse_count[2]);
    }
  
    MPI_Allreduce(&counters[3],&sum,1,zoltan_gno_mpi_type,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&counters[3],&min,1,zoltan_gno_mpi_type,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&counters[3],&max,1,zoltan_gno_mpi_type,MPI_MAX,zz->Communicator);
    ave = ((double) sum)/nprocs;
    if (proc == print_proc) 
      printf(" Max dots: ave = %g, min = " ZOLTAN_GNO_SPEC ", max = " ZOLTAN_GNO_SPEC "\n",ave,min,max);
    MPI_Barrier(zz->Communicator);
    if (stats > 1) 
      printf("    Proc %d max dots = " ZOLTAN_GNO_SPEC "\n",proc,counters[3]);
  
    MPI_Allreduce(&counters[4],&sum,1,zoltan_gno_mpi_type,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&counters[4],&min,1,zoltan_gno_mpi_type,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&counters[4],&max,1,zoltan_gno_mpi_type,MPI_MAX,zz->Communicator);
    ave = ((double) sum)/nprocs;
    if (proc == print_proc) 
      printf(" Max memory: ave = %g, min = " ZOLTAN_GNO_SPEC ", max = " ZOLTAN_GNO_SPEC "\n",ave,min,max);
    MPI_Barrier(zz->Communicator);
    if (stats > 1) 
      printf("    Proc %d max memory = " ZOLTAN_GNO_SPEC "\n",proc,counters[4]);
    
    if (reuse) {
      MPI_Allreduce(&counters[5],&sum,1,zoltan_gno_mpi_type,MPI_SUM,zz->Communicator);
      MPI_Allreduce(&counters[5],&min,1,zoltan_gno_mpi_type,MPI_MIN,zz->Communicator);
      MPI_Allreduce(&counters[5],&max,1,zoltan_gno_mpi_type,MPI_MAX,zz->Communicator);
      ave = ((double) sum)/nprocs;
      if (proc == print_proc) 
        printf(" # of Reuse: ave = %g, min = " ZOLTAN_GNO_SPEC ", max = " ZOLTAN_GNO_SPEC "\n",ave,min,max);
      MPI_Barrier(zz->Communicator);
      if (stats > 1) 
        printf("    Proc %d # of Reuse = " ZOLTAN_GNO_SPEC "\n",proc,counters[5]);
    }
  
    MPI_Allreduce(&counters[6],&sum,1,zoltan_gno_mpi_type,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&counters[6],&min,1,zoltan_gno_mpi_type,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&counters[6],&max,1,zoltan_gno_mpi_type,MPI_MAX,zz->Communicator);
    ave = ((double) sum)/nprocs;
    if (proc == print_proc) 
      printf(" # of OverAlloc: ave = %g, min = " ZOLTAN_GNO_SPEC ", max = " ZOLTAN_GNO_SPEC "\n",ave,min,max);
    MPI_Barrier(zz->Communicator);
    if (stats > 1) 
      printf("    Proc %d # of OverAlloc = " ZOLTAN_GNO_SPEC "\n",proc,counters[6]);

    par_median_print_counts(zz->Communicator, print_proc);
  }

  /* timer info */
    
  if (timetotal>0){

    MPI_Allreduce(&timers[0],&rsum,1,MPI_DOUBLE,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&timers[0],&rmin,1,MPI_DOUBLE,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&timers[0],&rmax,1,MPI_DOUBLE,MPI_MAX,zz->Communicator);
    ave = rsum/nprocs;
    if (proc == print_proc) {
      printf(" Start-up time (secs): ave = %g, min = %g, max = %g\n",
  	   ave,rmin,rmax);
      printf(" Start-up time (%%): ave = %g, min = %g, max = %g\n",
  	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
    }
    MPI_Barrier(zz->Communicator);
    if (stats > 1) 
      printf("    Proc %d start-up time = %g\n",proc,timers[0]);
    
    MPI_Allreduce(&timers[1],&rsum,1,MPI_DOUBLE,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&timers[1],&rmin,1,MPI_DOUBLE,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&timers[1],&rmax,1,MPI_DOUBLE,MPI_MAX,zz->Communicator);
    ave = rsum/nprocs;
    if (proc == print_proc) {
      printf(" Pre-median time (secs): ave = %g, min = %g, max = %g\n",
  	   ave,rmin,rmax);
      printf(" Pre-median time (%%): ave = %g, min = %g, max = %g\n",
  	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
    }
    MPI_Barrier(zz->Communicator);
    if (stats > 1) 
      printf("    Proc %d pre-median time = %g\n",proc,timers[1]);
    
    MPI_Allreduce(&timers[2],&rsum,1,MPI_DOUBLE,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&timers[2],&rmin,1,MPI_DOUBLE,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&timers[2],&rmax,1,MPI_DOUBLE,MPI_MAX,zz->Communicator);
    ave = rsum/nprocs;
    if (proc == print_proc) {
      printf(" Median time (secs): ave = %g, min = %g, max = %g\n",
  	   ave,rmin,rmax);
      printf(" Median time (%%): ave = %g, min = %g, max = %g\n",
  	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
    }
    MPI_Barrier(zz->Communicator);
    if (stats > 1) 
      printf("    Proc %d median time = %g\n",proc,timers[2]);
    
    MPI_Allreduce(&timers[3],&rsum,1,MPI_DOUBLE,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&timers[3],&rmin,1,MPI_DOUBLE,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&timers[3],&rmax,1,MPI_DOUBLE,MPI_MAX,zz->Communicator);
    ave = rsum/nprocs;
    if (proc == print_proc) {
      printf(" Comm time (secs): ave = %g, min = %g, max = %g\n",
  	   ave,rmin,rmax);
      printf(" Comm time (%%): ave = %g, min = %g, max = %g\n",
  	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
    }
    MPI_Barrier(zz->Communicator);
    if (stats > 1) 
      printf("    Proc %d comm time = %g\n",proc,timers[3]);
  
  }
  
  if (zz->LB.Method == RCB && stats > 1)  {
    /* RCB boxes for each proc */
    if (proc == print_proc) printf(" RCB sub-domain boxes:\n");
    for (i = 0; i < 3; i++) {
      MPI_Barrier(zz->Communicator);
      if (proc == print_proc) printf("    Dimension %d\n",i+1);
      MPI_Barrier(zz->Communicator);
      printf("      Proc = %d: Box = %g %g\n",
	     proc,rcbbox->lo[i],rcbbox->hi[i]);
    }
  }
  /* For comparison, display the values that are displayed by
   * graph and hypergraph partitioners when FINAL_OUTPUT=1
   */

  if (stats) {
    int tmpmax;
    for (i = 0, tmpmax=0; i < dotnum; i++) 
      if (dotpt->Part[i] > tmpmax) tmpmax = dotpt->Part[i];

    MPI_Allreduce(&tmpmax,&numParts,1,MPI_INT, MPI_MAX, zz->Communicator);
    numParts++;

    lpartWgt = (double *)ZOLTAN_CALLOC(numParts, sizeof(double));
    gpartWgt = (double *)ZOLTAN_MALLOC(numParts * sizeof(double));

    if (numParts && (!lpartWgt || !gpartWgt)){
      Zoltan_Multifree(__FILE__, __LINE__, 2, &lpartWgt, &gpartWgt);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return;
    }

    for (i = 0, move=0.0; i < dotnum; i++){
      if (wdim)
        lpartWgt[dotpt->Part[i]] += dotpt->Weight[i*wdim];
      else
        lpartWgt[dotpt->Part[i]] += dotpt->uniformWeight;
 
      if (dotpt->Input_Part[i] != dotpt->Part[i]) {
        if (dotpt->Size){
          move += (double)dotpt->Size[i];
        }
        else{
          move += 1;
        }
      }
    }

    MPI_Reduce(lpartWgt, gpartWgt, numParts, MPI_DOUBLE, MPI_SUM, 
               print_proc, zz->Communicator);
    MPI_Reduce(&move, &gmove, 1, MPI_DOUBLE, MPI_SUM, 
               print_proc, zz->Communicator);

    ZOLTAN_FREE(&lpartWgt);

    if (proc == print_proc) {
      static int nRuns=0;
      static double balsum, balmax, balmin;
      static double movesum, movemax, movemin;
      char *countType;

      max_imbal = 0.0;
   
      if (wttot) {
        for (i = 0; i < numParts; i++){
          if (part_sizes[i]) {
            ib= (gpartWgt[i]-part_sizes[i]*wttot)/(part_sizes[i]*wttot);
            if (ib>max_imbal)
              max_imbal = ib;
          }
        }
      }
  
      bal = 1.0 + max_imbal;
  
      if (nRuns){
        if (gmove > movemax) movemax = gmove;
        if (gmove < movemin) movemin = gmove;
        if (bal > balmax) balmax = bal;
        if (bal < balmin) balmin = bal;
        movesum += gmove;
        balsum += bal;
      }
      else{
        movemax = movemin = movesum = gmove;
        balmax = balmin = balsum = bal;
      }
  
      countType = "moveCnt";
      if (zz->Get_Obj_Size_Multi || zz->Get_Obj_Size) {
        countType = "moveVol";
      }
  
      nRuns++;

      printf(" STATS Runs %d  bal  CURRENT %f  MAX %f  MIN %f  AVG %f\n",
              nRuns, bal, balmax, balmin, balsum/nRuns);
      printf(" STATS Runs %d  %s CURRENT %f  MAX %f  MIN %f  AVG %f\n",
              nRuns, countType, gmove, movemax, movemin, movesum/nRuns);
    }

    {
      int cntmax = 0, cntmin = 0, cntsum = 0;
      MPI_Allreduce(&dotnum, &cntmin, 1, MPI_INT, MPI_MIN, zz->Communicator);
      MPI_Allreduce(&dotnum, &cntmax, 1, MPI_INT, MPI_MAX, zz->Communicator);
      MPI_Allreduce(&dotnum, &cntsum, 1, MPI_INT, MPI_SUM, zz->Communicator);

      if (proc == print_proc) {
        double avg = (double) cntsum / (double) nprocs;
        printf(" STATS DETAIL count:  min %d  max %d  avg %f  imbal %f\n",
               cntmin, cntmax, avg, (double) cntmax / avg);
      }
    }

    if (wdim) {

      double *wtarr = (double *) ZOLTAN_MALLOC(4 * wdim * sizeof(double));
      double *wtarrmin = wtarr    + wdim;
      double *wtarrmax = wtarrmin + wdim;
      double *wtarrsum = wtarrmax + wdim;
      int d;
      for (d = 0; d < 4*wdim; d++) wtarr[d] = 0.;

      for (i = 0; i < dotnum; i++) {
        for (d = 0; d < wdim; d++) {
          wtarr[d] += dotpt->Weight[i*wdim+d];
        }
      }
    
      MPI_Allreduce(wtarr, wtarrmin, wdim, MPI_DOUBLE,MPI_MIN,zz->Communicator);
      MPI_Allreduce(wtarr, wtarrmax, wdim, MPI_DOUBLE,MPI_MAX,zz->Communicator);
      MPI_Allreduce(wtarr, wtarrsum, wdim, MPI_DOUBLE,MPI_SUM,zz->Communicator);

      if (proc == print_proc) {
        for (d = 0; d < wdim; d++) {
          double avg = wtarrsum[d] / (double)nprocs;
          printf(" STATS DETAIL wgt %d:  min %f  max %f  avg %f  imbal %f\n",
                 d, wtarrmin[d], wtarrmax[d], avg, wtarrmax[d]/avg);
        }
      }
      ZOLTAN_FREE(&wtarr);
    }
    else {
      if (proc == print_proc) {
        printf(" STATS DETAIL wdim = %d; no detail available\n", wdim);
      }
    }

    ZOLTAN_FREE(&gpartWgt);
  }
  MPI_Barrier(zz->Communicator);
}

void Zoltan_Free_And_Reset_Dot_Structure(struct Dot_Struct *dots)
{
  ZOLTAN_FREE(&dots->X);
  ZOLTAN_FREE(&dots->Y);
  ZOLTAN_FREE(&dots->Z);
  ZOLTAN_FREE(&dots->Weight);
  ZOLTAN_FREE(&dots->Proc);
  ZOLTAN_FREE(&dots->Input_Part);
  ZOLTAN_FREE(&dots->Part);
  ZOLTAN_FREE(&dots->Size);

  dots->nWeights = dots->uniformWeight = 0;
}

static int reallocate_dot_structure(struct Dot_Struct *dots, int newsize)
{
  int fail = 0;
  int ierr = ZOLTAN_OK;;

  if ((dots->X = (double *) ZOLTAN_REALLOC(dots->X, newsize * sizeof(double))) == NULL)
    fail = 1;

  if (!fail && dots->Y){
    if ((dots->Y = (double *) ZOLTAN_REALLOC(dots->Y, newsize * sizeof(double))) == NULL)
      fail = 1;
    
    if (!fail && dots->Z){
      if ((dots->Z = (double *) ZOLTAN_REALLOC(dots->Z, newsize * sizeof(double))) == NULL)
        fail = 1;
    }
  }

  if (!fail && dots->nWeights){
    if ((dots->Weight = (double *) ZOLTAN_REALLOC(dots->Weight, newsize * dots->nWeights * sizeof(double))) == NULL){
      fail = 1;
    }
  }

  if (!fail)
    if ((dots->Proc = (int *) ZOLTAN_REALLOC(dots->Proc, newsize * sizeof(int))) == NULL)
      fail = 1;

  if (!fail)
    if ((dots->Input_Part = (int *) ZOLTAN_REALLOC(dots->Input_Part, newsize * sizeof(int))) == NULL)
      fail = 1;

  if (!fail)
    if ((dots->Part = (int *) ZOLTAN_REALLOC(dots->Part, newsize * sizeof(int))) == NULL)
      fail = 1;

  if (!fail && dots->Size)
    if ((dots->Size = (int *) ZOLTAN_REALLOC(dots->Size, newsize * sizeof(int))) == NULL)
      fail = 1;
  
  if (fail)
    ierr = ZOLTAN_MEMERR;

  return ierr;
}


static int send_receive_ids(ZOLTAN_ID_TYPE *c, int num_ids, int outgoing, 
                            int total, char *sendbuf, int firstStaying, 
                            int *reorder, int *dotmark, int set, 
                            ZOLTAN_COMM_OBJ *cobj, int message_tag)
{
  int i, j, next;
  ZOLTAN_ID_TYPE *idval, *from;

  if (outgoing > 0){
    idval = (ZOLTAN_ID_TYPE *)sendbuf;

    for (i=0; i < outgoing; i++){
      from = c + (reorder[i] * num_ids);
      for (j=0; j < num_ids; j++)
        *idval++ = *from++;
    }

    idval = c;
 
    for (i = firstStaying, next=0; i < total; i++) {
      if (dotmark[i] == set) {
        from = c + (i * num_ids);
        next++;
        for (j=0; j < num_ids; j++)
          *idval++ = from[j];  /* keepers */
      }
    }
  }
  else{
    next = total;
  }

  return Zoltan_Comm_Do(cobj, message_tag, sendbuf,
                        sizeof(ZOLTAN_ID_TYPE) * num_ids, 
                        (char *)(c + (next * num_ids)));
}

static int send_receive_ints(int *c, int outgoing, int total, char *sendbuf,
          int firstStaying, int *reorder, int *dotmark, int set, ZOLTAN_COMM_OBJ *cobj, int message_tag)
{
  int i, next;
  int *ival;

  if (outgoing > 0){
    ival = (int *)sendbuf;

    for (i=0; i < outgoing; i++){
      *ival++ = c[reorder[i]];
    }

    for (i = firstStaying, next=0; i < total; i++) {
      if (dotmark[i] == set) {
        c[next++] = c[i];  /* keepers */
      }
    }
  }
  else{
    next = total;
  }

  return Zoltan_Comm_Do(cobj, message_tag, sendbuf, sizeof(int),
                        (char *)(c + next));
}

static int send_receive_weights(double *c, int wgtdim, int outgoing, 
                                int total, char *sendbuf,
                                int firstStaying, int *reorder, int *dotmark, 
                                int set, ZOLTAN_COMM_OBJ *cobj, int message_tag)
{
  int i, j, next;
  double *dval;

  if (outgoing > 0){
    dval = (double *)sendbuf;

    for (i=0; i < outgoing; i++){
      for (j=0; j < wgtdim; j++){
        *dval++ = c[reorder[i]*wgtdim + j];
      }
    }

    for (i = firstStaying, next=0; i < total; i++) {
      if (dotmark[i] == set) {
        for (j=0; j < wgtdim; j++){
          c[next++] = c[i*wgtdim + j];  /* keepers */
        }
      }
    }
  }
  else{
    next = total*wgtdim;
  }

  return Zoltan_Comm_Do(cobj, message_tag, sendbuf, sizeof(double)*wgtdim,
                        (char *)(c + next));
}

static int send_receive_doubles(double *c, int outgoing, int total,
                                char *sendbuf, int firstStaying, int *reorder,
                                int *dotmark, int set, ZOLTAN_COMM_OBJ *cobj,
                                int message_tag)
{
  int i, next;
  double *dval;

  if (outgoing > 0){
    dval = (double *)sendbuf;

    for (i=0; i < outgoing; i++){
      *dval++ = c[reorder[i]];
    }

    for (i = firstStaying, next=0; i < total; i++) {
      if (dotmark[i] == set) {
        c[next++] = c[i];  /* keepers */
      }
    }
  }
  else{
    next = total;
  }

  return Zoltan_Comm_Do(cobj, message_tag, sendbuf, sizeof(double),
                        (char *)(c + next));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RB_Candidates_Copy_Input(
  ZZ *zz,                            /* Load-balancing structure */
  int dotnum,                        /* number of dots on this processor */
  ZOLTAN_ID_PTR gidpt,               /* pointer to array of global IDs. */
  ZOLTAN_ID_PTR lidpt,               /* pointer to array of local IDs. */
  struct Dot_Struct *dotpt,          /* pointer to array of Dots. */
  int *num_input,             /* output: number of objs on this proc. */
  ZOLTAN_ID_PTR *input_gids,  /* output: global IDs of objs on this proc. */
  ZOLTAN_ID_PTR *input_lids,  /* output: local IDs of objs on this proc. */
  int **input_procs,          /* output: NULL; not used for CANDIDATE_LISTS */
  int **input_to_part         /* output: NULL; not used for CANDIDATE_LISTS */
)
{
/*
 * For CANDIDATE_LISTS output, the GIDs as input through the callback
 * functions are returned in the import lists.  We copy them now.
 */
char *yo = "Zoltan_RB_Candidates_Copy_Input";
int ierr = ZOLTAN_OK;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;

  ZOLTAN_TRACE_ENTER(zz, yo);

  *num_input = dotnum;
  *input_gids = *input_lids = NULL;
  *input_procs = *input_to_part = NULL;

  if (*num_input && 
     (!Zoltan_Special_Malloc(zz, (void **)input_gids, dotnum,
                               ZOLTAN_SPECIAL_MALLOC_GID)
     || !Zoltan_Special_Malloc(zz, (void **)input_lids, dotnum,
                               ZOLTAN_SPECIAL_MALLOC_LID))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  memcpy(*input_gids, gidpt, dotnum*num_gid_entries*sizeof(ZOLTAN_ID_TYPE));
  memcpy(*input_lids, lidpt, dotnum*num_lid_entries*sizeof(ZOLTAN_ID_TYPE));

End:
  if (ierr < 0) {
    Zoltan_Special_Free(zz, (void **)input_gids, ZOLTAN_SPECIAL_MALLOC_GID);
    Zoltan_Special_Free(zz, (void **)input_lids, ZOLTAN_SPECIAL_MALLOC_LID);
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RB_Candidates_Output(
  ZZ *zz,                          /* Load-balancing structure */
  int dotnum,                      /* number of dots on this processor */
  int *dindx,                      /* Index that provides grouping of dots
                                      by part; non-NULL only if serial_rcb
                                      was invoked (i.e., if there are more
                                      than one parts on a proc). */
  ZOLTAN_ID_PTR gidpt,             /* pointer to array of global IDs currently
                                      on this proc. */
  ZOLTAN_ID_PTR lidpt,             /* pointer to array of local IDs currently
                                      on this proc. */
  struct Dot_Struct *dotpt,        /* pointer to array of Dots currently
                                      on this proc. */
  int num_input,                   /* number of global IDs input on this proc */
  ZOLTAN_ID_PTR input_gids,        /* array of global IDs input on this proc */
  int *num_candidate_gids,         /* output: # of objs input on this proc */
  ZOLTAN_ID_PTR *candidate_gids    /* output: global IDs of candidates for 
                                              each input obj. */
)
{
/*
 * Function to build the return arguments expected by Zoltan when
 * RETURN_LISTS=CANDIDATE_LISTS (used in matching only).
 * For each part, determine a candidate among all dots assigned to the part.
 * Then, for each input GID, return the candidate for its part in the 
 * candidate_gids array.
 */
char *yo = "Zoltan_RB_Candidates_Output";
int i;
int maxdotnum;
int ierr = ZOLTAN_OK;
int num_gid_entries = zz->Num_GID;
struct Zoltan_DD_Struct *dd = NULL;
ZOLTAN_ID_PTR dot_candidates = NULL;
ZOLTAN_ID_PTR current_candidate;

  ZOLTAN_TRACE_ENTER(zz, yo);

  *num_candidate_gids = num_input;
  *candidate_gids = NULL;

  /* Build data directory to find the candidates; RCB has dots distributed
   * according to RCB distribution, but need candidates according to the
   * input distribution
   */

  MPI_Allreduce(&dotnum, &maxdotnum, 1, MPI_INT, MPI_MAX, zz->Communicator);
  ierr = Zoltan_DD_Create(&dd, zz->Communicator, num_gid_entries, 
                          0, num_gid_entries*sizeof(ZOLTAN_ID_TYPE), 
                          maxdotnum, 0);

  /* Select candidates for each part; store them for each GID. */
  dot_candidates = (ZOLTAN_ID_PTR) 
                   ZOLTAN_MALLOC(dotnum*num_gid_entries*sizeof(ZOLTAN_ID_TYPE));
  if (dindx) {
    int prevpart = dotpt->Part[dindx[0]];
    current_candidate = &(gidpt[dindx[0]*num_gid_entries]);
    /* there is more than one part on this proc. */
    /* use dindx to access dots grouped by part */
    for (i = 0; i < dotnum; i++) {
      if (dotpt->Part[dindx[i]] != prevpart) {
        /* Assuming dots are grouped according to part number. */
        /* Then candidate is first dot in new part */
        current_candidate = &(gidpt[dindx[i]*num_gid_entries]);
        prevpart = dotpt->Part[dindx[i]];
      }
      ZOLTAN_SET_GID(zz, &(dot_candidates[dindx[i]*num_gid_entries]),
                         current_candidate);
    }
  }
  else {  
    /* dindx is NULL; all dots on this proc are in a single part */
    /* use first dot as candidate */
    ZOLTAN_ID_PTR current_dot_candidate = dot_candidates;
    current_candidate = &(gidpt[0]);
    for (i = 0; i < dotnum; i++) {
      ZOLTAN_SET_GID(zz, current_dot_candidate, current_candidate);
      current_dot_candidate += num_gid_entries;
    }
  }

  /* Update directory with candidates */
  ierr = Zoltan_DD_Update(dd, gidpt, NULL, (char*)dot_candidates, NULL, dotnum);

  /* Find candidates for input GIDs */
  if (num_input && 
      !Zoltan_Special_Malloc(zz, (void **)candidate_gids, num_input,
                               ZOLTAN_SPECIAL_MALLOC_GID)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
  }

  ierr = Zoltan_DD_Find(dd, input_gids, NULL, (char*)*candidate_gids, NULL,
                        num_input, NULL); 

End:

  Zoltan_DD_Destroy(&dd);
  ZOLTAN_FREE(&dot_candidates);

  if (ierr < 0) 
    Zoltan_Special_Free(zz, (void **)candidate_gids, ZOLTAN_SPECIAL_MALLOC_GID);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ierr);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
