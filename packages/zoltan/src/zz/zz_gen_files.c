/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002,2003, Sandia National Laboratories.          *
 * For more info, see the README file in the top-level Zoltan directory.     *   *****************************************************************************/
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
#define GLOBAL_GRAPH 2

/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines for generating output files
 *  from Zoltan that describe the data given to Zoltan. These 
 *  functions are all callable by the application, and can be
 *  invoked independent of load-balancing (Zoltan_LB_Balance).
 */
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Generate_Files(ZZ *zz, char *fname, int base_index)
{
/*
 *  Generate up to four output files:
 *   a) Current assignment of objects to processors (partitions?)
 *   b) Geometry of the objects.
 *   c) Graph if graph query functions are available.
 *   d) Hypergraph if hypergraph query functions are available.
 *
 *  The output is serialized, such that each processor
 *  will open and close each output file in turn.
 */

  int error=ZOLTAN_OK;
  ZOLTAN_ID_PTR local_ids = NULL;
  ZOLTAN_ID_PTR global_ids = NULL;
  FILE *fp;
  char full_fname[256];
  int *vtxdist, *xadj, *adjncy, *part;
  int *vwgt, *adjwgt;
  float *float_vwgt, *ewgts;
  double *xyz;
  int i, j, num_obj, num_geom, num_edges, glob_edges;
  char *yo = "Zoltan_Generate_Files";

  /* Initialize all local pointers to NULL. This is necessary
   * because we free all non-NULL pointers upon errors.
   */
  vtxdist = xadj = adjncy = vwgt = adjwgt = part = NULL;
  float_vwgt = ewgts = NULL;
  xyz = NULL;

  /* Zoltan_Get_Obj_List allocates memory for all return lists. */
  error = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids,
                              zz->Obj_Weight_Dim, &float_vwgt, &part);
  if (error){
    /* Return error */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Get_Obj_List returned error.");
    error = ZOLTAN_FATAL;
    goto End;
  }

  /* Build (ParMetis) graph data structures, or just get vtxdist. */
  error = Zoltan_Build_Graph(zz, GLOBAL_GRAPH, 1, num_obj,
         global_ids, local_ids, zz->Obj_Weight_Dim, zz->Edge_Weight_Dim,
         &vtxdist, &xadj, &adjncy, &ewgts);
  if (error != ZOLTAN_OK && error != ZOLTAN_WARN){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Zoltan_Build_Graph returned error.");
    error = ZOLTAN_FATAL;
    goto End;
  }
  /* Local number of edges. */
  if (xadj==NULL || xadj[num_obj]==0)
    num_edges = 0;
  else
    num_edges = xadj[num_obj];
  /* Compute global number of edges. */
  MPI_Reduce(&num_edges, &glob_edges, 1, MPI_INT, MPI_SUM, 0, 
      zz->Communicator);  
  glob_edges /= 2;

  /* Write to files, serialized. */
  /* Note: This will be slow (not scalable) for many procs. */
  Zoltan_Print_Sync_Start(zz->Communicator, 0); 

  /* Write object assignments to file. */
  /* Only write proc number, or partition number too?? */
  if (num_obj > 0){
    sprintf(full_fname, "%s.assign", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_WARN;
      goto End;
    }
    for (i=0; i<num_obj; i++)
      fprintf(fp, "%d\n", zz->Proc);
    fclose(fp);
  }

  /* Write geometry to file, if applicable. */
  if (zz->Get_Num_Geom != NULL && zz->Get_Geom != NULL) {
    num_geom = zz->Get_Num_Geom(zz->Get_Num_Geom_Data, &error);
    if (error != ZOLTAN_OK) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
             "Error returned from user function Get_Num_Geom.");
      goto End;
    }
    if (num_geom>0)
      xyz = (double *) ZOLTAN_MALLOC (num_geom * sizeof(double));

    sprintf(full_fname, "%s.coords", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_WARN;
      goto End;
    }
    for (i=0; i<num_obj; i++){
      zz->Get_Geom(zz->Get_Geom_Data, zz->Num_GID, zz->Num_LID,
              &(global_ids[i*zz->Num_GID]), &(local_ids[i*zz->Num_LID]),
              xyz, &error);
      if (error == ZOLTAN_FATAL || error == ZOLTAN_MEMERR) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo,
             "Error returned from user function Get_Geom.");
        goto End;
      }
      for (j=0; j<num_geom; j++)
        fprintf(fp, "%f ", xyz[j]);
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  /* Write graph to file, if applicable. */
  if (zz->Get_Num_Edges != NULL && zz->Get_Edge_List != NULL) {
    sprintf(full_fname, "%s.graph", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_WARN;
      goto End;
    }

    /* If proc 0, write first line. */
    if (zz->Proc == 0){
      fprintf(fp, "%d %d %d\n", vtxdist[zz->Num_Proc], glob_edges, 0);
    }

    /* Print edge list for each node (object). */
    for (i=0; i<num_obj; i++){
      for (j=xadj[i]; j<xadj[i+1]; j++){
        fprintf(fp, "%d ", adjncy[j]+base_index);
      }
      fprintf(fp, "\n");
    }
  }
  
  Zoltan_Print_Sync_End(zz->Communicator, 0); 

End:
  ZOLTAN_FREE(&xyz);
  ZOLTAN_FREE(&vtxdist);
  ZOLTAN_FREE(&xadj);
  ZOLTAN_FREE(&adjncy);
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);
  ZOLTAN_FREE(&float_vwgt);
  ZOLTAN_FREE(&part);
  return error;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

