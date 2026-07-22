// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include <stdio.h>
#include <stdlib.h>
#include "zoltan.h"

/****************************************************************************/
/* Test packing and unpacking of ZZ struct for RCB cut-tree communication   */
/* See issue #8476                                                          */
/****************************************************************************/


/*****************************************************************************
 * Toy mesh and callback functions on it 
 * Store global mesh on each processor.  
 * Each processor reports only a portion of the mesh to Zoltan in callbacks
 *****************************************************************************/
struct Mesh {
  int nCoords;        // Store global mesh on each proc
  double *x, *y, *z;  // Store global mesh on each proc
  int nMyCoords;      // Number of coordinates reported by this proc
  int myFirstCoord;   // First coordinate reported by this proc
};

void initMesh(struct Mesh *mesh, int me, int np, 
              int n_, double *x_, double *y_, double *z_) 
{
  mesh->nCoords = n_;
  mesh->x = x_;  mesh->y = y_;  mesh->z = z_;
  int coordsPerProc = mesh->nCoords/np;
  mesh->nMyCoords = (me == (np-1) ? mesh->nCoords - (np-1)*coordsPerProc 
                                  : coordsPerProc);
  mesh->myFirstCoord = me * coordsPerProc;
}

int nObj(void *data, int *ierr) {
  *ierr = ZOLTAN_OK; 
  return ((struct Mesh *) data)->nMyCoords;
}

void objMulti(void *data, int ngid, int nlid, 
              ZOLTAN_ID_PTR gid, ZOLTAN_ID_PTR lid, int wdim, float *wgt,
              int *ierr) 
{
  *ierr = ZOLTAN_OK;
  struct Mesh *mesh = (struct Mesh *) data;
  int i, j;
  for (i = 0; i < mesh->nMyCoords; i++) {
    lid[i*nlid] = i+mesh->myFirstCoord;
    gid[i*ngid] = i+mesh->myFirstCoord; 
    for (j = 0; j < wdim; j++) wgt[i*wdim+j] = 1.;
  }
}

int nGeom(void *data, int *ierr) { *ierr = ZOLTAN_OK; return 3; }

void geomMulti(void *data, int ngid, int nlid, int nobj, 
               ZOLTAN_ID_PTR gid, ZOLTAN_ID_PTR lid, int ndim,
               double *coords, int *ierr)
{
  struct Mesh *mesh = (struct Mesh *) data;
  int i;
  for (i = 0; i < nobj; i++) {
    coords[i*ndim]   = mesh->x[lid[i*nlid]];
    coords[i*ndim+1] = mesh->y[lid[i*nlid]];
    coords[i*ndim+2] = mesh->z[lid[i*nlid]];
  }
  *ierr = ZOLTAN_OK;
}

/**************************************************************************
 * Test packing and unpacking of ZZ struct for RCB cut-tree communication
 * Subcommunicator computes RCB decomposition and keeps cuts.
 * Some rank which is in both subcommunicator and MPI_COMM_WORLD (in this 
 * case, rank 0) packs the resulting Zoltan struct into a buffer and 
 * broadcasts it to all procs in MPI_COMM_WORLD. 
 * All procs unpack the buffer into a new Zoltan struct.
 * All procs compare LB_Point_Assign results using the new Zoltan struct 
 * to baselines using the original.
 **************************************************************************/
int run_test(       /* Test options relevant to serialization, such as... */
  int nparts,       /*      NUM_GLOBAL_PARTS                              */
  int ngid,         /*      NUM_GID_ENTRIES                               */
  int nlid,         /*      NUM_LID_ENTRIES                               */
  int usePartSizes, /*      Part_Info array                               */
  int useRemap      /*      Remap array                                   */
) 
{
  int ierr = ZOLTAN_OK;
  int me, np;
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  if (np < 3) printf("This test is more useful on three or more processors.\n");

  /* Coordinates to be partitioned */
  const int nOne = 27;
  double xOne[] = {0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2};
  double yOne[] = {0,0,0,1,1,1,2,2,2,0,0,0,1,1,1,2,2,2,0,0,0,1,1,1,2,2,2};
  double zOne[] = {0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2};
  struct Mesh meshOne;
  initMesh(&meshOne, me, np, nOne, xOne, yOne, zOne);

  /* Coordinates to use for testing */
  const int nTwo = 8;
  double xTwo[] = {0.1, 1.1, 2.1, 1.1, 0.1, 0.1, 2.1, 2.1};
  double yTwo[] = {0.1, 0.1, 0.1, 1.1, 1.1, 0.1, 1.1, 2.1};
  double zTwo[] = {0.1, 0.1, 0.1, 1.1, 1.1, 1.1, 1.1, 2.1};

  /* Create a subcommunicator in which to do partitioning. */
  /* For this test, we'll put rank 0 of MPI_COMM_WORLD in subComm */
  MPI_Comm subComm = MPI_COMM_NULL;
  MPI_Comm_split(MPI_COMM_WORLD, (me <= np / 2 ? 1 : MPI_UNDEFINED), 0,
                 &subComm);
  
  /* Compute the RCB partition on a subcommunicator subComm */
  struct Zoltan_Struct *zz = NULL;

  if (subComm != MPI_COMM_NULL) {

    zz = Zoltan_Create(subComm);

    Zoltan_Set_Num_Obj_Fn(zz, nObj, &meshOne);
    Zoltan_Set_Obj_List_Fn(zz, objMulti, &meshOne);
    Zoltan_Set_Num_Geom_Fn(zz, nGeom, &meshOne);
    Zoltan_Set_Geom_Multi_Fn(zz, geomMulti, &meshOne);

    Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
    Zoltan_Set_Param(zz, "KEEP_CUTS", "1");
    Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.02");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "PART");
    (useRemap ? Zoltan_Set_Param(zz, "REMAP", "1")
              : Zoltan_Set_Param(zz, "REMAP", "0"));

    char msg[22];
    sprintf(msg, "%d", nparts);
    Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", msg);
    sprintf(msg, "%d", ngid);
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", msg);
    sprintf(msg, "%d", nlid);
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", msg);

    if (usePartSizes) {
      int *parts = (int *) malloc(sizeof(int) * nparts);
      float *sizes = (float *) malloc(sizeof(float) * nparts);
      int i;
      for (i = 0; i < nparts; i++) {
        parts[i] = i;
        sizes[i] = i+1;
      }
      Zoltan_LB_Set_Part_Sizes(zz, 1, nparts, parts, NULL, sizes);
      free(parts);
      free(sizes);
    }

    int nChanges;
    int nGid, nLid;
    int nImp, nExp;
    ZOLTAN_ID_PTR iGid = NULL, iLid = NULL;
    ZOLTAN_ID_PTR eGid = NULL, eLid = NULL;
    int *iPart = NULL, *iProc = NULL;
    int *partOne = NULL, *eProc = NULL;

    Zoltan_LB_Partition(zz, &nChanges, &nGid, &nLid,
                        &nImp, &iGid, &iLid, &iProc, &iPart,
                        &nExp, &eGid, &eLid, &eProc, &partOne);

    Zoltan_LB_Free_Part(&iGid, &iLid, &iProc, &iPart);
    Zoltan_LB_Free_Part(&eGid, &eLid, &eProc, &partOne);
  }

  /* Pack the buffer; broadcast to all procs in MPI_COMM_WORLD. */
  /* Test assumes that rank 0 of MPI_COMM_WORLD is in subComm.  */

  /* First broadcast the buffer size; we make the user own the buffer */
  size_t bufSize;
  if (me == 0) bufSize = Zoltan_Serialize_Size(zz);

  MPI_Bcast((char *) &bufSize, sizeof(bufSize), MPI_CHAR, 0, MPI_COMM_WORLD);

  /* Then allocate and broadcast the buffer */
  char *buf = NULL;
  buf = (char *) malloc(bufSize * sizeof(char));

  if (me == 0) ierr = Zoltan_Serialize(zz, bufSize, buf);
  if (ierr != ZOLTAN_OK) {
    printf("TEST FAILED:  Error in Zoltan_Serialize\n"); fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  MPI_Bcast(buf, bufSize, MPI_CHAR, 0, MPI_COMM_WORLD);

  /* All processors unpack the buffer into a new ZZ struct */

  struct Zoltan_Struct *newZZ = Zoltan_Create(MPI_COMM_WORLD);
  ierr = Zoltan_Deserialize(newZZ, bufSize, buf);
  if (ierr != ZOLTAN_OK) {
    printf("TEST FAILED:  Error in Zoltan_Deserialize\n"); fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  free(buf); buf = NULL;
  
  /* Check the results */
  /* Compute and broadcast answer using original struct zz */

  int answer[nTwo];
  if (me == 0) {
    int i;
    for (i = 0; i < nTwo; i++) {
      double tmp[3]; tmp[0] = xTwo[i]; tmp[1] = yTwo[i]; tmp[2] = zTwo[i];
      Zoltan_LB_Point_PP_Assign(zz, tmp, NULL, &answer[i]);
      printf("Point (%f %f %f) on part %d\n", 
              xTwo[i], yTwo[i], zTwo[i], answer[i]);
    }
  }
  MPI_Bcast(answer, nTwo, MPI_INT, 0, MPI_COMM_WORLD);

  /* Each processor computes answer using new struct newZZ */

  int errCnt = 0;
  int i;
  for (i = 0; i < nTwo; i++) {
    int newAnswer;
    double tmp[3]; tmp[0] = xTwo[i]; tmp[1] = yTwo[i]; tmp[2] = zTwo[i];
    Zoltan_LB_Point_PP_Assign(newZZ, tmp, NULL, &newAnswer);
    printf("%d Point (%f %f %f) on part %d %d\n", 
           me, xTwo[i], yTwo[i], zTwo[i], answer[i], newAnswer);
    if (newAnswer != answer[i]) {
      errCnt++;
      printf("%d Error (%f %f %f):  part %d != new part %d\n", 
             me, xTwo[i], yTwo[i], zTwo[i], answer[i], newAnswer);
    }
  }

  if (zz) Zoltan_Destroy(&zz);
  if (newZZ) Zoltan_Destroy(&newZZ);
  if (subComm != MPI_COMM_NULL) MPI_Comm_free(&subComm);

  /* Gather global test result */

  int gErrCnt = 0;
  MPI_Allreduce(&errCnt, &gErrCnt, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (me == 0) 
    printf("TEST %s: gErrCnt=%d\n", (gErrCnt ? "FAILED" : "PASSED"), gErrCnt);

  return gErrCnt;
}

/**************************************************************************/
int main(int narg, char **arg) {

  MPI_Init(&narg, &arg);
  float ver;
  Zoltan_Initialize(narg, arg, &ver);
  int ierr = ZOLTAN_OK;
  int me, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  if (me == 0) {printf("TEST ONE\n"); fflush(stdout); }
  ierr += run_test(5, 1, 1, 1, 1);

  if (me == 0) {printf("TEST TWO\n"); fflush(stdout); }
  ierr += run_test(np+4, 1, 2, 0, 1);

  if (me == 0) {printf("TEST THREE\n"); fflush(stdout); }
  ierr += run_test((np-2>0 ? np-2 : np), 2, 3, 1, 0);

  if (me == 0) {printf("TEST FOUR\n"); fflush(stdout); }
  ierr += run_test(np, 3, 1, 0, 1);

  if (me == 0) {printf("TEST FIVE\n"); fflush(stdout); }
  ierr += run_test(np, 3, 1, 1, 0);

  if (me == 0) {printf("TEST SIX\n"); fflush(stdout); }
  ierr += run_test(np, 3, 1, 1, 1);

  MPI_Finalize();
  return ierr;
}
