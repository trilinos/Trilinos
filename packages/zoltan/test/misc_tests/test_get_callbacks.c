// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/* Test for the Zoltan_Get_Fn interface function */

#include <stdio.h>
#include <mpi.h>
#include "zoltan.h"

int numObjFn(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  printf("Greetings from numObjFn\n");
  return *((int *)data);
}

void objSizeMultiFn(void *data, int ngid, int nlid, int nids, 
                    ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, 
                    int *nbytes, int *ierr)
{
  int i;
  *ierr = ZOLTAN_OK;
  printf("Greetings from objSizeMulti\n");
  for (i = 0; i < nids; i++) nbytes[i] = ((int *)data)[i];
}

void packObjMultiFn(void *data, int ngid, int nlid, int nids,
                    ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int *dest,
                    int *size, int *index, char *buf, int *ierr)
{
  *ierr = ZOLTAN_OK;
  printf("Greetings from packObjMultiFn %d\n", *((int *)data));
  ((int *)buf)[0] = *((int *) data);
}

void unpackObjMultiFn(void *data, int ngid, int nids, ZOLTAN_ID_PTR gids, 
                      int *size, int *index, char *buf, int *ierr)
{
  *ierr = ZOLTAN_OK;
  printf("Greetings from unpackObjMultiFn %d\n", *((int *)data));
  ((int *)buf)[0] = *((int *) data);
}

/****************************************************************************/
/****************************************************************************/
int main (int narg, char **arg)
{
  float ver;
  struct Zoltan_Struct *zz = NULL;
  ZOLTAN_VOID_FN *fn = NULL;
  int *fndata = NULL;
  int me, np; 
  int i, ierr;
  int nobj = 123, sizes[123], fn_sizes[123], buf = -1;
  int fn_val = 0;
  int fndata_val = 0;
  int nerrs = 0;

  /* Initialize Zoltan */
  MPI_Init(&narg, &arg);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  /* Initialize some data */
  for (i = 0; i < nobj; i++) {
    sizes[i] = i * sizeof(int);
    fn_sizes[i] = 0;
  }

  /* Initialize Zoltan */
  Zoltan_Initialize(narg, arg, &ver);
  zz = Zoltan_Create(MPI_COMM_WORLD);

  /* Register callbacks */
  Zoltan_Set_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE,
                (ZOLTAN_VOID_FN *)numObjFn, (void *)(&nobj));
  Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE,
                (ZOLTAN_VOID_FN *)objSizeMultiFn, (void *)sizes);
  Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_MULTI_FN_TYPE,
                (ZOLTAN_VOID_FN *)packObjMultiFn, (void *)(&me));
  Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE,
                (ZOLTAN_VOID_FN *)unpackObjMultiFn, (void *)(&np));

  /* Retrieve callbacks and call them; test results */
  /* num_obj */
  ierr = Zoltan_Get_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE,
                       &fn, (void **)&fndata);
  fn_val = ((ZOLTAN_NUM_OBJ_FN *)fn)((void *)fndata, &ierr);
  fndata_val = *fndata;
  printf("%d of %d: ZOLTAN_NUM_OBJ_FN_TYPE nobj %d fn %d fndata %d\n", 
         me, np, nobj, fn_val, fndata_val);
  if (!((ierr == ZOLTAN_OK) && (fn_val == nobj) && (fndata_val == nobj))) {
    printf("FAIL\n");
    nerrs++;
  }

  /* obj_size */
  ierr = Zoltan_Get_Fn(zz, ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE,
                       &fn, (void **)&fndata);
  ((ZOLTAN_OBJ_SIZE_MULTI_FN *)fn)((void *)fndata, 1, 1, nobj, NULL, NULL, 
                                   fn_sizes, &ierr);
  for (i = 0; i < nobj; i++) {
    if (!((ierr == ZOLTAN_OK) && (fn_sizes[i] == sizes[i]) 
                              && (fndata[i] == sizes[i]))) {
      printf("%d of %d: ZOLTAN_OBJ_SIZE_FN_TYPE sizes %d "
             "fn_sizes %d fndata %d\n",
             me, np, sizes[i], fn_sizes[i], fndata[i]);
      printf("FAIL\n");
      nerrs++;
    }
  }

  /* pack */
  ierr = Zoltan_Get_Fn(zz, ZOLTAN_PACK_OBJ_MULTI_FN_TYPE,
                       &fn, (void **)&fndata);
  ((ZOLTAN_PACK_OBJ_MULTI_FN *)fn)((void *)fndata, 1, 1, nobj, 
                                   NULL, NULL, NULL, NULL, NULL, 
                                   (char *)&buf, &ierr);
  fndata_val = *fndata;
  printf("%d of %d: ZOLTAN_PACK_OBJ_MULTI_FN_TYPE fn %d fndata %d\n", 
         me, np, buf, fndata_val);
  if (!((ierr == ZOLTAN_OK) && (buf == me) && (fndata_val == me))) {
    printf("FAIL\n");
    nerrs++;
  }

  /* unpack */
  ierr = Zoltan_Get_Fn(zz, ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE,
                       &fn, (void **)&fndata);
  ((ZOLTAN_UNPACK_OBJ_MULTI_FN *)fn)((void *)fndata, 1, nobj,
                                     NULL, NULL, NULL, (char *)&buf, &ierr);
  fndata_val = *fndata;
  printf("%d of %d: ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE fn %d fndata %d\n", 
         me, np, buf, fndata_val);
  if (!((ierr == ZOLTAN_OK) && (buf == np) && (fndata_val == np))) {
    printf("FAIL\n");
    nerrs++;
  }

  /* Wrap up */
  if (nerrs == 0) printf("PASS\n");

  Zoltan_Destroy(&zz);
  MPI_Finalize();

  return 0;
}
