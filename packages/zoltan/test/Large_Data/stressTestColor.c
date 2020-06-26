/*
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine      kddevin@sandia.gov
 *                    Erik Boman        egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */
/************************************************************
* This is called a stress test because it builds an
* arbitrarily large graph.  It tests the HIER_ASSIST
* option to hierarchical partitioning.
*
* TODO:
* Create a function that performs communication where comm
* volume is proportional to the graph edge weight.  This is
* to test the value of partitioning to network hierarchy.
***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <signal.h>

#ifndef _MSC_VER
#include <getopt.h>
#endif /* _MSC_VER */

#include "zz_const.h"

#define NX 16
#define NY 16
#define NZ 16

static int myRank, numProcs;
static size_t numMyVertices;
ZOLTAN_ID_TYPE *vtxGID = NULL;
ZOLTAN_ID_TYPE *nborGID = NULL;
int *nborIndex = NULL;
int *nborProc = NULL;


static void check_error_status(int status, char *s)
{
int gstatus;

  MPI_Allreduce(&status, &gstatus, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if (gstatus > 0){
    if (myRank == 0){
      fprintf(stderr,"Error: %s\n",s);
    }
    MPI_Finalize();
    exit(1);
  }
}

static void free_graph()
{
  if (vtxGID) free(vtxGID);
  if (nborIndex) free(nborIndex);
  if (nborGID) free(nborGID);
  if (nborProc) free(nborProc);

  vtxGID = nborGID = NULL;
  nborIndex = nborProc = NULL;
}


static int create_a_graph(int nx, int ny, int nz)
{
  int i, sum, n, j;
  long gid, nbors[6], count[2];

  /* Divide z direction among procs */
  if (myRank == 0)
    printf("%d x %d x %d = %d\n", nx, ny, nz, nx*ny*nz);

  int myNz = nz / numProcs;
  numMyVertices = nx * ny * myNz;
  int myFirstVtx = myRank * numMyVertices;
  int myFirstZ = myRank * myNz;
  int nLayer = nx*ny;

  vtxGID = malloc(sizeof(ZOLTAN_ID_TYPE)*numMyVertices);
  nborGID = malloc(sizeof(ZOLTAN_ID_TYPE)*6*numMyVertices);  
  nborProc = malloc(sizeof(int)*6*numMyVertices);
  nborIndex = malloc(sizeof(int)*(numMyVertices+1));
  nborIndex[0] = 0;

  int vcnt = 0;
  int ecnt = 0;
  int x, y, z;
  for (z = 0; z < myNz; z++) {
    for (y = 0; y < ny; y++) {
      for (x = 0; x < nx; x++) {
        int id = vcnt + myFirstVtx;
        vtxGID[vcnt] = id;
        if (x != 0) { /* left */
          nborGID[ecnt] = id-1; 
          nborProc[ecnt] = myRank;
          ecnt++;
        }
        if (x != nx-1) { /* right */
          nborGID[ecnt] = id+1;
          nborProc[ecnt] = myRank;
          ecnt++;
        }
        if (y != 0) { /* down */
          nborGID[ecnt] = id-nx;
          nborProc[ecnt] = myRank;
          ecnt++;
        }
        if (y != ny-1) {  /* up */
          nborGID[ecnt] = id+nx;
          nborProc[ecnt] = myRank;
          ecnt++;
        }
        if (z+myFirstZ != 0) { /* front */
          nborGID[ecnt] = id-nLayer;
          nborProc[ecnt] = (z != 0 ? myRank : myRank-1);
          ecnt++;
        }
        if (z+myFirstZ != nz-1) { /* back */
          nborGID[ecnt] = id+nLayer;
          nborProc[ecnt] = (z != myNz-1 ? myRank : myRank+1);
          ecnt++;
        }
        vcnt++;
        nborIndex[vcnt] = ecnt;
      }
    }
  }

#if 0
  for (int i = 0; i < numMyVertices; i++) {
    printf("%d %d Vtx %d: ", myRank, i, vtxGID[i]);
    for (int j = nborIndex[i]; j < nborIndex[i+1]; j++) 
      printf("%d:%d ", nborGID[j], nborProc[j]);
    printf("\n");
  }
#endif
  
  return 0;
}


/* Zoltan query functions. */

static int get_number_of_vertices(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return numMyVertices;
}

static void get_vertex_list(void *data, int sizeGID, int sizeLID,
                  ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{
  int i;
  *ierr = ZOLTAN_OK;

  for (i=0; i < numMyVertices; i++){
    globalID[i] = vtxGID[i];
    localID[i] = i;
  }
}

static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int *numEdges, int *ierr)
{
int i;

  *ierr = ZOLTAN_OK;

  for (i=0; i < num_obj; i++){
    numEdges[i] = nborIndex[localID[i]+1] - nborIndex[localID[i]];
  }
}

static void get_edge_list(void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges,
        ZOLTAN_ID_PTR nbor, int *owner, int wgt_dim, float *ewgts, int *ierr)
{
int nextv, nextp, npins, p1, i, lid;

  *ierr = ZOLTAN_OK;

  for (nextv=0, nextp=0; nextv < num_obj; nextv++){

    lid = localID[nextv];
    p1 = nborIndex[lid];
    npins = nborIndex[lid+1] - p1;

    if (num_edges[nextv] != npins){
      fprintf(stderr,"num edges != num pins\n");
      *ierr = ZOLTAN_FATAL;
      return;
    }

    for (i=0; i <  npins; i++, nextp++){
      nbor[nextp] = nborGID[p1+i];
      owner[nextp] = nborProc[p1+i];
    }
  }
}

int main(int argc, char *argv[])
{
  int rc, status;
  float ver;
  struct Zoltan_Struct *zz;
  int ierr = ZOLTAN_OK;
  int *colors = NULL;

#ifndef _MSC_VER
  struct option opts[10];
#endif /* _MSC_VER */

  int nx = NX, ny = NY, nz = NZ;
  char *debug_level="1";

  status = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  Zoltan_Initialize(argc, argv, &ver);
  zz = Zoltan_Create(MPI_COMM_WORLD);

  /******************************************************************
  ** Check that this test makes sense.
  ******************************************************************/

  if (sizeof(long) < sizeof(ZOLTAN_ID_TYPE)){
    if (myRank == 0){
      printf("ERROR: This code assumes that a long is at least %d bytes\n",(int)sizeof(ZOLTAN_ID_TYPE));
    }
    status = 1;
  }

  check_error_status(status, "configuration error");

  /******************************************************************
  ** Initialize zoltan
  ******************************************************************/

#ifdef _MSC_VER
  if (myRank == 0) {
    printf("\n*** getopt not supported in Windows; ");
    printf("command-line arguments will be ignored ***\n\n");
  }
#else
  /* options for Unix runs; Windoze will use default values because it does
   * not have getopt
   */

  opts[0].name = "nx";
  opts[0].has_arg = 1;
  opts[0].flag = NULL;
  opts[0].val = 0;

  opts[1].name = "ny";
  opts[1].has_arg = 1;
  opts[1].flag = NULL;
  opts[1].val = 1;

  opts[2].name = "nz";
  opts[2].has_arg = 1;
  opts[2].flag = NULL;
  opts[2].val = 2;

  opts[3].name = "debug_level";
  opts[3].has_arg = 1;
  opts[3].flag = NULL;
  opts[3].val = 3;

  opts[4].name = 0;
  opts[4].has_arg = 0;
  opts[4].flag = NULL;
  opts[4].val = 0;

  status = 0;

  while (1){
    rc = getopt_long_only(argc, argv, "",  opts, NULL);

    if (rc == '?'){
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(0);
    }
    else if (rc == 0){
      nx = atoi(optarg);
    }
    else if (rc == 1){
      ny = atoi(optarg);
    }
    else if (rc == 2){
      nz = atoi(optarg);
    }
    else if (rc == 3){
      debug_level = optarg;
    }
    else if (rc <= 0){
      break;
    }
  }
#endif /* _MSC_VER */

  /* start */

  Zoltan_Memory_Debug(0);

  status = create_a_graph(nx, ny, nz);
  check_error_status(status, "creating the graph");

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", debug_level);
  Zoltan_Set_Param(zz, "GRAPH_BUILD_TYPE", "FAST_NO_DUP");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, NULL);
  Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, NULL);
  Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list,  NULL);
  Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list,  NULL);

  Zoltan_Set_Param(zz, "COLORING_PROBLEM", "DISTANCE-1");

  colors = (int *) malloc(sizeof(int) * numMyVertices);

  rc = Zoltan_Color(zz, 1, numMyVertices, vtxGID, colors);

  free_graph();
  free(colors);

  Zoltan_Destroy(&zz);

  MPI_Finalize();

  return status;
}


