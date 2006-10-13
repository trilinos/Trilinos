#ifndef SIMPLEGRAPH_H
#define SIMPLEGRAPH_H
/*
 * A very simple mesh or graph.
 *   Global IDs of vertices are 1 through 25.
 *   Vertex adjacencies are listed in the "edges" array.
 *
 *   Regarded as a mesh, vertex 1 is at location (0,0)
 *   and vertex 25 is at location (4,4);
 *
 *  21----22----23----24---25
 *  |     |     |     |    |
 *  16----17----18----19---20
 *  |     |     |     |    |
 *  11----12----13----14---15
 *  |     |     |     |    |
 *  6-----7-----8-----9----10
 *  |     |     |     |    |
 *  1-----2-----3-----4----5
 *
 */

static int simpleNumVertices=25;
static int simpleNumEdges[25] = {
2, 3, 3, 3, 2,
3, 4, 4, 4, 3,
3, 4, 4, 4, 3,
3, 4, 4, 4, 3,
2, 3, 3, 3, 2
};
static int simpleEdges[25][4]={
{2,6},       /* adjacent to vertex 1 */
{1,3,7},     /* adjacent to vertex 2 */
{2,8,4},
{3,9,5},
{4,10},
{1,7,11},
{6,2,8,12},
{7,3,9,13},
{8,4,10,14},
{9,5,15},
{6,12,16},
{11,7,13,17},
{12,8,14,18},
{13,9,15,19},
{14,10,20},
{11,17,21},
{16,12,18,22},
{17,13,19,23},
{18,14,20,24},
{19,15,25},
{16,22},
{21,17,23},
{22,18,24},
{23,19,25},
{24,20}      /* adjacent to vertex 25 */
};

#include <mpi.h>

/* All processes must call this, process 0 draws the picture */

static void draw_partitions(int nobj, int *gids, int wgt_dim, float *wgts)
{
int i, j, me, idx, nprocs, root=0;
int *incounts=NULL, *ingids=NULL, *indisp=NULL, *parts=NULL, *gid;
float *inwgts=NULL, *wgt_total=NULL, *wgt;

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  if (me == root){
    incounts = (int *)malloc(sizeof(int) * nprocs);
  }

  MPI_Gather(&nobj, 1, MPI_INT, incounts, 1, MPI_INT, root, MPI_COMM_WORLD);

  if (me == root){
    indisp = (int *)malloc(sizeof(int) * (nprocs + 1));
    indisp[0] = 0;
    for (i=1; i<=nprocs; i++){
      indisp[i] = indisp[i-1] + incounts[i-1];
    }
    ingids = (int *)malloc(sizeof(int) * indisp[nprocs]);
  }

  MPI_Gatherv(gids, nobj, MPI_INT, ingids, incounts, indisp, MPI_INT,
              root, MPI_COMM_WORLD);

  if (wgt_dim > 0){
    if (me == root){
      inwgts = (float *)malloc(sizeof(float) * wgt_dim * indisp[nprocs]);
      if (wgt_dim > 1){
        for (i=0; i < nprocs; i++){
          incounts[i] *= wgt_dim;
          indisp[i] *= wgt_dim;
        }
      }
    }
    MPI_Gatherv(wgts, nobj*wgt_dim, MPI_FLOAT, inwgts, incounts, indisp, 
                MPI_FLOAT, root, MPI_COMM_WORLD);

    if ((me == root) && (wgt_dim > 1)){
      for (i=0; i < nprocs; i++){
        incounts[i] /= wgt_dim;
        indisp[i] /= wgt_dim;
      }
    }
  }


  if (me == root){
    wgt = inwgts;
    gid = ingids;
    wgt_total = (float *)calloc(sizeof(float) , nprocs);
    parts = (int *)malloc(sizeof(int) * simpleNumVertices);


    for (i=0; i<nprocs; i++){
      for (j=0; j<incounts[i]; j++){
        idx = *gid++ - 1;
        parts[idx] = i;
        if (wgt_dim > 0){
          wgt_total[i] += *wgt;
          wgt += wgt_dim;
        }
      }
    }

    for (i=21; i>0; i-=5){
      for (j=0; j<5; j++){
        printf("%2d  ",parts[i+j-1]);
      }
      printf("\n\n");
    }

    if (wgt_dim > 0){
      for (i=0; i<nprocs; i++){
        printf("Total weight partition %d: %f\n",i, wgt_total[i]);   
      }
    }
  }
  
  if (incounts) free(incounts);
  if (indisp) free(indisp);
  if (ingids) free(ingids);
  if (inwgts) free(inwgts);
  if (parts) free(parts);
  if (wgt_total) free(wgt_total);
}

#endif
