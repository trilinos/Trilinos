/*
** $Id$
**
** Functions to support writing simple Zoltan examples.
**   Create a simple rectilinear mesh with global IDs and
**   divide it among the processes.  Define call backs
**   that return points in the mesh.
*/

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "exzoltan.h"

#ifdef __cplusplus
extern "C" {
#endif

static int Divisions = 10;
static float *Points=NULL;
static int *GlobalIds=NULL;
static int NumPoints=0;

/*
** "Divisions" is the number of Divisions in the x, y, and z directions
*/
void exSetDivisions(int div)
{
  if (div > 0) Divisions = div;
}

int exInitializePoints(float **retPts, int **retIds, int rank, int size)
{
  float dx, dy, dz;
  float xmin, ymin, zmin;
  float *pt;
  float x, y, z;
  float *pts;
  int *ids;
  int id, i, j, k;
  int *numPts;
  int ptsPerProc, ptsAssigned, mySize;
  MPI_Status stat;
  float *sendPts;
  int *sendIds;

  int npts = Divisions * Divisions * Divisions;

  if (rank == 0)
    {
    pts = (float *)malloc(npts * 3 * sizeof(float));
    ids = (int *)malloc(npts * sizeof(int));

    dx = .4;
    dy = .3;
    dz = .1;

    pt = pts;
    id = 0;
    xmin = -12.0;
    ymin = -25.0;
    zmin = -30.0;

    for (i=0; i<Divisions; i++)
      {
      x = xmin + (i * dx);
      for (j=0; j<Divisions; j++)
        {
        y = ymin + (j * dy);
        for (k=0; k<Divisions; k++)
          {
          z = zmin + (k * dz);
          *pt++ = x;
          *pt++ = y;
          *pt++ = z;
          ids[id] = id;
          id++;
          }
        }
      }
    }

  /* divide these to start */

  numPts = (int *)malloc(sizeof(int) * size);
  ptsPerProc = npts / size;
  ptsAssigned = 0;

  for (i=0; i<size-1; i++)
    {
    numPts[i] = ptsPerProc;
    ptsAssigned += ptsPerProc;
    }

  numPts[size-1] = npts - ptsAssigned;

  mySize = numPts[rank];
  if (rank == 0)
    {
    sendPts = pts + (3 * numPts[0]);
    sendIds = ids + numPts[0];

    for (i=1; i<size; i++)
      {
      MPI_Send(sendPts, 3 * numPts[i], MPI_FLOAT, i, 0x01,MPI_COMM_WORLD);
      MPI_Send(sendIds, numPts[i], MPI_INT, i, 0x03,MPI_COMM_WORLD);
      sendPts += (3 * numPts[i]);
      sendIds += numPts[i];
      }
    }
  else
    {
    pts = (float *)malloc(sizeof(float) * 3 * mySize);
    ids = (int *)malloc(sizeof(int) * mySize);
    MPI_Recv(pts, 3 * mySize, MPI_FLOAT, 0, 0x01, MPI_COMM_WORLD, &stat);
    MPI_Recv(ids, mySize, MPI_INT, 0, 0x03, MPI_COMM_WORLD, &stat);
    }     
          
  Points = *retPts = pts;
  GlobalIds = *retIds = ids;  
  NumPoints = mySize;
          
  return mySize;
}     
void exShowIds(int *ids, int n)
{
  int i;

  for (i=0; i<n; i++)
    {
    if (i && (i%10 == 0)) printf("\n");
    printf("  %d", ids[i]);
    }
  printf("\n");
}
void exShowPoints(float *pts, int n, int *ids)
{
  int i;
  float *pt = pts;
  for (i=0; i<n; i++)
    {
    printf("  %d:  %f, %f, %f\n",ids[i],pt[0],pt[1],pt[2]);
    pt += 3;
    }
}
/**********************
** call backs
**********************/

int exGetNumberOfAssignedObjects(void *userDefinedData, int *err)
{
  *err = 0;
  return NumPoints;
}
void exGetObjectList(void *userDefinedData, int numGlobalIds, int numLids,
  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int wgt_dim, float *obj_wgts,
  int *err)
{
  int i;
    
  for (i=0; i<NumPoints; i++)
    {
    gids[i] = GlobalIds[i];
    lids[i] = i;
    }
    
  *err = 0;
    
  return;
}
int exGetObjectSize(void *userDefinedData, int *err)
{
  *err = 0; 
  return 3;
} 
void exGetObject(void *userDefinedData, int numGlobalIds, int numLids, int numObjs,
  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int numDim, double *pts, int *err)
{ 
  int i, id, id3;
  int next = 0;
  
  if (numDim != 3)    
    {
    *err = 1;         
    return;
    }

  for (i=0; i<numObjs; i++)
    {
    id = lids[i];
  
    if ((id < 0) || (id >= NumPoints))
      {
      *err = 1;
      return;
      }

    id3 = lids[i] * 3;

    pts[next++] = (double)(Points[id3]);
    pts[next++] = (double)(Points[id3 + 1]);
    pts[next++] = (double)(Points[id3 + 2]);
    }
} 

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

