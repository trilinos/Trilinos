/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * tio.c
 *
 * This file contains routines related to I/O
 *
 * Started 10/19/94
 * George
 *
 * $Id$
 *
 */

#include <parmetis.h>


/*************************************************************************
* This function reads the CSR matrix
**************************************************************************/
void ReadTestGraph(GraphType *graph, char *filename, MPI_Comm comm)
{
  int i, j, k, l, npes, mype;
  int nvtxs, nedges, penum, snedges, snvtxs;
  idxtype *gxadj, *gadjncy;  
  idxtype *vtxdist, *sxadj, *ssize;
  MPI_Status status;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  vtxdist = graph->vtxdist = idxsmalloc(npes+1, 0, "ReadGraph: vtxdist");

  if (mype == 0) {
    ssize = idxsmalloc(npes, 0, "ReadGraph: ssize");

    ReadMetisGraph(filename, &nvtxs, &gxadj, &gadjncy);

    printf("Nvtxs: %d, Nedges: %d\n", nvtxs, gxadj[nvtxs]);

    /* Construct vtxdist and send it to all the processors */
    vtxdist[0] = 0;
    for (i=0,k=nvtxs; i<npes; i++) {
      l = k/(npes-i);
      vtxdist[i+1] = vtxdist[i]+l;
      k -= l;
    }
  }

  MPI_Bcast((void *)vtxdist, npes+1, IDX_DATATYPE, 0, comm);

  graph->gnvtxs = vtxdist[npes];
  graph->nvtxs = vtxdist[mype+1]-vtxdist[mype];
  graph->xadj = idxmalloc(graph->nvtxs+1, "ReadGraph: xadj");

  if (mype == 0) {
    for (penum=0; penum<npes; penum++) {
      snvtxs = vtxdist[penum+1]-vtxdist[penum];
      sxadj = idxmalloc(snvtxs+1, "ReadGraph: sxadj");

      idxcopy(snvtxs+1, gxadj+vtxdist[penum], sxadj);
      for (i=snvtxs; i>=0; i--)
        sxadj[i] -= sxadj[0];

      ssize[penum] = gxadj[vtxdist[penum+1]] - gxadj[vtxdist[penum]];

      if (penum == mype) 
        idxcopy(snvtxs+1, sxadj, graph->xadj);
      else
        MPI_Send((void *)sxadj, snvtxs+1, IDX_DATATYPE, penum, 1, comm); 

      free(sxadj);
    }
  }
  else 
    MPI_Recv((void *)graph->xadj, graph->nvtxs+1, IDX_DATATYPE, 0, 1, comm, &status);


  graph->nedges = graph->xadj[graph->nvtxs];
  graph->adjncy = idxmalloc(graph->nedges, "ReadGraph: graph->adjncy");

  if (mype == 0) {
    for (penum=0; penum<npes; penum++) {
      if (penum == mype) 
        idxcopy(ssize[penum], gadjncy+gxadj[vtxdist[penum]], graph->adjncy);
      else
        MPI_Send((void *)(gadjncy+gxadj[vtxdist[penum]]), ssize[penum], IDX_DATATYPE, penum, 1, comm); 
    }

    free(ssize);
  }
  else 
    MPI_Recv((void *)graph->adjncy, graph->nedges, IDX_DATATYPE, 0, 1, comm, &status);

  graph->vwgt = NULL;
  graph->adjwgt = NULL;

  if (mype == 0) 
    GKfree(&gxadj, &gadjncy, LTERM);

  MALLOC_CHECK(NULL);
}



/*************************************************************************
* This function reads the CSR matrix
**************************************************************************/
float *ReadTestCoordinates(GraphType *graph, char *filename, int ndims, MPI_Comm comm)
{
  int i, j, k, l, npes, mype;
  int nvtxs, nedges, penum, snedges, snvtxs;
  float *xyz, *txyz;
  FILE *fpin;
  idxtype *vtxdist, *sxadj, *sadjncy, *ssize, scratch;
  MPI_Status status;
  char xyzfile[256];

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  vtxdist = graph->vtxdist;

  xyz = fmalloc(graph->nvtxs*ndims, "io");

  if (mype == 0) {
    sprintf(xyzfile, "%s.xyz", filename);
    if ((fpin = fopen(xyzfile, "r")) == NULL) 
      errexit("Failed to open file %s\n", xyzfile);
  }

  if (mype == 0) {
    txyz = fmalloc(2*graph->nvtxs*ndims, "io");

    for (penum=0; penum<npes; penum++) {
      for (k=0, i=vtxdist[penum]; i<vtxdist[penum+1]; i++, k++) {
        for (j=0; j<ndims; j++)
          fscanf(fpin, "%e ", txyz+k*ndims+j);
      }

      if (penum == mype) 
        memcpy((void *)xyz, (void *)txyz, sizeof(float)*ndims*k);
      else {
        MPI_Send((void *)txyz, ndims*k, MPI_FLOAT, penum, 1, comm); 
      }
    }
    free(txyz);
    fclose(fpin);
  }
  else 
    MPI_Recv((void *)xyz, ndims*graph->nvtxs, MPI_FLOAT, 0, 1, comm, &status);

  return xyz;
}



/*************************************************************************
* This function reads the spd matrix
**************************************************************************/
void ReadMetisGraph(char *filename, int *r_nvtxs, idxtype **r_xadj, idxtype **r_adjncy)
{
  int i, j, k, l, edge, nvtxs, nedges;
  idxtype *xadj, *adjncy, *vwgt, *adjwgt;
  char *line, *oldstr, *newstr;
  FILE *fpin;

  line = (char *)malloc(sizeof(char)*(8192+1));

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fgets(line, 8192, fpin);
  sscanf(line, "%d %d", &nvtxs, &nedges);
  nedges *=2;

  xadj = idxmalloc(nvtxs+1, "ReadGraph: xadj");
  adjncy = idxmalloc(nedges, "ReadGraph: adjncy");

  /* Start reading the graph file */
  for (xadj[0]=0, k=0, i=0; i<nvtxs; i++) {
    fgets(line, 8192, fpin);
    oldstr = line;
    newstr = NULL;

    for (;;) {
      edge = (int)strtol(oldstr, &newstr, 10) -1;
      oldstr = newstr;

      if (edge < 0)
        break;

      adjncy[k++] = edge;
    } 
    xadj[i+1] = k;
  }

  fclose(fpin);

  free(line);

  *r_nvtxs = nvtxs;
  *r_xadj = xadj;
  *r_adjncy = adjncy;
}
