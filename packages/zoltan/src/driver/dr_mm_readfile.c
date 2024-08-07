// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include <stdlib.h>
#include <ctype.h>
#include "zoltan.h"

#include "dr_const.h"
#include "dr_hg_readfile.h"
#include "dr_util_const.h"
#include "dr_input_const.h" /* just for the matrix_obj define's */
#include "dr_mmio.h"
#include "ch_input_const.h"
#include "ch_init_dist_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* struct for indices i & j and a real value */
struct ijv
{
  int i, j;
  float v;
};

/* comparison functions. */
static int comp(const void *, const void *);

static int add_new_vals(int *newvals, int newCount,
       int **myVals, int *myCount, int *myMaxCount);

/*****************************************************************************/
/*****************************************************************************/
/* Read MatrixMarket sparse matrix. Weights currently not supported. */
/* Must be called by all processes.                                  */
/* Returns global number of pins.                                    */

int MM_readfile (
 int Proc,
 int Num_Proc,
 ZOLTAN_FILE* f,
 PARIO_INFO_PTR pio_info,
 int *nVtx, int *nEdge, int *nPins,
 int **index,   int **vertex,
 int *vwgt_dim, float **vwgt,
 int *ewgt_dim, float **ewgt,  /* hyperedge weights */
 int **ch_start, int **ch_adj,
 int *ch_ewgt_dim, float **ch_ewgts, /* graph weights */
 int *base,
 int *global_nPins)
{
int prev_edge;
const char *yo = "MM_readfile";
MM_typecode matcode;
int M, N, nz, gnz;
int j, k, tmp;
MPI_Status status;
int read_in_chunks = 0;
char line[128];
int *sendcount = NULL, *start = NULL, *outVals  = NULL;
int *myVals= NULL, *inVals= NULL, *inptr = NULL;
struct ijv *myIJV = NULL, *iptr = NULL;
int remaining, chunksize, amt, myInCount, myMaxCount, myCount; 
int rc, edge, vtx, idx,owner, num_lone_vertices;
short *assignments=NULL;
int error = 0;  /* flag to indicate status */

    /* Initialize */
    *nVtx  = *nEdge  = *nPins = *vwgt_dim = *ewgt_dim = 0;
    *index = *vertex = *ch_start = *ch_adj = NULL;
    *vwgt = *ewgt = NULL;
    *base = 0;   /* MatrixMarket is 1-based but I substract this 1 already */

/*
*  This code was adapted from NIST's MatrixMarket examples, see
*    http://math.nist.gov/MatrixMarket/formats.html
*
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/

    M=N=gnz=0;

    if (Proc == 0){
      ZOLTAN_FILE_rewind(f);
      if (mm_read_banner(f, &matcode) != 0) {
          fprintf(stderr,"%s Could not process Matrix Market banner.\n",yo);
      }
/*       else if (mm_is_complex(matcode) && mm_is_matrix(matcode) && */
/*               mm_is_sparse(matcode) ) { */
/*         /\*  This is how one can screen matrix types if their application *\/ */
/*         /\*  only supports a subset of the Matrix Market data types.      *\/ */

/*           fprintf(stderr,"%s Sorry, this application does not support ",yo); */
/*           fprintf(stderr,"%s Market Market type: [%s]\n", */
/*                   mm_typecode_to_str(matcode),yo); */
/*       } */
      else if (mm_read_mtx_crd_size(f, &M, &N, &gnz) !=0){
        /* find out size of sparse matrix .... */
          fprintf(stderr,"%s Can't read size of mtx file\n", yo);
          M=N=gnz=0;
      }
      else if ((pio_info->matrix_obj!=COLUMNS) && (pio_info->matrix_obj!=ROWS)){
        fprintf(stderr, "%s Invalid option for matrix objects.", yo);
        M=N=gnz=0;
      }
    }
    MPI_Bcast(&gnz, 1, MPI_INT, 0, zoltan_get_global_comm());
    MPI_Bcast(&M, 1, MPI_INT, 0, zoltan_get_global_comm());
    MPI_Bcast(&N, 1, MPI_INT, 0, zoltan_get_global_comm());

    if (pio_info->matrix_obj==COLUMNS){
      *nVtx = N;
      *nEdge = M;
    }
    else{
      *nVtx = M;
      *nEdge = N;
    }
    *global_nPins = gnz;

    if ((Num_Proc > 1) && pio_info->chunk_reader && (gnz > Num_Proc)){
      /* Process 0 will read the input file in chunks and send
       * pins to other processes before reading the next chunk.
       * This is so we can handle very large files.
       */
      read_in_chunks = 1;
    }

    if (!read_in_chunks && (Proc > 0)){
      return error;   /* global number of non-zeroes */
    }

    if (read_in_chunks){

      /* Process 0 reads in a portion of the file, and sends
       * pins to other processes using a scheme that ensures
       * edges are not divided across processes.  (We need that
       * so that we can build the graph in ch_adj later on.)
       */

      if (pio_info->init_dist_type == INITIAL_FILE) {
        /* It was just too complicated to get this working */
        fprintf(stderr,
     "fatal: Chaco assignments not supported for very large hypergraph files");
        error = 1;
        return error;
      }

      chunksize = 500*1024;   /* have a better idea? */

      /* Initialize which process owns which vertex, and global
       * variables used by graph callbacks. */

      ch_dist_init(Num_Proc, *nVtx, pio_info, &assignments, 0, zoltan_get_global_comm());

      if (Proc == 0){
        sendcount = (int *)malloc(Num_Proc * sizeof(int));
        start = (int *)malloc((Num_Proc+1) * sizeof(int));
        outVals = (int *)malloc(chunksize * 2 * sizeof(int));
        if (!sendcount || !start || !outVals){
          fprintf(stderr,"%s Memory allocation\n",yo);
          error = 1;
          return error;
        }
      }
      remaining = gnz;
      myCount=0;
      myMaxCount = chunksize;
      myVals = (int *)malloc(myMaxCount * 2 * sizeof(int));
      inVals = (int *)malloc(chunksize * 2 * sizeof(int));
      if (!myVals|| !inVals){
        fprintf(stderr,"%s Memory allocation\n",yo);
        error = 1;
        return error;
      }

      while (remaining > 0){
        amt = (remaining < chunksize) ? remaining : chunksize;
        if (Proc == 0){
          /* Read in a portion of disk file, and send pins to
           * other processes.
           */
          memset(sendcount, 0, sizeof(int) * Num_Proc);
          inptr = inVals;
          for (k=0; k<amt; k++){
            ZOLTAN_FILE_gets(line, 128, f);
            sscanf(line, "%d %d", &edge, &vtx);
            if (pio_info->matrix_obj==ROWS){
              tmp = edge;
              edge = vtx;
              vtx = tmp;
            }
            edge--;           /* so we can use it as array index */
	    vtx--;            /* Want the same numbering, start from 0 */
            *inptr++ = edge;
            *inptr++ = vtx;

            if (M == N){
              owner = ch_dist_proc(edge, NULL, 0);  /* yes I mean edge */
            }
            else{
              owner = edge % Num_Proc;
            }
            if ((owner < 0) || (owner >= Num_Proc)){
              printf("Confusion\n"); error = 1; return error;
            }
            sendcount[owner]++;
          }
          myInCount = sendcount[0];
          start[0] = 0;
          for (j=0; j<Num_Proc; j++){
            if (j > 0){
              MPI_Send(sendcount + j, 1, MPI_INT, j, 0x0101, zoltan_get_global_comm());
            }
            start[j+1] = start[j] + sendcount[j];
            sendcount[j] = 0;
          }
          inptr = inVals;
          for (k=0; k<amt; k++){
            edge = *inptr++;
            vtx = *inptr++;
            if (M == N){
              owner = ch_dist_proc(edge, NULL, 0);  /* yes I mean edge */
            }
            else{
              owner = edge % Num_Proc;
            }
            idx = start[owner] + sendcount[owner]++;
            outVals[idx*2] = edge;
            outVals[idx*2 + 1] = vtx;
          }
          for (j=1; j<Num_Proc; j++){
            MPI_Send(outVals + (2*start[j]), sendcount[j] * 2,
                     MPI_INT, j, 0x0102, zoltan_get_global_comm());
          }
          rc = add_new_vals(outVals, myInCount, &myVals, &myCount, &myMaxCount);
          if (rc) {
            fprintf(stderr,"Process %d out of memory\n",Proc);
            error = 1;
            return error;
          }
        }
        else{
          /* Await pins from process 0 and store them in a buffer
           */
          MPI_Recv(&myInCount, 1, MPI_INT, 0, 0x0101, zoltan_get_global_comm(),&status);
          MPI_Recv(inVals, myInCount*2, MPI_INT, 0, 0x0102, zoltan_get_global_comm(),&status);
          rc = add_new_vals(inVals, myInCount, &myVals, &myCount, &myMaxCount);
          if (rc) {
            fprintf(stderr,"Process %d out of memory\n",Proc);
            error = 1;
            return error;
          }
        }
        remaining -= amt;
      }

      if (Proc == 0){
        free(sendcount);
        free(start);
        free(outVals);
      }
      free(inVals);

      nz = myCount;
    }
    else{
      /* Not reading in chunks.
       * Process 0 will read in the entire file.
       */
      myIJV = (struct ijv *) malloc(gnz * sizeof(struct ijv));
      if (gnz && !myIJV){
        fprintf(stderr,"%s Memory allocation\n",yo);
        error = 1;
        return error;
      }

      for (k=0; k<gnz; k++){
        ZOLTAN_FILE_gets(line, 128, f);
        if (mm_is_pattern(matcode)) {
          sscanf(line, "%d %d", &(myIJV[k].i), &(myIJV[k].j));
          myIJV[k].v = 1.;
        }
        else {
          sscanf(line, "%d %d %f", &(myIJV[k].i), &(myIJV[k].j),
                 &(myIJV[k].v));
        }

        if (pio_info->matrix_obj==ROWS){
          tmp = myIJV[k].i;
          myIJV[k].i = myIJV[k].j;
          myIJV[k].j = tmp;
        }
        --(myIJV[k].i);
        --(myIJV[k].j); /* Want to start at 0 */
      }
      nz = gnz;
    }

    /************************/
    /* sort indices         */
    /************************/
    if (myIJV)
      qsort(myIJV, nz, sizeof(struct ijv), comp);
    else  /* chunky style: (i,j) pairs are packed in myVals */
      qsort(myVals, nz, 2*sizeof(int), comp);

    /************************/
    /* Debug: write out matrix */
    /************************/

    /* if ((Proc==0) && (Debug_Driver > 3)) { */
    if (0){ /* We don't know the debug level, so don't print. */
      printf("Debug: Sparse matrix in %s is:\n", yo);
      mm_write_banner(stdout, matcode);
      mm_write_mtx_crd_size(stdout, M, N, nz);
      for (k=0; k<nz; k++)
          printf("%d %d %f\n", myIJV[k].i+1, myIJV[k].j+1, myIJV[k].v);
    }

    /************************/
    /* Populate Zoltan hg data structs. */
    /************************/

    /* Default is no weights. For symmetric matrices, generate weights later! */
    *ewgt_dim = 0;
    *vwgt_dim = 0;

    *nPins = nz;

    /* allocate storage for hypergraph data */
    if (!(*index  = (int*) malloc ((*nEdge+1) * sizeof(int)))
     || (*nPins && !(*vertex = (int*) malloc  (*nPins   * sizeof(int))))) {
         fprintf(stderr, "%s Insufficient memory.", yo);
         error = 1;
         goto End;
    }

    /* Construct index and vertex arrays from matrix (myVals) */

    /* Initialize index array to -1. */

    for (j=0; j<(*nEdge+1); j++)
      (*index)[j] = -1;
    /* Loop over the nonzeros. For each hyperedge (same index i),
       make list of vertices (index j). Data are sorted by i index. */
    prev_edge = -1;
    num_lone_vertices=0;
    for (k=0; k<(*nPins); k++){
      if (myIJV)
        inptr=(int*)(myIJV+k);
      else
        inptr= myVals+2*k;
      /* next edge is mat[k].i */
      if (inptr[0] > prev_edge) {
        (*index)[inptr[0]] = k;
        num_lone_vertices++;
      }
      prev_edge = inptr[0];
      (*vertex)[k] = inptr[1];
    }
    /* Complete index array; overwrite -1s. */
    tmp = *nPins;
    for (j=(*nEdge); j; j--){
      if ((*index)[j] == -1)
        (*index)[j] = tmp;
      else
        tmp = (*index)[j];
    }
    (*index)[0] = 0;

    /*  Build Graph data structures, too, if matrix is square (N==M). */
    /*  Assume the input matrix is symmetric. */

    /*  If N==M and the graph is not symmetric we get errors later on
     *  in migration.  This happens with test directory hg_FluidDFT.
     */

    if (N == M){                   /* Square matrix */
      int sum;
      int *cnt = NULL;
      int *adj = NULL;     /* list of adjacent vertices */
      float *vwgts = NULL;     /* vertex weights for diagonal entries */
      float *ewgts = NULL;     /* edge weights for off-diagonals */
      start = NULL;   /* index of start of vertices' edge lists */

      /* Assume symmetric matrix. Always create weights, but they are
         used only if user sets Zoltan parameter obj_weight_dim=1. */
      *ch_ewgt_dim = 1;
      *vwgt_dim = 1;

      cnt = (int *) calloc(N, sizeof(int));
      start = (int *) malloc((N+1) * sizeof(int));
      if (*vwgt_dim)
        vwgts = (float *) calloc((*vwgt_dim) * N, sizeof(float));
      if (*ch_ewgt_dim)
        ewgts = (float *) calloc((*ch_ewgt_dim) * (*nPins), sizeof(float));
      if ((N && !cnt) || !start || ((*vwgt_dim && N) && !vwgts)
         || ((*ch_ewgt_dim && *nPins) && !ewgts)) {
        fprintf(stderr, "%s Insufficient memory.", yo);
        error = 1;
        goto End;
      }

      /* Count number of neighbors for each vertex */
      if (myIJV)
        for (k = 0,iptr=myIJV; k < nz; ++k, ++iptr) {
          if (iptr->i != iptr->j) { /* Don't include self-edges */
            cnt[iptr->i]++;
            if (cnt[iptr->i] == 1) num_lone_vertices--;
          }
        }
      else /* chunky style, use myVals */
        for (k = 0,inptr=myVals; k < nz; ++k, inptr+=2) {
          if (inptr[0] != inptr[1]) { /* Don't include self-edges */
            cnt[inptr[0]]++;
            if (cnt[inptr[0]] == 1) num_lone_vertices--;
          }
        }
      if (num_lone_vertices){
        printf("WARNING (Proc %d): in graph generated from hg, have %d vertices with no adjacencies\n",Proc,num_lone_vertices);
      }

      /* Compute start:  prefix sum of cnt */
      /* Also reset cnts to zero. */
      start[0] = 0;
      for (sum = 0, k = 1; k <= N; k++) {
        sum += cnt[k-1];
        start[k] = sum;
        cnt[k-1]=0;
      }

      /* Load the array of adjacencies */
      adj = (int *) malloc(sum * sizeof(int));
      if (sum && !adj) {
        fprintf(stderr, "%s Insufficient memory.", yo);
        error = 1;
        goto End;
      }
      if (myIJV)
        for (k = 0,iptr=myIJV; k < nz; ++k, ++iptr) {
          if (iptr->i == iptr->j) { /* Diagonal entry */
            if (*vwgt_dim)
              vwgts[iptr->i] = iptr->v;
          }
          else { /* Off-diagonal */
            idx = start[iptr->i]+cnt[iptr->i];
            adj[idx] = iptr->j;
            if (*ch_ewgt_dim) 
              ewgts[idx] = iptr->v;
            cnt[iptr->i]++;
          }
        }
      else /* chunky style, use myVals */
        for (k = 0,inptr=myVals; k < nz; ++k, inptr+=2) {
          if (inptr[0] == inptr[1]) { /* Diagonal entry */
            /* ignore self edges */
          }
          else { /* Off-diagonal */
            adj[start[inptr[0]]+cnt[inptr[0]]] = inptr[1];
            cnt[inptr[0]]++;
          }
      }
      safe_free((void **)(void *) &cnt);
      *ch_start = start;
      *ch_adj = adj;
      *ch_ewgts = ewgts;
      *vwgt = vwgts;
#if 0
      /* DEBUG */
      if (vwgts){
        printf("zdrive debug, vertex weights: ");
        for (k=0; k< 5; k++)
          printf("%f ", vwgts[k]);
        printf("\n");
        printf("zdrive debug, edge weights: ");
      }
      if (ewgts){
        for (k=0; k< 5; k++)
          printf("%f ", ewgts[k]);
        printf("\n");
      }
#endif
    }  /* N == M */

End:
    if (error){
      /* Set hypergraph to be empty. */
      *nVtx  = *nEdge  = *nPins = *vwgt_dim = *ewgt_dim = 0;
      *global_nPins = 0;
      safe_free((void **)(void *) index);
      safe_free((void **)(void *) vertex);
    }
    safe_free((void **)(void *) &myIJV);
    safe_free((void **)(void *) &myVals);
    return error;
}


/* sort by increasing i value */
/* secondary key is j value */
static int comp(const void *a, const void *b)
{
  int ai = ((int *)a)[0];
  int aj = ((int *)a)[1];
  int bi = ((int *)b)[0];
  int bj = ((int *)b)[1];

  if ((ai - bi) == 0)
    return (aj - bj);
  else
    return (ai - bi);
}
static int add_new_vals(int *newvals, int newCount, 
       int **myVals, int *myCount, int *myMaxCount)
{
int newsize;

  if (newCount == 0) return 0;

  if (*myCount + newCount > *myMaxCount){
    newsize = *myMaxCount + newCount*2;
    *myVals = (int *)realloc(*myVals, newsize * 2 * sizeof(int));

    if (! *myVals){
      return 1;
    }
    *myMaxCount = newsize;
  }

  memcpy(*myVals + (*myCount * 2), newvals, newCount * 2 * sizeof(int));

  *myCount += newCount;

  return 0;
}

#ifdef __cplusplus
}
#endif
