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

#include <stdlib.h>
#include <ctype.h>
#include "zoltan.h"

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

/* struct for indices i & j and a value */
struct ijv
{
  int i, j;
  /* double val; */
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
 FILE *f,
 PARIO_INFO_PTR pio_info,
 int *nVtx, int *nEdge, int *nPins,
 int **index,   int **vertex,
 int *vwgt_dim, float **vwgt,
 int *ewgt_dim, float **ewgt,
 int **ch_start, int **ch_adj,
 int *base)
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
int remaining, chunksize, amt, myInCount, myMaxCount, myCount; 
int rc, edge, vtx, idx,owner, num_lone_vertices;
short *assignments=NULL;

    /* Initialize */
    *nVtx  = *nEdge  = *nPins = *vwgt_dim = *ewgt_dim = 0;

    *base = 1;   /* MatrixMarket is 1-based. */

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
      rewind(f);
      if (mm_read_banner(f, &matcode) != 0) {
          fprintf(stderr,"%s Could not process Matrix Market banner.\n",yo);
      }
      else if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
              mm_is_sparse(matcode) ) {
        /*  This is how one can screen matrix types if their application */
        /*  only supports a subset of the Matrix Market data types.      */
      
          fprintf(stderr,"%s Sorry, this application does not support ",yo);
          fprintf(stderr,"%s Market Market type: [%s]\n", 
                  mm_typecode_to_str(matcode),yo);
      }
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
    MPI_Bcast(&gnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (gnz == 0){
      return gnz;      /* failure */
    }
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
    if (pio_info->matrix_obj==COLUMNS){
      *nVtx = N;
      *nEdge = M;
    }
    else{
      *nVtx = M;
      *nEdge = N;
    }
  
    if ((Num_Proc > 1) && pio_info->chunk_reader && (gnz > Num_Proc)){
      /* Process 0 will read the input file in chunks and send
       * pins to other processes before reading the next chunk.
       * This is so we can handle very large files.
       */
      read_in_chunks = 1;
    }
  
    if (!read_in_chunks && (Proc > 0)){
      return gnz;   /* global number of non-zeroes */
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
        return 0;
      }

      chunksize = 500*1024;   /* have a better idea? */

      /* Initialize which process owns which vertex, and global
       * variables used by graph callbacks. */

      ch_dist_init(Num_Proc, *nVtx, pio_info, &assignments, 0, MPI_COMM_WORLD);

      if (Proc == 0){
        sendcount = (int *)malloc(Num_Proc * sizeof(int));
        start = (int *)malloc((Num_Proc+1) * sizeof(int));
        outVals = (int *)malloc(chunksize * 2 * sizeof(int));
        if (!sendcount || !start || !outVals){
          fprintf(stderr,"%s Memory allocation\n",yo);
          return 0;
        }
      }
      remaining = gnz;
      myCount=0;
      myMaxCount = chunksize;
      myVals = (int *)malloc(myMaxCount * 2 * sizeof(int));
      inVals = (int *)malloc(chunksize * 2 * sizeof(int));
      if (!myVals|| !inVals){
        fprintf(stderr,"%s Memory allocation\n",yo);
        return 0;
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
            fgets(line, 128, f);
            sscanf(line, "%d %d", &edge, &vtx);
            if (pio_info->matrix_obj==ROWS){
              tmp = edge;
              edge = vtx;
              vtx = edge;
            }
            edge--;           /* so we can use it as array index */
            *inptr++ = edge;
            *inptr++ = vtx;

            if (M == N){
              owner = ch_dist_proc(edge, NULL, 0);  /* yes I mean edge */
            }
            else{
              owner = edge % Num_Proc;
            }
            if ((owner < 0) || (owner >= Num_Proc)){
              printf("Confusion\n"); return 0;
            }
            sendcount[owner]++;
          }
          myInCount = sendcount[0];
          start[0] = 0;
          for (j=0; j<Num_Proc; j++){
            if (j > 0){
              MPI_Send(sendcount + j, 1, MPI_INT, j, 0x0101, MPI_COMM_WORLD);
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
                     MPI_INT, j, 0x0102, MPI_COMM_WORLD);
          }
          rc = add_new_vals(outVals, myInCount, &myVals, &myCount, &myMaxCount);
          if (rc) {
            fprintf(stderr,"Process %d out of memory\n",Proc);
            return 0;  
          }
        }
        else{
          /* Await pins from process 0 and store them in a buffer
           */
          MPI_Recv(&myInCount, 1, MPI_INT, 0, 0x0101, MPI_COMM_WORLD,&status);
          MPI_Recv(inVals, myInCount*2, MPI_INT, 0, 0x0102, MPI_COMM_WORLD,&status);
          rc = add_new_vals(inVals, myInCount, &myVals, &myCount, &myMaxCount);
          if (rc) {
            fprintf(stderr,"Process %d out of memory\n",Proc);
            return 0;  
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
      myVals = (int *) malloc(gnz * sizeof(int) * 2);
      if (gnz && !myVals){
        fprintf(stderr,"%s Memory allocation\n",yo);
        return 0;
      }
      inptr = myVals;
      for (k=0; k<gnz; k++){
        fgets(line, 128, f);
        sscanf(line, "%d %d", inptr, inptr+1);
        
        if (pio_info->matrix_obj==ROWS){
          tmp = inptr[0];
          inptr[0] = inptr[1];
          inptr[1] = tmp;
        }
        inptr[0]--;
        inptr += 2;
      }
      nz = gnz;
    }

    /************************/
    /* sort indices         */
    /************************/
    qsort(myVals, nz, sizeof(int)*2, comp);

    /************************/
    /* Debug: write out matrix */
    /************************/

    /* if ((Proc==0) && (Debug_Driver > 3)) { */
    if (0){ /* We don't know the debug level, so don't print. */
      printf("Debug: Sparse matrix in %s is:\n", yo);
      mm_write_banner(stdout, matcode);
      mm_write_mtx_crd_size(stdout, M, N, nz);
      inptr = myVals;
      for (k=0; k<nz; k++,inptr+=2)
          printf("%d %d\n", inptr[0]+1, inptr[1]);
    }

    /************************/
    /* Populate Zoltan hg data structs. */
    /************************/

    /* Weights not supported. */
    *vwgt_dim = 0;
    *ewgt_dim = 0;

    *nPins = nz;

    /* allocate storage for hypergraph data */    
    if (!(*index  = (int*) malloc ((*nEdge+1) * sizeof(int)))
     || !(*vertex = (int*) malloc  (*nPins   * sizeof(int)))) {
         fprintf(stderr, "%s Insufficient memory.", yo);
         gnz = 0;
         goto End;
    }

    /* Construct index and vertex arrays from mat (ijv) array */

    /* Initialize index array to -1. */

    for (j=0; j<(*nEdge+1); j++)
      (*index)[j] = -1;
    /* Loop over the nonzeros. For each hyperedge (same index i), 
       make list of vertices (index j). Data are sorted by i index. */
    prev_edge = -1;
    num_lone_vertices=0;
    for (k=0,inptr=myVals; k<(*nPins); k++,inptr+=2){
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
      start = NULL;   /* index of start of vertices' edge lists */
      int *adj = NULL;     /* list of adjacent vertices */
     
      cnt = (int *) calloc(N, sizeof(int));
      start = (int *) malloc((N+1) * sizeof(int));
      if (N && (!cnt || !start)) {
        fprintf(stderr, "%s Insufficient memory.", yo);
        gnz = 0;
        goto End;
      }
      
      /* Count number of neighbors for each vertex */

      for (k = 0,inptr=myVals; k < nz; k++,inptr+=2) {
        if (inptr[0]+1 != inptr[1]) { /* Don't include self-edges */
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
        gnz = 0;
        goto End;
      }
      for (k = 0,inptr=myVals; k < nz; k++,inptr+=2) {
        if (inptr[0]+1 != inptr[1]) { /* Don't include self-edges */
          adj[start[inptr[0]]+cnt[inptr[0]]] = inptr[1];
          cnt[inptr[0]]++;
        }
      }
      safe_free((void **) &cnt);
      *ch_start = start;
      *ch_adj = adj;
    }  /* N == M */

End:
    if (gnz == 0){
      /* Set hypergraph to be empty. */
      *nVtx  = *nEdge  = *nPins = *vwgt_dim = *ewgt_dim = 0;
      safe_free((void **) &index);
      safe_free((void **) &vertex);
    }
    safe_free((void **) &myVals);
    return gnz;
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
    *myVals = realloc(*myVals, newsize * 2 * sizeof(int));

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
