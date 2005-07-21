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

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* struct for indices i & j and a value */
struct ijv
{
  int i, j;
};

/* comparison functions. */
static int comp(const void *, const void *);


/*****************************************************************************/
/*****************************************************************************/
/* Read MatrixMarket sparse matrix. Weights currently not supported. */

int MM_readfile (
 int Proc,
 FILE *f,
 int *nVtx, int *nEdge, int *nPins,
 int **index,   int **vertex,
 int *vwgt_dim, float **vwgt,
 int *ewgt_dim, float **ewgt,
 int *base, int matrix_obj)
{
int err = ZOLTAN_OK;
int prev_edge;
const char *yo = "MM_readfile";
MM_typecode matcode;
int M, N, nz;   
int j, k, tmp;
struct ijv *mat;

    /* Initialize */
    *nVtx  = *nEdge  = *nPins = *vwgt_dim = *ewgt_dim = 0;

    *base = 1;   /* MatrixMarket is 1-based. */
    rewind(f);   /* need to read first line again! */
          

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



    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if (mm_read_mtx_crd_size(f, &M, &N, &nz) !=0)
        exit(1);


    /* reserve memory for matrices */

    mat = (struct ijv *) malloc(nz * sizeof(struct ijv));

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (k=0; k<nz; k++) {
        char line[1000];
        fgets(line, 1000, f);
        sscanf(line, "%d %d", &mat[k].i, &mat[k].j);
        if (matrix_obj==ROWS){
          /* swap (i,j) to transpose matrix and switch vertices/hyperedges */
          tmp = mat[k].i;
          mat[k].i = mat[k].j;
          mat[k].j = tmp;
        }
        (mat[k].i)--;  
    }

    /* if (f !=stdin) fclose(f); */

    /************************/
    /* sort indices         */
    /************************/
    qsort(mat, nz, sizeof(struct ijv), comp);

    /************************/
    /* Debug: write out matrix */
    /************************/

    /* if ((Proc==0) && (Debug_Driver > 3)) { */
    if (0){ /* We don't know the debug level, so don't print. */
      printf("Debug: Sparse matrix in %s is:\n", yo);
      mm_write_banner(stdout, matcode);
      mm_write_mtx_crd_size(stdout, M, N, nz);
      for (k=0; k<nz; k++)
          printf("%d %d\n", mat[k].i+1, mat[k].j);
    }

    /************************/
    /* Populate Zoltan hg data structs. */
    /************************/

    /* Weights not supported. */
    *vwgt_dim = 0;
    *ewgt_dim = 0;

    /* Rows or columns are hedges? */
    if (matrix_obj==COLUMNS){
      *nVtx = N;
      *nEdge = M;
    }
    else if (matrix_obj==ROWS){
      *nVtx = M;
      *nEdge = N;
    }
    else {
       fprintf(stderr, "%s Invalid option for matrix objects.", yo);
       err = ZOLTAN_FATAL;
       goto End;
    }
    *nPins = nz;


    /* allocate storage for hypergraph data */    
    if (!(*index  = (int*) malloc ((*nEdge+1) * sizeof(int)))
     || !(*vertex = (int*) malloc  (*nPins   * sizeof(int)))) {
         fprintf(stderr, "%s Insufficient memory.", yo);
         err = ZOLTAN_MEMERR;
         goto End;
    }

    /* Construct index and vertex arrays from mat (ijv) array */

    /* Initialize index array to -1. */
    for (j=0; j<(*nEdge+1); j++)
      (*index)[j] = -1;
    /* Loop over the nonzeros. For each hyperedge (same index i), 
       make list of vertices (index j). Data are sorted by i index. */
    prev_edge = -1;
    for (k=0; k<(*nPins); k++){
      /* next edge is mat[k].i */
      if (mat[k].i > prev_edge) 
        (*index)[mat[k].i] = k;
      prev_edge = mat[k].i;
      (*vertex)[k] = mat[k].j;
    }
    /* Complete index array; overwrite -1s. */
    tmp = *nPins;
    for (j=(*nEdge); j; j--){
      if ((*index)[j] == -1)
        (*index)[j] = tmp;
      else
        tmp = (*index)[j];
    }

End:
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
      /* Set hypergraph to be empty. */
      *nVtx  = *nEdge  = *nPins = *vwgt_dim = *ewgt_dim = 0;
      safe_free((void **) &index);
      safe_free((void **) &vertex);
    }
    safe_free((void **) &mat);
    return err;
}


/* sort by increasing i value */
/* secondary key is j value */
static int comp(const void *a, const void *b)
{
  if ((((struct ijv *)a)->i - ((struct ijv *)b)->i) == 0)
    return ((struct ijv *)a)->j - ((struct ijv *)b)->j;
  else
    return ((struct ijv *)a)->i - ((struct ijv *)b)->i;
}

#ifdef __cplusplus
}
#endif
