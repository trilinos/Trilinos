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

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdlib.h>
#include <ctype.h>
#include "zoltan.h"

#include "dr_hg_readfile.h"
#include "dr_util_const.h"
#include "dr_mmio.h"


#define BUF_LEN 1000
#define ERROR(proc,yo,msg,err) \
 {printf("Proc(%d) From(%s): %s\n",proc,yo,msg);return(err);}

/* struct for indices i & j and a value */
struct ijv
{
  int i, j;
  double val;
};

/* comparison functions. */
static int comp(const void *, const void *);


static int   readfile  (int, FILE*, int*, int*, int*, int**, int**, int*,
 float**, int*, float**, int*);
  
static int MM_readfile (int, FILE*, int*, int*, int*, int**, int**, int*,
 float**, int*, float**, int*);

static int nextstr (FILE *f, char *string);




/*****************************************************************************/

int Zoltan_HG_Readfile (
  int Proc,
  FILE *f,
  int *nVtx, int *nEdge, int *nPins,
  int **hindex,  int **hvertex,
  int *vwgt_dim, float **vwgt,
  int *ewgt_dim, float **ewgt,
  int *base)
{
  char string[BUF_LEN];
  char *yo = "Zoltan_HG_Readfile";

  /* Initialize return values in case of error. */
  *nVtx  = *nEdge  = *nPins = *vwgt_dim = *ewgt_dim = *base = 0;
  *hindex = *hvertex = NULL;
  *vwgt = *ewgt = NULL;
  
  while (1)
     {
     fgets (string, BUF_LEN-1, f);
     if (strncmp (string, "%%Matrix", 8) == 0)
        return MM_readfile (Proc, f, nVtx, nEdge, nPins, hindex, hvertex,
         vwgt_dim, vwgt, ewgt_dim, ewgt, base);     
     
     if (*string != '%')
        break;
     }
     
  /* Extended format: nVtx nEdge nPins NM where both N&M are [0,9]. N is the
   * number of dimension for vertex weights and M is number of dimensions for
   * the edge weights. */   
  if (sscanf(string, "%d %d %d %1d%1d", nVtx,nEdge,nPins,vwgt_dim,ewgt_dim) != 5)
    ERROR (Proc, yo, "Unrecognized file format.", ZOLTAN_FATAL);

  if (*nVtx <= 1)                    /* wrong guess if nVtx really is one! */
     {
     sscanf (string, "%d %d %d %d %d", base, nVtx, nEdge, nPins, vwgt_dim);
     return readfile (Proc, f, nVtx, nEdge, nPins, hindex, hvertex, vwgt_dim, 
      vwgt, ewgt_dim, ewgt, base);                         /* patoh format */
     }
  else
     {
     *base = 1;
     return readfile (Proc, f, nVtx, nEdge, nPins, hindex, hvertex, vwgt_dim,
      vwgt, ewgt_dim, ewgt, base);                         /* IBM format */
     }
}


/*****************************************************************************/
/* Read IBM, Patoh, & Chaco like hypergraph formats */

static int readfile (
 int Proc,
 FILE *f,
 int *nVtx, int *nEdge, int *nPins,
 int **index,   int **vertex,
 int *vwgt_dim, float **vwgt,
 int *ewgt_dim, float **ewgt,
 int *base)
{
  int i, j, pin;
  char string[BUF_LEN];
  float *p;
  char *yo = "readfile";

  if (!(*index  = (int*) malloc ((*nEdge+1) * sizeof(int)))
   || !(*vertex = (int*) malloc ((*nPins    * sizeof(int)))))
      ERROR (Proc, yo, "Insufficient memory", ZOLTAN_MEMERR);
      
  if (*vwgt_dim > 0 
   && !(*vwgt = (float*) malloc (*nVtx * *vwgt_dim * sizeof (float))))
      ERROR (Proc, yo, "Insufficient memory", ZOLTAN_MEMERR)  
  
  if (*ewgt_dim > 0 
   && !(*ewgt = (float*) malloc (*nEdge * *ewgt_dim * sizeof (float))))
      ERROR (Proc, yo, "Insufficient memory", ZOLTAN_MEMERR)
             
  pin = 0;
  p = *ewgt;  
  for (i = 0; i < *nEdge; i++)
     {
     (*index)[i] = pin;
     for (j = 0; j < *ewgt_dim; j++)
        {
        if (nextstr (f, string) != 1)
           ERROR (Proc, yo, "Wrong number of edge weights", ZOLTAN_FATAL)
        *p++ = atof (string);      
        }      
     while (nextstr (f,string) == 1)
        (*vertex)[pin++] = atoi(string);
     }
  (*index)[i] = pin; 
  if (pin != *nPins)
     ERROR (Proc, yo, "Assertion Error in readfile.", ZOLTAN_FATAL);

  p = *vwgt;     
  if (*vwgt_dim > 0)
     for (i = 0; i < *nVtx; i++)
        {
        for (j = 0; nextstr (f, string) == 1 && j < *vwgt_dim; j++)
           *p++ = atof (string);
        if (j != *vwgt_dim)
           ERROR (Proc, yo, "Wrong number of vertex weights", ZOLTAN_FATAL)    
        }
        
  if (nextstr (f, string) == 1)
     ERROR (Proc, yo, "Input file is longer than expected", ZOLTAN_FATAL)      
  return ZOLTAN_OK;
}       



/*****************************************************************************/
/* Read MatrixMarket sparse matrix. Weights currently not supported. */

static int MM_readfile (
 int Proc,
 FILE *f,
 int *nVtx, int *nEdge, int *nInput,
 int **index,   int **vertex,
 int *vwgt_dim, float **vwgt,
 int *ewgt_dim, float **ewgt,
 int *base)
{
int err = ZOLTAN_OK;
int prev_edge;
int rowhedges=1; /* default is row hyperedge model */
char *yo = "MM_readfile";
int ret_code;
MM_typecode matcode;
int M, N, nz;   
int j, k, tmp;
struct ijv *mat;

    *base = 0;   /* MatrixMarket is 1-based, but we convert to 0-based. */
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

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    /* reserve memory for matrices */

    mat = (struct ijv *) malloc(nz * sizeof(struct ijv));

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (k=0; k<nz; k++)
    {
        fscanf(f, "%d %d %lg\n", &mat[k].i, &mat[k].j, &mat[k].val);
        (mat[k].i)--;  /* adjust from 1-based to 0-based */
        (mat[k].j)--;
        if (!rowhedges){
          /* swap (i,j) to transpose matrix and make rows hyperedges. */
          tmp = mat[k].i;
          mat[k].i = mat[k].j;
          mat[k].j = tmp;
        }
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
          printf("%d %d %20.19g\n", mat[k].i+1, mat[k].j+1, mat[k].val);
    }

    /************************/
    /* Populate Zoltan hg data structs. */
    /************************/

    /* Weights not supported. */
    *vwgt_dim = 0;
    *ewgt_dim = 0;

    /* Rows or columns are hedges? */
    if (rowhedges){
      *nVtx = N;
      *nEdge = M;
    }
    else{
      *nVtx = M;
      *nEdge = N;
    }
    *nInput = nz;


    /* allocate storage for hypergraph data */    
    if (!(*index  = (int*) malloc ((*nEdge+1) * sizeof(int)))
     || !(*vertex = (int*) malloc  (*nInput   * sizeof(int)))) {
         fprintf(stderr, "%s Insufficient memory.", yo);
         err = ZOLTAN_MEMERR;
         goto End;
    }

    /* Construct index and vertex arrays from mat (ijv) array */

    /* Initialize index array to -1. */
    for (j=0; j<(*nEdge+1); j++)
      (*index)[j] = -1;
    /* Loop over the nonzeros. For each hyperedge, make list of vertices. */
    prev_edge = -1;
    for (k=0; k<(*nInput); k++){
      /* next edge is mat[k].i */
      if (mat[k].i > prev_edge) 
        (*index)[mat[k].i] = k;
      prev_edge = mat[k].i;
      (*vertex)[k] = mat[k].j;
    }
    /* Complete index array; overwrite -1s. */
    tmp = *nInput;
    for (j=(*nEdge); j; j--){
      if ((*index)[j] == -1)
        (*index)[j] = tmp;
      else
        tmp = (*index)[j];
    }

End:
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
       *nVtx  = *nEdge  = *nInput = *vwgt_dim = *ewgt_dim = 0;
       Zoltan_Multifree(__FILE__, __LINE__, 4, index, vertex, ewgt, vwgt);
    }
    return err;
}


static int comp(const void *a, const void *b)
{
  return ((struct ijv *)a)->i - ((struct ijv *)b)->i;
}



/*****************************************************************************/
/* Uses character reads to read arbitrary long hyperedge definition lines.   */
static int nextstr (FILE *f, char *string)
{
char ch;
  
  *string = 0;
  while (1) {
    ch = fgetc(f);    
    if (ch == EOF)          return 0;
    if (ch == '\n')         return -1;       
    if (isspace(ch) == 0)   break;
  }
    
  while (1)  {
    if (ch == EOF || isspace(ch)) {
      *string = 0;
      ungetc (ch, f);  /* ??? */
      return 1;
    }
    *string++ = ch; 
    ch = fgetc(f);
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
