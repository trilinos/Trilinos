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
#ifdef HGEXEC
#include "hypergraph.h"
#else
#include "zoltan.h"
#endif

#include "dr_hg_readfile.h"
#include "dr_util_const.h"
#include "dr_mmio.h"

#define BUF_LEN 1000000

/* struct for indices i & j and a value */
struct ijv
{
  int i, j;
  double val;
};

/* comparison functions. */
static int comp(const void *, const void *);


static int old_readfile  (int, FILE*, int*, int*, int*, int**, int**, int*,
 float**, int*, float**, int*);
static int patoh_readfile (int, FILE*, int*, int*, int*, int**, int**, int*,
 float**, int*, float**, int*);
static int MM_readfile (int, FILE*, int*, int*, int*, int**, int**, int*,
 float**, int*, float**, int*);

static int nextstr (FILE *f, char *string);

/*****************************************************************************/

int Zoltan_HG_Readfile (
  int Proc,
  FILE *f,
  int *nVtx, 
  int *nEdge, 
  int *nInput,
  int **hindex, 
  int **hvertex,
  int *vwgt_dim, 
  float **vwgt,
  int *ewgt_dim, 
  float **ewgt,
  int *base
)
{
char string[81], *s;
int err = ZOLTAN_OK;
char *yo = "Zoltan_HG_Readfile";
int MM = 0;

  /* Initialize return values in case of error. */
  *nVtx   = *nEdge   = *nInput = *vwgt_dim = *ewgt_dim = 0;
  *hindex = *hvertex = NULL;
  *vwgt   = *ewgt    = NULL;
  *base   = 0;

  do {
    if (!fgets (string, 80, f)) {
      fprintf(stderr, "%s ERROR...not able to read input file\n", yo);
      err = ZOLTAN_FATAL;
      goto End;
    }
    s = strtok(string, " \t\n");
    if (!MM) 
      MM = !strncmp(s, "%%Matrix", 8); /* MM !=0 if match */

  } while (*s == '%');         /* Skip leading comment lines. */

  if (atoi(s) < 2) /* Guess Patoh -- Not correct for files 
                                     with only one vertex. */
    err = patoh_readfile (Proc, f, nVtx, nEdge, nInput, hindex, hvertex,
                          vwgt_dim, vwgt, ewgt_dim, ewgt,base);
  else if (atoi(s) > 1) {
    if (MM) /* MatrixMarket format */
      err = MM_readfile (Proc, f, nVtx, nEdge, nInput, hindex, hvertex,
                         vwgt_dim, vwgt, ewgt_dim, ewgt, base);
    else /* Assume our HG format */
      err = old_readfile (Proc, f, nVtx, nEdge, nInput, hindex, hvertex,
                          vwgt_dim, vwgt, ewgt_dim, ewgt, base);
  }

End:
  return  err;
}

/*****************************************************************************/



static int old_readfile (
 int Proc,
 FILE *f,
 int *nVtx, int *nEdge, int *nInput,
 int **index,   int **vertex,
 int *vwgt_dim, float **vwgt,
 int *ewgt_dim, float **ewgt,
 int *base)
{
int count, check, err = ZOLTAN_OK;
char errstr[200];
int Hedge=0, code=0, pin, i;
char string[BUF_LEN];
int vdim=0, edim=0;
char *yo = "old_readfile";
          
    *base = 1;   /* IBM-format files are assumed to be 1-based. */
    rewind(f);       /* need to read first line again! */
    do {
       if (!fgets (string, 80, f)) {
          fprintf(stderr, "%s ERROR...not able to read input file\n", yo);
          err = ZOLTAN_FATAL;
          goto End;
       }
    } while (string[0] == '%');  /* Skip leading comment lines. */

    /* read input command line, note final variables are optional */
    if (sscanf (string, "%d %d %d %d %d %d", nVtx, nEdge, nInput, &code, &vdim,
                &edim) < 4) {
       fprintf(stderr, "%s ERROR first line of file must be: |V| |E| |I| (code)"
        "\n", yo);
       err = ZOLTAN_FATAL;
       goto End;
    }        /* consider removing the restriction requiring nInput later */

    /* allocate storage for hyperedge weights */
    if (code == 1 || code == 11)  {
       *ewgt_dim = (edim > 0) ? edim : 1;
       *ewgt = (float*) malloc (*nEdge * (*ewgt_dim) * sizeof(float));     
       if (*ewgt == NULL) {
          fprintf(stderr, "%s Insufficient memory for ewgt.", yo);
          err = ZOLTAN_MEMERR;
          goto End;
       }    
    }     
          
    /* allocate storage for hypergraph data */    
    if (!(*index  = (int*) malloc ((*nEdge+1) * sizeof(int)))
     || !(*vertex = (int*) malloc  (*nInput   * sizeof(int)))) {
         fprintf(stderr, "%s Insufficient memory.", yo);
         err = ZOLTAN_MEMERR;
         goto End;
    }

    /* read list of vertices in each hyperedge */
    check = 0;
    Hedge = 0;
    for (i = 0; i < *nEdge; i++) {
       (*index)[i] = Hedge;

       /* read optional hyperedge weights first */
       if (code == 1 || code == 11) {
          for (count = 0; count < *ewgt_dim; count++)
             if (nextstr (f, string) != 1)  {
                sprintf(errstr, "%s ERROR...reading hyperedge weight %d\n", 
                        yo, i);
                fprintf(stderr,  errstr);
                err = ZOLTAN_FATAL;
                goto End;
             }             
          (*ewgt)[*ewgt_dim * i + count] = (float) atoi(string);          
       }
       
       /* now read vertices in current hyperedge */
       while (nextstr(f,string) == 1) {
          ++check;
          pin = atoi(string);
          if (pin <= 0 || pin > *nVtx) {
              sprintf(errstr,  "%s ERROR pin %d of hyperedge %d is out of range"
               " [%d,%d]!\n", yo, pin, i, 1, *nVtx);
              fprintf(stderr, errstr);
              err = ZOLTAN_FATAL;
              goto End;
          }
          (*vertex)[Hedge++] = pin;
       }
    }
    (*index)[*nEdge] = Hedge;
    if (check != *nInput) {
       sprintf(errstr, "%s ERROR...PIN sanity check failed\n", yo);
       fprintf(stderr,  errstr);
       err = ZOLTAN_FATAL;
       goto End;
    }                 

    /* conditionally read nVtx vertex weight vectors */
    if (code == 10 || code == 11) {
       *vwgt_dim = (vdim > 0) ? vdim : 1;
       *vwgt = (float*) malloc (*nVtx * (*vwgt_dim) * sizeof(float)); 
       if (*vwgt == NULL) {
          fprintf(stderr, "%s Insufficient memory for vwgt.", yo);
          err = ZOLTAN_MEMERR;
          goto End;
       }
       for (i = 0; i < *nVtx; i++) {
          for (count = 0; count < *vwgt_dim; count++) {
             if (nextstr(f, string) != 1) {
                sprintf(errstr, "%s ERROR...reading vertex weight %d\n", yo, i);
                fprintf(stderr,  errstr);
                err = ZOLTAN_FATAL;
                goto End;
             }             
             (*vwgt)[*vwgt_dim * i + count] = (float) atoi(string);
          }
          if (nextstr(f,string) != -1)  {
             sprintf(errstr, "%s ERROR...reading vertex weight %d\n",yo, i);
             fprintf(stderr,  errstr);
             err = ZOLTAN_FATAL;
             goto End;
          }                          
       }
    }

    if (nextstr (f, string) == 1) {
       fprintf(stderr, "%s Input file is longer than expected!\n",yo);
       err = ZOLTAN_WARN;
    }

End:
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
       *nVtx  = *nEdge  = *nInput = *vwgt_dim = *ewgt_dim = 0;
       Zoltan_Multifree(__FILE__, __LINE__, 4, index, vertex, ewgt, vwgt);
    }
    return err;
}


/*****************************************************************************/



/* Read PaToH formatted hypergraph file. */
static int patoh_readfile (
int Proc,
FILE *f,
int *nVtx, int *nEdge, int *nInput,
int **index,   int **vertex,
int *vwgt_dim, float **vwgt,
int *ewgt_dim, float **ewgt,
int *base)
{
int i, j;
int count, err = ZOLTAN_OK;
char errstr[200];
int Hedge=0, code=0, dims=1, pin;
    /* KDD -- should dims be initialized to zero (for no vwgts) instead of 1? */
char string[BUF_LEN], *s;
char *yo = "patoh_readfile";

  /* TODO: edge weights, multiple edge/vertex weights */

  rewind(f);
  do {
    if (!fgets (string, 80, f)) {
      fprintf(stderr, "%s ERROR...not able to read input file\n", yo);
      err = ZOLTAN_FATAL;
      goto End;
    }
  } while (string[0] == '%');  /* Skip leading comment lines. */

  count = sscanf (string, "%d %d %d %d %d %d",  base, nVtx, nEdge, nInput,
                  &code, &dims);
  if (count <  4) {  /* code and dims are optional */
    fprintf(stderr, 
    "%s Control line of file must be: base |V| |E| |I| (code) (dims)\n", yo);
    err = ZOLTAN_FATAL;
    goto End;
  }

  /* nEdge HYPEREDGE LINES */
  /* KDD -- This logic is wrong if no pins are specified. */
  if (!(*index  = (int*) malloc ((*nEdge+1) * sizeof(int)))
   || !(*vertex = (int*) malloc  (*nInput   * sizeof(int)))) {
    fprintf(stderr, "%s Insufficient memory.", yo);
    err = ZOLTAN_MEMERR;
    goto End;
  }
  if (code == 2 || code == 3) {
    /* KDD -- This logic is wrong if no edges are specified. */
    *ewgt = (float*) malloc (*nEdge * sizeof(float));
    if (*ewgt == NULL) {
      fprintf(stderr, "%s Insufficient memory.", yo);
       err = ZOLTAN_MEMERR;
       goto End;
    }
  }
             
  /* nEdge HYPEREDGE LINES */
  /* KDD -- This logic is wrong if no pins are specified. */
  if (!(*index  = (int*) malloc ((*nEdge+1) * sizeof(int)))
   || !(*vertex = (int*) malloc  (*nInput   * sizeof(int)))) {
    fprintf(stderr, "%s Insufficient memory.", yo);
    err = ZOLTAN_MEMERR;
    goto End;
  }
  if (code == 10 || code == 11) {
    /* KDD -- This logic is wrong if no edges are specified. */
    *ewgt = (float*) malloc (*nEdge * sizeof(float));
    if (*ewgt == NULL) {
      fprintf(stderr, "%s Insufficient memory.", yo);
      err = ZOLTAN_MEMERR;
      goto End;
    }
  }

  Hedge = 0;
  j = 0;
  for (i = 0; i < *nEdge; i++) {
    (*index)[i] = Hedge;

    if (!(fgets (string, BUF_LEN, f))) {
      sprintf(errstr, "%s ERROR... read hvertex %d\n",yo, i);
      fprintf(stderr, errstr);
      err = ZOLTAN_FATAL;
      goto End;
    }

    s = strtok(string," \n");
    if (s == NULL)
      continue;
    if (*s == '%') {
      i--;       /* don't count comment line in data file */
      continue;
    }
    if (code == 10 || code == 11) {
      (*ewgt)[j++] = (float) atoi(s);
      s = strtok (NULL, " \t\n");
    }
    while (s != NULL) {
      pin = atoi(s);
      if (pin < 0 + *base || pin >= *nVtx + *base) {
        sprintf(errstr, "%s ERROR...pin %d of edge %d is out of range "
                "[%d,%d]!\n",yo, pin,i,0+*base, *nVtx-1+*base);
        fprintf(stderr, errstr);
        err = ZOLTAN_FATAL;
        goto End;
      }
      (*vertex)[Hedge++] = pin;
      s = strtok(NULL," \n\t");
    }
  }
  (*index)[*nEdge] = Hedge;

  if (code == 0  || code == 10)
    goto End;

  /* nVtx vertex weights */
  /* KDD -- shouldn't this code use dims in some way? */
  *vwgt_dim = 1;
  if (!((*vwgt) = (float *) malloc (*nVtx * sizeof(float)))) {
    fprintf(stderr, "%s Insufficient memory for vwgt.", yo);
    err = ZOLTAN_MEMERR;
    goto End;
  }
  i = 0;
  while (fgets (string, BUF_LEN, f) != NULL && i < *nVtx) {
    s = strtok (string, " \n");
    if (*s == '%')
      continue;
    while (s != NULL && i < *nVtx) {
      (*vwgt)[i++] = (float) atoi(s);
      s = strtok (NULL, " \t\n");
    }
  }
End:
  if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
    *nVtx  = *nEdge  = *nInput = *vwgt_dim = *ewgt_dim = 0;
    safe_free((void **) index);
    safe_free((void **) vertex);
    safe_free((void **) ewgt);
    safe_free((void **) vwgt);
  }
  return err;
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
char errstr[200];
int prev_edge;
int rowhedges=1; /* default is row hyperedge model */
char *yo = "MM_readfile";

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


    int ret_code;
    MM_typecode matcode;
    int M, N, nz;   
    int i, j, k, tmp;
    struct ijv *mat;

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

    if (f !=stdin) fclose(f);

    /************************/
    /* sort indices         */
    /************************/
    qsort(mat, nz, sizeof(struct ijv), comp);

    /************************/
    /* Debug: write out matrix */
    /************************/

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    for (k=0; k<nz; k++)
        fprintf(stdout, "%d %d %20.19g\n", mat[k].i+1, mat[k].j+1, mat[k].val);

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


/*
#include <ctype.h>
*/

static int nextstr (FILE *f, char *string)
{
int i = 0;
char ch;
  
  string[0] = 0;
  while (1) {
    ch = fgetc(f);    
    if (ch == EOF)
      return 0;
    if (ch == '\n')
      return -1;       
    if (isspace(ch) == 0)
      break;
  }
    
  while (1)  {
    if (ch == EOF || isspace(ch)) {
      string[i] = 0;
      ungetc (ch, f);
      return 1;
    }
    string[i++] = ch; 
    ch = fgetc(f);
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
