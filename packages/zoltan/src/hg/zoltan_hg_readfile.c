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


#include "hypergraph.h"

#define BUF_LEN 1000000


static int old_readfile  (int, FILE*, int*, int*, int*, int**, int**, int*, 
                          float**, int*, float**, int *) ;
static int patoh_readfile (int, FILE*, int*, int*, int*, int**, int**, int*, 
                           float**, int*, float**, int *) ;


/*****************************************************************************/

int Zoltan_HG_Readfile (
 int Proc,
 FILE *f,
 int *nVtx, int *nEdge, int *nPin,
 int **hindex,   int **hvertex,
 int *vwgt_dim, float **vwgt,
 int *ewgt_dim, float **ewgt,
 int *base)
    {
    int err = ZOLTAN_OK ;
    char string[81], *s;
    char *yo = "Zoltan_HG_Readfile" ;

    /* Initialize return values in case of error. */
    *nVtx   = *nEdge   = *nPin = *vwgt_dim = *ewgt_dim = 0;
    *hindex = *hvertex = NULL;
    *vwgt   = *ewgt    = NULL;
    *base   = 0;

    do {
       if (!fgets (string, 80, f))
          {
          ZOLTAN_PRINT_ERROR(Proc, yo, "ERROR...not able to read input file\n");
          err = ZOLTAN_FATAL;
          goto End;
          }
       s = strtok(string, " \t\n") ;
    } while (*s == '%');  /* Skip leading comment lines. */

    if (atoi(s) < 2) /* Note -- this logic is not correct for files 
                        with only one vertex. */
       err = patoh_readfile (Proc,f,nVtx,nEdge,nPin,hindex,hvertex,
                             vwgt_dim,vwgt,ewgt_dim,ewgt,base) ;
    else if (atoi(s) > 1) 
       err = old_readfile   (Proc,f,nVtx,nEdge,nPin,hindex,hvertex,
                             vwgt_dim,vwgt,ewgt_dim,ewgt,base) ;

End:
    return  err ;
    }

/*****************************************************************************/


static int old_readfile (int Proc,
 FILE *f,
 int *nVtx, int *nEdge, int *nPin,
 int **index,   int **vertex,
 int *vwgt_dim, float **vwgt,
 int *ewgt_dim, float **ewgt,
 int *base)
    {
    int count, ierr = ZOLTAN_OK;
    char errstr[200] ;
    int Hedge=0, code=0, pin, i;
    char string[BUF_LEN], *s;
    char *yo = "old_readfile" ;

    /* TODO: edge weights, multiple edge/vertex weights */

    /* IBM-format files are assumed to be 1-based. */
    *base = 1;

    rewind(f) ;

    do {
       if (!fgets (string, 80, f))
          {
          ZOLTAN_PRINT_ERROR(Proc, yo, "ERROR...not able to read input file\n");
          ierr = ZOLTAN_FATAL;
          goto End;
          }
    } while (string[0] == '%');  /* Skip leading comment lines. */

    count = sscanf (string, "%d %d %d %d", nVtx, nEdge, nPin, &code) ;
    if (count <  3)
       {
       ZOLTAN_PRINT_ERROR (Proc, yo, 
                   "ERROR, first line of file must be: |V| |E| |P| (code)\n") ;
       ierr = ZOLTAN_FATAL;
       goto End;
       }

    /* nEdge HYPEREDGE LINES */
    /* KDD -- This logic is wrong if no pins are specified. */
    if (!((*index)  = (int *) ZOLTAN_MALLOC (sizeof (int) * (*nEdge+1))) ||
        !((*vertex) = (int *) ZOLTAN_MALLOC (sizeof (int) *  *nPin)))
           {
           ZOLTAN_PRINT_ERROR(Proc, yo, "Insufficient memory.");
           ierr = ZOLTAN_MEMERR;
           goto End;
           }

    Hedge = 0 ;
    for (i = 0 ; i < *nEdge; i++)
       {
       (*index)[i] = Hedge;
       if (!(fgets (string, BUF_LEN, f)))
          {
          sprintf(errstr, "ERROR... read hvertex %d\n",i);
          ZOLTAN_PRINT_ERROR (Proc, yo, errstr) ;
          ierr = ZOLTAN_FATAL;
          goto End;
          }
       s = strtok(string," \n");
       while (s)
          {
          pin = atoi(s);
          if (pin <= 0 || pin > *nVtx)
              {
              sprintf(errstr, 
                      "ERROR...pin %d of vertex %d is out of range [%d,%d]!\n",
                      pin,i, 1, *nVtx);
              ZOLTAN_PRINT_ERROR (Proc, yo, errstr) ;
              ierr = ZOLTAN_FATAL;
              goto End;
              }
          (*vertex)[Hedge++] = pin;
          s = strtok(NULL," \n");
          }
       }
    (*index)[*nEdge] = Hedge;

    /* nVtx vertex weights */
    if ((code/10)%10 == 1)
       {
       if (!((*vwgt) = (float *) ZOLTAN_MALLOC (sizeof (float) * *nVtx)))
          {
          ZOLTAN_PRINT_ERROR (Proc, yo, "Insufficient memory for vwgt.");
          ierr = ZOLTAN_MEMERR;
          goto End;
          }
       for (i = 0 ; i < *nVtx; i++)
          {
          if (!(fgets (string, BUF_LEN, f)))
             {
             sprintf(errstr, "ERROR... reading weight %d\n",i);
             ZOLTAN_PRINT_ERROR (Proc, yo, errstr) ;
             ierr = ZOLTAN_FATAL;
             goto End;
             }
          s = strtok (string, " \n");
          (*vwgt)[i] = (float)(atoi(s));
          }
       }

    if (fscanf(f, "%d", &i) != EOF)
       {
       ZOLTAN_PRINT_WARN (Proc, yo, "Input file is longer than expected!\n") ;
       ierr = ZOLTAN_WARN;
       }

End:
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) 
       {
       *nVtx  = *nEdge  = *nPin = *vwgt_dim = *ewgt_dim = 0;
       Zoltan_Multifree(__FILE__, __LINE__, 4, index, vertex, ewgt, vwgt);
       }
    return ierr;
    }


/*****************************************************************************/

static int patoh_readfile (int Proc,
 FILE *f,
 int *nVtx, int *nEdge, int *nPin,
 int **index,   int **vertex,
 int *vwgt_dim, float **vwgt,
 int *ewgt_dim, float **ewgt,
 int *base)
    {
    int i, j ;
    int count, ierr = ZOLTAN_OK;
    char errstr[200];
    int Hedge=0, code=0, dims=1, pin;
    /* KDD -- should dims be initialized to zero (for no vwgts) instead of 1? */
    char string[BUF_LEN], *s;
    char *yo = "patoh_readfile" ;


    /* TODO: edge weights, multiple edge/vertex weights */

    /* Read PaToH file. */
    rewind(f) ;
    do {
       if (!fgets (string, 80, f))
          {
          ZOLTAN_PRINT_ERROR(Proc, yo, "ERROR...not able to read input file\n");
          ierr = ZOLTAN_FATAL;
          goto End;
          }
    } while (string[0] == '%');  /* Skip leading comment lines. */

    count = sscanf (string, "%d %d %d %d %d %d", 
                    base, nVtx, nEdge, nPin, &code, &dims) ;
    if (count <  4)  /* code and dims are optional */ 
       {
       ZOLTAN_PRINT_ERROR (Proc, yo, 
         "Control line of file must be: base |V| |E| |P| (code) (dims)\n");
       ierr = ZOLTAN_FATAL ;
       goto End;
       }

    /* nEdge HYPEREDGE LINES */
    /* KDD -- This logic is wrong if no pins are specified. */
    if (!((*index)  = (int *) ZOLTAN_MALLOC (sizeof (int) * (*nEdge+1))) ||
        !((*vertex) = (int *) ZOLTAN_MALLOC (sizeof (int) *  *nPin)))
           {
           ZOLTAN_PRINT_ERROR(Proc, yo, "Insufficient memory.");
           ierr = ZOLTAN_MEMERR;
           goto End;
           }
    if (code==2 || code==3)
       {
       /* KDD -- This logic is wrong if no edges are specified. */
       *ewgt = (float *) ZOLTAN_MALLOC (sizeof (float) * *nEdge) ;
       if (*ewgt == NULL)
          {
          ZOLTAN_PRINT_ERROR(Proc, yo, "Insufficient memory.");
          ierr = ZOLTAN_MEMERR;
          goto End;
          }
       }

    Hedge = 0 ;
    j = 0 ;
    for (i = 0 ; i < *nEdge; i++)
       {
       (*index)[i] = Hedge;

       if (!(fgets (string, BUF_LEN, f)))
          {
          sprintf(errstr, "ERROR... read hvertex %d\n",i);
          ZOLTAN_PRINT_ERROR (Proc, yo, errstr) ;
          ierr = ZOLTAN_FATAL;
          goto End;
          }

       s = strtok(string," \n");
       if (s == NULL)
          continue ;
       if (*s == '%')
          {
          i-- ;       /* don't count comment line in data file */
          continue ;
          }
       if (code == 2 || code == 3)
          {
          (*ewgt)[j++] = (float) atoi(s) ;
          s = strtok (NULL, " \t\n") ;
          }
       while (s != NULL)
          {
          pin = atoi(s) ;
          if (pin < 0+*base || pin >= *nVtx+*base)
              {
              sprintf(errstr, 
                      "ERROR...pin %d of edge %d is out of range [%d,%d]!\n",
                      pin,i,0+*base, *nVtx-1+*base);
              ZOLTAN_PRINT_ERROR (Proc, yo, errstr) ;
              ierr = ZOLTAN_FATAL;
              goto End;
              }
          (*vertex)[Hedge++] = pin ;
          s = strtok(NULL," \n\t");
          }
       }
    (*index)[*nEdge] = Hedge;

    if (code == 0  || code == 2)
       goto End;

    /* nVtx vertex weights */
    /* KDD -- shouldn't this code use dims in some way? */
    if (!((*vwgt) = (float *) ZOLTAN_MALLOC (sizeof (float) * *nVtx)))
        {
        ZOLTAN_PRINT_ERROR (Proc, yo, "Insufficient memory for vwgt.");
        ierr = ZOLTAN_MEMERR;
        goto End;
        }
    i = 0 ;
    while (fgets (string, BUF_LEN, f) != NULL && i < *nVtx)
        {
        s = strtok (string, " \n") ;
        if (*s == '%')
            continue ;
        while (s != NULL && i < *nVtx)
           {
           (*vwgt)[i++] = (float)(atoi(s));
           s = strtok (NULL, " \t\n") ;
           }
        }
End:
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) 
       {
       *nVtx  = *nEdge  = *nPin = *vwgt_dim = *ewgt_dim = 0;
       Zoltan_Multifree(__FILE__, __LINE__, 4, index, vertex, ewgt, vwgt);
       }
    return ierr ;

    }
/*****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
