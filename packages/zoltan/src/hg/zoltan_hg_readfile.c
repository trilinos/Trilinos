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

int Zoltan_HG_Readfile   (ZZ*, FILE*, int*, int*, int*, int**, int**, int*, float**, int*, float**) ;
static int ibm_readfile  (ZZ*, FILE*, int*, int*, int*, int**, int**, int*, float**, int*, float**) ;
static int umit_readfile (ZZ*, FILE*, int*, int*, int*, int**, int**, int*, float**, int*, float**) ;


int Zoltan_HG_Readfile (ZZ *zz,
 FILE *f,
 int *nVtx, int *nEdge, int *nPin,
 int **index,   int **vertex,
 int *vwgt_dim, float **vwgt,
 int *ewgt_dim, float **ewgt)
    {
    int count ;
    char errstr[200] ;
    int Hedge=0, code=0, pin, i;
    char string[BUF_LEN], *s;
    char *yo = "ibm_readfile" ;

    /* TODO: edge weights, multiple edge/vertex weights */

    rewind(f) ;
    if (!fgets (string, 80, f))
       {
       sprintf(errstr, "ERROR...not able to read input file\n");
       ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
       return ZOLTAN_WARN;
       }

    count = sscanf (string, "%d %d %d %d", nVtx, nEdge, nPin, &code) ;
    if (count <  4)
       {
       sprintf (errstr, "ERROR, first line of file must be: |V| |E| |P| (code)\n");
       ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
       }

    /* nEdge HYPEREDGE LINES */
    if (!((*index)  = (int *) ZOLTAN_MALLOC (sizeof (int) * (*nEdge+1))) ||
        !((*vertex) = (int *) ZOLTAN_MALLOC (sizeof (int) *  *nPin)))
           {
           ZOLTAN_FREE ((void **) index) ;
           ZOLTAN_FREE ((void **) vertex) ;
           ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
           return ZOLTAN_MEMERR;
           }

    Hedge = 0 ;
    for (i = 0 ; i < *nEdge; i++)
       {
       (*index)[i] = Hedge;
       if (!(fgets (string, BUF_LEN, f)))
          {
          sprintf(errstr, "ERROR... read hvertex %d\n",i);
          ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
          return ZOLTAN_WARN;
          }
       s = strtok(string," \n");
       while (s)
          {
          pin = atoi(s);
          s = strtok(NULL," \n");
          if (pin <= 0 || pin > *nVtx)
              {
              sprintf(errstr, "ERROR...pin %d of vertex %d is out of range [%d,%d]!\n",
               pin,i+1,1, *nVtx);
              ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
              return ZOLTAN_WARN;
              }
          (*vertex)[Hedge++] = pin-1;
          }
       }
    (*index)[*nEdge] = Hedge;

    /* nVtx vertex weights */
    if ((code/10)%10 == 1)
       {
       if (!((*vwgt) = (float *) ZOLTAN_MALLOC (sizeof (float) * *nVtx)))
          {
          ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory for vwgt.");
          return ZOLTAN_MEMERR;
          }
       for (i = 0 ; i < *nVtx; i++)
          {
          if (!(fgets (string, BUF_LEN, f)))
             {
             sprintf(errstr, "ERROR... reading weight %d\n",i);
             ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
             return ZOLTAN_WARN;
             }
          s = strtok (string, " \n");
          (*vwgt)[i] = (float)(atoi(s));
          }
       }

    if (fscanf(f, "%d", &i) != EOF)
       {
       sprintf(errstr, "ERROR... file is too long!\n");
       ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
       return ZOLTAN_WARN;
       }
printf ("friday returning from ibm read routine\n") ;

    return ZOLTAN_OK;
    }



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
