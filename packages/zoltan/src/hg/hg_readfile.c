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

int Zoltan_HG_Readfile (
  ZZ *zz, 
  HGraph *hg,
  char *hgraphfile
)
{ 
  FILE  *f;
  int   i, Hedge=0, pin, code=0;
  char  *s, string[BUF_LEN];
  char errstr[200] ;
  char *yo = "hgraph_load" ;

/* TODO: edge weights, multiple edge/vertex weights */

  hg->info = 0;
  hg->vindex = hg->vedge = NULL;
  hg->vwgt = hg->ewgt = NULL;
  hg->coor = NULL;

  if (!(f = fopen( hgraphfile, "r" )))
  { sprintf(errstr, "ERROR...not able to read file %s!\n",hgraphfile);
    ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
    return ZOLTAN_WARN;
  }
  if (!(fgets (string, BUF_LEN, f)))
  { sprintf(errstr, "ERROR... read line 1 of file %s!\n",hgraphfile);
    ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
    return ZOLTAN_WARN;
  }

/* HEAD LINE */
  hg->nVtx = atoi (strtok(string," \n\t"));
  if ((s=strtok(0," \n\t")))
    hg->nEdge = atoi(s);
  else
  { sprintf(errstr, "ERROR...missing number of hyperedges in head line");
    ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
    return ZOLTAN_WARN;
  }
  if ((s=strtok(0," \n\t")))
    hg->nPin = atoi(s);
  else
  { sprintf(errstr,"ERROR...missing number of pins in head line");
    ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
    return ZOLTAN_WARN;
  }
  if ((s=strtok(0," \n\t")))
    code = atoi (s);

/* nEdge HYPEREDGE LINES */
  hg->hindex  = (int *) ZOLTAN_MALLOC (sizeof (int) * (hg->nEdge+1));
  hg->hvertex = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nPin);
  if (hg->hindex == NULL || hg->hvertex == NULL)
     {
     ZOLTAN_FREE ((void **) &hg->hindex) ;
     ZOLTAN_FREE ((void **) &hg->hvertex) ;
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i=0; i<hg->nEdge; i++)
  { hg->hindex[i] = Hedge;
    if (!(fgets (string, BUF_LEN, f)))
    { sprintf(errstr, "ERROR... read hvertex %d from %s!\n",i,hgraphfile);
      ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
      return ZOLTAN_WARN;
    }
    s = strtok(string," \n");

    while (s)
    { pin = atoi(s);
      s = strtok(NULL," \n");
      if (pin<=0 || pin>hg->nVtx)
      { sprintf(errstr, "ERROR...pin %d of vertex %d is out of range [%d,%d]!\n",
         pin,i+1,1,hg->nVtx);
        ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
        return ZOLTAN_WARN;
      }
      hg->hvertex[Hedge++] = pin-1;
  } }
  hg->hindex[hg->nEdge] = Hedge;

/* nVtx vertex weights */
  if ((code/10)%10 == 1)
  { hg->vwgt = (float *) ZOLTAN_MALLOC (sizeof (float) * hg->nVtx);
    if (hg->vwgt == NULL)
       {
       ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
       return ZOLTAN_MEMERR;
       }
    for (i=0; i<hg->nVtx; i++)
    { if (!(fgets (string, BUF_LEN, f)))
        { sprintf(errstr, "ERROR... read weight %d from %s!\n",i,hgraphfile);
          ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
          return ZOLTAN_WARN;
        }
      s = strtok(string," \n");
      hg->vwgt[i] = (float)(atoi(s));
  } }

  if (fscanf(f, "%d", &i) != EOF)
  { sprintf(errstr, "ERROR... file %s is too long!\n", hgraphfile);
    ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
    return ZOLTAN_WARN;
  }

  if (Zoltan_HG_Create_Mirror (zz, hg))
    return ZOLTAN_WARN;
  if (fclose(f))
    return ZOLTAN_WARN;
  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
