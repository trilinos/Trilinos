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

#include "zz_const.h"
#include "zz_util_const.h"
#include "hg_hypergraph.h"
    
#define MEMORY_ERROR { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
  ierr = ZOLTAN_MEMERR; \
  goto End; \
}


int Zoltan_HG_Hypergraph_Callbacks(
  ZZ *zz,
  ZHG *zhg,                /* Input:   Pointer to Zoltan's structure with
                                       GIDs, LIDs.
                              Output:  Removed edge fields of zhg are changed
                                       if dense edges are removed. */
  int gnVtx,               /* Input:   Global number of vertices in hgraph */
  float esize_threshold,   /* Input:   %age of gnVtx considered a dense edge */
  int return_removed,      /* Input:   flag indicating whether to return
                                       removed edge info */
  int *nedges,             /* Output:  Number of hyperedges on this processor */
  ZOLTAN_ID_PTR *egids,    /* Output:  GIDs of hyperedges on this processor */
  ZOLTAN_ID_PTR *elids,    /* Output:  LIDs of hyperedges on this processor */
  int **esizes,            /* Output:  # of vertices for each hyperedge on
                                       this processor */
  float **ewgts,           /* Output:  edge weights for each hyperedge on
                                       this processor */
  int *npins,              /* Output:  # of pins on this processor = 
                                       sum esizes[i], i = 0..nedges-1. */
  ZOLTAN_ID_PTR *pins,     /* Output:  vertex GIDs of pins */
  int **pin_procs          /* Output:  processors owning pin vertices */
)
{
/* Function to call the Zoltan Hypergraph callback functions */
/* Can be called from serial or parallel hypergraph partitioners */
static char *yo = "Zoltan_HG_Hypergraph_Callbacks";
int ierr = ZOLTAN_OK;
int i, j;
int ewgtdim = zz->Edge_Weight_Dim;
int nwgt;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;
float gesize_threshold;   /* Edges with more vertices than gesize_threshold
                             are considered to be dense. */
int nkeep;                /* Number of edges below gesize_threshold. */
int nremove;              /* Number of edges to be removed; i.e., number of
                             edges above gesize_threshold. */
ZOLTAN_ID_PTR remove_egids = NULL;  /* Edge GIDs for removed edges */
ZOLTAN_ID_PTR remove_elids = NULL;  /* Edge LIDs for removed edges */
int *remove_esizes = NULL;          /* Edge sizes (# pins) for removed edges */
float *remove_ewgts = NULL;         /* Edge weights for removed edges */

  *nedges = zz->Get_Num_HG_Edges(zz->Get_Num_HG_Edges_Data, &ierr);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_Num_HG_Edges");
    goto End;
  }
  
  if (*nedges > 0) {

    /* Get info about the edges:  GIDs, LIDs, sizes, edge weights */
    *egids = ZOLTAN_MALLOC_GID_ARRAY(zz, *nedges);
    *elids = ZOLTAN_MALLOC_LID_ARRAY(zz, *nedges);
    *esizes = (int *) ZOLTAN_MALLOC(*nedges * sizeof(int));
    nwgt = *nedges * ewgtdim;
    if (nwgt) 
      *ewgts = (float *) ZOLTAN_MALLOC(nwgt * sizeof(float));
    if (!*esizes || !*egids || (num_lid_entries && !*elids) ||
        (nwgt && !*ewgts)) MEMORY_ERROR;

    ierr = zz->Get_HG_Edge_Info(zz->Get_HG_Edge_Info_Data,
                                     num_gid_entries, num_lid_entries,
                                     *nedges, ewgtdim, 
                                     *egids, *elids, *esizes, *ewgts); 
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_HG_Edge_Info");
      goto End;
    }
                     
    /* Remove dense edges from input list */
    gesize_threshold = esize_threshold * gnVtx;
    nremove = 0;
    for (i = 0; i < *nedges; i++) 
      if ((*esizes)[i] > gesize_threshold)  {
        nremove++;
      }

    if (nremove) {
      if (return_removed) {
        /* Keep a record of removed edges so we can get their edge lists
         * later if needed (e.g., to evaluate total partition quality) */
        zhg->nRemove = nremove;
        zhg->Remove_EGIDs = remove_egids = ZOLTAN_MALLOC_GID_ARRAY(zz, nremove);
        zhg->Remove_ELIDs = remove_elids = ZOLTAN_MALLOC_LID_ARRAY(zz, nremove);
        zhg->Remove_Esize = remove_esizes 
                          = (int *) ZOLTAN_MALLOC(nremove * sizeof(int));
        if (ewgtdim)
          zhg->Remove_Ewgt = remove_ewgts 
                           = (float *) ZOLTAN_MALLOC(nremove * ewgtdim
                                                              * sizeof(float));
  
        if (!remove_egids || (num_lid_entries && !remove_elids) 
            || !remove_esizes || (ewgtdim && !remove_ewgts)) MEMORY_ERROR;
      }
      
      nremove = nkeep = 0;
      for (i = 0; i < *nedges; i++)
        if ((*esizes)[i] <= gesize_threshold) {
          /* Keep the edge in egids/elids to obtain its pins. */
          if (nkeep != i) {
            ZOLTAN_SET_GID(zz, &((*egids)[nkeep*num_gid_entries]),
                               &((*egids)[i*num_gid_entries]));
            if (num_lid_entries)
              ZOLTAN_SET_LID(zz, &((*elids)[nkeep*num_lid_entries]),
                                 &((*elids)[i*num_lid_entries]));
            (*esizes)[nkeep] = (*esizes)[i];
            for (j = 0; j < ewgtdim; j++)
              (*ewgts)[nkeep * ewgtdim + j] = (*ewgts)[i*ewgtdim + j];
          }
          nkeep++;
        }
        else if (return_removed) {
          /* Remove the edges from egids/elids; don't want to have to 
             allocate memory for its pins */
          ZOLTAN_SET_GID(zz, &(remove_egids[nremove*num_gid_entries]),
                             &((*egids)[i*num_gid_entries]));
          if (num_lid_entries)
            ZOLTAN_SET_LID(zz, &(remove_elids[nremove*num_lid_entries]),
                               &((*elids)[i*num_lid_entries]));
          remove_esizes[nremove] = (*esizes)[i];
          for (j = 0; j < ewgtdim; j++)
            remove_ewgts[nremove * ewgtdim + j] = (*ewgts)[i*ewgtdim + j];
          (nremove)++;
        }
      *nedges = nkeep;
    }
    
    /* Now get lists of vertices that are in hyperedges */

    *npins = 0;
    for (i = 0; i < *nedges; i++)
      *npins += (*esizes)[i];

    if (*npins) {
      *pins = ZOLTAN_MALLOC_GID_ARRAY(zz, *npins);
      *pin_procs = (int *) ZOLTAN_MALLOC(*npins * sizeof(int));
      if (!*pins || !*pin_procs) MEMORY_ERROR;

      ierr = zz->Get_HG_Edge_List(zz->Get_HG_Edge_List_Data, 
                                  num_gid_entries, num_lid_entries, *nedges, 
                                  *egids, *elids,
                                  *esizes, *pins, 
                                  *pin_procs);

      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo,"Error returned from Get_HG_Edge_List");
        goto End;
      }
    }
  }
  
  /* 
   * KDDKDD -- Assuming hyperedges are given to Zoltan by one processor only.
   * KDDKDD -- Eventually, will loosen that constraint and remove duplicates.
   * KDDKDD -- Or the code might work (although with extra communication)
   * KDDKDD -- with the duplicates.
   * KDDKDD -- Might be easier once we have edge GIDs.
   */

End:
  
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
