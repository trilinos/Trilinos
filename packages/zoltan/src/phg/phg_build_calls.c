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

#include "phg.h"
#include "zz_const.h"
#include "zz_util_const.h"
    
#define MEMORY_ERROR { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
  ierr = ZOLTAN_MEMERR; \
  goto End; \
}


int Zoltan_PHG_Hypergraph_Callbacks(
  ZZ *zz,
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
static char *yo = "Zoltan_PHG_Hypergraph_Callbacks";
int ierr = ZOLTAN_OK;
int i;
int nwgt;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;

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
    nwgt = *nedges * zz->Edge_Weight_Dim;
    if (nwgt) 
      *ewgts = (float *) ZOLTAN_MALLOC(nwgt * sizeof(float));
    if (!*esizes || !*egids || (num_lid_entries && !*elids) ||
        (nwgt && !*ewgts)) MEMORY_ERROR;

    ierr = zz->Get_HG_Edge_Info(zz->Get_HG_Edge_Info_Data,
                                     num_gid_entries, num_lid_entries,
                                     *nedges, zz->Edge_Weight_Dim, 
                                     *egids, *elids, *esizes, *ewgts); 
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_HG_Edge_Info");
      goto End;
    }
                     
    /* KDDKDD ADD HERE:  DENSE EDGE REMOVAL FROM egids/elids */

    *npins = 0;
    for (i = 0; i < *nedges; i++)
      *npins += (*esizes)[i];

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
