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


#include <stdio.h>
#include <memory.h>
#include "zz_const.h"
#include "zz_rand.h"
#include "params_const.h"
#include "all_allo_const.h"

/* Random partitioning, does not attempt balance! */
/* Each processor selects a random subset of objects (default is all)
   to give to random partitions (and processors). */

/*****************************************************************************/
/*  Parameters structure for Random method. */

/*
No parameters for now; uncomment to add!
static PARAM_VARS Random_params[] = {
                  { NULL, NULL, NULL, 0 } };
*/
/*****************************************************************************/

int Zoltan_Random2d(
  ZZ *zz                       /* The Zoltan structure.                     */
)
{
  int ierr = ZOLTAN_OK;
  int i, have_pins, num_obj;
  ZOLTAN_ID_PTR rows = NULL;
  ZOLTAN_ID_PTR cols = NULL; 
  int *parts = NULL;
  static char *yo = "Zoltan_Random2d";
  int nl, np, format;
  int *cptr = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Require hypergraph queries.
     TODO: Replace with matrix queries.
   */
  if (!zz->Get_HG_Size_CS || !zz->Get_HG_CS){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Hypergraph query functions undefined");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_FATAL;
  }

  /* Get size and type of compressed pin storage */
  zz->Get_HG_Size_CS(zz->Get_HG_Size_CS_Data, &nl, &np, &format, &ierr);

  if ((format != ZOLTAN_COMPRESSED_EDGE)&&(format != ZOLTAN_COMPRESSED_VERTEX)){    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
      "Invalid matrix compression format returned in Get_HG_Size_CS");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_FATAL;
  }
  ZOLTAN_TRACE_DETAIL(zz, yo, "done with Get_HG_Size_CS");

  have_pins = ((nl > 0) && (np > 0));
  if (format != ZOLTAN_COMPRESSED_EDGE){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "ZOLTAN_COMPRESSED_EDGE is required.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_FATAL;
  }

  /* Get the hypergraph pins in compressed storage format */

  if (have_pins){

    rows = ZOLTAN_MALLOC_GID_ARRAY(zz, nl);
    cptr = (int *)ZOLTAN_MALLOC(nl * sizeof(int));
    cols = ZOLTAN_MALLOC_GID_ARRAY(zz, np);

    if (!rows|| !cptr || !cols){
      Zoltan_Multifree(__FILE__, __LINE__, 3, &rows, &cptr, &cols);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "memory allocation");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_FATAL;
    }
    zz->Get_HG_CS(zz->Get_HG_CS_Data, zz->Num_GID,
            nl, np, format, rows, cptr, cols, &ierr);

    /* Randomly assign partitions to pins. */

    ZOLTAN_TRACE_DETAIL(zz, yo, "done with Get_HG_CS");

  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

