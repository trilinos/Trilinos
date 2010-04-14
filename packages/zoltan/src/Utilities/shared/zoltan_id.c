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
#include "zoltan_types.h"
#include "zoltan_id.h"
#include "zoltan_util.h"
#include "zoltan_mem.h"

/*****************************************************************************/
/*
 *  This file contains routines for manipulating 
 *  the global and local IDs used by Zoltan.
 *
 *  Some manipulations are performed via macros.  In particular, macros
 *  specifying whether global or local IDs are to be manipulated are
 *  provided.  Use of these macros is recommended over use of these
 *  basic functions.
 *  See zoltan_id.h for definitions of these macros.  
 */
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Routines for allocating and initializing IDs.
 */

ZOLTAN_ID_PTR ZOLTAN_Malloc_ID(int n, char *file, int line)
{
/*
 * Allocates an array of size n of ZOLTAN_ID_TYPEs and initializes them.
 */

ZOLTAN_ID_PTR tmp;
char *yo = "ZOLTAN_Malloc_ID";

  /* 
   * Don't use ZOLTAN_MALLOC macro here; prefer to pass file and line 
   * where ZOLTAN_Malloc_ID was called.
   */
  tmp = (ZOLTAN_ID_PTR) Zoltan_Malloc(n * sizeof(ZOLTAN_ID_TYPE), file, line);

  if (tmp != NULL) {
    ZOLTAN_INIT_ID(n,tmp);
  }
  else if (n > 0) {
    char msg[256];
    sprintf(msg, "NULL pointer returned; malloc called from %s, line %d.",
            file, line);
    ZOLTAN_PRINT_ERROR(-1, yo, msg);
  }

  return tmp;
}

/*****************************************************************************/
/*****************************************************************************/
/*
 *  Routines for printing IDs.
 */

void ZOLTAN_PRINT_ID(int n, ZOLTAN_ID_PTR a)
{
/* Prints a single ID. */
int i;

  printf("(");
  for (i = 0; i < n; i++)
    printf("%u ",a[i]);
  printf(") ");
}
  
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Routine to compare Global IDs.  Return 1 if equal.
 *  For performance reasons, this function has been replaced by a macro.
 */
/*****************************************************************************/

#ifdef USE_ZOLTAN_EQ_ID_FUNC
int ZOLTAN_EQ_ID(int n, ZOLTAN_ID_PTR a, ZOLTAN_ID_PTR b)
{
/* 
 * Returns 1 if a == b; 0 otherwise.
 * a == b if for all i, a[i] == b[i].
 */
int i;
  for (i = 0; i < n; i++)
    if (a[i] != b[i])
      return(0);
  return(1);
}
#endif

/*****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
