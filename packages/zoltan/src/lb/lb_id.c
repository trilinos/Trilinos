/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lb_const.h"
#include "lb_id_const.h"

/*****************************************************************************/
/*
 *  This file contains routines for manipulating 
 *  the global and local IDs used by Zoltan.
 *
 *  Some manipulations are performed via macros.  In particular, macros
 *  specifying whether global or local IDs are to be manipulated are
 *  provided.  Use of these macros is recommended over use of these
 *  basic functions.
 *  See lb_id_const.h for definitions of these macros.  
 */
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Routines for allocating and initializing IDs.
 */

LB_ID_PTR LB_Malloc_ID(int n, char *file, int line)
{
/*
 * Allocates an array of size n of LB_ID_TYPEs and initializes them.
 */

LB_ID_PTR tmp;
char *yo = "LB_Malloc_ID";

  tmp = (LB_ID_PTR) LB_Malloc(n * sizeof(LB_ID_TYPE), file, line);

  if (tmp != NULL) {
    LB_INIT_ID(n,tmp);
  }
  else if (n > 0) {
    char msg[256];
    sprintf(msg, "NULL pointer returned; malloc called from %s, line %d.",
            file, line);
    LB_PRINT_ERROR(-1, yo, msg);
  }

  return tmp;
}

/*****************************************************************************/
/*****************************************************************************/
/*
 *  Routines for printing IDs.
 */

void LB_PRINT_ID(int n, LB_ID_PTR a)
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
 *  Routines to compare Global IDs.  
 *  Functions are provided to test whether two LB_GIDs are equal (EQ),
 *  less than (LT), and greater than (GT).
 *  The negation operator can be used to test whether two LB_GIDs are 
 *  not equal (!LB_EQ_GID(lb,a,b)), less than or equal (!LB_GT_GID(lb,a,b))
 *  or greater than or equal (!LB_LT_GID(lb,a,b)).
 *  Comparison functions are not needed for LB_LIDs as LB_LIDs are not used
 *  within the load-balancing routines.
 */
/*****************************************************************************/

int LB_EQ_ID(int n, LB_ID_PTR a, LB_ID_PTR b)
{
/* 
 * Returns TRUE if a == b; FALSE otherwise.
 * a == b if for all i, a[i] == b[i].
 */
int i;
  for (i = 0; i < n; i++)
    if (a[i] != b[i])
      return(FALSE);
  return(TRUE);
}

/*****************************************************************************/

int LB_LT_ID(int n, LB_ID_PTR a, LB_ID_PTR b)
{
/* 
 * Returns TRUE if a < b; FALSE otherwise.
 * a < b if for some i, a[i] < b[i] and a[j] == b[j] for all j < i.
 */
int i;
  for (i = 0; i < n; i++)
    if (a[i] == b[i])
      continue;
    else if (a[i] > b[i])
      return(FALSE);
    else /* a[i] < b[i] */
      return(TRUE);

  return(FALSE); /* because a == b */
}

/*****************************************************************************/

int LB_GT_ID(int n, LB_ID_PTR a, LB_ID_PTR b)
{
/* 
 * Returns TRUE if a < b; FALSE otherwise.
 * a > b if for some i, a[i] > b[i] and a[j] == b[j] for all j < i.
 */
int i;
  for (i = 0; i < n; i++)
    if (a[i] == b[i])
      continue;
    else if (a[i] < b[i])
      return(FALSE);
    else /* a[i] > b[i] */
      return(TRUE);

  return(FALSE); /* because a == b */
}
