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

#ifndef __PARAMS_CONST_H
#define __PARAMS_CONST_H

#include "lb_const.h"

#ifndef HAVE_PROTOTYPES
#   if defined(__STDC__) || defined(__GNUC__) || defined(__cplusplus) || defined(c_plusplus)
#       define	HAVE_PROTOTYPES
#   endif
#endif

#undef PROTO
#ifdef HAVE_PROTOTYPES
#   define	PROTO(x)	x
#else
#   define	PROTO(x)	()
#endif

/*
 * Type used to pass list of available parameters to lower routines.
 */

typedef struct Param_Vars {
  char *name;			/* Parameter variable name (all CAPS) */
  void *ptr;			/* Pointer to parameter variable */
  char *type;			/* type of parameter: */
				/* INT, DOUBLE, LONG, STRING, or CHAR */
} PARAM_VARS;

/* string length limit for param val. Allocate space of this + 1 */
#define MAX_PARAM_STRING_LEN 50

/* Universal parameter value struct */
typedef struct Param_Utype {
  int def;   /* default flag */
  int ival;
  double dval;
  long lval;
  char sval[MAX_PARAM_STRING_LEN + 1];
  char cval;
} PARAM_UTYPE;

/* API for general parameter setting functions */
typedef int ZOLTAN_SET_PARAM_FN(char *, char *); 

/* function declarations for parameter modification routines */

extern int Zoltan_Assign_Param_Vals(ZOLTAN_PARAM *, PARAM_VARS *, int, int, int);
extern int Zoltan_Bind_Param(PARAM_VARS *, char *, void *);
extern int Zoltan_Set_Key_Param(ZZ *, char *, char *);
extern void Zoltan_Print_Key_Params(ZZ *);
extern void Zoltan_Print_Params(ZOLTAN_PARAM *ptr);
extern int Zoltan_Check_Param(char *, char *, PARAM_VARS *,
    PARAM_UTYPE *, int *);
extern void Zoltan_Free_Params(ZOLTAN_PARAM **);

#endif
