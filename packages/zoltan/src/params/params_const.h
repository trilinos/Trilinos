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


#ifndef __PARAMS_CONST_H
#define __PARAMS_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


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
				/* INT, FLOAT, DOUBLE, LONG, STRING, or CHAR */
  int length;			/* length of vector; 0 if scalar */
} PARAM_VARS;

/* string length limit for param val. Allocate space of this + 1 */
#define MAX_PARAM_STRING_LEN 50

/* Universal parameter value struct. (This can be a C union to save space.) */
typedef struct Param_Utype {
  int def;   /* default flag */
  int ival;
  float fval;
  double dval;
  long lval;
  char sval[MAX_PARAM_STRING_LEN + 1];
  char cval;
} PARAM_UTYPE;

/*
 * Type used to store linked list of new values for parameters.
 */

typedef struct Param_List {
  char *name;
  int index;
  char *new_val;
  struct Param_List *next;
} PARAM_LIST;

/* API for general parameter setting functions */
typedef int ZOLTAN_SET_PARAM_FN(char *, char *); 
typedef int ZOLTAN_SET_PARAM_VEC_FN(char *, int, char *);

/* function declarations for parameter modification routines */

extern int Zoltan_Assign_Param_Vals(PARAM_LIST *, PARAM_VARS *, int, int, int);
extern int Zoltan_Bind_Param(PARAM_VARS *, char *, void *);
extern int Zoltan_Bind_Param_Vec(PARAM_VARS *, char *, void *, int);
extern void Zoltan_Print_Params(PARAM_LIST *ptr);
extern int Zoltan_Check_Param(char *, char *, PARAM_VARS *,
    PARAM_UTYPE *, int *);
extern void Zoltan_Free_Params(PARAM_LIST **);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
