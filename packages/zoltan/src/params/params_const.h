/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


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
  char sval[MAX_PARAM_STRING_LEN];
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
extern int Zoltan_Check_Param(const char *, const char *, PARAM_VARS *,
    PARAM_UTYPE *, int *);
extern void Zoltan_Free_Params(PARAM_LIST **);
extern int Zoltan_Copy_Params(PARAM_LIST **to, PARAM_LIST const *from);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
