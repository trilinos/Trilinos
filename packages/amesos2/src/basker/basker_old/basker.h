// @HEADER
// *****************************************************************************
//                   Basker: A Direct Linear Solver package
//
// Copyright 2011 NTESS and the Basker contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

/* =================== basker.h =============================================== */
/* basker.c decribes these methods */
int basker_basker
(
	int Ap [],
	int Ai [],
	double Ax [],
	int anrow,
	int ancol,
	int ws [],
	double X [],
	int *Lp,
	int **Li,
	double **Lx,
	int *Up,
	int **Ui,
	double **Ux,
	int *lnnz,
	int *unnz,
	int *pinv
) ;

long basker_basker_l
(
	long Ap [],
	long Ai [],
	double Ax [],
	long anrow,
	long ancol,
	long ws [],
	double X [],
	long *Lp,
	long **Li,
	double **Lx,
	long *Up,
	long **Ui,
	double **Ux,
	long *lnnz,
	long *unnz,
	long *pinv
) ;
/* Uncomment next line to turn the debugging on. */
/*#define DEBUG*/


#ifdef MATLAB_MEX_FILE
#include "matrix.h"
#include "mex.h"

#ifdef DEBUG
#define BASKERASSERT(a) mxAssert(a, "")
#else
#define BASKERASSERT(a)
#endif

#define BASKERREALLOC(ptr, size) mxRealloc(ptr, size)

#else /* MATLAB_MEX_FILE */

#include <stdio.h>
#include <assert.h>

#ifdef DEBUG
#define BASKERASSERT(a) assert(a)
#else
#define BASKERASSERT(a)
#endif

#define BASKERREALLOC(ptr, size) realloc(ptr, size)

#endif /* MATLAB_MEX_FILE */

#ifdef DEBUG  /* DEBUG */

#define PRINT(params) printf params

#else /* DEBUG */

#define PRINT(params)

#endif /* DEBUG */
