/* ========================================================================== */
/* === Cholesky/cholmod_lsolve ============================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Cholesky version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Solve Lx=b with unit or non-unit diagonal, or solve LDx=b.
 * This file is not compiled separately.  It is #include'd in cholmod_solve.c
 * instead.
 *
 * workspace: none
 *
 * Macros for handling different types of factors:
 *
 * LSOLVE(k)	defines the name of a routine for an n-by-k right-hand-side.
 * LNZ(j)	returns the number of nonzeros in column L, incl. the diagonal.
 * PEND		finds the start (p) and end (pend) of column j, and the #
 *		nonzeros in column j, including the diagonal (lnz).
 * GET_LNZ	gets a copy of L->nz pointer, for unpacked or dynamic factors
 */

/* undefine all prior definitions */
#undef FORM_NAME
#undef LSOLVE
#undef LNZ
#undef PEND
#undef GET_LNZ
#undef DIAG

/* -------------------------------------------------------------------------- */
/* define the method */
/* -------------------------------------------------------------------------- */

#if defined (LL)
/* LL' packed, solve Lx=b with non-unit diagonal */
#define FORM_NAME(rank) ll_lsolve_ ## rank
#define DIAG
#define PACKED

#elif defined (LD) && !defined (PACKED)
/* LDL' unpacked or dynamic, solve LDx=b */
#define FORM_NAME(rank) ldl_unpacked_ldsolve_ ## rank
#define DIAG

#elif defined (LD) && defined (PACKED)
/* LDL' packed, solve LDx=b */
#define FORM_NAME(rank) ldl_packed_ldsolve_ ## rank
#define DIAG

#elif defined (PACKED)
/* LDL' packed, solve Lx=b with unit diagonal */
#define FORM_NAME(rank) ldl_packed_lsolve_ ## rank

#else
/* LDL' unpacked or dynamic, solve Lx=b with unit diagonal */
#define FORM_NAME(rank) ldl_unpacked_lsolve_ ## rank

#endif

/* -------------------------------------------------------------------------- */
/* define the name of the routine */
/* -------------------------------------------------------------------------- */

#define LSOLVE(rank) FORM_NAME(rank)

/* -------------------------------------------------------------------------- */
/* how to access a column of L: packed or unpacked */
/* -------------------------------------------------------------------------- */

#ifdef PACKED
/* L is packed, stored in Li/Lx [Lp [j] ... Lp [j+1]-1] */
#define LNZ(j) (Lp [j+1] - Lp [j])
#define PEND(j,p,pend,lnz) { p = Lp [j] ; pend = Lp [j+1] ; lnz = pend - p ; }
#define GET_LNZ
#else
/* L is unpacked, stored in Li/Lx [Lp [j] ... Lp [j] + Lnz [j]] */
#define LNZ(j) (Lnz [j])
#define PEND(j,p,pend,lnz) { p = Lp [j] ; lnz = Lnz [j] ; pend = p + lnz ; }
#define GET_LNZ int *Lnz = L->nz
#endif

/* ========================================================================== */
/* === LSOLVE (1) =========================================================== */
/* ========================================================================== */

static void LSOLVE (1)
(
    cholmod_factor *L,
    double X [ ]		    /* n-by-1 in row form */
)
{
    double *Lx = L->x ;
    int *Li = L->i ;
    int *Lp = L->p ;
    int j, n = L->n ;
    GET_LNZ ;

    for (j = 0 ; j < n ; )
    {
	/* get the start, end, and length of column j */
	int p, pend, lnz ;
	PEND (j, p, pend, lnz) ;

	/* find a chain of supernodes (up to j, j+1, and j+2) */
	if (lnz < 4 || lnz != LNZ (j+1) + 1 || Li [p+1] != j+1)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y = X [j] ;
#ifdef LL
	    y /= Lx [p] ;
	    X [j] = y ;
#elif defined (LD)
	    X [j] = y / Lx [p] ;
#endif

	    for (p++ ; p < pend ; p++)
	    {
		X [Li [p]] -= Lx [p] * y ;
	    }
	    j++ ;

	}
	else if (lnz != LNZ (j+2) + 2 || Li [p+2] != j+2)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2] ;
	    int q = Lp [j+1] ;

#ifdef LL
	    y [0] = X [j] / Lx [p] ;
	    y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
	    X [j  ] = y [0] ;
	    X [j+1] = y [1] ;
#elif defined (LD)
	    y [0] = X [j] ;
	    y [1] = X [j+1] - Lx [p+1] * y [0] ;
	    X [j  ] = y [0] / Lx [p] ;
	    X [j+1] = y [1] / Lx [q] ;
#else
	    y [0] = X [j] ;
	    y [1] = X [j+1] - Lx [p+1] * y [0] ;
	    X [j+1] = y [1] ;
#endif

	    for (p += 2, q++ ; p < pend ; p++, q++)
	    {
		X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] ;
	    }
	    j += 2 ;

	}
	else
	{


	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3] ;
	    int q = Lp [j+1] ;
	    int r = Lp [j+2] ;
#ifdef LL
	    y [0] = X [j] / Lx [p] ;
	    y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
	    y [2] = (X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1]) / Lx [r] ;
	    X [j  ] = y [0] ;
	    X [j+1] = y [1] ;
	    X [j+2] = y [2] ;
#elif defined (LD)
	    y [0] = X [j] ;
	    y [1] = X [j+1] - Lx [p+1] * y [0] ;
	    y [2] = X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1] ;
	    X [j  ] = y [0] / Lx [p] ;
	    X [j+1] = y [1] / Lx [q] ;
	    X [j+2] = y [2] / Lx [r] ;
#else
	    y [0] = X [j] ;
	    y [1] = X [j+1] - Lx [p+1] * y [0] ;
	    y [2] = X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1] ;
	    X [j+1] = y [1] ;
	    X [j+2] = y [2] ;
#endif

	    for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
	    {
		X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] + Lx [r] * y [2] ;
	    }
	    j += 3 ;
	}
    }
}


/* ========================================================================== */
/* === LSOLVE (2) =========================================================== */
/* ========================================================================== */

/* Solve Lx=b, where b has 2 columns */

static void LSOLVE (2)
(
    cholmod_factor *L,
    double X [ ]		    /* n-by-2 in row form */
)
{
    double *Lx = L->x ;
    int *Li = L->i ;
    int *Lp = L->p ;
    int j, n = L->n ;
    GET_LNZ ;

    for (j = 0 ; j < n ; )
    {
	/* get the start, end, and length of column j */
	int p, pend, lnz ;
	PEND (j, p, pend, lnz) ;

	/* find a chain of supernodes (up to j, j+1, and j+2) */
	if (lnz < 4 || lnz != LNZ (j+1) + 1 || Li [p+1] != j+1)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y [2] ;
	    y [0] = X [2*j  ] ;
	    y [1] = X [2*j+1] ;
#ifdef LL
	    y [0] /= Lx [p] ;
	    y [1] /= Lx [p] ;
	    X [2*j  ] = y [0] ;
	    X [2*j+1] = y [1] ;
#elif defined (LD)
	    X [2*j  ] = y [0] / Lx [p] ;
	    X [2*j+1] = y [1] / Lx [p] ;
#endif

	    for (p++ ; p < pend ; p++)
	    {
		int i = 2 * Li [p] ;
		X [i  ] -= Lx [p] * y [0] ;
		X [i+1] -= Lx [p] * y [1] ;
	    }
	    j++ ;

	}
	else if (lnz != LNZ (j+2) + 2 || Li [p+2] != j+2)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2][2] ;
	    int q = Lp [j+1] ;
	    y [0][0] = X [2*j  ] ;
	    y [0][1] = X [2*j+1] ;
#ifdef LL
	    y [0][0] /= Lx [p] ;
	    y [0][1] /= Lx [p] ;
	    y [1][0] = (X [2*j+2] - Lx [p+1] * y [0][0]) / Lx [q] ;
	    y [1][1] = (X [2*j+3] - Lx [p+1] * y [0][1]) / Lx [q] ;
	    X [2*j  ] = y [0][0] ;
	    X [2*j+1] = y [0][1] ;
	    X [2*j+2] = y [1][0] ;
	    X [2*j+3] = y [1][1] ;
#elif defined (LD)
	    y [1][0] = X [2*j+2] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [2*j+3] - Lx [p+1] * y [0][1] ;
	    X [2*j  ] = y [0][0] / Lx [p] ;
	    X [2*j+1] = y [0][1] / Lx [p] ;
	    X [2*j+2] = y [1][0] / Lx [q] ;
	    X [2*j+3] = y [1][1] / Lx [q] ;
#else
	    y [1][0] = X [2*j+2] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [2*j+3] - Lx [p+1] * y [0][1] ;
	    X [2*j+2] = y [1][0] ;
	    X [2*j+3] = y [1][1] ;
#endif

	    for (p += 2, q++ ; p < pend ; p++, q++)
	    {
		int i = 2 * Li [p] ;
		X [i  ] -= Lx [p] * y [0][0] + Lx [q] * y [1][0] ;
		X [i+1] -= Lx [p] * y [0][1] + Lx [q] * y [1][1] ;
	    }
	    j += 2 ;

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3][2] ;
	    int q = Lp [j+1] ;
	    int r = Lp [j+2] ;
	    y [0][0] = X [2*j  ] ;
	    y [0][1] = X [2*j+1] ;
#ifdef LL
	    y [0][0] /= Lx [p] ;
	    y [0][1] /= Lx [p] ;
	    y [1][0] = (X[2*j+2] - Lx[p+1] * y[0][0]) / Lx [q] ;
	    y [1][1] = (X[2*j+3] - Lx[p+1] * y[0][1]) / Lx [q] ;
	    y [2][0] = (X[2*j+4] - Lx[p+2] * y[0][0] - Lx[q+1] * y[1][0])/Lx[r];
	    y [2][1] = (X[2*j+5] - Lx[p+2] * y[0][1] - Lx[q+1] * y[1][1])/Lx[r];
	    X [2*j  ] = y [0][0] ;
	    X [2*j+1] = y [0][1] ;
	    X [2*j+2] = y [1][0] ;
	    X [2*j+3] = y [1][1] ;
	    X [2*j+4] = y [2][0] ;
	    X [2*j+5] = y [2][1] ;
#elif defined (LD)
	    y [1][0] = X [2*j+2] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [2*j+3] - Lx [p+1] * y [0][1] ;
	    y [2][0] = X [2*j+4] - Lx [p+2] * y [0][0] - Lx [q+1] * y [1][0] ;
	    y [2][1] = X [2*j+5] - Lx [p+2] * y [0][1] - Lx [q+1] * y [1][1] ;
	    X [2*j  ] = y [0][0] / Lx [p] ;
	    X [2*j+1] = y [0][1] / Lx [p] ;
	    X [2*j+2] = y [1][0] / Lx [q] ;
	    X [2*j+3] = y [1][1] / Lx [q] ;
	    X [2*j+4] = y [2][0] / Lx [r] ;
	    X [2*j+5] = y [2][1] / Lx [r] ;
#else
	    y [1][0] = X [2*j+2] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [2*j+3] - Lx [p+1] * y [0][1] ;
	    y [2][0] = X [2*j+4] - Lx [p+2] * y [0][0] - Lx [q+1] * y [1][0] ;
	    y [2][1] = X [2*j+5] - Lx [p+2] * y [0][1] - Lx [q+1] * y [1][1] ;
	    X [2*j+2] = y [1][0] ;
	    X [2*j+3] = y [1][1] ;
	    X [2*j+4] = y [2][0] ;
	    X [2*j+5] = y [2][1] ;
#endif

	    for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
	    {
		int i = 2 * Li [p] ;
		X[i  ] -= Lx[p] * y[0][0] + Lx[q] * y[1][0] + Lx[r] * y[2][0] ;
		X[i+1] -= Lx[p] * y[0][1] + Lx[q] * y[1][1] + Lx[r] * y[2][1] ;
	    }
	    j += 3 ;
	}
    }
}


/* ========================================================================== */
/* === LSOLVE (3) =========================================================== */
/* ========================================================================== */

/* Solve Lx=b, where b has 3 columns */


static void LSOLVE (3)
(
    cholmod_factor *L,
    double X [ ]		    /* n-by-3 in row form */
)
{
    double *Lx = L->x ;
    int *Li = L->i ;
    int *Lp = L->p ;
    int j, n = L->n ;
    GET_LNZ ;

    for (j = 0 ; j < n ; )
    {
	/* get the start, end, and length of column j */
	int p, pend, lnz ;
	PEND (j, p, pend, lnz) ;

	/* find a chain of supernodes (up to j, j+1, and j+2) */
	if (lnz < 4 || lnz != LNZ (j+1) + 1 || Li [p+1] != j+1)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y [3] ;
	    y [0] = X [3*j  ] ;
	    y [1] = X [3*j+1] ;
	    y [2] = X [3*j+2] ;
#ifdef LL
	    y [0] /= Lx [p] ;
	    y [1] /= Lx [p] ;
	    y [2] /= Lx [p] ;
	    X [3*j  ] = y [0] ;
	    X [3*j+1] = y [1] ;
	    X [3*j+2] = y [2] ;
#elif defined (LD)
	    X [3*j  ] = y [0] / Lx [p] ;
	    X [3*j+1] = y [1] / Lx [p] ;
	    X [3*j+2] = y [2] / Lx [p] ;
#endif

	    for (p++ ; p < pend ; p++)
	    {
		int i = 3 * Li [p] ;
		double lx = Lx [p] ;
		X [i  ] -= lx * y [0] ;
		X [i+1] -= lx * y [1] ;
		X [i+2] -= lx * y [2] ;
	    }
	    j++ ;

	}
	else if (lnz != LNZ (j+2) + 2 || Li [p+2] != j+2)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2][3] ;
	    int q = Lp [j+1] ;
	    y [0][0] = X [3*j  ] ;
	    y [0][1] = X [3*j+1] ;
	    y [0][2] = X [3*j+2] ;
#ifdef LL
	    y [0][0] /= Lx [p] ;
	    y [0][1] /= Lx [p] ;
	    y [0][2] /= Lx [p] ;
	    y [1][0] = (X [3*j+3] - Lx [p+1] * y [0][0]) / Lx [q] ;
	    y [1][1] = (X [3*j+4] - Lx [p+1] * y [0][1]) / Lx [q] ;
	    y [1][2] = (X [3*j+5] - Lx [p+1] * y [0][2]) / Lx [q] ;
	    X [3*j  ] = y [0][0] ;
	    X [3*j+1] = y [0][1] ;
	    X [3*j+2] = y [0][2] ;
	    X [3*j+3] = y [1][0] ;
	    X [3*j+4] = y [1][1] ;
	    X [3*j+5] = y [1][2] ;
#elif defined (LD)
	    y [1][0] = X [3*j+3] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [3*j+4] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [3*j+5] - Lx [p+1] * y [0][2] ;
	    X [3*j  ] = y [0][0] / Lx [p] ;
	    X [3*j+1] = y [0][1] / Lx [p] ;
	    X [3*j+2] = y [0][2] / Lx [p] ;
	    X [3*j+3] = y [1][0] / Lx [q] ;
	    X [3*j+4] = y [1][1] / Lx [q] ;
	    X [3*j+5] = y [1][2] / Lx [q] ;
#else
	    y [1][0] = X [3*j+3] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [3*j+4] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [3*j+5] - Lx [p+1] * y [0][2] ;
	    X [3*j+3] = y [1][0] ;
	    X [3*j+4] = y [1][1] ;
	    X [3*j+5] = y [1][2] ;
#endif

	    for (p += 2, q++ ; p < pend ; p++, q++)
	    {
		int i = 3 * Li [p] ;
		double lx [2] ;
		lx [0] = Lx [p] ;
		lx [1] = Lx [q] ;
		X [i  ] -= lx [0] * y [0][0] + lx [1] * y [1][0] ;
		X [i+1] -= lx [0] * y [0][1] + lx [1] * y [1][1] ;
		X [i+2] -= lx [0] * y [0][2] + lx [1] * y [1][2] ;
	    }
	    j += 2 ;

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3][3] ;
	    int q = Lp [j+1] ;
	    int r = Lp [j+2] ;
	    y [0][0] = X [3*j  ] ;
	    y [0][1] = X [3*j+1] ;
	    y [0][2] = X [3*j+2] ;
#ifdef LL
	    y [0][0] /= Lx [p] ;
	    y [0][1] /= Lx [p] ;
	    y [0][2] /= Lx [p] ;
	    y [1][0] = (X[3*j+3] - Lx[p+1] * y[0][0]) / Lx [q] ;
	    y [1][1] = (X[3*j+4] - Lx[p+1] * y[0][1]) / Lx [q] ;
	    y [1][2] = (X[3*j+5] - Lx[p+1] * y[0][2]) / Lx [q] ;
	    y [2][0] = (X[3*j+6] - Lx[p+2] * y[0][0] - Lx[q+1] * y[1][0])/Lx[r];
	    y [2][1] = (X[3*j+7] - Lx[p+2] * y[0][1] - Lx[q+1] * y[1][1])/Lx[r];
	    y [2][2] = (X[3*j+8] - Lx[p+2] * y[0][2] - Lx[q+1] * y[1][2])/Lx[r];
	    X [3*j  ] = y [0][0] ;
	    X [3*j+1] = y [0][1] ;
	    X [3*j+2] = y [0][2] ;
	    X [3*j+3] = y [1][0] ;
	    X [3*j+4] = y [1][1] ;
	    X [3*j+5] = y [1][2] ;
	    X [3*j+6] = y [2][0] ;
	    X [3*j+7] = y [2][1] ;
	    X [3*j+8] = y [2][2] ;
#elif defined (LD)
	    y [1][0] = X [3*j+3] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [3*j+4] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [3*j+5] - Lx [p+1] * y [0][2] ;
	    y [2][0] = X [3*j+6] - Lx [p+2] * y [0][0] - Lx [q+1] * y [1][0] ;
	    y [2][1] = X [3*j+7] - Lx [p+2] * y [0][1] - Lx [q+1] * y [1][1] ;
	    y [2][2] = X [3*j+8] - Lx [p+2] * y [0][2] - Lx [q+1] * y [1][2] ;
	    X [3*j  ] = y [0][0] / Lx [p] ;
	    X [3*j+1] = y [0][1] / Lx [p] ;
	    X [3*j+2] = y [0][2] / Lx [p] ;
	    X [3*j+3] = y [1][0] / Lx [q] ;
	    X [3*j+4] = y [1][1] / Lx [q] ;
	    X [3*j+5] = y [1][2] / Lx [q] ;
	    X [3*j+6] = y [2][0] / Lx [r] ;
	    X [3*j+7] = y [2][1] / Lx [r] ;
	    X [3*j+8] = y [2][2] / Lx [r] ;
#else
	    y [1][0] = X [3*j+3] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [3*j+4] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [3*j+5] - Lx [p+1] * y [0][2] ;
	    y [2][0] = X [3*j+6] - Lx [p+2] * y [0][0] - Lx [q+1] * y [1][0] ;
	    y [2][1] = X [3*j+7] - Lx [p+2] * y [0][1] - Lx [q+1] * y [1][1] ;
	    y [2][2] = X [3*j+8] - Lx [p+2] * y [0][2] - Lx [q+1] * y [1][2] ;
	    X [3*j+3] = y [1][0] ;
	    X [3*j+4] = y [1][1] ;
	    X [3*j+5] = y [1][2] ;
	    X [3*j+6] = y [2][0] ;
	    X [3*j+7] = y [2][1] ;
	    X [3*j+8] = y [2][2] ;
#endif

	    for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
	    {
		int i = 3 * Li [p] ;
		double lx [3] ;
		lx [0] = Lx [p] ;
		lx [1] = Lx [q] ;
		lx [2] = Lx [r] ;
		X [i  ] -= lx[0] * y[0][0] + lx[1] * y[1][0] + lx[2] * y[2][0] ;
		X [i+1] -= lx[0] * y[0][1] + lx[1] * y[1][1] + lx[2] * y[2][1] ;
		X [i+2] -= lx[0] * y[0][2] + lx[1] * y[1][2] + lx[2] * y[2][2] ;
	    }
	    j += 3 ;
	}
    }
}


/* ========================================================================== */
/* === LSOLVE (4) =========================================================== */
/* ========================================================================== */

/* Solve Lx=b, where b has 4 columns */

static void LSOLVE (4)
(
    cholmod_factor *L,
    double X [ ]		    /* n-by-4 in row form */
)
{
    double *Lx = L->x ;
    int *Li = L->i ;
    int *Lp = L->p ;
    int j, n = L->n ;
    GET_LNZ ;

    for (j = 0 ; j < n ; )
    {
	/* get the start, end, and length of column j */
	int p, pend, lnz ;
	PEND (j, p, pend, lnz) ;

	/* find a chain of supernodes (up to j, j+1, and j+2) */
	if (lnz < 4 || lnz != LNZ (j+1) + 1 || Li [p+1] != j+1)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y [4] ;
	    y [0] = X [4*j  ] ;
	    y [1] = X [4*j+1] ;
	    y [2] = X [4*j+2] ;
	    y [3] = X [4*j+3] ;
#ifdef LL
	    y [0] /= Lx [p] ;
	    y [1] /= Lx [p] ;
	    y [2] /= Lx [p] ;
	    y [3] /= Lx [p] ;
	    X [4*j  ] = y [0] ;
	    X [4*j+1] = y [1] ;
	    X [4*j+2] = y [2] ;
	    X [4*j+3] = y [3] ;
#elif defined (LD)
	    X [4*j  ] = y [0] / Lx [p] ;
	    X [4*j+1] = y [1] / Lx [p] ;
	    X [4*j+2] = y [2] / Lx [p] ;
	    X [4*j+3] = y [3] / Lx [p] ;
#endif

	    for (p++ ; p < pend ; p++)
	    {
		int i = 4 * Li [p] ;
		double lx = Lx [p] ;
		X [i  ] -= lx * y [0] ;
		X [i+1] -= lx * y [1] ;
		X [i+2] -= lx * y [2] ;
		X [i+3] -= lx * y [3] ;
	    }
	    j++ ;

	}
	else if (lnz != LNZ (j+2) + 2 || Li [p+2] != j+2)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2][4] ;
	    int q = Lp [j+1] ;
	    y [0][0] = X [4*j  ] ;
	    y [0][1] = X [4*j+1] ;
	    y [0][2] = X [4*j+2] ;
	    y [0][3] = X [4*j+3] ;
#ifdef LL
	    y [0][0] /= Lx [p] ;
	    y [0][1] /= Lx [p] ;
	    y [0][2] /= Lx [p] ;
	    y [0][3] /= Lx [p] ;
	    y [1][0] = (X [4*j+4] - Lx [p+1] * y [0][0]) / Lx [q] ;
	    y [1][1] = (X [4*j+5] - Lx [p+1] * y [0][1]) / Lx [q] ;
	    y [1][2] = (X [4*j+6] - Lx [p+1] * y [0][2]) / Lx [q] ;
	    y [1][3] = (X [4*j+7] - Lx [p+1] * y [0][3]) / Lx [q] ;
	    X [4*j  ] = y [0][0] ;
	    X [4*j+1] = y [0][1] ;
	    X [4*j+2] = y [0][2] ;
	    X [4*j+3] = y [0][3] ;
	    X [4*j+4] = y [1][0] ;
	    X [4*j+5] = y [1][1] ;
	    X [4*j+6] = y [1][2] ;
	    X [4*j+7] = y [1][3] ;
#elif defined (LD)
	    y [1][0] = X [4*j+4] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [4*j+5] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [4*j+6] - Lx [p+1] * y [0][2] ;
	    y [1][3] = X [4*j+7] - Lx [p+1] * y [0][3] ;
	    X [4*j  ] = y [0][0] / Lx [p] ;
	    X [4*j+1] = y [0][1] / Lx [p] ;
	    X [4*j+2] = y [0][2] / Lx [p] ;
	    X [4*j+3] = y [0][3] / Lx [p] ;
	    X [4*j+4] = y [1][0] / Lx [q] ;
	    X [4*j+5] = y [1][1] / Lx [q] ;
	    X [4*j+6] = y [1][2] / Lx [q] ;
	    X [4*j+7] = y [1][3] / Lx [q] ;
#else
	    y [1][0] = X [4*j+4] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [4*j+5] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [4*j+6] - Lx [p+1] * y [0][2] ;
	    y [1][3] = X [4*j+7] - Lx [p+1] * y [0][3] ;
	    X [4*j+4] = y [1][0] ;
	    X [4*j+5] = y [1][1] ;
	    X [4*j+6] = y [1][2] ;
	    X [4*j+7] = y [1][3] ;
#endif

	    for (p += 2, q++ ; p < pend ; p++, q++)
	    {
		int i = 4 * Li [p] ;
		double lx [2] ;
		lx [0] = Lx [p] ;
		lx [1] = Lx [q] ;
		X [i  ] -= lx [0] * y [0][0] + lx [1] * y [1][0] ;
		X [i+1] -= lx [0] * y [0][1] + lx [1] * y [1][1] ;
		X [i+2] -= lx [0] * y [0][2] + lx [1] * y [1][2] ;
		X [i+3] -= lx [0] * y [0][3] + lx [1] * y [1][3] ;
	    }
	    j += 2 ;

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3][4] ;
	    int q = Lp [j+1] ;
	    int r = Lp [j+2] ;
	    y [0][0] = X [4*j  ] ;
	    y [0][1] = X [4*j+1] ;
	    y [0][2] = X [4*j+2] ;
	    y [0][3] = X [4*j+3] ;
#ifdef LL
	    y [0][0] /= Lx [p] ;
	    y [0][1] /= Lx [p] ;
	    y [0][2] /= Lx [p] ;
	    y [0][3] /= Lx [p] ;
	    y [1][0] = (X[4*j+4] - Lx[p+1] * y[0][0]) / Lx [q] ;
	    y [1][1] = (X[4*j+5] - Lx[p+1] * y[0][1]) / Lx [q] ;
	    y [1][2] = (X[4*j+6] - Lx[p+1] * y[0][2]) / Lx [q] ;
	    y [1][3] = (X[4*j+7] - Lx[p+1] * y[0][3]) / Lx [q] ;
	    y [2][0] = (X[4*j+8] - Lx[p+2] * y[0][0] - Lx[q+1] * y[1][0])/Lx[r];
	    y [2][1] = (X[4*j+9] - Lx[p+2] * y[0][1] - Lx[q+1] * y[1][1])/Lx[r];
	    y [2][2] = (X[4*j+10]- Lx[p+2] * y[0][2] - Lx[q+1] * y[1][2])/Lx[r];
	    y [2][3] = (X[4*j+11]- Lx[p+2] * y[0][3] - Lx[q+1] * y[1][3])/Lx[r];
	    X [4*j  ] = y [0][0] ;
	    X [4*j+1] = y [0][1] ;
	    X [4*j+2] = y [0][2] ;
	    X [4*j+3] = y [0][3] ;
	    X [4*j+4] = y [1][0] ;
	    X [4*j+5] = y [1][1] ;
	    X [4*j+6] = y [1][2] ;
	    X [4*j+7] = y [1][3] ;
	    X [4*j+8] = y [2][0] ;
	    X [4*j+9] = y [2][1] ;
	    X [4*j+10]= y [2][2] ;
	    X [4*j+11]= y [2][3] ;
#elif defined (LD)
	    y [1][0] = X [4*j+4] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [4*j+5] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [4*j+6] - Lx [p+1] * y [0][2] ;
	    y [1][3] = X [4*j+7] - Lx [p+1] * y [0][3] ;
	    y [2][0] = X [4*j+8] - Lx [p+2] * y [0][0] - Lx [q+1] * y [1][0] ;
	    y [2][1] = X [4*j+9] - Lx [p+2] * y [0][1] - Lx [q+1] * y [1][1] ;
	    y [2][2] = X [4*j+10]- Lx [p+2] * y [0][2] - Lx [q+1] * y [1][2] ;
	    y [2][3] = X [4*j+11]- Lx [p+2] * y [0][3] - Lx [q+1] * y [1][3] ;
	    X [4*j  ] = y [0][0] / Lx [p] ;
	    X [4*j+1] = y [0][1] / Lx [p] ;
	    X [4*j+2] = y [0][2] / Lx [p] ;
	    X [4*j+3] = y [0][3] / Lx [p] ;
	    X [4*j+4] = y [1][0] / Lx [q] ;
	    X [4*j+5] = y [1][1] / Lx [q] ;
	    X [4*j+6] = y [1][2] / Lx [q] ;
	    X [4*j+7] = y [1][3] / Lx [q] ;
	    X [4*j+8] = y [2][0] / Lx [r] ;
	    X [4*j+9] = y [2][1] / Lx [r] ;
	    X [4*j+10]= y [2][2] / Lx [r] ;
	    X [4*j+11]= y [2][3] / Lx [r] ;
#else
	    y [1][0] = X [4*j+4] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [4*j+5] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [4*j+6] - Lx [p+1] * y [0][2] ;
	    y [1][3] = X [4*j+7] - Lx [p+1] * y [0][3] ;
	    y [2][0] = X [4*j+8] - Lx [p+2] * y [0][0] - Lx [q+1] * y [1][0] ;
	    y [2][1] = X [4*j+9] - Lx [p+2] * y [0][1] - Lx [q+1] * y [1][1] ;
	    y [2][2] = X [4*j+10]- Lx [p+2] * y [0][2] - Lx [q+1] * y [1][2] ;
	    y [2][3] = X [4*j+11]- Lx [p+2] * y [0][3] - Lx [q+1] * y [1][3] ;
	    X [4*j+4] = y [1][0] ;
	    X [4*j+5] = y [1][1] ;
	    X [4*j+6] = y [1][2] ;
	    X [4*j+7] = y [1][3] ;
	    X [4*j+8] = y [2][0] ;
	    X [4*j+9] = y [2][1] ;
	    X [4*j+10]= y [2][2] ;
	    X [4*j+11]= y [2][3] ;
#endif

	    for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
	    {
		int i = 4 * Li [p] ;
		double lx [3] ;
		lx [0] = Lx [p] ;
		lx [1] = Lx [q] ;
		lx [2] = Lx [r] ;
		X [i  ] -= lx[0] * y[0][0] + lx[1] * y[1][0] + lx[2] * y[2][0] ;
		X [i+1] -= lx[0] * y[0][1] + lx[1] * y[1][1] + lx[2] * y[2][1] ;
		X [i+2] -= lx[0] * y[0][2] + lx[1] * y[1][2] + lx[2] * y[2][2] ;
		X [i+3] -= lx[0] * y[0][3] + lx[1] * y[1][3] + lx[2] * y[2][3] ;
	    }
	    j += 3 ;
	}
    }
}


/* ========================================================================== */
/* === LSOLVE (k) =========================================================== */
/* ========================================================================== */

static void LSOLVE (k)
(
    cholmod_factor *L,
    double X [ ],		    /* n-by-(1,2,3, or 4) in row form */
    int nr
)
{
    switch (nr)
    {
	case 1: LSOLVE (1) (L, X) ; break ;
	case 2: LSOLVE (2) (L, X) ; break ;
	case 3: LSOLVE (3) (L, X) ; break ;
	case 4: LSOLVE (4) (L, X) ; break ;
    }
}

/* prepare for the next inclusion of this file in cholmod_solve.c */
#undef LL
#undef LD
#undef PACKED
