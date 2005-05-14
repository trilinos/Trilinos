/* ========================================================================== */
/* === Cholesky/ltsolve ===================================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Cholesky version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Solve L'x=b with unit or non-unit diagonal, or solve DL'x=b.
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
#define FORM_NAME(rank) ll_ltsolve_ ## rank
#define DIAG
#define PACKED

#elif defined (LD) && !defined (PACKED)
/* LDL' unpacked or dynamic, solve LDx=b */
#define FORM_NAME(rank) ldl_unpacked_dltsolve_ ## rank
#define DIAG

#elif defined (LD) && defined (PACKED)
/* LDL' packed, solve LDx=b */
#define FORM_NAME(rank) ldl_packed_dltsolve_ ## rank
#define DIAG

#elif defined (PACKED)
/* LDL' packed, solve Lx=b with unit diagonal */
#define FORM_NAME(rank) ldl_packed_ltsolve_ ## rank

#else
/* LDL' unpacked or dynamic, solve Lx=b with unit diagonal */
#define FORM_NAME(rank) ldl_unpacked_ltsolve_ ## rank

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

    for (j = n-1 ; j >= 0 ; )
    {
	/* get the start, end, and length of column j */
	int p, pend, lnz ;
	PEND (j, p, pend, lnz) ;

	/* find a chain of supernodes (up to j, j-1, and j-2) */
	if (j == 0 || lnz != LNZ (j-1) - 1 || Li [Lp [j-1]+1] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y = X [j] ;
#ifdef DIAG
	    double d = Lx [p] ;
#endif
#ifdef LD
	    y /= d ;
#endif

	    for (p++ ; p < pend ; p++)
	    {
		y -= Lx [p] * X [Li [p]] ;
	    }

#ifdef LL
	    X [j] = y / d ;
#else
	    X [j] = y ;
#endif

	    j-- ;

	}
	else if (j == 1 || lnz != LNZ (j-2)-2 || Li [Lp [j-2]+2] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2], t ;
	    int q = Lp [j-1] ;
#ifdef DIAG
	    double d [2] ;
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
#endif
	    t = Lx [q+1] ;
#ifdef LD
	    y [0] = X [j  ] / d [0] ;
	    y [1] = X [j-1] / d [1] ;
#else
	    y [0] = X [j  ] ;
	    y [1] = X [j-1] ;
#endif

	    for (p++, q += 2 ; p < pend ; p++, q++)
	    {
		int i = Li [p] ;
		y [0] -= Lx [p] * X [i] ;
		y [1] -= Lx [q] * X [i] ;
	    }

#ifdef LL
	    y [0] /= d [0] ;
	    y [1] = (y [1] - t * y [0]) / d [1] ;
#else
	    y [1] -= t * y [0] ;
#endif

	    X [j  ] = y [0] ;
	    X [j-1] = y [1] ;

	    j -= 2 ;

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3], t [3] ;
	    int q = Lp [j-1] ;
	    int r = Lp [j-2] ;
#ifdef DIAG
	    double d [3] ;
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
	    d [2] = Lx [r] ;
#endif
	    t [0] = Lx [q+1] ;
	    t [1] = Lx [r+1] ;
	    t [2] = Lx [r+2] ;
#ifdef LD
	    y [0] = X [j]   / d [0] ;
	    y [1] = X [j-1] / d [1] ;
	    y [2] = X [j-2] / d [2] ;
#else
	    y [0] = X [j] ;
	    y [1] = X [j-1] ;
	    y [2] = X [j-2] ;
#endif

	    for (p++, q += 2, r += 3 ; p < pend ; p++, q++, r++)
	    {
		int i = Li [p] ;
		y [0] -= Lx [p] * X [i] ;
		y [1] -= Lx [q] * X [i] ;
		y [2] -= Lx [r] * X [i] ;
	    }

	    q = Lp [j-1]  ;
	    r = Lp [j-2]  ;

#ifdef LL
	    y [0] /= d [0] ;
	    y [1] = (y [1] - t [0] * y [0]) / d [1] ;
	    y [2] = (y [2] - t [2] * y [0] - t [1] * y [1]) / d [2] ;
#else
	    y [1] -= t [0] * y [0] ;
	    y [2] -= t [2] * y [0] + t [1] * y [1] ;
#endif

	    X [j-2] = y [2] ;
	    X [j-1] = y [1] ;
	    X [j  ] = y [0] ;

	    j -= 3 ;
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

    for (j = n-1 ; j >= 0 ; )
    {
	/* get the start, end, and length of column j */
	int p, pend, lnz ;
	PEND (j, p, pend, lnz) ;

	/* find a chain of supernodes (up to j, j+1, and j+2) */
	if (j == 0 || lnz != LNZ (j-1) - 1 || Li [Lp [j-1]+1] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y [2] ;
#ifdef DIAG
	    double d = Lx [p] ;
#endif
#ifdef LD
	    y [0] = X [2*j  ] / d ;
	    y [1] = X [2*j+1] / d ;
#else
	    y [0] = X [2*j  ] ;
	    y [1] = X [2*j+1] ;
#endif

	    for (p++ ; p < pend ; p++)
	    {
		int i = 2 * Li [p] ;
		y [0] -= Lx [p] * X [i  ] ;
		y [1] -= Lx [p] * X [i+1] ;
	    }

#ifdef LL
	    X [2*j  ] = y [0] / d ;
	    X [2*j+1] = y [1] / d ;
#else
	    X [2*j  ] = y [0] ;
	    X [2*j+1] = y [1] ;
#endif

	    j-- ;

	}
	else if (j == 1 || lnz != LNZ (j-2)-2 || Li [Lp [j-2]+2] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2][2], t ;
	    int q = Lp [j-1] ;
#ifdef DIAG
	    double d [2] ;
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
#endif
	    t = Lx [q+1] ;
#ifdef LD
	    y [0][0] = X [2*j  ] / d [0] ;
	    y [0][1] = X [2*j+1] / d [0] ;
	    y [1][0] = X [2*j-2] / d [1] ;
	    y [1][1] = X [2*j-1] / d [1] ;
#else
	    y [0][0] = X [2*j  ] ;
	    y [0][1] = X [2*j+1] ;
	    y [1][0] = X [2*j-2] ;
	    y [1][1] = X [2*j-1] ;
#endif

	    for (p++, q += 2 ; p < pend ; p++, q++)
	    {
		int i = 2 * Li [p] ;
		y [0][0] -= Lx [p] * X [i] ;
		y [0][1] -= Lx [p] * X [i+1] ;
		y [1][0] -= Lx [q] * X [i] ;
		y [1][1] -= Lx [q] * X [i+1] ;
	    }

#ifdef LL
	    y [0][0] /= d [0] ;
	    y [0][1] /= d [0] ;
	    y [1][0] = (y [1][0] - t * y [0][0]) / d [1] ;
	    y [1][1] = (y [1][1] - t * y [0][1]) / d [1] ;
#else
	    y [1][0] -= t * y [0][0] ;
	    y [1][1] -= t * y [0][1] ;
#endif
	    X [2*j  ] = y [0][0] ;
	    X [2*j+1] = y [0][1] ;
	    X [2*j-2] = y [1][0] ;
	    X [2*j-1] = y [1][1] ;

	    j -= 2 ;

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3][2], t [3] ;
	    int q = Lp [j-1] ;
	    int r = Lp [j-2] ;
#ifdef DIAG
	    double d [3] ; 
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
	    d [2] = Lx [r] ;
#endif
	    t [0] = Lx [q+1] ;
	    t [1] = Lx [r+1] ;
	    t [2] = Lx [r+2] ;
#ifdef LD
	    y [0][0] = X [2*j  ] / d [0] ;
	    y [0][1] = X [2*j+1] / d [0] ;
	    y [1][0] = X [2*j-2] / d [1] ;
	    y [1][1] = X [2*j-1] / d [1] ;
	    y [2][0] = X [2*j-4] / d [2] ;
	    y [2][1] = X [2*j-3] / d [2] ;
#else
	    y [0][0] = X [2*j  ] ;
	    y [0][1] = X [2*j+1] ;
	    y [1][0] = X [2*j-2] ;
	    y [1][1] = X [2*j-1] ;
	    y [2][0] = X [2*j-4] ;
	    y [2][1] = X [2*j-3] ;
#endif

	    for (p++, q += 2, r += 3 ; p < pend ; p++, q++, r++)
	    {
		int i = 2 * Li [p] ;
		y [0][0] -= Lx [p] * X [i] ; y [0][1] -= Lx [p] * X [i+1] ;
		y [1][0] -= Lx [q] * X [i] ; y [1][1] -= Lx [q] * X [i+1] ;
		y [2][0] -= Lx [r] * X [i] ; y [2][1] -= Lx [r] * X [i+1] ;
	    }

#ifdef LL
	    y [0][0] /= d [0] ;
	    y [0][1] /= d [0] ;
	    y [1][0] = (y [1][0] - t [0] * y [0][0]) / d [1] ;
	    y [1][1] = (y [1][1] - t [0] * y [0][1]) / d [1] ;
	    y [2][0] = (y [2][0] - t [2] * y [0][0] - t [1] * y [1][0]) / d [2];
	    y [2][1] = (y [2][1] - t [2] * y [0][1] - t [1] * y [1][1]) / d [2];
#else
	    y [1][0] -= t [0] * y [0][0] ;
	    y [1][1] -= t [0] * y [0][1] ;
	    y [2][0] -= t [2] * y [0][0] + t [1] * y [1][0] ;
	    y [2][1] -= t [2] * y [0][1] + t [1] * y [1][1] ;
#endif

	    X [2*j  ] = y [0][0] ;
	    X [2*j+1] = y [0][1] ;
	    X [2*j-2] = y [1][0] ;
	    X [2*j-1] = y [1][1] ;
	    X [2*j-4] = y [2][0] ;
	    X [2*j-3] = y [2][1] ;

	    j -= 3 ;

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

    for (j = n-1 ; j >= 0 ; )
    {
	/* get the start, end, and length of column j */
	int p, pend, lnz ;
	PEND (j, p, pend, lnz) ;

	/* find a chain of supernodes (up to j, j+1, and j+2) */
	if (j == 0 || lnz != LNZ (j-1) - 1 || Li [Lp [j-1]+1] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y [3] ;
#ifdef DIAG
	    double d = Lx [p] ;
#endif
#ifdef LD
	    y [0] = X [3*j  ] / d ;
	    y [1] = X [3*j+1] / d ;
	    y [2] = X [3*j+2] / d ;
#else
	    y [0] = X [3*j  ] ;
	    y [1] = X [3*j+1] ;
	    y [2] = X [3*j+2] ;
#endif

	    for (p++ ; p < pend ; p++)
	    {
		int i = 3 * Li [p] ;
		y [0] -= Lx [p] * X [i  ] ;
		y [1] -= Lx [p] * X [i+1] ;
		y [2] -= Lx [p] * X [i+2] ;
	    }

#ifdef LL
	    X [3*j  ] = y [0] / d ;
	    X [3*j+1] = y [1] / d ;
	    X [3*j+2] = y [2] / d ;
#else
	    X [3*j  ] = y [0] ;
	    X [3*j+1] = y [1] ;
	    X [3*j+2] = y [2] ;
#endif

	    j-- ;

	}
	else if (j == 1 || lnz != LNZ (j-2)-2 || Li [Lp [j-2]+2] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2][3], t ;
	    int q = Lp [j-1] ;
#ifdef DIAG
	    double d [2] ;
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
#endif
	    t = Lx [q+1] ;
#ifdef LD
	    y [0][0] = X [3*j  ] / d [0] ;
	    y [0][1] = X [3*j+1] / d [0] ;
	    y [0][2] = X [3*j+2] / d [0] ;
	    y [1][0] = X [3*j-3] / d [1] ;
	    y [1][1] = X [3*j-2] / d [1] ;
	    y [1][2] = X [3*j-1] / d [1] ;
#else
	    y [0][0] = X [3*j  ] ;
	    y [0][1] = X [3*j+1] ;
	    y [0][2] = X [3*j+2] ;
	    y [1][0] = X [3*j-3] ;
	    y [1][1] = X [3*j-2] ;
	    y [1][2] = X [3*j-1] ;
#endif

	    for (p++, q += 2 ; p < pend ; p++, q++)
	    {
		int i = 3 * Li [p] ;
		y [0][0] -= Lx [p] * X [i] ;
		y [0][1] -= Lx [p] * X [i+1] ;
		y [0][2] -= Lx [p] * X [i+2] ;
		y [1][0] -= Lx [q] * X [i] ;
		y [1][1] -= Lx [q] * X [i+1] ;
		y [1][2] -= Lx [q] * X [i+2] ;
	    }

	    q = Lp [j-1] ;

#ifdef LL
	    y [0][0] /= d [0] ;
	    y [0][1] /= d [0] ;
	    y [0][2] /= d [0] ;
	    y [1][0] = (y [1][0] - t * y [0][0]) / d [1] ;
	    y [1][1] = (y [1][1] - t * y [0][1]) / d [1] ;
	    y [1][2] = (y [1][2] - t * y [0][2]) / d [1] ;
#else
	    y [1][0] -= t * y [0][0] ;
	    y [1][1] -= t * y [0][1] ;
	    y [1][2] -= t * y [0][2] ;
#endif
	    X [3*j  ] = y [0][0] ;
	    X [3*j+1] = y [0][1] ;
	    X [3*j+2] = y [0][2] ;
	    X [3*j-3] = y [1][0] ;
	    X [3*j-2] = y [1][1] ;
	    X [3*j-1] = y [1][2] ;

	    j -= 2 ;

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3][3], t [3] ;
	    int q = Lp [j-1] ;
	    int r = Lp [j-2] ;
#ifdef DIAG
	    double d [3] ;
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
	    d [2] = Lx [r] ;
#endif
	    t [0] = Lx [q+1] ;
	    t [1] = Lx [r+1] ;
	    t [2] = Lx [r+2] ;
#ifdef LD
	    y [0][0] = X [3*j  ] / d [0] ;
	    y [0][1] = X [3*j+1] / d [0] ;
	    y [0][2] = X [3*j+2] / d [0] ;
	    y [1][0] = X [3*j-3] / d [1] ;
	    y [1][1] = X [3*j-2] / d [1] ;
	    y [1][2] = X [3*j-1] / d [1] ;
	    y [2][0] = X [3*j-6] / d [2] ;
	    y [2][1] = X [3*j-5] / d [2] ;
	    y [2][2] = X [3*j-4] / d [2] ;
#else
	    y [0][0] = X [3*j  ] ;
	    y [0][1] = X [3*j+1] ;
	    y [0][2] = X [3*j+2] ;
	    y [1][0] = X [3*j-3] ;
	    y [1][1] = X [3*j-2] ;
	    y [1][2] = X [3*j-1] ;
	    y [2][0] = X [3*j-6] ;
	    y [2][1] = X [3*j-5] ;
	    y [2][2] = X [3*j-4] ;
#endif

	    for (p++, q += 2, r += 3 ; p < pend ; p++, q++, r++)
	    {
		int i = 3 * Li [p] ;
		y [0][0] -= Lx [p] * X [i] ;
		y [0][1] -= Lx [p] * X [i+1] ;
		y [0][2] -= Lx [p] * X [i+2] ;
		y [1][0] -= Lx [q] * X [i] ;
		y [1][1] -= Lx [q] * X [i+1] ;
		y [1][2] -= Lx [q] * X [i+2] ;
		y [2][0] -= Lx [r] * X [i] ;
		y [2][1] -= Lx [r] * X [i+1] ;
		y [2][2] -= Lx [r] * X [i+2] ;
	    }

#ifdef LL
	    y [0][0] /= d [0] ;
	    y [0][1] /= d [0] ;
	    y [0][2] /= d [0] ;
	    y [1][0] = (y [1][0] - t [0] * y [0][0]) / d [1] ;
	    y [1][1] = (y [1][1] - t [0] * y [0][1]) / d [1] ;
	    y [1][2] = (y [1][2] - t [0] * y [0][2]) / d [1] ;
	    y [2][0] = (y [2][0] - t [2] * y [0][0] - t [1] * y [1][0]) / d [2];
	    y [2][1] = (y [2][1] - t [2] * y [0][1] - t [1] * y [1][1]) / d [2];
	    y [2][2] = (y [2][2] - t [2] * y [0][2] - t [1] * y [1][2]) / d [2];
#else
	    y [1][0] -= t [0] * y [0][0] ;
	    y [1][1] -= t [0] * y [0][1] ;
	    y [1][2] -= t [0] * y [0][2] ;
	    y [2][0] -= t [2] * y [0][0] + t [1] * y [1][0] ;
	    y [2][1] -= t [2] * y [0][1] + t [1] * y [1][1] ;
	    y [2][2] -= t [2] * y [0][2] + t [1] * y [1][2] ;
#endif

	    X [3*j  ] = y [0][0] ;
	    X [3*j+1] = y [0][1] ;
	    X [3*j+2] = y [0][2] ;
	    X [3*j-3] = y [1][0] ;
	    X [3*j-2] = y [1][1] ;
	    X [3*j-1] = y [1][2] ;
	    X [3*j-6] = y [2][0] ;
	    X [3*j-5] = y [2][1] ;
	    X [3*j-4] = y [2][2] ;

	    j -= 3 ;
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

    for (j = n-1 ; j >= 0 ; )
    {
	/* get the start, end, and length of column j */
	int p, pend, lnz ;
	PEND (j, p, pend, lnz) ;

	/* find a chain of supernodes (up to j, j+1, and j+2) */
	if (j == 0 || lnz != LNZ (j-1) - 1 || Li [Lp [j-1]+1] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y [4] ;
#ifdef DIAG
	    double d = Lx [p] ;
#endif
#ifdef LD
	    y [0] = X [4*j  ] / d ;
	    y [1] = X [4*j+1] / d ;
	    y [2] = X [4*j+2] / d ;
	    y [3] = X [4*j+3] / d ;
#else
	    y [0] = X [4*j  ] ;
	    y [1] = X [4*j+1] ;
	    y [2] = X [4*j+2] ;
	    y [3] = X [4*j+3] ;
#endif

	    for (p++ ; p < pend ; p++)
	    {
		int i = 4 * Li [p] ;
		y [0] -= Lx [p] * X [i  ] ;
		y [1] -= Lx [p] * X [i+1] ;
		y [2] -= Lx [p] * X [i+2] ;
		y [3] -= Lx [p] * X [i+3] ;
	    }

#ifdef LL
	    X [4*j  ] = y [0] / d ;
	    X [4*j+1] = y [1] / d ;
	    X [4*j+2] = y [2] / d ;
	    X [4*j+3] = y [3] / d ;
#else
	    X [4*j  ] = y [0] ;
	    X [4*j+1] = y [1] ;
	    X [4*j+2] = y [2] ;
	    X [4*j+3] = y [3] ;
#endif

	    j-- ;

	}
	else if (j == 1 || lnz != LNZ (j-2)-2 || Li [Lp [j-2]+2] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2][4], t ;
	    int q = Lp [j-1] ;
#ifdef DIAG
	    double d [2] ;
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
#endif
	    t = Lx [q+1] ;
#ifdef LD
	    y [0][0] = X [4*j  ] / d [0] ;
	    y [0][1] = X [4*j+1] / d [0] ;
	    y [0][2] = X [4*j+2] / d [0] ;
	    y [0][3] = X [4*j+3] / d [0] ;
	    y [1][0] = X [4*j-4] / d [1] ;
	    y [1][1] = X [4*j-3] / d [1] ;
	    y [1][2] = X [4*j-2] / d [1] ;
	    y [1][3] = X [4*j-1] / d [1] ;
#else
	    y [0][0] = X [4*j  ] ;
	    y [0][1] = X [4*j+1] ;
	    y [0][2] = X [4*j+2] ;
	    y [0][3] = X [4*j+3] ;
	    y [1][0] = X [4*j-4] ;
	    y [1][1] = X [4*j-3] ;
	    y [1][2] = X [4*j-2] ;
	    y [1][3] = X [4*j-1] ;
#endif

	    for (p++, q += 2 ; p < pend ; p++, q++)
	    {
		int i = 4 * Li [p] ;
		y [0][0] -= Lx [p] * X [i] ;
		y [0][1] -= Lx [p] * X [i+1] ;
		y [0][2] -= Lx [p] * X [i+2] ;
		y [0][3] -= Lx [p] * X [i+3] ;
		y [1][0] -= Lx [q] * X [i] ;
		y [1][1] -= Lx [q] * X [i+1] ;
		y [1][2] -= Lx [q] * X [i+2] ;
		y [1][3] -= Lx [q] * X [i+3] ;
	    }

#ifdef LL
	    y [0][0] /= d [0] ;
	    y [0][1] /= d [0] ;
	    y [0][2] /= d [0] ;
	    y [0][3] /= d [0] ;
	    y [1][0] = (y [1][0] - t * y [0][0]) / d [1] ;
	    y [1][1] = (y [1][1] - t * y [0][1]) / d [1] ;
	    y [1][2] = (y [1][2] - t * y [0][2]) / d [1] ;
	    y [1][3] = (y [1][3] - t * y [0][3]) / d [1] ;
#else
	    y [1][0] -= t * y [0][0] ;
	    y [1][1] -= t * y [0][1] ;
	    y [1][2] -= t * y [0][2] ;
	    y [1][3] -= t * y [0][3] ;
#endif
	    X [4*j  ] = y [0][0] ;
	    X [4*j+1] = y [0][1] ;
	    X [4*j+2] = y [0][2] ;
	    X [4*j+3] = y [0][3] ;
	    X [4*j-4] = y [1][0] ;
	    X [4*j-3] = y [1][1] ;
	    X [4*j-2] = y [1][2] ;
	    X [4*j-1] = y [1][3] ;

	    j -= 2 ;

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3][4], t [3] ;
	    int q = Lp [j-1] ;
	    int r = Lp [j-2] ;
#ifdef DIAG
	    double d [3] ;
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
	    d [2] = Lx [r] ;
#endif
	    t [0] = Lx [q+1] ;
	    t [1] = Lx [r+1] ;
	    t [2] = Lx [r+2] ;
#ifdef LD
	    y [0][0] = X [4*j  ] / d [0] ;
	    y [0][1] = X [4*j+1] / d [0] ;
	    y [0][2] = X [4*j+2] / d [0] ;
	    y [0][3] = X [4*j+3] / d [0] ;
	    y [1][0] = X [4*j-4] / d [1] ;
	    y [1][1] = X [4*j-3] / d [1] ;
	    y [1][2] = X [4*j-2] / d [1] ;
	    y [1][3] = X [4*j-1] / d [1] ;
	    y [2][0] = X [4*j-8] / d [2] ;
	    y [2][1] = X [4*j-7] / d [2] ;
	    y [2][2] = X [4*j-6] / d [2] ;
	    y [2][3] = X [4*j-5] / d [2] ;
#else
	    y [0][0] = X [4*j  ] ;
	    y [0][1] = X [4*j+1] ;
	    y [0][2] = X [4*j+2] ;
	    y [0][3] = X [4*j+3] ;
	    y [1][0] = X [4*j-4] ;
	    y [1][1] = X [4*j-3] ;
	    y [1][2] = X [4*j-2] ;
	    y [1][3] = X [4*j-1] ;
	    y [2][0] = X [4*j-8] ;
	    y [2][1] = X [4*j-7] ;
	    y [2][2] = X [4*j-6] ;
	    y [2][3] = X [4*j-5] ;
#endif

	    for (p++, q += 2, r += 3 ; p < pend ; p++, q++, r++)
	    {
		int i = 4 * Li [p] ;
		y [0][0] -= Lx [p] * X [i] ;
		y [0][1] -= Lx [p] * X [i+1] ;
		y [0][2] -= Lx [p] * X [i+2] ;
		y [0][3] -= Lx [p] * X [i+3] ;
		y [1][0] -= Lx [q] * X [i] ;
		y [1][1] -= Lx [q] * X [i+1] ;
		y [1][2] -= Lx [q] * X [i+2] ;
		y [1][3] -= Lx [q] * X [i+3] ;
		y [2][0] -= Lx [r] * X [i] ;
		y [2][1] -= Lx [r] * X [i+1] ;
		y [2][2] -= Lx [r] * X [i+2] ;
		y [2][3] -= Lx [r] * X [i+3] ;
	    }

#ifdef LL
	    y [0][0] /= d [0] ;
	    y [0][1] /= d [0] ;
	    y [0][2] /= d [0] ;
	    y [0][3] /= d [0] ;
	    y [1][0] = (y [1][0] - t [0] * y [0][0]) / d [1] ;
	    y [1][1] = (y [1][1] - t [0] * y [0][1]) / d [1] ;
	    y [1][2] = (y [1][2] - t [0] * y [0][2]) / d [1] ;
	    y [1][3] = (y [1][3] - t [0] * y [0][3]) / d [1] ;
	    y [2][0] = (y [2][0] - t [2] * y [0][0] - t [1] * y [1][0]) / d [2];
	    y [2][1] = (y [2][1] - t [2] * y [0][1] - t [1] * y [1][1]) / d [2];
	    y [2][2] = (y [2][2] - t [2] * y [0][2] - t [1] * y [1][2]) / d [2];
	    y [2][3] = (y [2][3] - t [2] * y [0][3] - t [1] * y [1][3]) / d [2];
#else
	    y [1][0] -= t [0] * y [0][0] ;
	    y [1][1] -= t [0] * y [0][1] ;
	    y [1][2] -= t [0] * y [0][2] ;
	    y [1][3] -= t [0] * y [0][3] ;
	    y [2][0] -= t [2] * y [0][0] + t [1] * y [1][0] ;
	    y [2][1] -= t [2] * y [0][1] + t [1] * y [1][1] ;
	    y [2][2] -= t [2] * y [0][2] + t [1] * y [1][2] ;
	    y [2][3] -= t [2] * y [0][3] + t [1] * y [1][3] ;
#endif

	    X [4*j  ] = y [0][0] ;
	    X [4*j+1] = y [0][1] ;
	    X [4*j+2] = y [0][2] ;
	    X [4*j+3] = y [0][3] ;
	    X [4*j-4] = y [1][0] ;
	    X [4*j-3] = y [1][1] ;
	    X [4*j-2] = y [1][2] ;
	    X [4*j-1] = y [1][3] ;
	    X [4*j-8] = y [2][0] ;
	    X [4*j-7] = y [2][1] ;
	    X [4*j-6] = y [2][2] ;
	    X [4*j-5] = y [2][3] ;

	    j -= 3 ;
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
