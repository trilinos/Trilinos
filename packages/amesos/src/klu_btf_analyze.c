/* ========================================================================== */
/* === klu_btf_analyze ====================================================== */
/* ========================================================================== */

/* Order the matrix using BTF, and then AMD on the blocks */

#include "klu_btf.h"
#include "klu_kernel.h"
#include "cbtf.h"
#include "klu_dump.h"

/* ========================================================================== */
/* === klu_btf_analyze2 ===================================================== */
/* ========================================================================== */

/* klu_btf_analyze2 is not-user callable.  See klu_btf_analyze below */

static int klu_btf_analyze2	/* returns max nz in any block or < 0 if error*/
(
    /* inputs, not modified */
    int n,		/* A is n-by-n */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    int nblocks,	/* # of blocks */
    int Pbtf [ ],	/* BTF row permutation */
    int Qbtf [ ],	/* BTF col permutation */
    int R [ ],		/* size n+1, but only Rbtf [0..nblocks] is used */

    /* output only, not defined on input */
    int P [ ],		/* size n */
    int Q [ ],		/* size n */
    double Lnz [ ],	/* size n, but only Lnz [0..nblocks-1] is used */
    int *p_maxnz,
    int *p_nzoff,
    double *p_lnz,

    /* workspace, not defined on input or output */
    int Pamd [ ],	/* size maxblock */
    int Cp [ ],		/* size maxblock+1 */
    int Ci [ ],		/* size maxnz */
    int Ep [ ],		/* size maxblock+1 */
    int Ei [ ],		/* size maxnz */
    int Pinv [ ]	/* size maxblock */
)
{
    int *RowCount, *W, k1, k2, nk, k, block, oldcol, pend, row, oldrow, newcol,
	result, pc, p, newrow, col, maxnz, nzoff ;
    double Info [AMD_INFO], lnz, lnz1 ;

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    /* use Pamd as workspace for RowCount and W */
    RowCount = Pamd ;
    W = Pamd ;

    /* compute the inverse of Pbtf */
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++)
    {
	P [k] = EMPTY ;
	Q [k] = EMPTY ;
	Pinv [k] = EMPTY ;
    }
#endif
    for (k = 0 ; k < n ; k++)
    {
	ASSERT (Pbtf [k] >= 0 && Pbtf [k] < n) ;
	Pinv [Pbtf [k]] = k ;
    }
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++) ASSERT (Pinv [k] != EMPTY) ;
#endif
    nzoff = 0 ;
    lnz = 0 ;
    maxnz = 1 ;

    /* ---------------------------------------------------------------------- */
    /* order each block using AMD */
    /* ---------------------------------------------------------------------- */

    for (block = 0 ; block < nblocks ; block++)
    {

	/* ------------------------------------------------------------------ */
	/* the block is from rows/columns k1 to k2-1 */
	/* ------------------------------------------------------------------ */

	k1 = R [block] ;
	k2 = R [block+1] ;
	nk = k2 - k1 ;
	PRINTF (("BLOCK %d, k1 %d k2-1 %d nk %d\n", block, k1, k2-1, nk)) ;

	if (nk == 1)
	{

	    /* -------------------------------------------------------------- */
	    /* singleton case */
	    /* -------------------------------------------------------------- */

	    Lnz [block] = 1 ;
	    P [k1] = Pbtf [k1] ;
	    Q [k1] = Qbtf [k1] ;
	    oldcol = Q [k1] ;
	    nzoff += ((Ap [oldcol+1] - Ap [oldcol]) - 1) ;
	    PRINTF (("nzoff so far %d\n", nzoff)) ;

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* construct the kth block, C */
	    /* -------------------------------------------------------------- */

	    for (row = 0 ; row < nk ; row++)
	    {
		RowCount [row] = 0 ;
	    }
	    pc = 0 ;
	    for (k = k1 ; k < k2 ; k++)
	    {
		newcol = k-k1 ;
		Cp [newcol] = pc ;
		oldcol = Qbtf [k] ;
		pend = Ap [oldcol+1] ;
		for (p = Ap [oldcol] ; p < pend ; p++)
		{
		    oldrow = Ai [p] ;
		    newrow = Pinv [oldrow] ;
		    if (newrow < k1)
		    {
			nzoff++ ;
		    }
		    else
		    {
			/* (newrow,newcol) is an entry in the block */
			ASSERT (newrow < k2) ;
			newrow -= k1 ;
			Ci [pc++] = newrow ;
			RowCount [newrow]++ ;
		    }
		}
	    }
	    Cp [nk] = pc ;
	    maxnz = MAX (maxnz, pc) ;
	    ASSERT (klu_valid (nk, Cp, Ci, (double *) NULL)) ;

	    /* since AMD requires a sorted input, compute E = C' */
	    Ep [0] = 0 ;
	    for (row = 0 ; row < nk ; row++)
	    {
		Ep [row+1] = Ep [row] + RowCount [row] ;
	    }
	    ASSERT (Ep [nk] == pc) ;
	    for (row = 0 ; row < nk ; row++)
	    {
		W [row] = Ep [row] ;
	    }
	    for (col = 0 ; col < nk ; col++)
	    {
		pend = Cp [col+1] ;
		for (p = Cp [col] ; p < pend ; p++)
		{
		    ASSERT (W [Ci [p]] >= 0 && W [Ci [p]] < pc) ;
		    Ei [W [Ci [p]]++] = col ;
		}
	    }
	    ASSERT (klu_valid (nk, Ep, Ei, (double *) NULL)) ;

	    /* -------------------------------------------------------------- */
	    /* order the block (AMD computes E+E'), using default parameters */
	    /* -------------------------------------------------------------- */

	    PRINTF (("calling AMD\n")) ;
	    result = amd_order (nk, Ep, Ei, Pamd, (double *) NULL, Info) ;
	    PRINTF (("AMD done\n")) ;
	    if (result == AMD_OUT_OF_MEMORY)
	    {
		return (KLU_OUT_OF_MEMORY) ;
	    }
	    else if (result == AMD_INVALID || ((int) (Info [AMD_NZDIAG]) != nk))
	    {
		/* fatal error - something is corrupted in the matrix E */
		PRINTF (("AMD invalid!\n")) ;
		return (-3) ;
	    }

	    /* get the ordering statistics from AMD */
	    lnz1 = Info [AMD_LNZ] + nk ;
	    Lnz [block] = lnz1 ;
	    lnz += lnz1 ;

	    /* -------------------------------------------------------------- */
	    /* combine the AMD ordering with the BTF ordering */
	    /* -------------------------------------------------------------- */

	    PRINTF (("Pamd, 1-based:\n")) ;
	    for (k = 0 ; k < nk ; k++)
	    {
		ASSERT (k + k1 < n) ;
		ASSERT (Pamd [k] + k1 < n) ;
		Q [k + k1] = Qbtf [Pamd [k] + k1] ;
	    }
	    for (k = 0 ; k < nk ; k++)
	    {
		ASSERT (k + k1 < n) ;
		ASSERT (Pamd [k] + k1 < n) ;
		P [k + k1] = Pbtf [Pamd [k] + k1] ;
	    }
	}
    }

    ASSERT (nzoff >= 0 && nzoff <= Ap [n] - n) ;

    /* return estimates of # of nonzeros in L including diagonal */
    *p_lnz = lnz ;

    *p_maxnz = maxnz ;
    *p_nzoff = nzoff ;

    return (KLU_OK) ;
}

/* ========================================================================== */
/* === klu_btf_analyze ====================================================== */
/* ========================================================================== */

klu_symbolic *klu_btf_analyze	/* returns NULL if error, or a valid
				   klu_symbolic object if successful */
(
    /* inputs, not modified */
    int n,		/* A is n-by-n */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ]		/* size nz, row indices */
)
{
    int nblocks, *Qbtf, nz, *Cp, *Ci, *Qamd, *Pinv, *Pamd, *Ep, *Ei, *Pbtf, ok,
	block, maxblock, k1, k2, nk, maxnz, *P, *Q, *R, nzoff, result ;
    klu_symbolic *Symbolic ;
    double *Lnz, lnz ;

    nz = Ap [n] ;
    if (nz <= 0)
    {
	PRINTF (("klu_btf_analyze: singular\n")) ;
	return ((klu_symbolic *) NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate the Symbolic object */
    /* ---------------------------------------------------------------------- */

    Symbolic = (klu_symbolic *) ALLOCATE (sizeof (klu_symbolic)) ;
    if (Symbolic == (klu_symbolic *) NULL)
    {
	return ((klu_symbolic *) NULL) ;
    }

    P = (int *) ALLOCATE (n * sizeof (int)) ;
    Q = (int *) ALLOCATE (n * sizeof (int)) ;
    R = (int *) ALLOCATE ((n+1) * sizeof (int)) ;
    Lnz = (double *) ALLOCATE (n * sizeof (double)) ;

    Symbolic->n = n ;
    Symbolic->P = P ;
    Symbolic->Q = Q ;
    Symbolic->R = R ;
    Symbolic->Lnz = Lnz ;

    if ((P == (int *) NULL) || (Q == (int *) NULL) || (R == (int *) NULL) ||
	(Lnz == (double *) NULL))
    {
	klu_btf_free_symbolic (&Symbolic) ;
	return ((klu_symbolic *) NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace for BTF permutation */
    /* ---------------------------------------------------------------------- */

    Pbtf = (int *) ALLOCATE (n * sizeof (int)) ;
    Qbtf = (int *) ALLOCATE (n * sizeof (int)) ;
    if ((Pbtf == (int *) NULL) || (Qbtf == (int *) NULL))
    {
	FREE (Pbtf, int) ; 
	FREE (Qbtf, int) ; 
	klu_btf_free_symbolic (&Symbolic) ;
	return ((klu_symbolic *) NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* find the block triangular form */
    /* ---------------------------------------------------------------------- */

    PRINTF (("calling cbtf\n")) ;
    nblocks = cbtf (n, Ap, Ai, Pbtf, Qbtf, R) ;
    PRINTF (("cbtf nblocks %d\n", nblocks)) ;
    if (nblocks == 0)
    {
	FREE (Pbtf, int) ; 
	FREE (Qbtf, int) ; 
	klu_btf_free_symbolic (&Symbolic) ;
	return ((klu_symbolic *) NULL) ;
    }
    else if (nblocks < 0)
    {
	FREE (Pbtf, int) ; 
	FREE (Qbtf, int) ; 
	klu_btf_free_symbolic (&Symbolic) ;
	return ((klu_symbolic *) NULL) ;
    }
    Symbolic->nblocks = nblocks ;

    /* ---------------------------------------------------------------------- */
    /* find the size of the largest block */
    /* ---------------------------------------------------------------------- */

    maxblock = 1 ;
    for (block = 0 ; block < nblocks ; block++)
    {
	k1 = R [block] ;
	k2 = R [block+1] ;
	nk = k2 - k1 ;
	PRINTF (("block %d size %d\n", block, nk)) ;
	maxblock = MAX (maxblock, nk) ;
    }
    /* TODO: Ci and Ei could be of size max (nnz of blocks of A) */
    maxnz = nz ;
    PRINTF (("maxblock size %d maxnz %d\n", maxblock, maxnz)) ;
    Symbolic->maxblock = maxblock ;

    /* TODO: merge adjacent 1-by-1 blocks into an upper triangular block */

    /* ---------------------------------------------------------------------- */
    /* allocate more workspace, for klu_btf_analyze2 */
    /* ---------------------------------------------------------------------- */

    Pamd = (int *) ALLOCATE (maxblock * sizeof (int)) ;
    Cp = (int *) ALLOCATE ((maxblock+1) * sizeof (int)) ;
    Ep = (int *) ALLOCATE ((maxblock+1) * sizeof (int)) ;
    Ci = (int *) ALLOCATE ((maxnz+1) * sizeof (int)) ;
    Ei = (int *) ALLOCATE ((maxnz+1) * sizeof (int)) ;
    Pinv = (int *) ALLOCATE (n * sizeof (int)) ;

    ok = (Pamd != (int *) NULL) &&
	(Cp != (int *) NULL) &&
	(Ci != (int *) NULL) &&
	(Ep != (int *) NULL) &&
	(Ei != (int *) NULL) &&
	(Pinv != (int *) NULL) ;

    /* ---------------------------------------------------------------------- */
    /* order each block of the BTF matrix using AMD */
    /* ---------------------------------------------------------------------- */

    if (ok)
    {
	PRINTF (("calling klu_btf_analyze2\n")) ;
	result = klu_btf_analyze2 (n, Ap, Ai, nblocks, Pbtf, Qbtf,
	    R, P, Q, Lnz, &maxnz, &nzoff, &lnz, Pamd, Cp, Ci, Ep, Ei, Pinv) ;
	PRINTF (("klu_btf_analyze2 done\n")) ;
	ok = (result == KLU_OK) ;
    }
    Symbolic->lnz = lnz ;
    Symbolic->unz = lnz ;	/* estimated fill-in is symmetric, currently */
    Symbolic->maxnz = maxnz ;
    Symbolic->nzoff = nzoff ;

    /* ---------------------------------------------------------------------- */
    /* free all workspace */
    /* ---------------------------------------------------------------------- */

    FREE (Pamd, int) ;
    FREE (Cp, int) ;
    FREE (Ci, int) ;
    FREE (Ep, int) ;
    FREE (Ei, int) ;
    FREE (Pinv, int) ;
    FREE (Pbtf, int) ;
    FREE (Qbtf, int) ;

    PRINTF (("klu_btf_analyze done\n")) ;

    /* ---------------------------------------------------------------------- */
    /* return the symbolic object */
    /* ---------------------------------------------------------------------- */

    if (!ok)
    {
	klu_btf_free_symbolic (&Symbolic) ;
	return ((klu_symbolic *) NULL) ;
    }
    else
    {
	return (Symbolic) ;
    }
}
