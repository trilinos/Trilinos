/* ========================================================================== */
/* === TRILINOS_KLU_sort ============================================================= */
/* ========================================================================== */

/* sorts the columns of L and U so that the row indices appear in strictly
 * increasing order.
 */

#include "trilinos_klu_internal.h"

/* ========================================================================== */
/* === sort ================================================================= */
/* ========================================================================== */

/* Sort L or U using a double-transpose */

static void sort (Int n, Int *Xip, Int *Xlen, Unit *LU, Int *Tp, Int *Tj,
    Entry *Tx, Int *W)
{
    Int *Xi ;
    Entry *Xx ;
    Int p, i, j, len, nz, tp, xlen, pend ;

    ASSERT (TRILINOS_KLU_valid_LU (n, FALSE, Xip, Xlen, LU)) ;

    /* count the number of entries in each row of L or U */ 
    for (i = 0 ; i < n ; i++)
    {
	W [i] = 0 ;
    }
    for (j = 0 ; j < n ; j++)
    {
	GET_POINTER (LU, Xip, Xlen, Xi, Xx, j, len) ;
	for (p = 0 ; p < len ; p++)
	{
	    W [Xi [p]]++ ;
	}
    }

    /* construct the row pointers for T */
    nz = 0 ;
    for (i = 0 ; i < n ; i++)
    {
	Tp [i] = nz ;
	nz += W [i] ;
    }
    Tp [n] = nz ;
    for (i = 0 ; i < n ; i++)
    {
	W [i] = Tp [i] ;
    }

    /* transpose the matrix into Tp, Ti, Tx */
    for (j = 0 ; j < n ; j++)
    {
	GET_POINTER (LU, Xip, Xlen, Xi, Xx, j, len) ;
	for (p = 0 ; p < len ; p++)
	{
	    tp = W [Xi [p]]++ ;
	    Tj [tp] = j ;
	    Tx [tp] = Xx [p] ;
	}
    }

    /* transpose the matrix back into Xip, Xlen, Xi, Xx */
    for (j = 0 ; j < n ; j++)
    {
	W [j] = 0 ;
    }
    for (i = 0 ; i < n ; i++)
    {
	pend = Tp [i+1] ;
	for (p = Tp [i] ; p < pend ; p++)
	{
	    j = Tj [p] ;
	    GET_POINTER (LU, Xip, Xlen, Xi, Xx, j, len) ;
	    xlen = W [j]++ ;
	    Xi [xlen] = i ;
	    Xx [xlen] = Tx [p] ;
	}
    }

    ASSERT (TRILINOS_KLU_valid_LU (n, FALSE, Xip, Xlen, LU)) ;
}


/* ========================================================================== */
/* === TRILINOS_KLU_sort ============================================================= */
/* ========================================================================== */

Int TRILINOS_KLU_sort
(
    TRILINOS_KLU_symbolic *Symbolic,
    TRILINOS_KLU_numeric *Numeric,
    TRILINOS_KLU_common *Common
)
{
    Int *R, *W, *Tp, *Ti, *Lip, *Uip, *Llen, *Ulen ;
    Entry *Tx ;
    Unit **LUbx ;
    Int n, nk, nz, block, nblocks, maxblock, k1 ;
    size_t m1 ;

    if (Common == NULL)
    {
	return (FALSE) ;
    }
    Common->status = TRILINOS_KLU_OK ;

    n = Symbolic->n ;
    R = Symbolic->R ;
    nblocks = Symbolic->nblocks ;
    maxblock = Symbolic->maxblock ;

    Lip  = Numeric->Lip ;
    Llen = Numeric->Llen ;
    Uip  = Numeric->Uip ;
    Ulen = Numeric->Ulen ;
    LUbx = (Unit **) Numeric->LUbx ;

    m1 = ((size_t) maxblock) + 1 ;

    /* allocate workspace */
    nz = MAX (Numeric->max_lnz_block, Numeric->max_unz_block) ;
    W  = (Int*) TRILINOS_KLU_malloc (maxblock, sizeof (Int), Common) ;
    Tp = (Int*) TRILINOS_KLU_malloc (m1, sizeof (Int), Common) ;
    Ti = (Int*) TRILINOS_KLU_malloc (nz, sizeof (Int), Common) ;
    Tx = (Entry*) TRILINOS_KLU_malloc (nz, sizeof (Entry), Common) ;

    PRINTF (("\n======================= Start sort:\n")) ;

    if (Common->status == TRILINOS_KLU_OK)
    {
	/* sort each block of L and U */
	for (block = 0 ; block < nblocks ; block++)
	{
	    k1 = R [block] ;
	    nk = R [block+1] - k1 ;
	    if (nk > 1)
	    {
		PRINTF (("\n-------------------block: %d nk %d\n", block, nk)) ;
		sort (nk, Lip + k1, Llen + k1, LUbx [block], Tp, Ti, Tx, W) ;
		sort (nk, Uip + k1, Ulen + k1, LUbx [block], Tp, Ti, Tx, W) ;
	    }
	}
    }

    PRINTF (("\n======================= sort done.\n")) ;

    /* free workspace */
    TRILINOS_KLU_free (W, maxblock, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Tp, m1, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Ti, nz, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Tx, nz, sizeof (Entry), Common) ;
    return (Common->status == TRILINOS_KLU_OK) ;
}
